""" es_runners for conformer
"""

import shutil
import time
import random
import subprocess
import os
import numpy

import automol
import elstruct
import autofile
from autofile import fs
from mechanalyzer.inf import thy as tinfo
from mechlib import filesys
from mechlib.amech_io.printer import info_message, warning_message
from mechlib.amech_io.printer import debug_message, error_message, obj
from mechlib.amech_io.printer import existing_path, bad_conformer, checking
from mechroutines.es import runner as es_runner
from mechroutines.es._routines import _util as util
from mechroutines.es._routines._geom import remove_imag

from automol.extern import Ring_Reconstruction as RR
from phydat import phycon

# Initial conformer
def initial_conformer(spc_dct_i, spc_info, ini_method_dct, method_dct,
                      ini_cnf_save_fs, cnf_run_fs, cnf_save_fs,
                      es_keyword_dct):
    """ Assess if a conformer layer with a geometry exists in the save
        filesys for the given species.
        If not, attempt to generate some guess structure using InChI strings
        or input geom from user.
        and optimize
        it with input method. Then assess if the optimized structure
        corresponds to genuine minimum on the PES via a frequency calculation.
        If a minimum is found, save the conformer geometry, zmatrix, energy,
        and torsions to the save filesys.
        Also, the function assessess if the species is unstable and will
        save the appropriate information.
    """
    ini_thy_info = tinfo.from_dct(ini_method_dct)
    thy_info = tinfo.from_dct(method_dct)
    mod_thy_info = tinfo.modify_orb_label(
        thy_info, spc_info)
    mod_ini_thy_info = tinfo.modify_orb_label(
        ini_thy_info, spc_info)
    [kickoff_size, kickoff_backward] = spc_dct_i['kickoff']

    _, cnf_path = filesys.mincnf.min_energy_conformer_locators(
        cnf_save_fs, mod_thy_info)
    overwrite = es_keyword_dct['overwrite']
    if not cnf_path:
        info_message(
            'No conformer found in save filesys. Checking for running jobs...')
        if _init_geom_is_running(cnf_run_fs) and not overwrite:
            _run = False
        else:
            info_message(
                'No conformers are running in run filesys.' +
                'Proceeding with optimization...')
            _run = True
    elif overwrite:
        info_message(
            'User specified to overwrite energy with new run...')
        _run = True
    else:
        _run = False

    if _run:
        info_message('Obtaining some initial guess geometry.')
        geo_init = _obtain_ini_geom(spc_dct_i, ini_cnf_save_fs,
                                    mod_ini_thy_info,
                                    overwrite)

        if geo_init is not None:
            info_message(
                'Assessing if there are any functional groups',
                'that cause instability')

            zma_init = automol.geom.zmatrix(geo_init)

            rid = autofile.schema.generate_new_ring_id()
            cid = autofile.schema.generate_new_conformer_id()

            # Determine if there is an instability, if so return prods
            instab_zmas = automol.reac.instability_product_zmas(zma_init)
            if not instab_zmas:

                # Build a cid and a run fs
                cnf_run_fs[-1].create((rid, cid))
                run_fs = autofile.fs.run(cnf_run_fs[-1].path((rid, cid)))

                if not automol.geom.is_atom(geo_init):
                    geo_found = _optimize_molecule(
                        spc_info, zma_init,
                        method_dct,
                        cnf_save_fs, (rid, cid),
                        run_fs,
                        overwrite,
                        kickoff_size=kickoff_size,
                        kickoff_backward=kickoff_backward)
                else:
                    geo_found = _optimize_atom(
                        spc_info, zma_init,
                        method_dct,
                        cnf_save_fs, (rid, cid),
                        run_fs,
                        overwrite)
            else:
                info_message(
                    'Found functional groups that cause instabilities')
                filesys.save.instability(
                    zma_init, instab_zmas, cnf_save_fs,
                    rng_locs=(rid,), tors_locs=(cid,), zma_locs=(0,))
                geo_found = True
        else:
            geo_found = False
            warning_message(
                'Unable to obtain an initial guess geometry')
    else:
        existing_path('Initial geometry', cnf_path)
        geo_found = True

    return geo_found


def _obtain_ini_geom(spc_dct_i, ini_cnf_save_fs,
                     mod_ini_thy_info, overwrite):
    """ Obtain an initial geometry to be optimized. Checks a hieratchy
        of places to obtain the initial geom.
            (1) Geom dict which is the input from the user
            (2) Geom from inchi
    """

    geo_init = None
    # Obtain geom from thy fs or remove the conformer filesystem if needed
    if not overwrite:
        ini_min_locs, ini_path = filesys.mincnf.min_energy_conformer_locators(
            ini_cnf_save_fs, mod_ini_thy_info)
        if ini_path:
            geo_init = ini_cnf_save_fs[-1].file.geometry.read(ini_min_locs)
            path = ini_cnf_save_fs[-1].path(ini_min_locs)
            info_message(
                f'Getting inital geometry from inplvl at path {path}')
    else:
        debug_message(
            'Removing original conformer save data for instability')
        for locs in ini_cnf_save_fs[-1].existing():
            cnf_save_path = ini_cnf_save_fs[-1].path(locs)
            debug_message(f'Removing {cnf_save_path}')
            shutil.rmtree(cnf_save_path)

    if geo_init is None:
        if 'geo' in spc_dct_i:
            geo_init = spc_dct_i['geo']
            info_message(
                'Getting initial geometry from geom dictionary')

    if geo_init is None:
        geo_init = automol.chi.geometry(spc_dct_i['canon_enant_ich'])
        info_message('Getting initial geometry from inchi')

    # Check if the init geometry is connected
    if geo_init is not None:
        if not automol.geom.is_connected(geo_init):
            geo_init = None

    return geo_init


def _optimize_atom(spc_info, zma_init,
                   method_dct,
                   cnf_save_fs, locs,
                   run_fs,
                   overwrite):
    """ Deal with an atom separately
    """

    thy_info = tinfo.from_dct(method_dct)
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
    script_str, kwargs = es_runner.qchem_params(
        method_dct, spc_info=spc_info,
        geo=automol.zmat.geometry(zma_init))

    # Call the electronic structure optimizer
    success, ret = es_runner.execute_job(
        job=elstruct.Job.ENERGY,
        script_str=script_str,
        run_fs=run_fs,
        geo=zma_init,
        spc_info=spc_info,
        thy_info=mod_thy_info,
        zrxn=None,
        overwrite=overwrite,
        **kwargs
    )

    if success:
        info_message('Succesful reference geometry optimization')
        filesys.save.atom(
            ret, cnf_save_fs, mod_thy_info[1:], zma_init,
            rng_locs=(locs[0],), tors_locs=(locs[1],))

    return success


def _optimize_molecule(spc_info, zma_init,
                       method_dct,
                       cnf_save_fs, locs,
                       run_fs,
                       overwrite,
                       kickoff_size=0.1, kickoff_backward=False):
    """ Optimize a proper geometry
    """
    thy_info = tinfo.from_dct(method_dct)
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
    script_str, kwargs = es_runner.qchem_params(
        method_dct,
        spc_info=spc_info,
        geo=automol.zmat.geometry(zma_init),
        job=elstruct.Job.OPTIMIZATION)

    # Call the electronic structure optimizer
    success, ret = es_runner.execute_job(
        job=elstruct.Job.OPTIMIZATION,
        script_str=script_str,
        run_fs=run_fs,
        geo=zma_init,
        spc_info=spc_info,
        thy_info=mod_thy_info,
        zrxn=None,
        overwrite=overwrite,
        **kwargs
    )

    # read the geometry
    geo_conn = False
    if success:
        inf_obj, _, out_str = ret
        geo = elstruct.reader.opt_geometry(inf_obj.prog, out_str)
        zma = elstruct.reader.opt_zmatrix(inf_obj.prog, out_str)
        if zma is None:
            zma = automol.geom.zmatrix(geo)
        geo_conn = bool(automol.geom.is_connected(geo))

    # If connected, check for imaginary modes and fix them if possible
    conf_found = False
    if geo_conn and success:

        # Remove the imaginary mode
        geo, ret = remove_imag(
            geo, ret, spc_info, method_dct, run_fs,
            kickoff_size=kickoff_size,
            kickoff_backward=kickoff_backward,
            kickoff_mode=0)

        # Recheck connectivity for imag-checked geometry
        if geo is not None:
            conf_found = True
            conn = automol.geom.is_connected(geo)
            proper_stereo = _inchi_are_same(spc_info[0], geo)
            if conn and proper_stereo:
                info_message(
                    'Saving structure as the first conformer...', newline=1)
                filesys.save.conformer(
                    ret, None, cnf_save_fs, mod_thy_info[1:],
                    rng_locs=(locs[0],), tors_locs=(locs[1],),
                    init_zma=zma)
            else:
                if not conn:
                    info_message('Saving disconnected species...')
                    filesys.save.instability(
                        zma_init, zma, cnf_save_fs,
                        rng_locs=(locs[0],), tors_locs=(locs[1],),
                        zma_locs=(0,))
        else:
            warning_message('No geom found...', newline=1)
            conf_found = False
    elif success:
        info_message('Saving disconnected species...')
        conf_found = False
        filesys.save.instability(
            zma_init, zma, cnf_save_fs,
            rng_locs=(locs[0],), tors_locs=(locs[1],), zma_locs=(0,))

    return conf_found


def single_conformer(zma, spc_info, mod_thy_info,
                     cnf_run_fs, cnf_save_fs,
                     script_str, overwrite,
                     retryfail=True, zrxn=None,
                     use_locs=None, resave=False,
                     **kwargs):
    """ generate single optimized geometry to be saved into a
        filesystem
    """
    skip_job = False

    if resave:
        _presamp_save(
            spc_info, cnf_run_fs, cnf_save_fs,
            mod_thy_info, zrxn=zrxn, rid=None, ref_zma=zma)
        if use_locs is None:
            print('getting rid')
            rid = rng_loc_for_geo(
                automol.zmat.geometry(zma), cnf_save_fs)
            if rid is not None:
                cid = autofile.schema.generate_new_conformer_id()
                locs = (rid, cid)

    if this_conformer_is_running(zma, cnf_run_fs):
        skip_job = True
    elif this_conformer_was_run_in_save(zma, cnf_save_fs):
        skip_job = True
    if not skip_job:
        run_in_run, _ = filesys.mincnf.this_conformer_was_run_in_run(
            zma, cnf_run_fs, cnf_save_fs, mod_thy_info)
        if run_in_run:
            skip_job = True

    if not skip_job:
        # Build the filesystem
        if use_locs is None:
            rid = autofile.schema.generate_new_ring_id()
            cid = autofile.schema.generate_new_conformer_id()
            locs = (rid, cid)
        else:
            locs = use_locs
            rid = locs[0]

        cnf_run_fs[-1].create(locs)
        cnf_run_path = cnf_run_fs[-1].path(locs)
        run_fs = autofile.fs.run(cnf_run_path)

        # Run the optimization
        info_message('Optimizing a single conformer...')
        success, ret = es_runner.execute_job(
            job=elstruct.Job.OPTIMIZATION,
            script_str=script_str,
            run_fs=run_fs,
            geo=zma,
            spc_info=spc_info,
            thy_info=mod_thy_info,
            zrxn=zrxn,
            overwrite=overwrite,
            frozen_coordinates=(),
            saddle=bool(zrxn is not None),
            retryfail=retryfail,
            **kwargs
        )

        if success:
            inf_obj, _, out_str = ret
            prog = inf_obj.prog
            method = inf_obj.method
            ene = elstruct.reader.energy(prog, method, out_str)
            geo = elstruct.reader.opt_geometry(prog, out_str)
            # zma = elstruct.reader.opt_zmatrix(prog, out_str)
            saved_locs, saved_geos, saved_enes = _saved_cnf_info(
                cnf_save_fs, mod_thy_info)

            opt_zma = None
            if zma is not None:
                opt_zma = filesys.save.read_zma_from_geo(zma, geo)
            if zma is None:
                opt_zma = filesys.save.read_job_zma(ret, init_zma=zma)
            viable = _geo_connected(geo, zrxn)
            if viable:
                if zrxn:
                    viable = _ts_geo_viable(
                        opt_zma, zrxn, cnf_save_fs, mod_thy_info, ref_zma=zma)
                else:
                    viable = _inchi_are_same(spc_info[0], geo)

            if viable:
                if _geo_unique(geo, ene, saved_geos, saved_enes, zrxn=zrxn):
                    sym_id = _sym_unique(
                        geo, ene, saved_geos, saved_enes)
                    if sym_id is None:
                        if cnf_save_fs[0].file.info.exists():
                            debug_message(
                                'inf_obj path', cnf_save_fs[0].path())
                            rinf_obj = cnf_save_fs[0].file.info.read()
                            rinf = rinf_obj
                            debug_message(
                                'inf_obj for r', rinf)
                            rnsampd = rinf_obj.nsamp
                            rnsampd += 1
                            rinf.nsamp = rnsampd
                        else:
                            rinf = autofile.schema.info_objects.conformer_trunk(0)
                            rinf.nsamp = 1
                        if cnf_save_fs[1].file.info.exists([locs[0]]):
                            cinf_obj_s = cnf_save_fs[1].file.info.read([locs[0]])
                            cinf = cinf_obj_s
                            cnsampd = cinf_obj_s.nsamp
                            cnsampd += 1
                            cinf.nsamp = cnsampd
                        else:
                            cinf = autofile.schema.info_objects.conformer_branch(0)
                            cinf.nsamp = 1
                        cnf_save_fs[1].create([locs[0]])
                        cnf_save_fs[0].file.info.write(rinf)
                        cnf_save_fs[1].file.info.write(cinf, [locs[0]])
                        filesys.save.conformer(
                            ret, None, cnf_save_fs, mod_thy_info[1:],
                            zrxn=zrxn, init_zma=zma,
                            rng_locs=(locs[0],), tors_locs=(locs[1],))
                        saved_geos.append(geo)
                        saved_enes.append(ene)
                        saved_locs.append(locs)

                        # Update the conformer trajectory file
                        obj('vspace')
                        filesys.mincnf.traj_sort(
                            cnf_save_fs, mod_thy_info)
                        filesys.mincnf.traj_sort(
                            cnf_save_fs, mod_thy_info, locs[0])
                    else:
                        sym_locs = saved_locs[sym_id]
                        filesys.save.sym_indistinct_conformer(
                            geo, cnf_save_fs, locs, sym_locs)
                        if cnf_save_fs[-1].exists(locs):
                            cnf_save_path = cnf_save_fs[-1].path(locs)
                        if cnf_run_fs[-1].exists(locs):
                            cnf_run_path = cnf_run_fs[-1].path(locs)


def conformer_sampling(zma, spc_info, thy_info,
                       cnf_run_fs, cnf_save_fs, rid,
                       script_str, overwrite,
                       nsamp_par=(False, 3, 3, 1, 50, 50),
                       tors_names=(),
                       zrxn=None, two_stage=False,
                       retryfail=False, resave=False,
                       repulsion_thresh=40.0, print_debug=True,
                       **kwargs):
    """ run sampling algorithm to find conformers
    """

    # Check if any saving needs to be done before hand
    ref_rid = rid
    cnf_run_fs[1].create([rid])
    if resave:
        _presamp_save(
            spc_info, cnf_run_fs, cnf_save_fs, thy_info, zrxn=zrxn, rid=rid, ref_zma=zma)

    # Build filesys
    cnf_save_fs[1].create([rid])
    inf_obj = autofile.schema.info_objects.conformer_branch(0)

    # Set the samples
    nsamp, tors_range_dct = util.calc_nsamp(
        tors_names, nsamp_par, zma, zrxn=zrxn)
    nsamp0 = nsamp
    nsampd = util.calc_nsampd(cnf_save_fs, cnf_run_fs, rid)

    tot_samp = nsamp - nsampd
    brk_tot_samp = nsamp * 5

    info_message(
        ' - Number of samples that have been currently run:', nsampd)
    info_message(' - Number of samples requested:', nsamp)

    if nsamp-nsampd > 0:
        info_message(
            f'Running {nsamp-nsampd} samples...', newline=1)

    # Generate all of the conformers, as needed
    samp_idx = 1
    samp_attempt_idx = 1
    while True:
        nsamp = nsamp0 - nsampd
        rid = ref_rid
        # Break the while loop if enough sampls completed
        if nsamp <= 0:
            info_message(
                'Requested number of samples have been completed.',
                'Conformer search complete.')
            break
        if samp_attempt_idx == brk_tot_samp:
            info_message(
                f'Max sample num: 5*{nsamp} attempted, ending search',
                'Run again if more samples desired.')
            break

        # Run the conformer sampling
        if nsampd > 0:
            samp_zma, = automol.zmat.samples(zma, 1, tors_range_dct)
        else:
            samp_zma = zma

        info_message(
            'Generating sample Z-Matrix that does not have',
            'high intramolecular repulsion...')
        bad_geo_cnt = 0
        while not automol.zmat.has_low_relative_repulsion_energy(samp_zma, zma) and bad_geo_cnt < 1000:
            if print_debug:
                warning_message('Structure has high repulsion.')
                warning_message(
                    'Generating new sample Z-Matrix')
            samp_zma, = automol.zmat.samples(zma, 1, tors_range_dct)
            bad_geo_cnt += 1

        cid = autofile.schema.generate_new_conformer_id()
        locs = [rid, cid]

        cnf_run_fs[-1].create(locs)
        cnf_run_path = cnf_run_fs[-1].path(locs)
        run_fs = autofile.fs.run(cnf_run_path)

        info_message(f"Run {samp_idx}/{tot_samp}")
        tors_names = tuple(tors_range_dct.keys())
        print('two_stage test:', two_stage, tors_names)
        if two_stage and tors_names:
            frozen_coords_lst = (tors_names, ())
            success, ret = es_runner.multi_stage_optimization(
                script_str=script_str,
                run_fs=run_fs,
                geo=samp_zma,
                spc_info=spc_info,
                thy_info=thy_info,
                frozen_coords_lst=frozen_coords_lst,
                zrxn=zrxn,
                overwrite=overwrite,
                saddle=bool(zrxn is not None),
                retryfail=retryfail,
                **kwargs
            )
        else:
            success, ret = es_runner.execute_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                run_fs=run_fs,
                geo=samp_zma,
                spc_info=spc_info,
                thy_info=thy_info,
                zrxn=zrxn,
                overwrite=overwrite,
                saddle=bool(zrxn is not None),
                retryfail=retryfail,
                **kwargs
            )

        # save function added here
        if success:
            inf_obj_temp, _, out_str = ret
            prog = inf_obj_temp.prog
            samp_geo = elstruct.reader.opt_geometry(prog, out_str)
            # Determine ring state and update rid
            rid = rng_loc_for_geo(samp_geo, cnf_save_fs)
            if rid is None:
                rid = autofile.schema.generate_new_ring_id()
            locs = [rid, cid]
            save_conformer(
                ret, cnf_run_fs, cnf_save_fs, locs, thy_info,
                zrxn=zrxn, orig_ich=spc_info[0], rid_traj=True,
                init_zma=samp_zma, ref_zma=samp_zma)

            nsampd = util.calc_nsampd(cnf_save_fs, cnf_run_fs, rid)
            nsampd += 1
            samp_idx += 1
            inf_obj.nsamp = nsampd
            cnf_save_fs[1].file.info.write(inf_obj, [rid])
            cnf_run_fs[1].file.info.write(inf_obj, [rid])

        # Increment attempt counter
        samp_attempt_idx += 1


def ring_conformer_sampling(
        zma, spc_info, thy_info,
        cnf_run_fs, cnf_save_fs,
        script_str, overwrite,
        algorithm='torsions',
        thresholds='default',
        eps=0.2,
        checks=1,
        skip=False,
        nsamp_par=(False, 3, 1, 3, 50, 50),
        ring_tors_dct=None,
        zrxn=None, two_stage=False, retryfail=False,
        **kwargs):
    """ run sampling algorithm to find conformers
    """

    if skip: 
        print("\nRing puckering: Skipping...\n")
        return
# I can use resave in conf_samp
    # _presamp_save(
    #         spc_info, cnf_run_fs, cnf_save_fs,
    #         thy_info, zrxn=zrxn, rid=None, ref_zma=zma)
    
    import rdkit

    # Build filesys
    cnf_save_fs[0].create()
    inf_obj = autofile.schema.info_objects.conformer_trunk(0)
    # Set up torsions
    geo = automol.zmat.geometry(zma)
    tors_dcts = ring_tors_dct.items() if ring_tors_dct else {}
    rings_atoms = []
    for ring_atoms, samp_range_dct in tors_dcts:
        rings_atoms.append([int(idx)-1 for idx in ring_atoms.split('-')])
    gra = automol.geom.graph(geo)
    ngbs = automol.graph.atoms_sorted_neighbor_atom_keys(gra)
    
    # Exclude from puckering sampling the ring atoms after a higher order bond
    rot_bset = set(automol.graph.rotational_bond_keys(gra,with_rings_rotors=True))
    rings_bonds = {frozenset([ ring[i-1], atom ]
                   ) for ring in rings_atoms for i, atom in enumerate(ring)}
    rings_planar_bonds = rings_bonds.difference(rot_bset)
    rings_planar_atoms = set()
    for ring_atoms,_ in tors_dcts:
        ring_atoms = [int(idx)-1 for idx in ring_atoms.split('-')]
        for planar_bond in rings_planar_bonds:
            if all(atom in ring_atoms for atom in planar_bond):
                ring_list_index = max(planar_bond) + 1
                if ring_list_index < len(ring_atoms): 
                    rings_planar_atoms.add(ring_atoms[ring_list_index])
    # Find dihedrals of atoms in ring_planar_atoms
    coos = automol.zmat.coordinates(zma)
    planar_dih, dih_remover, remove_ring = [],[],[]
    for cord,dih_atoms in {key: value for key,value in coos.items(
                                ) if key.startswith("D")}.items():
        if dih_atoms[0][0] in rings_planar_atoms: 
            planar_dih.append(cord)
    # Check rings_tors_dct and remove dihedrals of ring_planar_atoms if present
    for ring_atoms,ring_dihs in tors_dcts:
        dih_remover.extend([(ring_atoms,dih) for dih in ring_dihs if dih in planar_dih])
    for ring_atoms,dih in dih_remover:
        del ring_tors_dct[ring_atoms][dih]
    # Check if any ring now has no dihedrals to sample, and remove it
    for ring_atoms,ring_dihs in tors_dcts:
        if ring_dihs == {}: remove_ring.append(ring_atoms) 
    for ring_atoms in remove_ring:  
        del ring_tors_dct[ring_atoms]

    # Check whether all dihedrals were removed and eventually return from functions
    if len(tors_dcts)==0: 
        print("\nRing puckering: No dihedrals to sample\n")
        return

    # Collect information on total ring atoms and create a unique dic for sampled dihs
    # Create lists for check only on ring atoms not connected in Z matrices
    all_ring_atoms, all_samp_range_dct, all_unconnected_lst = [], {}, []
    bonds_from_zma = [frozenset(value[0]) for key,value in coos.items(
                                ) if key.startswith("R")]
    for ring_atoms, samp_range_dct in tors_dcts:
        ring_atoms = [int(idx)-1 for idx in ring_atoms.split('-')]
        all_ring_atoms.extend(ring_atoms)
        all_samp_range_dct.update(samp_range_dct)
        # Find unconnected ats comparing bonds from zmat to rings_bonds
        # Compare the bonds of the ring to those of the Z-Matrix, get those present
        # only in in the former
        ring_bonds = {frozenset([ring_atoms[i-1], atom]) for i, atom in enumerate(
                        ring_atoms)}
        unconnected_ats = [list(el) for el in ring_bonds.difference(bonds_from_zma)]
        for unconnected_bond in unconnected_ats: 
            #note that there is one redundancy for fused rings
            all_unconnected_lst.append((unconnected_bond, 
                                    automol.zmat.ring_distances(zma, unconnected_bond)))
    all_ring_atoms = list(set(all_ring_atoms))
    print("All Unconnected rings atoms: ",all_unconnected_lst)

    ### ALGORITHM CHOICE ###
    # Initialize variables for CREST and PUCKER algorithms
    samp_zmas,samp_zmas_crest,samp_zmas_pucker,samp_zmas_torsions = {}, [], [], []
    vma =  automol.zmat.vmatrix(zma)

    nsamp = util.ring_samp_zmas(all_ring_atoms, nsamp_par, len(rings_atoms))
    if algorithm in ["torsions","pucker","robust"]:
        print("nsamp for pucker or torsions algorithms: ",nsamp)


    if algorithm == 'crest' or algorithm == 'robust': 
        ring_puckering_with_crest(geo, zrxn, spc_info, vma, samp_zmas_crest)
        samp_zmas["crest"] = samp_zmas_crest


    if algorithm == 'pucker' or algorithm == 'robust': 
        ring_puckering_with_cremerpople(geo, vma, tors_dcts, ngbs, nsamp, 
                                        all_ring_atoms, coos, samp_zmas_pucker,
                                        bonds_from_zma)
    
        # Now sample with dihedrals also 4-membered rings
        for ring_atoms, samp_range_dct in tors_dcts:
            ring_atoms = [int(idx)-1 for idx in ring_atoms.split('-')]
            if len(ring_atoms) == 4:
                print("Found a 4 membered ring")
                print(samp_range_dct)
                for i,samp_zmai in enumerate(samp_zmas_pucker):
                    samp_zmas_pucker[i] = automol.zmat.samples(
                                            samp_zmai, 1, samp_range_dct)[0]
        samp_zmas["pucker"] = samp_zmas_pucker


    if algorithm == 'torsions' or algorithm == 'robust':
        samp_zmas_torsions = list(automol.zmat.samples(zma, nsamp, all_samp_range_dct))
        all_ring_atoms_list_all, new_coord_rings_all = {}, {}
        
        for key_dct,_ in tors_dcts:
            ring_atoms = [int(idx)-1 for idx in key_dct.split('-')]
            all_ring_atoms_list_all[key_dct] = ring_atoms
        first_ring_at_sub_dct, last_ring_at_sub_dct = subs_analysis(all_ring_atoms,
                                                        all_ring_atoms_list_all, ngbs, geo)

        for i,samp_zmai in enumerate(samp_zmas_torsions):
            samp_geoi = automol.zmat.geometry(samp_zmai)
            for key_dct,ring_atoms in all_ring_atoms_list_all.items():
                new_coord_ring_all = automol.geom.coordinates(
                                    samp_geoi, tuple(ring_atoms),angstrom=True)
                new_coord_rings_all[key_dct] = [tuple(xyz) for xyz in new_coord_ring_all]
            samp_zmas_torsions[i] = fixings_subs_positions(
                                samp_zmai, all_ring_atoms_list_all, geo, coos,
                                first_ring_at_sub_dct, last_ring_at_sub_dct,
                                new_coord_rings_all)         
        samp_zmas["torsions"] = samp_zmas_torsions


    if algorithm == 'torsions2' or algorithm == 'robust':
        average_dih = {}
        for key_dct, samp_range_dct in tors_dcts:
            ring_atoms = [int(idx)-1 for idx in key_dct.split('-')] 
            dihs = []
            print(ring_atoms)
            for i in range(len(ring_atoms)):
                atm_idxs = [ring_atoms[j%len(ring_atoms)] for j in range(i,i+4)]
                print(atm_idxs)
                dih = automol.geom.dihedral_angle(geo, *atm_idxs)
                if dih > numpy.pi: 
                    dih -= 2*numpy.pi
                elif dih < -numpy.pi: 
                    dih += 2*numpy.pi
                dihs.append(dih)
            avg_dih = 0.
            for dih in dihs:
                avg_dih += abs(dih)
            avg_dih /= len(ring_atoms)
            average_dih[key_dct] = avg_dih

        samp_zmas_torsions_two = list(automol.zmat.samples_avg_dih(zma, geo, 
                                  tors_dcts, average_dih, ring_tors_dct,dih_remover))
        
        all_ring_atoms_list_all, new_coord_rings_all = {}, {}       
        for key_dct,_ in tors_dcts:
            ring_atoms = [int(idx)-1 for idx in key_dct.split('-')]
            all_ring_atoms_list_all[key_dct] = ring_atoms
        first_ring_at_sub_dct, last_ring_at_sub_dct = subs_analysis(all_ring_atoms,
                                                        all_ring_atoms_list_all, ngbs, geo)

        for i,samp_zmai in enumerate(samp_zmas_torsions_two):
            samp_geoi = automol.zmat.geometry(samp_zmai)
            for key_dct,ring_atoms in all_ring_atoms_list_all.items():
                new_coord_ring_all = automol.geom.coordinates(
                                        samp_geoi, tuple(ring_atoms),angstrom=True)
                new_coord_rings_all[key_dct] = [tuple(xyz) for xyz in new_coord_ring_all]
            samp_zmas_torsions_two[i] = fixings_subs_positions(
                                samp_zmai, all_ring_atoms_list_all, geo, coos,
                                first_ring_at_sub_dct, last_ring_at_sub_dct,
                                new_coord_rings_all)         
        samp_zmas["torsions2"] = samp_zmas_torsions_two

    with open("allsamples.xyz","w") as f:
        for algo,s_zmas in samp_zmas.items():
            for zmai in s_zmas:
                new_geo = automol.zmat.geometry(zmai)
                geo_string = automol.geom.xyz_string(new_geo, comment="")
                f.write(geo_string+"\n")        


    print("\n\nENTERING CHECKS LOOP\n\n")
    # Control dictionaries
    check_dct = {
        'dist': 3.5e-1,
        'coulomb': 1.5e-2,
    }
    _, saved_geos, _ = _saved_cnf_info(
        cnf_save_fs, thy_info)
    frag_saved_geos = [automol.geom.ring_fragments_geometry(
                       geoi, rings_atoms, ngbs) for geoi in saved_geos]

    relax_thresh = {
                    "dist":1.,
                    "angles":1.,
                    "potential":1.,
                    }
    
    if thresholds == 'relaxed': 
        relax_thresh = {
                "dist":1.7,
                "angles":0.8,
                "potential":1.5,
                }
        
    if algorithm == "robust":
        algorithm = ["crest","torsions2","torsions","pucker"]
    else:
        algorithm = [algorithm]

    ### CHECKS LOOP ###
    unique_zmas = []
    for algo in algorithm:
        print(f"Working on algorithm {algo}")
        unique_zmas.extend(ring_checks_loops(
                                            checks, samp_zmas[algo], all_unconnected_lst,
                                            relax_thresh, algo, check_dct, rings_atoms,
                                            ngbs, frag_saved_geos, cnf_run_fs, cnf_save_fs,
                                            thy_info, spc_info, vma, geo
                                            ))
        print(f"Valid samples after checks: {len(unique_zmas)}")
    # Set up the DBSCAN clustering
    unique_geos = [automol.zmat.geometry(zmai) for zmai in unique_zmas]

    if len(unique_geos) > 0:
        unique_geos = automol.geom.dbscan(
                            unique_geos,
                            rings_atoms, 
                            eps, 
                            min_samples=1
                            )
    else:
        print("No valid samples! Try changing the protocol...")
    
    unique_zmas = [automol.zmat.base.from_geometry(vma, geoi) for geoi in unique_geos]
    print(f"Valid samples after clustering: {len(unique_zmas)}")

    with open("uniques.xyz","w") as f:  
        for geoi in unique_geos:
            geo_string = automol.geom.xyz_string(geoi, comment="")
            f.write(geo_string+"\n")  


    # Set the samples
    nsamp = len(unique_zmas)
    nsamp0 = nsamp
    nsampd = util.calc_nsampd(cnf_save_fs, cnf_run_fs)

    tot_samp = nsamp + nsampd
    info_message(' - Number of samples that have been currently run:', nsampd)
    info_message(' - Number of new samples requested:', nsamp)
    info_message(' - Number of total samples:', tot_samp)

    if nsamp > 0:
        info_message(
            f'Running {nsamp} samples...', newline=1)
    

    # Create list of saved geos; initialize with saved geos
    num_saved = len(saved_geos)
    print("Initial len saved geos: ", len(saved_geos))
    
    frag_saved_geos = [automol.geom.ring_fragments_geometry(
                        geoi,rings_atoms) for geoi in saved_geos]
    
    rings_geos_strings = [automol.geom.xyz_string(geoi
                        ) for geoi in frag_saved_geos]
    
    mols = [rdkit.Chem.rdmolfiles.MolFromXYZBlock(geoi
            ) for geoi in rings_geos_strings]
        
    for samp_idx,samp_zma in enumerate(unique_zmas):
        nsamp = tot_samp - nsampd

        info_message(f"\nSample {samp_idx+1}/{nsamp0}")      
        # Build the filesystem
        rid = autofile.schema.generate_new_ring_id()
        cid = autofile.schema.generate_new_conformer_id()  
        locs = (rid, cid)  

        cnf_run_fs[-1].create(locs) 
        cnf_run_path = cnf_run_fs[-1].path(locs)
        run_fs = autofile.fs.run(cnf_run_path)

        # For first step of opt, fix all dihs of the rings and relax subs
        full_ring_tors = automol.zmat.all_rings_dihedrals(zma,rings_atoms)
        tors_names = tuple(set(names
                            for tors_dct in full_ring_tors
                            for names in tors_dct.keys()))
        print("tors_names",tors_names)
        # tors_names = tuple(set(names
        #                     for tors_dct in ring_tors_dct.values()
        #                     for names in tors_dct.keys()))
        
        if two_stage and tors_names:
            frozen_coords_lst = (tors_names, ())
            success, ret = es_runner.multi_stage_optimization(
                script_str=script_str,
                run_fs=run_fs,
                geo=samp_zma,
                spc_info=spc_info,
                thy_info=thy_info,
                frozen_coords_lst=frozen_coords_lst,
                zrxn=zrxn,
                overwrite=overwrite,
                saddle=bool(zrxn is not None),
                retryfail=retryfail,
                **kwargs
            )
        else:
            success, ret = es_runner.execute_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                run_fs=run_fs,
                geo=samp_zma,
                spc_info=spc_info,
                thy_info=thy_info,
                zrxn=zrxn,
                overwrite=overwrite,
                saddle=bool(zrxn is not None),
                retryfail=retryfail,
                **kwargs
            )
        print("adl - Finished opt. Success? ", success)
        if success:
            print("I am in if, I should show an rmsd next (lenmols= ", len(mols))
            # Get ring-subgeom
            inf_obj_temp, _, out_str = ret
            prog = inf_obj_temp.prog
            good_geo = elstruct.reader.opt_geometry(prog, out_str)
            good_ring_geo = automol.geom.ring_fragments_geometry(good_geo,rings_atoms)
            string_ring_geo = automol.geom.xyz_string(good_ring_geo)
            new_mol = rdkit.Chem.rdmolfiles.MolFromXYZBlock(string_ring_geo)

            # Compare to previously saved ring geos
            for i in range(len(mols)):
                rdkitrmsd = rdkit.Chem.AllChem.GetBestRMS(mols[i], new_mol)
                print(f"rmsd new - {i}: {rdkitrmsd}")
                if rdkitrmsd < 0.01:
                    print("Ring state previously saved in filsys")
                    print("Moving to the next conformer...")
                    break
            # If bestrmsd is > 0.01 always then save
            # Typical values of ca. 0.005-8 obseved for equivalent geos
            else:
                save_conformer(
                    ret, cnf_run_fs, cnf_save_fs, locs, thy_info,
                    zrxn=zrxn, orig_ich=spc_info[0], rid_traj=False,
                    init_zma=samp_zma)
                nsampd = util.calc_nsampd(cnf_save_fs, cnf_run_fs)
                nsampd += 1
                inf_obj.nsamp = nsampd
                cnf_save_fs[0].file.info.write(inf_obj)
                cnf_run_fs[0].file.info.write(inf_obj)
             
                _,saved_geos,_ = _saved_cnf_info(
                                    cnf_save_fs, thy_info)
                print("len saved geos now: ", len(saved_geos))
                
                if num_saved < len(saved_geos):
                    # rings_geos_strings.append(string_ring_geo) #Useful for printing
                    # mols.append(new_mol)
                    num_saved = len(saved_geos)

    with open("final-rings-stru.xyz","w") as f:
        for geo_string in rings_geos_strings:
            f.write(geo_string+"\n")


########### CHECKS LOOPS #############
def ring_checks_loops(
                        checks, samp_zmas, all_unconnected_lst,
                        relax_thresh, algorithm, check_dct,
                        rings_atoms, ngbs, frag_saved_geos,
                        cnf_run_fs, cnf_save_fs, thy_info,
                        spc_info, vma, geo
                        ):
    '''
    adl Two implementations of loops of checks for the sampled geometries.
    1) ring closure - repulsive potential - unique ring structure - not already run
    2) ring closure - CREST topology and energy checks
    Both followed by DBSCAN clustering
    '''
    def ring_closure_check():
        samp_check = True
        for unconnected_ats,unconnected_dist_value_dct in all_unconnected_lst:
            if not automol.zmat.ring_distances_reasonable( 
                            samp_zma, unconnected_ats, unconnected_dist_value_dct,
                            0.3 * relax_thresh["dist"]
                            ): 
                samp_check = False
        return samp_check
    
    def unique_check():
        samp_check = True
        try:
            frag_samp_geo = automol.geom.ring_fragments_geometry(
                samp_geo, rings_atoms, ngbs)
            samp_check = automol.geom.is_unique(frag_samp_geo,
                                frag_saved_geos, check_dct) and samp_check
            samp_check = automol.geom.is_unique(frag_samp_geo,
                                unique_frag_geos, check_dct) and samp_check
        except:
            # Accounts for weird unconnected geos that can be skipped
            samp_check = False                         
        return samp_check, frag_samp_geo
            
    
    unique_frag_geos, unique_geos = [], []
    ref_geo_for_rep_pot = geo

    if checks == 1:
        for i,samp_zma in enumerate(samp_zmas):
            print("samp",i+1)
            # if zrxn != None: # Use it for TS??
            #     if not automol.reac.similar_saddle_point_structure(samp_zma, zma, zrxn): 
            #       continue
            #     print('   -0.5) reasonable saddle point structure')
            if ring_closure_check() == False: 
                continue
            samp_geo = automol.zmat.geometry(samp_zma)
            print('   -1) reasonable distances')

            if len(unique_frag_geos) == 0:
                if algorithm in ["crest","torsions2"]:
                    ref_geo_for_rep_pot = samp_geo
            if automol.geom.has_low_relative_repulsion_energy(
                                    samp_geo, ref_geo_for_rep_pot, 
                                    thresh = 40. * relax_thresh['potential']):
                print('   -2) reasonable rel repuls')

                samp_check, frag_samp_geo = unique_check()
                if samp_check == False: 
                    continue
                print('   -3) unique check')

                run_in_run, _ = filesys.mincnf.this_conformer_was_run_in_run(
                                samp_zma, cnf_run_fs, cnf_save_fs, thy_info)               
                if run_in_run: 
                    continue
                else:
                    print('   -4) not run in run, ok')
                    unique_geos.append(samp_geo)
                    unique_frag_geos.append(frag_samp_geo)
                    if algorithm in ["crest","torsions2"]:
                        if len(unique_frag_geos) == 1:
                            ref_geo_for_rep_pot = samp_geo

    elif checks == 2:
        for i,samp_zma in enumerate(samp_zmas):
            print("samp",i+1)

            if ring_closure_check() == False: 
                continue
            samp_geo = automol.zmat.geometry(samp_zma)
            unique_frag_geos.append(samp_geo)
            print('   -1) reasonable distances')
            
        # Get xyz of all unique sampling points found, after original geo
        with open("pucker_checks.xyz","w") as f:
            geo_string = automol.geom.xyz_string(geo, comment=" ")
            f.write(geo_string+"\n")
            for samp_geo in unique_frag_geos:
                geo_string = automol.geom.xyz_string(samp_geo, comment=" ")
                f.write(geo_string+"\n")

        crest_geos = automol.geom.checks_with_crest(
                                    "pucker_checks.xyz",
                                    spc_info,
                                    )
        unique_frag_geos = []
        for samp_geo in crest_geos:
            samp_check, frag_samp_geo = unique_check()
            if samp_check == False: 
                continue
            print('   -3) unique check')

            samp_zma = automol.zmat.base.from_geometry(vma, samp_geo)
            run_in_run, _ = filesys.mincnf.this_conformer_was_run_in_run(
                            samp_zma, cnf_run_fs, cnf_save_fs, thy_info)               
            if run_in_run: 
                continue
            else:
                print('   -4) not run in run, ok')
                unique_geos.append(samp_geo)
                unique_frag_geos.append(frag_samp_geo)

    
    return [automol.zmat.base.from_geometry(vma, geoi) for geoi in unique_geos]


########### RING PUCKERING WITH CREST #############
def ring_puckering_with_crest(geo, zrxn, spc_info, vma, samp_zmas_crest):
    """
    adl Setup CREST calculation:
    - convert geo to string
    - write string to xyz file
    - call crest and wait that execution ends
    - read various xyzs from rotamers
    - convert them in zmatrix and proceed with opt
    - Add constraints file! 
    For species, do not add constraints, for TS constrain distances of reacting atoms
    """
    constrained_atoms, crest_constrain = {}, " "
    geo_string = automol.geom.xyz_string(geo, comment="")
    filename = "crest_inp_pucker.xyz"
    with open(filename,"w") as f:
        f.write(geo_string)

    if zrxn is not None:
        ts_gra = automol.reac.ts_graph(zrxn)
        ts_bond_keys = automol.graph.ts.reacting_bond_keys(ts_gra)
        for bond in ts_bond_keys:
            constrained_atoms[','.join(map(str,[el+1 for el in bond])
                            )] = automol.geom.distance(geo, *list(bond), angstrom=True)

    # Setup crest subfolder
    crest_dir_prefix = "crest_calc"
    dirs_lst = [dir for dir in os.listdir() if crest_dir_prefix in dir]
    folder_nums = []
    if not dirs_lst: crest_dir = crest_dir_prefix+"_1"
    else:
        for direc in dirs_lst:
            folder_nums.append(int(direc.split("_")[2])) 
        crest_dir = f"{crest_dir_prefix}_{max(folder_nums) + 1}"
    os.system(f"mkdir -p {crest_dir}") 
    print(f"\n####\nWorking in {crest_dir}\n####\n")

    # Create constraints file if constrained atoms are present
    if constrained_atoms: #Not empty dict
        crest_constrain = f"crest {filename} --constrain 1" # Generate coord.ref
        p = subprocess.Popen(crest_constrain, stdout=subprocess.PIPE, shell=True)
        output, err = p.communicate()  
        #This makes the wait possible
        p_status = p.wait()
        # Only bonds of unconnected atoms are constrained
        # For TSs I also need to fix distance of reacting atoms!
        with open(".xcontrol.sample","w") as f:
            f.write("$constrain\n")
            f.write("force constant = 0.10\n")
            for constrained_ats,constrained_dist in constrained_atoms.items():
                f.write(f"distance: {constrained_ats}, {constrained_dist:.4f}\n")
            f.write("$end\n")
        crest_constrain = " --cinp .xcontrol.sample " # Deafult name of constraints file
        os.system(f"cp coord.ref .xcontrol.sample {crest_dir}")
    crest_sampl = f'''cp {filename} {crest_dir}
                    echo {spc_info[-2]} > {crest_dir}/.CHRG
                    echo {int(spc_info[-1])-1} > {crest_dir}/.UHF
                    cd {crest_dir}
                    crest {filename} --gfn2 --mrest 10 --noreftopo --ewin 200.{crest_constrain}-T 8 > crest_ouput.out
                    '''
    with subprocess.Popen(crest_sampl, stdout=subprocess.PIPE, shell=True) as p:
        p.communicate()  
        p.wait()
    with open(f"{crest_dir}/crest_conformers.xyz","r") as f:
        samp_geos = (geoi for geoi,_ in automol.geom.from_xyz_trajectory_string(f.read()))
    for geoi in samp_geos:
        samp_zma = automol.zmat.from_geometry(vma, geoi)
        samp_zmas_crest.append(samp_zma)


########### RING PUCKERING WITH CREMER POPLE PARAMS #############
def ring_puckering_with_cremerpople(geo, vma_adl, tors_dcts, ngbs, nsamp, all_ring_atoms,
                                    coos, samp_zmas_pucker, dist_thresh=1.1):                               
    """
    Valid for 5 to 16 membered rings
    """
    all_ring_atoms_list = {}
    all_ring_bonds_list, all_bond_lengths = {}, {}
    all_ring_angles_list, all_angles_lengths = {}, {}
    geo_string = automol.geom.string(geo, angstrom=True)
    geo_list = [ [float(x) for x in line.split()[1:]] for line in geo_string.split('\n') ]
    for key_dct,_ in tors_dcts:
        # list of ring atoms
        ring_atoms = [int(idx)-1 for idx in key_dct.split('-')]
        N_atoms = len(ring_atoms)
        if N_atoms < 5: continue # We treat 4 membered rings separately
        # list of lists of ring atoms for each ring
        all_ring_atoms_list[key_dct] = ring_atoms
        # list of tuples of ring bonds
        ring_bonds_list = [(atom, ring_atoms[(i+1)%N_atoms] 
                            ) for i, atom in enumerate(ring_atoms)]
        # list of lists of bonds of each ring
        all_ring_bonds_list[key_dct] = ring_bonds_list
        # list of bond lengths for current ring
        bond_lengths = [automol.geom.distance(geo,i,j,angstrom=True) for i,j in ring_bonds_list]
        # list of lists of bond lengths for each ring
        all_bond_lengths[key_dct] = bond_lengths
        # list of tuples of ring angles
        ring_angle_list = [(atom, ring_atoms[(i+1)%N_atoms], ring_atoms[(i+2)%N_atoms], 
                            ) for i, atom in enumerate(ring_atoms)]
        # list of lists of angles of each ring
        all_ring_angles_list[key_dct] = ring_angle_list
        # list of angle amplitudes for current ring
        angle_amplit = [automol.geom.central_angle(geo,i,j,k) for i,j,k in ring_angle_list]
        # list of lists of angle amplitudes for each ring
        all_angles_lengths[key_dct] = angle_amplit    

    # Reorder all_ring_atoms_list dictionary, ascending atom ordering
    all_ring_atoms_list = dict(sorted(all_ring_atoms_list.items(), key=lambda item: item[1]))
    print("reordered all atoms rings list", all_ring_atoms_list)

    # Compute puckering parameters for each ring present in molecular structure
    rings_puckering_params,initZ = {}, {}
    for key_dct, ring_atoms in all_ring_atoms_list.items():
        coord = [xyz for i,xyz in enumerate(geo_list) if i in ring_atoms]
        rings_puckering_params[key_dct],initZ[key_dct] = automol.geom.cremer_pople_params(coord)
    print("Initial puckering_params", rings_puckering_params)

    first_ring_at_sub_dct, last_ring_at_sub_dct = subs_analysis(
                 all_ring_atoms,all_ring_atoms_list, ngbs, geo)

    # Generate random values for cremer pople parameters
    rng = numpy.random.default_rng() # seed=int or array[ints] for reproducibility
    puck_combos = [[] for i in range(nsamp)]
    for nsa in range(nsamp):
        for key_dct,params in rings_puckering_params.items():
            n_atoms_ring = len(all_ring_atoms_list[key_dct])
            pucker_q_range = 1 + 0.1*(n_atoms_ring%6) if n_atoms_ring>5 else 1               
            pucker_q = rng.random(len(params[0])-1)*pucker_q_range if len(params[0])>1 else 0
            if n_atoms_ring%2 == 0: 
                pucker_qlast = -pucker_q_range + rng.random(1) * 2*pucker_q_range
            else: 
                pucker_qlast = rng.random(1)*pucker_q_range
            pucker_phi = -numpy.pi + rng.random(len(params[1])) * 2 * numpy.pi
            if n_atoms_ring > 5: 
                curr_ring_sampled_params = (pucker_q.tolist(),
                                            pucker_qlast.tolist(),pucker_phi.tolist())
            elif n_atoms_ring == 5: 
                curr_ring_sampled_params = ([],pucker_qlast.tolist(),pucker_phi.tolist())
            puck_combos[nsa].append(curr_ring_sampled_params)

    # Cycle through different values of puckering
    for combo in puck_combos:
        new_coord_rings,new_Z_rings = {}, {} 
         # Initialize the zmatrix of the sampling points
        samp_zma = automol.zmat.from_geometry(vma_adl, geo)  
        for key_dct,ring_params in zip(all_ring_atoms_list,combo):
            pucker_q = ring_params[0] + ring_params[1]
            pucker_phi = list(ring_params[2])
            new_coord_ring,newZ = RR.SetRingPuckerCoords(all_ring_atoms_list[key_dct], 
                                                    pucker_q,
                                                    pucker_phi, 
                                                    all_bond_lengths[key_dct], 
                                                    all_angles_lengths[key_dct])

            new_coord_rings[key_dct] = [tuple(xyz) for xyz in new_coord_ring]
            new_Z_rings[key_dct] = newZ

        # Cycle through the rings and update torsions relative to the atoms of each ring
        for key_dct,ring_atoms in all_ring_atoms_list.items():
            # Create geo data structure for ring atoms only (needed to use dihderal_angle function)
            at_syms = automol.geom.symbols(geo, ring_atoms) 
            ring_geo = automol.geom.from_data(''.join(at_syms),new_coord_rings[key_dct],angstrom=True)
            # Update Z Matrix
            new_key_dct = {}
            for name, cord in coos.items():
                atm_idxs = cord[0]
                if len(atm_idxs) == 2:
                    new_key_dct[name] = automol.geom.distance(geo, *atm_idxs, angstrom=True)
                elif len(atm_idxs) == 3:
                    new_key_dct[name] = automol.geom.central_angle(geo, *atm_idxs, degree=True)
                elif len(set(atm_idxs) & set(ring_atoms)) == 4 : 
                    # From fourth ring atom to end of ring (also accounts for fused rings this way)
                    # Refer to indexes on ring only geometry
                    indexes = [i for i,at in enumerate(ring_atoms) if at in atm_idxs]
                    indexes = sorted(indexes, reverse=True)
                    new_key_dct[name] =  automol.geom.dihedral_angle(ring_geo, *indexes, degree=True)
                    # Added part for fused rings!
                    # Fixes case of DH in ring that in zmatrix is defined with 
                    # respect to an atom not in the ring
                elif (len(set(atm_idxs[:3]) & set(ring_atoms)) == 3
                      ) and (ring_atoms.index(atm_idxs[0]) > 2):
                    indexes = [i for i,at in enumerate(ring_atoms) if at in atm_idxs]
                    indexes = sorted(indexes, reverse=True)
                    indexes.append(indexes[-1]-1)
                
                    for name2, cord2 in coos.items(): 
                        if cord2[0][0] == atm_idxs[-1] and len(cord2[0])==4:
                            name_of_dh_of_last_of_at_idxs = name2
                            break
                    new_key_dct[name] =  (automol.geom.dihedral_angle(ring_geo, *indexes, degree=True
                        ) - new_key_dct[name_of_dh_of_last_of_at_idxs]) % 360.
                    # if new_key_dct[name] > 180.: new_key_dct[name] -= 360.
                    # if new_key_dct[name] < 45. and new_key_dct[name] > 0: new_key_dct[name] += 45. 
                    # elif new_key_dct[name] > -45. and new_key_dct[name] < 0: new_key_dct[name] -= 45.

                elif len(atm_idxs) == 4: # All other atoms keep original value
                    new_key_dct[name] = automol.zmat.value(samp_zma, name, degree=True) 

            samp_zma = automol.zmat.set_values_by_name(samp_zma, new_key_dct)
     
        samp_zma = fixings_subs_positions(samp_zma, all_ring_atoms_list, geo, coos,
                           first_ring_at_sub_dct, last_ring_at_sub_dct,
                           new_coord_rings, dist_thresh)

        samp_zmas_pucker.append(samp_zma)


########### GET INFORMATION ON SUBS ON FIRST AND LAST RING ATOM #############
def subs_analysis(all_ring_atoms,all_ring_atoms_list, ngbs, geo):
    last_ring_at_sub_dct,first_ring_at_sub_dct = {}, {}

    for key_dct, ring_atoms in all_ring_atoms_list.items():

        # Gets position of first substituent on first and last atom of ring (which usually messes up)
        first_ring_at,last_ring_at = min(ring_atoms),max(ring_atoms)
        first_ring_at_ngbs,last_ring_at_ngbs = set(ngbs[first_ring_at]), set(ngbs[last_ring_at])
        first_ring_at_subs = first_ring_at_ngbs.difference(set(all_ring_atoms))
        last_ring_at_subs = last_ring_at_ngbs.difference(set(all_ring_atoms))
        # Get xyz of ring and sub
        coord_ring = [xyz for i,(_,xyz) in enumerate(geo) if i in ring_atoms]
        first_coord_sub = [xyz for i,(_,xyz) in enumerate(geo) if i in first_ring_at_subs]
        last_coord_sub = [xyz for i,(_,xyz) in enumerate(geo) if i in last_ring_at_subs]
        # Call alpha beta calculator
        first_sub_params = [tuple(RR.GetRingSubstituentPosition(coord_ring,coord_sub,-1)
                                  ) for coord_sub in first_coord_sub]
        last_sub_params = [tuple(RR.GetRingSubstituentPosition(coord_ring,coord_sub,-1)
                                 ) for coord_sub in last_coord_sub]
        # Update sub dictionaries
        last_ring_at_sub_dct[key_dct] = {key:value for key,value in zip(
                                        last_ring_at_subs,last_sub_params)}
        first_ring_at_sub_dct[key_dct] = {key:value for key,value in zip(
                                        first_ring_at_subs,first_sub_params)}

    return first_ring_at_sub_dct,last_ring_at_sub_dct


########### FIX POSITIONS OF SUBS OF FIRST AND LAST RING ATOM #############
def fixings_subs_positions(samp_zma, all_ring_atoms_list, geo, coos,
                           first_ring_at_sub_dct, last_ring_at_sub_dct,
                           new_coord_rings, dist_thresh=1.2):
    
    samp_geo = automol.zmat.geometry(samp_zma)
    # I need to perform the substituents check here AFTER I have built the samp ZMat!
    for key_dct,ring_atoms in all_ring_atoms_list.items():

        at_syms = automol.geom.symbols(geo, ring_atoms) 
        new_key_dct = {}
        for name, cord in coos.items():
            atm_idxs = cord[0]
            if len(atm_idxs) == 2:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, angstrom=True)
            elif len(atm_idxs) == 3:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, degree=True)
            elif len(atm_idxs) == 4:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, degree=True)
                # First sub of first ring atom
                # Check that I am working on the third ring atom
                if atm_idxs[0] == ring_atoms[2]:
                    change_dh = False
                    for atm in first_ring_at_sub_dct[key_dct]:
                        dist_sub_to_last = automol.geom.distance(
                                samp_geo, atm,max(ring_atoms), angstrom=True)
                        if dist_sub_to_last < dist_thresh: change_dh = True

                        for atm2 in last_ring_at_sub_dct[key_dct]:
                            dist_sub_to_sub = automol.geom.distance(
                                samp_geo, atm,atm2, angstrom=True)
                            if dist_sub_to_sub < dist_thresh: change_dh = True
                    if change_dh: 
                        new_key_dct[name] -= 60.

        samp_zma = automol.zmat.set_values_by_name(samp_zma, new_key_dct)
    samp_geo = automol.zmat.geometry(samp_zma)

    # I need to perform the substituents check here AFTER AFTER I have built the samp ZMat!
    for key_dct,ring_atoms in all_ring_atoms_list.items():
        new_key_dct = {}
        for name, cord in coos.items():
            atm_idxs = cord[0]
            if len(atm_idxs) == 2:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, angstrom=True)
            elif len(atm_idxs) == 3:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, degree=True)
            elif len(atm_idxs) == 4:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, degree=True)
                # First sub of first ring atom
                # If previous -60 rotation didn't help, rotate 60 on the other direction
                if atm_idxs[0] == ring_atoms[2]:
                    change_dh = False
                    for atm in first_ring_at_sub_dct[key_dct]:
                        dist_sub_to_last = automol.geom.distance(
                                samp_geo, atm,max(ring_atoms), angstrom=True)
                        if dist_sub_to_last < dist_thresh: change_dh = True

                        for atm2 in last_ring_at_sub_dct[key_dct]:
                            dist_sub_to_sub = automol.geom.distance(
                                samp_geo, atm,atm2, angstrom=True)
                            if dist_sub_to_sub < dist_thresh: change_dh = True
                    if change_dh: 
                        new_key_dct[name] += 120.
                           
        samp_zma = automol.zmat.set_values_by_name(samp_zma, new_key_dct)
    samp_geo = automol.zmat.geometry(samp_zma)

    # I need to perform the substituents check here AFTER AFTER AFTER I have built the samp ZMat!
    for key_dct,ring_atoms in all_ring_atoms_list.items():
        new_key_dct = {}
        for name, cord in coos.items():
            atm_idxs = cord[0]
            if len(atm_idxs) == 2:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, angstrom=True)
            elif len(atm_idxs) == 3:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, degree=True)
            elif len(atm_idxs) == 4:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, degree=True)

                # Last ring atom first sub. 
                if atm_idxs[0] in last_ring_at_sub_dct[key_dct]:
                    change_dh = False
                    if atm_idxs[0] == min(last_ring_at_sub_dct[key_dct]):
                        for atm in last_ring_at_sub_dct[key_dct]:
                            dist_sub_to_first = automol.geom.distance(
                                                samp_geo, atm,min(ring_atoms), angstrom=True)
                            if dist_sub_to_first < dist_thresh: change_dh = True
                            for atm2 in first_ring_at_sub_dct[key_dct]:
                                dist_sub_to_sub = automol.geom.distance(
                                    samp_geo, atm,atm2, angstrom=True)
                                if dist_sub_to_sub < dist_thresh: change_dh = True
                    if change_dh:  
                        new_key_dct[name] -= 60.     
        samp_zma = automol.zmat.set_values_by_name(samp_zma, new_key_dct)
    samp_geo = automol.zmat.geometry(samp_zma)

    for key_dct,ring_atoms in all_ring_atoms_list.items():
        new_key_dct = {}
        for name, cord in coos.items():
            atm_idxs = cord[0]
            if len(atm_idxs) == 2:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, angstrom=True)
            elif len(atm_idxs) == 3:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, degree=True)
            elif len(atm_idxs) == 4:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, degree=True)

                # Last ring atom first sub. 
                if atm_idxs[0] in last_ring_at_sub_dct[key_dct]: 
                    change_dh = False
                    if atm_idxs[0] == min(last_ring_at_sub_dct[key_dct]):
                        for atm in last_ring_at_sub_dct[key_dct]:
                            dist_sub_to_first = automol.geom.distance(
                                                samp_geo, atm,min(ring_atoms), angstrom=True)
                            if dist_sub_to_first < dist_thresh: change_dh = True
                            for atm2 in first_ring_at_sub_dct[key_dct]:
                                dist_sub_to_sub = automol.geom.distance(
                                    samp_geo, atm,atm2, angstrom=True)
                                if dist_sub_to_sub < dist_thresh: change_dh = True
                    if change_dh:  
                        new_key_dct[name] += 120.      
        samp_zma = automol.zmat.set_values_by_name(samp_zma, new_key_dct)

    return samp_zma


def _presamp_save(spc_info, cnf_run_fs, cnf_save_fs,
                  thy_info, zrxn=None, rid=None, ref_zma=None):
    """ Loop over the RUN filesys and save conformers
    """

    job = elstruct.Job.OPTIMIZATION

    if not cnf_run_fs[0].exists():
        print(" - No conformers in RUN filesys to save.")
    else:
        print(" - Found conformers in RUN filesys to save.\n")
        for locs in cnf_run_fs[-1].existing():
            cnf_run_path = cnf_run_fs[-1].path(locs)
            run_fs = autofile.fs.run(cnf_run_path)
            if run_fs[-1].file.info.exists([job]):
                inf_obj = run_fs[-1].file.info.read([job])
                if inf_obj.status == autofile.schema.RunStatus.SUCCESS:
                    print(f"\nReading from conformer run at {cnf_run_path}")

                    # Read the electronic structure optimization job
                    success, ret = es_runner.read_job(
                        job=job, run_fs=run_fs)

                    if success:
                        if run_fs[-1].file.zmatrix.exists([job]):
                            init_zma = run_fs[-1].file.zmatrix.read([job])
                        else:
                            init_zma = None
                        save_conformer(
                            ret, cnf_run_fs, cnf_save_fs, locs, thy_info,
                            zrxn=zrxn, orig_ich=spc_info[0],
                            init_zma=init_zma, ref_zma=ref_zma)

        # Update the conformer trajectory file
        print('')
        for rloc in cnf_save_fs[-2].existing():
            filesys.mincnf.traj_sort(cnf_save_fs, thy_info, rid=rloc[0])


def save_conformer(ret, cnf_run_fs, cnf_save_fs, locs, thy_info, zrxn=None,
                   orig_ich='', rid_traj=False, init_zma=None, ref_zma=None):
    """ save the conformers that have been found so far
          # Only go through save procedure if conf not in save
          # may need to get geo, ene, etc; maybe make function
    """

    saved_locs, saved_geos, saved_enes = _saved_cnf_info(
        cnf_save_fs, thy_info, locs)

    inf_obj, _, out_str = ret
    prog = inf_obj.prog
    method = inf_obj.method
    ene = elstruct.reader.energy(prog, method, out_str)
    geo = elstruct.reader.opt_geometry(prog, out_str)
    zma = None
    if init_zma is not None:
        zma = filesys.save.read_zma_from_geo(init_zma, geo)
    if zma is None:
        zma = filesys.save.read_job_zma(ret, init_zma=init_zma)

    # Assess if geometry is properly connected
    viable = _geo_connected(geo, zrxn)
    if viable:
        if zrxn:
            viable = _ts_geo_viable(
                zma, zrxn, cnf_save_fs, thy_info, ref_zma=ref_zma)
        else:
            viable = _inchi_are_same(orig_ich, geo)

    # Determine uniqueness of conformer, save if needed
    if viable:
        if _geo_unique(geo, ene, saved_geos, saved_enes, zrxn):
            sym_id = _sym_unique(
                geo, ene, saved_geos, saved_enes)
            print('save_conformer locs:', locs, sym_id)
            if sym_id is None:
                filesys.save.conformer(
                    ret, None, cnf_save_fs, thy_info[1:],
                    init_zma=init_zma,  zrxn=zrxn,
                    rng_locs=(locs[0],), tors_locs=(locs[1],))
            # else:
            #     sym_locs = saved_locs[sym_id]
            #     filesys.save.sym_indistinct_conformer(
            #         geo, cnf_save_fs, locs, sym_locs)
            #     if cnf_save_fs[-1].exists(locs):
            #         cnf_save_path = cnf_save_fs[-1].path(locs)
            #         shutil.rmtree(cnf_save_path)
            #     if cnf_run_fs[-1].exists(locs):
            #         cnf_run_path = cnf_run_fs[-1].path(locs)
            #         shutil.rmtree(cnf_run_path)
        else:
            if cnf_save_fs[-1].exists(locs):
                cnf_save_path = cnf_save_fs[-1].path(locs)
                shutil.rmtree(cnf_save_path)
            if cnf_run_fs[-1].exists(locs):
                cnf_run_path = cnf_run_fs[-1].path(locs)
                shutil.rmtree(cnf_run_path)

        # Update the conformer trajectory file
        obj('vspace')
        rid = None
        if rid_traj:
            rid = locs[0]
        filesys.mincnf.traj_sort(cnf_save_fs, thy_info, rid=rid)


def _saved_cnf_info(cnf_save_fs, mod_thy_info, orig_locs=None):
    """ get the locs, geos and enes for saved conformers
    """

    saved_locs = list(cnf_save_fs[-1].existing())
    saved_locs = [locs for locs in saved_locs if not locs == orig_locs]
    saved_geos = [cnf_save_fs[-1].file.geometry.read(locs)
                  for locs in saved_locs]
    found_saved_locs = []
    found_saved_geos = []
    found_saved_enes = []
    for idx, locs in enumerate(saved_locs):
        path = cnf_save_fs[-1].path(locs)
        sp_save_fs = autofile.fs.single_point(path)
        if sp_save_fs[-1].file.energy.exists(mod_thy_info[1:4]):
            found_saved_enes.append(sp_save_fs[-1].file.energy.read(
                mod_thy_info[1:4]))
            found_saved_locs.append(saved_locs[idx])
            found_saved_geos.append(saved_geos[idx])
        else:
            info_message(
                f'No energy saved in single point directory for {path}')
            # geo_inf_obj = cnf_save_fs[-1].file.geometry_info.read(
            #     mod_thy_info[1:4])
            geo_inf_obj = cnf_save_fs[-1].file.geometry_info.read(
                locs)
            geo_end_time = geo_inf_obj.utc_end_time
            current_time = autofile.schema.utc_time()
            if (current_time - geo_end_time).total_seconds() < 120:
                last_time = (current_time - geo_end_time).total_seconds()
                wait_time = 120 - last_time
                info_message(
                    f'Geo was saved in the last {last_time:3.2f} seconds, '
                    f'waiting for {wait_time:3.2f} seconds')
                time.sleep(wait_time)
                if sp_save_fs[-1].file.energy.exists(mod_thy_info[1:4]):
                    found_saved_enes.append(sp_save_fs[-1].file.energy.read(
                        mod_thy_info[1:4]))
                    found_saved_locs.append(saved_locs[idx])
                    found_saved_geos.append(saved_geos[idx])
                    info_message('the energy is now found')
                else:
                    info_message('waiting helped nothing')

    return found_saved_locs, found_saved_geos, found_saved_enes


def _init_geom_is_running(cnf_run_fs):
    """ Check the RUN filesystem for currently running initial geometry submissions
    """
    running = False
    jobs = [elstruct.Job.OPTIMIZATION, elstruct.Job.HESSIAN]
    locs = cnf_run_fs[-1].existing()
    if not locs:
        wait_time = random.randint(10, 60)
        print('lets wait a bit', wait_time)
        time.sleep(wait_time)
        locs = cnf_run_fs[-1].existing()
    print('im going to check locs',  cnf_run_fs[-1].existing())
    for locs in cnf_run_fs[-1].existing():
        cnf_run_path = cnf_run_fs[-1].path(locs)
        run_fs = autofile.fs.run(cnf_run_path)
        print('im going to check here',  cnf_run_path)
        for job in jobs:
            if not run_fs[-1].file.info.exists([job]):
                continue
            print('im going to check it job',  job)
            inf_obj = run_fs[-1].file.info.read([job])
            status = inf_obj.status
            print('its job status is', status)
            if status == autofile.schema.RunStatus.RUNNING:
                start_time = inf_obj.utc_start_time
                current_time = autofile.schema.utc_time()
                _time = (current_time - start_time).total_seconds()
                print('its jtime is', _time)
                if _time < 3000000:
                    path = cnf_run_fs[-1].path(locs)
                    info_message(
                        'init_geom was started in the last '
                        f'{_time/3600:3.4f} hours in {path}.')
                    running = True
                    break
    return running


def this_conformer_was_run_in_save(zma, cnf_fs):
    """ Assess if a conformer was run in save
    """
    running = False
    for locs in cnf_fs[-1].existing(ignore_bad_formats=True):
        cnf_path = cnf_fs[-1].path(locs)
        if cnf_fs[-1].file.geometry_input.exists(locs):
            inp_str = cnf_fs[-1].file.geometry_input.read(locs)
            inp_str = inp_str.replace('=', '')
            inf_obj = cnf_fs[-1].file.geometry_info.read(locs)
            prog = inf_obj.prog
            inp_zma = elstruct.reader.inp_zmatrix(prog, inp_str)
            if inp_zma is not None:
                if automol.zmat.almost_equal(inp_zma, zma,
                                             dist_rtol=0.018, ang_atol=.2):
                    info_message(
                        f'This conformer was already run in {cnf_path}.')
                    running = True
                    break
            else:
                info_message(f'Program {prog} lacks inp ZMA reader for check')
            sym_fs = autofile.fs.symmetry(cnf_path)
            for sym_locs in sym_fs[-1].existing(ignore_bad_formats=True):
                if sym_fs[-1].file.geometry.exists(sym_locs):
                    sym_geo = sym_fs[-1].file.geometry.read(sym_locs)
                    try:
                        sym_zma = filesys.save.rebuild_zma_from_opt_geo(
                            inp_zma, sym_geo)
                        if automol.zmat.almost_equal(
                                sym_zma, zma,
                                dist_rtol=0.018, ang_atol=.2):
                            info_message(
                                f'This conformer was already run in sym of {cnf_path}.')
                            running = True
                        break
                    except:
                        print('zmat in different format')
    return running


def this_conformer_is_running(zma, cnf_run_fs):
    """ Check the RUN filesystem for similar geometry
        submissions that are currently running
    """

    running = False
    job = elstruct.Job.OPTIMIZATION
    for locs in cnf_run_fs[-1].existing(ignore_bad_formats=True):
        cnf_run_path = cnf_run_fs[-1].path(locs)
        run_fs = autofile.fs.run(cnf_run_path)
        run_path = run_fs[-1].path([job])
        if run_fs[-1].file.info.exists([job]):
            inf_obj = run_fs[-1].file.info.read([job])
            status = inf_obj.status
            if status == autofile.schema.RunStatus.RUNNING:
                start_time = inf_obj.utc_start_time
                current_time = autofile.schema.utc_time()
                if (current_time - start_time).total_seconds() < 3000000:
                    subrun_fs = autofile.fs.subrun(run_path)
                    inp_str = subrun_fs[0].file.input.read([0, 0])
                    inp_str = inp_str.replace('=', '')
                    prog = inf_obj.prog
                    if 'molpro' in prog:
                        print('Warning: Since using Molpro, check for running '
                              'conformer is disabled!')
                        running = False
                    else:
                        inp_zma = elstruct.reader.inp_zmatrix(prog, inp_str)
                        if automol.zmat.almost_equal(
                            inp_zma, zma, dist_rtol=0.018, ang_atol=.2):
                            _hr = (current_time - start_time).total_seconds()/3600.
                            info_message(
                                'This conformer was started in the last ' +
                                f'{_hr:3.4f} hours in {run_path}.')
                            running = True
                            break
    return running


def _geo_connected(geo, rxn):
    """ Assess if geometry is connected. Right now only works for
        minima
    """

    # Determine connectivity (only for minima)
    if rxn is None:
        gra = automol.geom.graph(geo)
        conns = automol.graph.connected_components(gra)
        lconns = len(conns)
    else:
        lconns = 1

    # Check connectivity
    if lconns == 1:
        connected = True
    else:
        bad_conformer('disconnected')
        connected = False

    return connected


def _geo_unique(geo, ene, seen_geos, seen_enes, zrxn=None):
    """ Assess if a geometry is unique to saved geos
        Need to pass the torsions
    """

    if zrxn is None:
        check_dct = {'dist': 0.3, 'tors': None}
    else:
        check_dct = {'dist': 0.3}

    no_similar_energies = True
    for sene in seen_enes:
        if abs(sene - ene) < 1e-5:
            no_similar_energies = False

    if no_similar_energies:
        unique = True
    else:
        unique, _ = automol.geom.is_unique(
            geo, seen_geos, check_dct=check_dct)

    if not unique:
        bad_conformer('not unique')

    return unique


def _inchi_are_same(orig_ich, geo):
    """ Assess if a geometry has the same connectivity to
     saved geos evaluated in temrs of inchi
    """
    same = False
    ich = automol.geom.chi(geo)
    assert automol.chi.is_complete(orig_ich), (
        f'the inchi {orig_ich} orig_ich is not complete')
    if ich == orig_ich:
        same = True
    if not same:
        warning_message(
            f" - new inchi {ich} not the same as old {orig_ich}")

    return same


def _check_old_inchi(orig_ich, seen_geos, saved_locs, cnf_save_fs):
    """
    This assumes you already have bad geos in your save
    """
    for i, geoi in enumerate(seen_geos):
        if not orig_ich == automol.geom.chi(geoi):
            smi = automol.geom.smiles(geoi)
            path = cnf_save_fs[-1].path(saved_locs[i])
            error_message(
                f'inchi do not match for {smi} at {path}')


def _sym_unique(geo, ene, saved_geos, saved_enes, ethresh=1.0e-5):
    """ Check if a conformer is symmetrically distinct from the
        existing conformers in the filesystem
    """

    sym_idx = None
    new_saved_geos = []
    idx_dct = {}
    for i, (sene, sgeo) in enumerate(zip(saved_enes, saved_geos)):
        if abs(ene - sene) < ethresh:
            idx_dct[len(new_saved_geos)] = i
            new_saved_geos.append(sgeo)
    if new_saved_geos:
        _, sym_idx = automol.geom.is_unique(
            geo, new_saved_geos, check_dct={'coulomb': 1e-2})

    if sym_idx is not None:
        print(' - Structure is not symmetrically unique.')
        sym_idx = idx_dct[sym_idx]

    return sym_idx


def _ts_geo_viable(zma, zrxn, cnf_save_fs, mod_thy_info, zma_locs=(0,), ref_zma=None):
    """ Perform a series of checks to assess the viability
        of a transition state geometry prior to saving
    """

    # Obtain the min-ene zma and bond keys
    viable = False
    if ref_zma is not None:
        sens = 2.
    else:
        _, cnf_save_path = filesys.mincnf.min_energy_conformer_locators(
            cnf_save_fs, mod_thy_info)
        if cnf_save_path:
            zma_save_fs = fs.zmatrix(cnf_save_path)
            ref_zma = zma_save_fs[-1].file.zmatrix.read(zma_locs)
            sens = 10.
    if ref_zma is not None:
        viable =  automol.reac.similar_saddle_point_structure(zma, ref_zma, zrxn, sens)
    else:
        print('nothing to compares zmat to')
    return viable


def unique_fs_ring_confs(
        cnf_save_fs, cnf_locs_lst,
        ini_cnf_save_fs, ini_cnf_locs_lst):
    """ Assess which structures from the cnf_save_fs currently exist
        within the ini_cnf_save_fs. Generate a lst of unique structures
        in the ini_cnf_save_fs.
    """
    uni_ini_rng_locs = []
    uni_ini_cnf_locs = []
    rng_dct = {}

    for ini_locs in ini_cnf_locs_lst:
        ini_rid, _ = ini_locs
        if ini_rid in [locs[0] for locs in uni_ini_rng_locs]:
            uni_ini_rng_locs.append(ini_locs)
            continue
        ini_cnf_save_path = ini_cnf_save_fs[-1].path(ini_locs)
        ini_zma_fs = autofile.fs.zmatrix(ini_cnf_save_path)
        inizma = ini_zma_fs[-1].file.zmatrix.read((0,))
        inigeo = automol.zmat.geometry(inizma)
        checking('structures', ini_cnf_save_path)
        # Check to see if a similar ring pucker is in the runlvl filesystem
        found_rid = None
        if ini_rid in rng_dct:
            found_rid = rng_dct[ini_rid]
        else:
            frag_ini_geo = automol.geom.ring_fragments_geometry(inigeo)
            if frag_ini_geo is not None:
                frag_ini_zma = automol.geom.zmatrix(frag_ini_geo)
            skip_trid = []
            for tlocs in cnf_save_fs[-1].existing():
                trid, _ = tlocs
                if frag_ini_geo is None:
                    found_rid = trid
                    rng_dct[ini_rid] = trid
                    break
                if trid in skip_trid:
                    continue
                cnf_save_path = cnf_save_fs[-1].path(tlocs)
                zma_fs = autofile.fs.zmatrix(cnf_save_path)
                zma = zma_fs[-1].file.zmatrix.read((0,))
                geo = automol.zmat.geometry(zma)
                # geo = cnf_save_fs[-1].file.geometry.read(tlocs)
                frag_geo = automol.geom.ring_fragments_geometry(geo)
                if frag_geo is None:
                    found_rid = trid
                    rng_dct[ini_rid] = trid
                    break
                frag_zma = automol.geom.zmatrix(frag_geo)
                if automol.zmat.almost_equal(frag_ini_zma, frag_zma,
                                             dist_rtol=0.1, ang_atol=.4):
                    rng_dct[ini_rid] = trid
                    found_rid = trid
                    break
                skip_trid.append(trid)
        # If no similar runlvl ring pucker, then add it to unique rings
        if found_rid is None:
            uni_ini_rng_locs.append(ini_locs)
            continue

        # If similar ring is found,
        # check actual conformers under that ring orientation
        found = False
        for locs in cnf_locs_lst:
            rid, _ = locs
            if rid != found_rid:
                continue
            cnf_save_path = cnf_save_fs[-1].path(locs)
            zma_fs = autofile.fs.zmatrix(cnf_save_path)
            zma = zma_fs[-1].file.zmatrix.read((0,))
            geo = automol.zmat.geometry(zma)
            if automol.zmat.almost_equal(inizma, zma,
                                         dist_rtol=0.1, ang_atol=.4):
                info_message(
                    f'- Similar structure found at {cnf_save_path}')
                found = True
                break

        # If no match was found, add to unique locs lst
        if not found:
            uni_ini_cnf_locs.append((ini_locs, found_rid))

    return uni_ini_rng_locs, uni_ini_cnf_locs


def unique_fs_confs(cnf_save_fs, cnf_save_locs_lst,
                    ini_cnf_save_fs, ini_cnf_save_locs_lst):
    """ Assess which structures from the cnf_save_fs currently exist
        within the ini_cnf_save_fs. Generate a lst of unique structures
        in the ini_cnf_save_fs.
    """

    uni_ini_cnf_save_locs = []
    for ini_locs in ini_cnf_save_locs_lst:

        # Initialize variable to see if initial struct is found
        found = False

        # Loop over structs in cnf_save, see if they match the current struct
        inigeo = ini_cnf_save_fs[-1].file.geometry.read(ini_locs)
        inizma = automol.geom.zmatrix(inigeo)
        # inizma =  ini_cnf_save_fs[-1].file.zmatrix.read(ini_locs)
        ini_cnf_save_path = ini_cnf_save_fs[-1].path(ini_locs)
        checking('structures', ini_cnf_save_path)
        for locs in cnf_save_locs_lst:
            geo = cnf_save_fs[-1].file.geometry.read(locs)
            zma = automol.geom.zmatrix(geo)
            # zma =  cnf_save_fs[-1].file.zmatrix.read(locs)
            if automol.zmat.almost_equal(inizma, zma,
                                         dist_rtol=0.1, ang_atol=.4):
                cnf_save_path = cnf_save_fs[-1].path(locs)
                info_message(
                    f'- Similar structure found at {cnf_save_path}')
                found = True
                break

        # If no match was found, add to unique locs lst
        if not found:
            uni_ini_cnf_save_locs.append(ini_locs)

    return uni_ini_cnf_save_locs


def rng_loc_for_geo(geo, cnf_save_fs):
    """ Find the ring-conf locators for a given geometry in the
        conformamer save filesystem
    """

    rid = None
    frag_geo = automol.geom.ring_fragments_geometry(geo)
    if frag_geo is not None:
        frag_zma = automol.geom.zmatrix(frag_geo)
    checked_rids = []
    for locs in cnf_save_fs[-1].existing():
        current_rid, _ = locs
        if current_rid in checked_rids:
            continue
        checked_rids.append(current_rid)
        locs_geo = cnf_save_fs[-1].file.geometry.read(locs)
        frag_locs_geo = automol.geom.ring_fragments_geometry(locs_geo)
        if frag_locs_geo is None or frag_geo is None:
            rid = locs[0]
            break
        frag_locs_zma = automol.geom.zmatrix(frag_locs_geo)
        # for now: set tolerances to include all ring puckering
        # previous tolerances: dist_rtol=0.15, ang_atol=.45):
        if automol.zmat.almost_equal(frag_locs_zma, frag_zma,
                                     dist_rtol=150., ang_atol=45.):
            rid = locs[0]
            break

    return rid

