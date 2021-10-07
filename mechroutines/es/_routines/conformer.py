""" es_runners for conformer
"""

import shutil
import time
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
            cnf_path = ini_cnf_save_fs[-1].path(locs)
            debug_message(f'Removing {cnf_path}')
            shutil.rmtree(cnf_path)

    if geo_init is None:
        if 'geo_inp' in spc_dct_i:
            geo_init = spc_dct_i['geo_inp']
            info_message(
                'Getting initial geometry from geom dictionary')

    if geo_init is None:
        geo_init = automol.inchi.geometry(spc_dct_i['inchi'])
        info_message('Getting initial geometry from inchi')

    # Check if the init geometry is connected
    if geo_init is not None:
        if not automol.geom.connected(geo_init):
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
        method_dct, job=elstruct.Job.OPTIMIZATION)

    # Call the electronic structure optimizer
    success, ret = es_runner.execute_job(
        job=elstruct.Job.ENERGY,
        script_str=script_str,
        run_fs=run_fs,
        geo=zma_init,
        spc_info=spc_info,
        thy_info=mod_thy_info,
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
        method_dct, job=elstruct.Job.OPTIMIZATION)

    # Call the electronic structure optimizer
    success, ret = es_runner.execute_job(
        job=elstruct.Job.OPTIMIZATION,
        script_str=script_str,
        run_fs=run_fs,
        geo=zma_init,
        spc_info=spc_info,
        thy_info=mod_thy_info,
        overwrite=overwrite,
        **kwargs
    )

    # read the geometry
    if success:
        inf_obj, _, out_str = ret
        geo = elstruct.reader.opt_geometry(inf_obj.prog, out_str)
        zma = elstruct.reader.opt_zmatrix(inf_obj.prog, out_str)
        if zma is None:
            zma = automol.geom.zmatrix(geo)
        geo_conn = bool(automol.geom.connected(geo))

    # If connected, check for imaginary modes and fix them if possible
    if geo_conn:

        # Remove the imaginary mode
        geo, ret = remove_imag(
            geo, ret, spc_info, method_dct,
            run_fs, kickoff_size, kickoff_backward, kickoff_mode=0,
            overwrite=overwrite)

        # Recheck connectivity for imag-checked geometry
        if geo is not None:
            conf_found = True
            conn = automol.geom.connected(geo)
            proper_stereo = _inchi_are_same(spc_info[0], geo)
            if conn and proper_stereo:
                info_message(
                    'Saving structure as the first conformer...', newline=1)
                filesys.save.conformer(
                    ret, None, cnf_save_fs, mod_thy_info[1:],
                    rng_locs=(locs[0],), tors_locs=(locs[1],))
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
    else:
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
                     use_locs=None,
                     **kwargs):
    """ generate single optimized geometry to be saved into a
        filesystem
    """
    skip_job = False
    if this_conformer_is_running(zma, cnf_run_fs):
        skip_job = True
    elif this_conformer_was_run_in_save(zma, cnf_save_fs):
        skip_job = True
    if not skip_job:
        run_in_run, _ = filesys.mincnf.this_conformer_was_run_in_run(
            zma, cnf_run_fs)
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
        cnf_run_fs[-1].create(locs)
        cnf_run_path = cnf_run_fs[-1].path(locs)
        run_fs = autofile.fs.run(cnf_run_path)

        # Run the optimization
        info_message('Optimizing a single conformer', zrxn)
        success, ret = es_runner.execute_job(
            job=elstruct.Job.OPTIMIZATION,
            script_str=script_str,
            run_fs=run_fs,
            geo=zma,
            spc_info=spc_info,
            thy_info=mod_thy_info,
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
            zma = elstruct.reader.opt_zmatrix(prog, out_str)
            saved_locs, saved_geos, saved_enes = _saved_cnf_info(
                cnf_save_fs, mod_thy_info)

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
                        ret, None, cnf_save_fs, mod_thy_info[1:], zrxn=zrxn,
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
    cnf_run_fs[1].create([rid])
    if resave:
        _presamp_save(
            spc_info, cnf_run_fs, cnf_save_fs, thy_info, zrxn=zrxn, rid=rid)

    # Build filesys
    cnf_save_fs[1].create([rid])
    inf_obj = autofile.schema.info_objects.conformer_branch(0)

    # Set the samples
    nsamp, tors_range_dct = _calc_nsamp(tors_names, nsamp_par, zma, zrxn=zrxn)
    nsamp0 = nsamp
    nsampd = _calc_nsampd(cnf_save_fs, cnf_run_fs, rid)

    tot_samp = nsamp - nsampd
    brk_tot_samp = nsamp * 5

    info_message(
        ' - Number of samples that have been currently run:', nsampd)
    info_message(' - Number of samples requested:', nsamp)

    if nsamp-nsampd > 0:
        info_message(
            f'Running {nsamp-nsampd} samples...', newline=1)
    samp_idx = 1
    samp_attempt_idx = 1
    while True:
        nsamp = nsamp0 - nsampd
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
        ref_pot = automol.pot.intramol_interaction_potential_sum(
            automol.zmat.geometry(zma))
        samp_pot = automol.pot.intramol_interaction_potential_sum(
            automol.zmat.geometry(samp_zma))
        while samp_pot-ref_pot > repulsion_thresh and bad_geo_cnt < 1000:
            if print_debug:
                warning_message('Structure has high repulsion.')
                warning_message(
                    'Sums of intramol LJ potential interactions [kcal/mol]:',
                    f'Ref:{ref_pot:.2f}, Test:{samp_pot:.2f}, '
                    f'Diff:{samp_pot-ref_pot:.2f}')
                warning_message(
                    'Generating new sample Z-Matrix')
            samp_zma, = automol.zmat.samples(zma, 1, tors_range_dct)
            samp_pot = automol.pot.intramol_interaction_potential_sum(
                automol.zmat.geometry(samp_zma))
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
                overwrite=overwrite,
                saddle=bool(zrxn is not None),
                retryfail=retryfail,
                **kwargs
            )

        # save function added here
        if success:
            _save_conformer(
                ret, cnf_save_fs, locs, thy_info,
                zrxn=zrxn, orig_ich=spc_info[0], rid_traj=True,
                init_zma=samp_zma)

            nsampd = _calc_nsampd(cnf_save_fs, cnf_run_fs, rid)
            nsampd += 1
            samp_idx += 1
            inf_obj.nsamp = nsampd
            cnf_save_fs[1].file.info.write(inf_obj, [rid])
            cnf_run_fs[1].file.info.write(inf_obj, [rid])

        # Increment attempt counter
        samp_attempt_idx += 1


def _num_samp_zmas(ring_atoms, nsamp_par):
    """ choose starting number of sample zmas
    """
    ntors = len(ring_atoms) - 3
    apar, bpar, cpar = nsamp_par[0:3]
    return 50 * (apar + bpar * cpar**ntors)


def ring_conformer_sampling(
        zma, spc_info, thy_info,
        cnf_run_fs, cnf_save_fs,
        script_str, overwrite,
        nsamp_par=(False, 3, 1, 3, 50, 50),
        ring_tors_dct=None,
        zrxn=None, two_stage=False, retryfail=False,
        **kwargs):
    """ run sampling algorithm to find conformers
    """

    # Build filesys
    cnf_save_fs[0].create()
    inf_obj = autofile.schema.info_objects.conformer_trunk(0)

    # Set up torsions
    geo = automol.zmat.geometry(zma)
    tors_dcts = ring_tors_dct.items() if ring_tors_dct is not None else {}
    rings_atoms = []
    for ring_atoms, samp_range_dct in tors_dcts:
        rings_atoms.append([int(idx)-1 for idx in ring_atoms.split('-')])
    gra = automol.geom.graph(geo)
    ngbs = automol.graph.atoms_sorted_neighbor_atom_keys(gra)

    check_dct = {
        'dist': 3.5e-1,
        'coulomb': 1.5e-2,
    }
    _, saved_geos, _ = _saved_cnf_info(
        cnf_save_fs, thy_info)
    frag_saved_geos = []
    for geoi in saved_geos:
        frag_saved_geos.append(
            automol.geom.ring_fragments_geometry(geoi, rings_atoms, ngbs))

    # Make sample zmas
    unique_geos, unique_frag_geos, unique_zmas = [], [], []
    for ring_atoms, samp_range_dct in tors_dcts:
        ring_atoms = [int(idx)-1 for idx in ring_atoms.split('-')]
        dist_value_dct = automol.zmat.ring_distances(zma, ring_atoms)
        nsamp = _num_samp_zmas(ring_atoms, nsamp_par)
        samp_zmas = automol.zmat.samples(zma, nsamp, samp_range_dct)
        for samp_zma in samp_zmas:
            if automol.zmat.ring_distances_reasonable(
                    samp_zma, ring_atoms, dist_value_dct):
                samp_geo = automol.zmat.geometry(samp_zma)
                frag_samp_geo = automol.geom.ring_fragments_geometry(
                    samp_geo, rings_atoms, ngbs)
                if automol.geom.ring_angles_reasonable(samp_geo, ring_atoms):
                    if not automol.pot.low_repulsion_struct(geo, samp_geo):
                        frag_samp_unique = automol.geom.is_unique(
                            frag_samp_geo, frag_saved_geos, check_dct)
                        samp_unique = automol.geom.is_unique(
                            frag_samp_geo, unique_frag_geos, check_dct)
                        if frag_samp_unique:
                            if samp_unique:
                                unique_zmas.append(samp_zma)
                                unique_geos.append(samp_geo)
                                unique_frag_geos.append(frag_samp_geo)

    # Set the samples
    nsamp = len(unique_zmas)
    nsamp0 = nsamp
    nsampd = _calc_nsampd(cnf_save_fs, cnf_run_fs)

    tot_samp = nsamp - nsampd
    info_message(
        ' - Number of samples that have been currently run:', nsampd)
    info_message(' - Number of samples requested:', nsamp)

    if nsamp-nsampd > 0:
        info_message(
            f'Running {nsamp-nsampd} samples...', newline=1)
    samp_idx = 1

    for samp_zma in unique_zmas:
        nsamp = nsamp0 - nsampd
        # Break the while loop if enough sampls completed
        if nsamp <= 0:
            info_message(
                'Requested number of samples have been completed.',
                'Conformer search complete.')
            break

        # Run the conformer sampling
        samp_geo = automol.zmat.geometry(samp_zma)
        rid = autofile.schema.generate_new_ring_id()
        cid = autofile.schema.generate_new_conformer_id()
        locs = (rid, cid)

        cnf_run_fs[-1].create(locs)
        cnf_run_path = cnf_run_fs[-1].path(locs)
        run_fs = autofile.fs.run(cnf_run_path)

        info_message(f"\nSample {samp_idx}/{tot_samp}")
        tors_names = tuple(set(names
                               for tors_dct in ring_tors_dct.values()
                               for names in tors_dct.keys()))
        if two_stage and tors_names:
            frozen_coords_lst = (tors_names, ())
            success, ret = es_runner.multi_stage_optimization(
                script_str=script_str,
                run_fs=run_fs,
                geo=samp_zma,
                spc_info=spc_info,
                thy_info=thy_info,
                frozen_coords_lst=frozen_coords_lst,
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
                overwrite=overwrite,
                saddle=bool(zrxn is not None),
                retryfail=retryfail,
                **kwargs
            )

        # save function added here
        if success:
            _save_conformer(
                ret, cnf_save_fs, locs, thy_info,
                zrxn=zrxn, orig_ich=spc_info[0], rid_traj=False,
                init_zma=samp_zma)

            nsampd = _calc_nsampd(cnf_save_fs, cnf_run_fs)
            nsampd += 1
            samp_idx += 1
            inf_obj.nsamp = nsampd
            cnf_save_fs[0].file.info.write(inf_obj)
            cnf_run_fs[0].file.info.write(inf_obj)


def _calc_nsamp(tors_names, nsamp_par, zma, zrxn=None):
    """ Determine the number of samples to od
    """

    if not any(tors_names):
        info_message(
            " - No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1
        tors_range_dct = None
    else:
        if zrxn is None:
            gra = automol.zmat.graph(zma)
            ntaudof = len(automol.graph.rotational_bond_keys(
               gra, with_h_rotors=False))
        else:
            ntaudof = len(automol.reac.rotational_bond_keys(
                zrxn, zma, with_h_rotors=False))
            # ntaudof = len(tors_names)
        nsamp = util.nsamp_init(nsamp_par, ntaudof)

        tors_ranges = tuple((0, 2*numpy.pi) for tors in tors_names)
        tors_range_dct = dict(zip(tors_names, tors_ranges))

    return nsamp, tors_range_dct


def _calc_nsampd(cnf_save_fs, cnf_run_fs, rid=None):
    """ Determine the number of samples completed
    """

    if rid is None:
        cnf_save_fs[0].create()
        if cnf_save_fs[0].file.info.exists():
            inf_obj_s = cnf_save_fs[0].file.info.read()
            nsampd = inf_obj_s.nsamp
        elif cnf_run_fs[0].file.info.exists():
            inf_obj_r = cnf_run_fs[0].file.info.read()
            nsampd = inf_obj_r.nsamp
        else:
            nsampd = 0
    else:
        cnf_save_fs[1].create([rid])
        if cnf_save_fs[1].file.info.exists([rid]):
            inf_obj_s = cnf_save_fs[1].file.info.read([rid])
            nsampd = inf_obj_s.nsamp
        elif cnf_run_fs[1].file.info.exists([rid]):
            inf_obj_r = cnf_run_fs[1].file.info.read([rid])
            nsampd = inf_obj_r.nsamp
        else:
            nsampd = 0

    return nsampd


def _presamp_save(spc_info, cnf_run_fs, cnf_save_fs,
                  thy_info, zrxn=None, rid=None):
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
                        _save_conformer(
                            ret, cnf_save_fs, locs, thy_info,
                            zrxn=zrxn, orig_ich=spc_info[0],
                            init_zma=init_zma)

        # Update the conformer trajectory file
        print('')
        filesys.mincnf.traj_sort(cnf_save_fs, thy_info, rid=rid)


def _save_conformer(ret, cnf_save_fs, locs, thy_info, zrxn=None,
                    orig_ich='', rid_traj=False, init_zma=None):
    """ save the conformers that have been found so far
          # Only go through save procedure if conf not in save
          # may need to get geo, ene, etc; maybe make function
    """

    saved_locs, saved_geos, saved_enes = _saved_cnf_info(
        cnf_save_fs, thy_info)

    inf_obj, _, out_str = ret
    prog = inf_obj.prog
    method = inf_obj.method
    ene = elstruct.reader.energy(prog, method, out_str)
    geo = elstruct.reader.opt_geometry(prog, out_str)
    zma = filesys.save.read_job_zma(ret, init_zma=init_zma)

    # Assess if geometry is properly connected
    viable = _geo_connected(geo, zrxn)
    if viable:
        if zrxn:
            viable = _ts_geo_viable(
                zma, zrxn, cnf_save_fs, thy_info)
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
                    ret, None, cnf_save_fs, thy_info[1:], zrxn=zrxn,
                    rng_locs=(locs[0],), tors_locs=(locs[1],))
            else:
                sym_locs = saved_locs[sym_id]
                filesys.save.sym_indistinct_conformer(
                    geo, cnf_save_fs, locs, sym_locs)

        # Update the conformer trajectory file
        obj('vspace')
        rid = None
        if rid_traj:
            rid = locs[0]
        filesys.mincnf.traj_sort(cnf_save_fs, thy_info, rid=rid)


def _saved_cnf_info(cnf_save_fs, mod_thy_info):
    """ get the locs, geos and enes for saved conformers
    """

    saved_locs = list(cnf_save_fs[-1].existing())
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
    job = elstruct.Job.OPTIMIZATION
    for locs in cnf_run_fs[-1].existing():
        cnf_run_path = cnf_run_fs[-1].path(locs)
        run_fs = autofile.fs.run(cnf_run_path)
        inf_obj = run_fs[-1].file.info.read([job])
        status = inf_obj.status
        if status == autofile.schema.RunStatus.RUNNING:
            start_time = inf_obj.utc_start_time
            current_time = autofile.schema.utc_time()
            _time = (current_time - start_time).total_seconds()
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
            print('checking input at ', cnf_path)
            inp_str = cnf_fs[-1].file.geometry_input.read(locs)
            inp_str = inp_str.replace('=', '')
            inf_obj = cnf_fs[-1].file.geometry_info.read(locs)
            prog = inf_obj.prog
            inp_zma = elstruct.reader.inp_zmatrix(prog, inp_str)
            if automol.zmat.almost_equal(inp_zma, zma,
                                         dist_rtol=0.018, ang_atol=.2):
                info_message(
                    'This conformer was already run in {cnf_path}.')
                running = True
                break
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
                    inp_zma = elstruct.reader.inp_zmatrix(prog, inp_str)
                    if automol.zmat.almost_equal(inp_zma, zma,
                                                 dist_rtol=0.018, ang_atol=.2):
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
    ich = automol.geom.inchi(geo)
    assert automol.inchi.is_complete(orig_ich), (
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
        if not orig_ich == automol.geom.inchi(geoi):
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


def _ts_geo_viable(zma, zrxn, cnf_save_fs, mod_thy_info, zma_locs=(0,)):
    """ Perform a series of checks to assess the viability
        of a transition state geometry prior to saving
    """

    # Obtain the min-ene zma and bond keys
    _, cnf_save_path = filesys.mincnf.min_energy_conformer_locators(
        cnf_save_fs, mod_thy_info)
    zma_save_fs = fs.zmatrix(cnf_save_path)
    ref_zma = zma_save_fs[-1].file.zmatrix.read(zma_locs)

    return automol.reac.similar_saddle_point_structure(zma, ref_zma, zrxn)


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
        if frag_locs_geo is None:
            rid = locs[0]
            break
        frag_locs_zma = automol.geom.zmatrix(frag_locs_geo)
        if automol.zmat.almost_equal(frag_locs_zma, frag_zma,
                                     dist_rtol=150., ang_atol=45.):
                                     # for now set to include all ring puckering
                                     # dist_rtol=0.15, ang_atol=.45):
            rid = locs[0]
            break

    return rid
