""" es_runners for conformer
"""

import numpy
import automol
import elstruct
import autofile
from autofile import fs
from routines.es._routines import _util as util
from routines.es import runner as es_runner
from lib import filesys
from lib.structure import geom as geomprep
from lib.phydat import bnd


def conformer_sampling(zma, spc_info,
                       mod_thy_info, thy_save_fs,
                       cnf_run_fs, cnf_save_fs,
                       script_str, overwrite,
                       saddle=False, nsamp_par=(False, 3, 3, 1, 50, 50),
                       tors_names='', dist_info=(),
                       two_stage=False, retryfail=True,
                       rxn_class='', **kwargs):
    """ Find the minimum energy conformer by optimizing from nsamp random
    initial torsional states
    """

    ich = spc_info[0]
    coo_names = []

    # Read the geometry and zma from the ini file system
    if saddle:
        coo_names.append(tors_names)

    tors_ranges = tuple((0, 2*numpy.pi) for tors in tors_names)
    tors_range_dct = dict(zip(tors_names, tors_ranges))
    if not saddle:
        gra = automol.inchi.graph(ich)
        ntaudof = len(
            automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
        nsamp = util.nsamp_init(nsamp_par, ntaudof)
    else:
        ntaudof = len(tors_names)
        nsamp = util.nsamp_init(nsamp_par, ntaudof)

    print('\nSaving any conformers in run filesys...')
    save_conformers(
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
        thy_info=mod_thy_info,
        saddle=saddle,
        dist_info=dist_info,
        rxn_class=rxn_class
    )

    print('\nSampling for more conformers if needed...')
    run_conformers(
        zma=zma,
        spc_info=spc_info,
        thy_info=mod_thy_info,
        nsamp=nsamp,
        tors_range_dct=tors_range_dct,
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
        script_str=script_str,
        overwrite=overwrite,
        saddle=saddle,
        two_stage=two_stage,
        retryfail=retryfail,
        **kwargs,
    )

    print('\nSaving any newly found conformers in run filesys...')
    save_conformers(
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
        thy_info=mod_thy_info,
        saddle=saddle,
        dist_info=dist_info,
        rxn_class=rxn_class
    )

    # Save information about the minimum energy conformer in top directory
    min_cnf_locs = filesys.mincnf.min_energy_conformer_locators(cnf_save_fs)
    if min_cnf_locs:
        geo = cnf_save_fs[-1].file.geometry.read(min_cnf_locs)
        if not saddle:
            assert automol.zmatrix.almost_equal(zma, automol.geom.zmatrix(geo))
            thy_save_fs[-1].file.geometry.write(geo, mod_thy_info[1:4])
        else:
            thy_save_fs[0].file.geometry.write(geo)


def single_conformer(zma, spc_info, thy_info,
                     thy_save_fs, cnf_run_fs, cnf_save_fs,
                     overwrite, saddle=False, dist_info=()):
    """ generate single optimized geometry for
        randomly sampled initial torsional angles
    """
    opt_script_str, _, kwargs, _ = es_runner.par.run_qchem_par(*thy_info[0:2])
    conformer_sampling(
        zma=zma,
        spc_info=spc_info,
        mod_thy_info=thy_info,
        thy_save_fs=thy_save_fs,
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
        script_str=opt_script_str,
        overwrite=overwrite,
        nsamp_par=[False, 0, 0, 0, 0, 1],
        saddle=saddle,
        dist_info=dist_info,
        two_stage=saddle,
        retryfail=False,
        **kwargs,
    )


def run_conformers(
        zma, spc_info, thy_info, nsamp, tors_range_dct,
        cnf_run_fs, cnf_save_fs, script_str, overwrite,
        saddle, two_stage, retryfail,
        **kwargs):
    """ run sampling algorithm to find conformers
    """
    if not tors_range_dct:
        print(" - No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1

    cnf_save_fs[0].create()
    vma = automol.zmatrix.var_(zma)
    if cnf_save_fs[0].file.vmatrix.exists():
        existing_vma = cnf_save_fs[0].file.vmatrix.read()
        assert vma == existing_vma
    cnf_save_fs[0].file.vmatrix.write(vma)
    nsamp0 = nsamp
    inf_obj = autofile.schema.info_objects.conformer_trunk(0, tors_range_dct)
    if cnf_save_fs[0].file.info.exists():
        inf_obj_s = cnf_save_fs[0].file.info.read()
        nsampd = inf_obj_s.nsamp
    elif cnf_run_fs[0].file.info.exists():
        inf_obj_r = cnf_run_fs[0].file.info.read()
        nsampd = inf_obj_r.nsamp
    else:
        nsampd = 0

    tot_samp = nsamp - nsampd
    print(' - Number of samples that have been currently run:', nsampd)
    print(' - Number of samples requested:', nsamp)

    if nsamp-nsampd > 0:
        print('\nRunning {} samples...'.format(nsamp-nsampd))
    samp_idx = 1
    while True:
        nsamp = nsamp0 - nsampd
        # Break the while loop if enough sampls completed
        if nsamp <= 0:
            print('Requested number of samples have been completed. '
                  'Conformer search complete.')
            break

        # Run the conformer sampling
        if nsampd > 0:
            samp_zma, = automol.zmatrix.samples(zma, 1, tors_range_dct)
        else:
            samp_zma = zma

        cid = autofile.schema.generate_new_conformer_id()
        locs = [cid]

        cnf_run_fs[-1].create(locs)
        cnf_run_path = cnf_run_fs[-1].path(locs)
        run_fs = autofile.fs.run(cnf_run_path)

        print("Run {}/{}".format(samp_idx, tot_samp))
        tors_names = list(tors_range_dct.keys())
        if two_stage and tors_names:
            print('Stage one beginning, holding the coordinates constant',
                  tors_names)
            es_runner.run_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                run_fs=run_fs,
                geom=samp_zma,
                spc_info=spc_info,
                thy_info=thy_info,
                overwrite=overwrite,
                frozen_coordinates=[tors_names],
                saddle=saddle,
                retryfail=retryfail,
                **kwargs
            )
            print('Stage one success, reading for stage 2')
            ret = es_runner.read_job(
                job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if ret:
                sinf_obj, _, out_str = ret
                prog = sinf_obj.prog
                samp_zma = elstruct.reader.opt_zmatrix(prog, out_str)
                print('Stage one success beginning stage two on', samp_zma)
                es_runner.run_job(
                    job=elstruct.Job.OPTIMIZATION,
                    script_str=script_str,
                    run_fs=run_fs,
                    geom=samp_zma,
                    spc_info=spc_info,
                    thy_info=thy_info,
                    overwrite=overwrite,
                    saddle=saddle,
                    retryfail=False,
                    **kwargs
                )
        else:
            es_runner.run_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                run_fs=run_fs,
                geom=samp_zma,
                spc_info=spc_info,
                thy_info=thy_info,
                overwrite=overwrite,
                saddle=saddle,
                retryfail=retryfail,
                **kwargs
            )

        if cnf_save_fs[0].file.info.exists():
            inf_obj_s = cnf_save_fs[0].file.info.read()
            nsampd = inf_obj_s.nsamp
        elif cnf_run_fs[0].file.info.exists():
            inf_obj_r = cnf_run_fs[0].file.info.read()
            nsampd = inf_obj_r.nsamp
        nsampd += 1
        samp_idx += 1
        inf_obj.nsamp = nsampd
        cnf_save_fs[0].file.info.write(inf_obj)
        cnf_run_fs[0].file.info.write(inf_obj)


def save_conformers(cnf_run_fs, cnf_save_fs, thy_info, saddle=False,
                    dist_info=(), rxn_class=''):
    """ save the conformers that have been found so far
        # Only go through save procedure if conf not in save
        # may need to get geo, ene, etc; maybe make function
    """

    locs_lst = cnf_save_fs[-1].existing()
    seen_geos = [cnf_save_fs[-1].file.geometry.read(locs)
                 for locs in locs_lst]
    seen_enes = [cnf_save_fs[-1].file.energy.read(locs)
                 for locs in locs_lst]

    if not cnf_run_fs[0].exists():
        print(" - No conformers in run filesys to save.")
    else:
        print(" - Found conformers in run filesys to save.\n")
        for locs in cnf_run_fs[-1].existing():
            cnf_run_path = cnf_run_fs[-1].path(locs)
            run_fs = autofile.fs.run(cnf_run_path)
            print("Reading from conformer run at {}".format(cnf_run_path))

            # Read the electronic structure optimization job
            ret = es_runner.read_job(
                job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)

            # Assess the geometry and save it if so
            if ret:
                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)
                geo = elstruct.reader.opt_geometry(prog, out_str)
                # zma = elstruct.reader.opt_zmatrix(prog, out_str)
                # zma = automol.geom.zmatrix(geo)

                # Assess if geometry is properly connected
                if _geo_connected(geo, saddle):

                    # Assess viability of transition state conformer
                    if saddle:
                        if not _ts_geo_viable():
                            continue

                    # Determine uniqueness of conformer, save if needed
                    if _geo_unique(geo, ene, seen_geos, seen_enes, saddle):
                        if _is_proper_isomer():
                            _save_conformer(cnf_save_fs, ret)
                    else:
                        if not _sym_distinct(geo, zma, cnf_save_path, cid):
                            _save_sym_indistinct_conformer(
                                geo, zma, cnf_save_path, cid)

                    # Append to list of all geos and enes in filesys
                    seen_geos.append(geo)
                    seen_enes.append(ene)

        # Update the conformer trajectory file
        print('')
        filesys.mincnf.traj_sort(cnf_save_fs)


def _geo_connected(geo, saddle):
    """ Assess if geometry is connected. Right now only works for 
        minima
    """

    # Determine connectivity (only for minima)
    if not saddle:
        gra = automol.geom.graph(geo)
        conns = automol.graph.connected_components(gra)
        lconns = len(conns)
    else:
        lconns = 1

    # Check connectivity
    if lconns == 1:
        connected = True
    else:
        print(" - Geometry is disconnected. Conformer will not be saved.")
        connected = False

    return connected


def _geo_unique(geo, ene, seen_geos, seen_enes, saddle):
    """ Assess if a geometry is unique to saved geos
    """

    unique = geomprep.is_unique_tors_dist_mat_energy(
        geo, ene, seen_geos, seen_enes, saddle)
    if not unique:
        print(" - Geometry is not unique. Conformer will not be saved.")

    return unique


def _sym_unique(geo, ene, saved_geos, sved_enes, ethresh=1.0e-5):
    """ Check if a conformer is symmetrically distinct from the
        existing conformers in the filesystem
    """

    # geo_sym, geo_sym2, ene_sym = [], [], []
    for geoi, enei in zip(saved_geos, saved_enes):
        if enei - enes[0] < ethresh:
            unique = is_unique_coulomb_energy(geo, ene, [geoi], [enei])
            if not unique:
                geo_sim.append(geoi)
                ene_sim.append(enei)


def _is_proper_isomer():
    """ Check if geom is the same isomer as those in the filesys
    """
    if cnf_save_fs[0].file.vmatrix.exists():
        exist_vma = cnf_save_fs[0].file.vmatrix.read()
        if vma != exist_vma:
            print(" - Isomer is not the same as starting isomer. Skipping...")
            proper_isomer = False
        else:
            proper_isomer = True
    else:
        proper_isomer = False

    return proper_isomer


def _ts_geo_viable():
    """ Perform a series of checks to assess the viability
        of a transition state geometry prior to saving
    """

    # ts_class, ts_original_zma, ts_tors_names,
    # ts_dist_info
    # geo, zma, final_dist = check_filesys_for_ts(
    #     ts_dct, ts_zma, cnf_save_fs, overwrite,
    #     typ, dist_info, dist_name, bkp_ts_class_data)
    # zma = cnf_save_fs[-1].file.zmatrix.read(
    # cnf_save_locs)

    # # Add an angle check which is added
    # to spc dct for TS (crap code...)
    # vals = automol.zmatrix.values(zma)
    # final_dist = vals[dist_name]
    # dist_info[1] = final_dist
    # angle = ts.chk.check_angle(
    #     ts_dct['zma'],
    #     ts_dct['dist_info'],
    #     ts_dct['class'])
    # ts_dct['dist_info'][1] = final_dist
    # ts_dct['dist_info'].append(angle)
    zma = elstruct.reader.opt_zmatrix(prog, out_str)
    dist_name = dist_info[0]
    dist_len = dist_info[1]
    ts_bnd = automol.zmatrix.bond_idxs(zma, dist_name)
    ts_bnd1 = min(ts_bnd)
    ts_bnd2 = max(ts_bnd)
    conf_dist_len = automol.zmatrix.values(zma)[dist_name]
    brk_name = dist_info[3]
    cent_atm = None
    ldist = len(dist_info)
    # print('zma test:/n', automol.zmatrix.string(zma))
    # print('ldist test:', ldist, dist_name, brk_name)
    if dist_name and brk_name and ldist > 4:
        angle = dist_info[4]
        brk_bnd = automol.zmatrix.bond_idxs(zma, brk_name)
        ang_atms = [0, 0, 0]
        cent_atm = list(set(brk_bnd) & set(ts_bnd))
        if cent_atm:
            ang_atms[1] = cent_atm[0]
            for idx in brk_bnd:
                if idx != ang_atms[1]:
                    ang_atms[0] = idx
            for idx in ts_bnd:
                if idx != ang_atms[1]:
                    ang_atms[2] = idx
            geom = automol.zmatrix.geometry(zma)
            conf_ang = automol.geom.central_angle(
                geom, *ang_atms)
    max_disp = 0.6
    if 'addition' in rxn_class:
        max_disp = 0.8
    if 'abstraction' in rxn_class:
        max_disp = 1.4

    # Check forming bond angle similar to ini config
    if cent_atm and 'elimination' not in rxn_class:
        if abs(conf_ang - angle) > .44:
            print(" - Transition State conformer has",
                  "diverged from original structure of",
                  "angle {:.3f} with angle {:.3f}".format(
                      angle, conf_ang))
            viable = False

    # Check if radical atom is closer to some atom other than the bonding atom
    if 'add' in rxn_class or 'abst' in rxn_class:
        cls = geomprep.is_atom_closest_to_bond_atom(
            zma, ts_bnd2, conf_dist_len)
        if not cls:
            print(" - Transition State conformer has",
                  "diverged from original structure of",
                  "dist {:.3f} with dist {:.3f}".format(
                      dist_len, conf_dist_len))
            print(' - Radical atom now has a new nearest neighbor')
            viable = False
        if abs(conf_dist_len - dist_len) > max_disp:
            print(" - Transition State conformer has",
                  "diverged from original structure of",
                  "dist {:.3f} with dist {:.3f}".format(
                      dist_len, conf_dist_len))
            viable = False

        # Check distances
        symbols = automol.zmatrix.symbols(zma)
        symb1, symb2 = symbols[ts_bnd1], symbols[ts_bnd2]
        if (symb1, symb2) in bnd.LEN_DCT:
            equi_bnd = bnd.LEN_DCT[(symb1, symb2)]
        elif (symb2, symb1) in bnd.LEN_DCT:
            equi_bnd = bnd.LEN_DCT[(symb2, symb1)]
        else:
            equi_bnd = 0.0
        displace_from_equi = conf_dist_len - equi_bnd
        dchk1 = abs(conf_dist_len - dist_len) > 0.2
        dchk2 = displace_from_equi < 0.2
        if dchk1 and dchk2:
            print(" - Transition State conformer has",
                  "converged to an",
                  "equilibrium structure with dist",
                  " {:.3f} comp with equil {:.3f}".format(
                      conf_dist_len, equi_bnd))
            viable = False
    else:
        if abs(conf_dist_len - dist_len) > 0.4:
            print(" - Transition State conformer has",
                  "diverged from original structure of",
                  "dist {:.3f} with dist {:.3f}".format(
                      dist_len, conf_dist_len))
            viable = False
            
    return viable


def _save_unique_conformer(cnf_save_fs, ret):
    """ Save the conformer in the filesystem 
    """
    
    # Set the path to the conformer save filesystem
    save_path = cnf_save_fs[-1].path(locs)

    # Unpack the ret object
    inf_obj, inp_str, out_str = ret
    prog = inf_obj.prog
    method = inf_obj.method
    ene = elstruct.reader.energy(prog, method, out_str)
    geo = elstruct.reader.opt_geometry(prog, out_str)

    # Build the conformer filesystem and save the structural info
    print(" - Geometry is unique. Saving...")
    print(" - Save path: {}".format(save_path))
    cnf_save_fs[-1].create(locs)
    cnf_save_fs[-1].file.geometry_info.write(inf_obj, locs)
    cnf_save_fs[-1].file.geometry_input.write(inp_str, locs)
    cnf_save_fs[-1].file.energy.write(ene, locs)
    cnf_save_fs[-1].file.geometry.write(geo, locs)

    # Build the zma filesystem and save the z-matrix
    zma_save_fs = fs.manager(save_path, 'ZMATRIX')
    zma_save_fs[-1].create([0])
    zma_save_fs[-1].file.zmatrix.write(zma, [0])

    # Saving the energy to a SP filesystem
    print(" - Saving energy of unique geometry...")
    sp_save_fs = autofile.fs.single_point(save_path)
    sp_save_fs[-1].create(thy_info[1:4])
    sp_save_fs[-1].file.input.write(inp_str, thy_info[1:4])
    sp_save_fs[-1].file.info.write(inf_obj, thy_info[1:4])
    sp_save_fs[-1].file.energy.write(ene, thy_info[1:4])


def _save_sym_indistinct_conformer(geo, zma, cnf_save_path, cid):
    """ Save a structure into the SYM directory of a conformer
    """

    # Build the sym file sys
    print(" - Saving structure in a sym directory...")
    sym_save_fs = fs.manager(cnf_save_path, 'SYMMETRIC')
    sym_save_fs[-1].create(cid)
    sym_save_fs[-1].file.geometry.write(geo, cid)
    sym_save_fs[-1].file.zmatrix.write(zma, cid)
