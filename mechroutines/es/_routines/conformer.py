""" es_runners for conformer
"""

import shutil
import numpy
import automol
import elstruct
import autofile
from autofile import fs
from routines.es._routines import _util as util
from routines.es import runner as es_runner
from lib import filesys
from lib.structure import geom as geomprep
from lib.structure import ts as tsprep
from lib.phydat import bnd


# Initial conformer
def initial_conformer(spc_dct_i, spc_info,
                      mod_thy_info,
                      thy_run_fs, thy_save_fs,
                      cnf_save_fs,
                      ini_thy_save_path, mod_ini_thy_info,
                      run_fs,
                      opt_script_str, overwrite,
                      kickoff_size=0.1, kickoff_backward=False,
                      **opt_kwargs):
    """ determine what to use as the reference geometry for all future runs
    If ini_thy_info refers to geometry dictionary then use that,
    otherwise values are from a hierarchy of:
    running level of theory, input level of theory, inchis.
    From the hierarchy an optimization is performed followed by a check for
    an imaginary frequency and then a conformer file system is set up.
    """

    # Initialize empty return
    ret = None

    if run_fs[0].file.info.exists([]):
        inf_obj = run_fs[0].file.info.read([])
        if inf_obj.status == autofile.schema.RunStatus.RUNNING:
            print('Reference geometry already running',
                  'in {}'.format(run_fs[0].path([])))
            return ret
    else:
        [prog, method, basis, _] = mod_thy_info
        status = autofile.schema.RunStatus.RUNNING
        inf_obj = autofile.schema.info_objects.run(
            job='', prog=prog, version='version', method=method, basis=basis,
            status=status)
        run_fs[0].file.info.write(inf_obj, [])

    exists = thy_save_fs[-1].file.geometry.exists(mod_thy_info[1:4])
    if not exists:
        print('No energy found in save filesys. Running energy...')
        _run = True
    elif overwrite:
        print('User specified to overwrite energy with new run...')
        _run = True
    else:
        _run = False

    if _run:
        print('Obtaining some initial guess geometry.')
        geo_init = _obtain_ini_geom(spc_dct_i,
                                    ini_thy_save_path, mod_ini_thy_info,
                                    overwrite)

        if geo_init is not None:
            print('Assessing if there are any functional groups',
                  'that cause instability')
            print('geo str\n', automol.geom.string(geo_init))
            if _functional_groups_stable(geo_init, thy_save_fs, mod_thy_info):
                zma_init = automol.geom.zmatrix(geo_init)
                if not automol.geom.is_atom(geo_init):
                    geo_found = _optimize_molecule(
                        spc_info, zma_init,
                        mod_thy_info, thy_run_fs, thy_save_fs,
                        cnf_save_fs,
                        run_fs,
                        opt_script_str, overwrite,
                        kickoff_size=kickoff_size,
                        kickoff_backward=kickoff_backward,
                        **opt_kwargs)
                else:
                    geo_found = _optimize_atom(
                        spc_info, zma_init,
                        mod_thy_info, thy_run_fs,
                        cnf_save_fs, run_fs,
                        overwrite, opt_script_str, **opt_kwargs)
            else:
                geo_found = True
                print('Found functional groups that cause instabilities')
        else:
            geo_found = False
            print('Unable to obtain an initial guess geometry')
    else:
        geo_found = True
        thy_path = thy_save_fs[-1].path(mod_thy_info[1:4])
        print('Initial geometry found and saved previously at {}'.format(
            thy_path))

    # Write the job status into the run filesystem
    if geo_found:
        inf_obj.status = autofile.schema.RunStatus.SUCCESS
        run_fs[0].file.info.write(inf_obj, [])
    else:
        inf_obj.status = autofile.schema.RunStatus.FAILURE
        run_fs[0].file.info.write(inf_obj, [])

    return geo_found


def _obtain_ini_geom(spc_dct_i, ini_thy_save_path,
                     mod_ini_thy_info, overwrite):
    """ Obtain an initial geometry to be optimized. Checks a hieratchy
        of places to obtain the initial geom.
            (1) Geom dict which is the input from the user
            (2) Geom from inchi
    """

    geo_init = None

    # Obtain geom from thy fs or remove the conformer filesystem if needed
    if not overwrite:
        ini_cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)
        ini_min_cnf_locs, _ = filesys.mincnf.min_energy_conformer_locators(
            ini_cnf_save_fs, mod_ini_thy_info)
        if ini_min_cnf_locs:
            geo_init = ini_cnf_save_fs[-1].file.geometry.read(ini_min_cnf_locs)
            print('Getting inital geometry from inplvl at path',
                  '{}'.format(ini_cnf_save_fs[-1].path(ini_min_cnf_locs)))
    else:
        print('Removing original conformer save data for instability')
        ini_cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)
        for locs in ini_cnf_save_fs[-1].existing():
            cnf_path = ini_cnf_save_fs[-1].path(locs)
            print('Removing {}'.format(cnf_path))
            shutil.rmtree(cnf_path)

    if geo_init is None:
        if 'geo_inp' in spc_dct_i:
            geo_init = spc_dct_i['geo_inp']
            print('Getting initial geometry from geom dictionary')

    if geo_init is None:
        geo_init = automol.inchi.geometry(spc_dct_i['inchi'])
        print('Getting initial geometry from inchi')

    # Check if the init geometry is connected
    if geo_init is not None:
        if not automol.geom.connected(geo_init):
            geo_init = None

    return geo_init


def _optimize_atom(spc_info, zma_init,
                   mod_thy_info, thy_run_fs,
                   cnf_save_fs, run_fs,
                   overwrite, opt_script_str, **opt_kwargs):
    """ Deal with an atom separately
    """

    geo, zma = _init_geom_opt(zma_init, spc_info, mod_thy_info,
                              run_fs, thy_run_fs,
                              opt_script_str, overwrite, **opt_kwargs)

    if geo is not None and zma is not None:
        locs = [autofile.schema.generate_new_conformer_id()]
        job = elstruct.Job.OPTIMIZATION
        filesys.save_struct(run_fs, cnf_save_fs, locs, job, mod_thy_info,
                            zma_locs=(0,), in_zma_fs=False)
        conf_found = True
    else:
        conf_found = False

    return conf_found


def _optimize_molecule(spc_info, zma_init,
                       mod_thy_info, thy_run_fs, thy_save_fs,
                       cnf_save_fs,
                       run_fs,
                       opt_script_str, overwrite,
                       kickoff_size=0.1, kickoff_backward=False,
                       **opt_kwargs):
    """ Optimize a proper geometry
    """

    # Optimize the initial geometry
    geo, zma = _init_geom_opt(
        zma_init, spc_info, mod_thy_info, run_fs, thy_run_fs,
        opt_script_str, overwrite, **opt_kwargs)

    # If connected, check for imaginary modes and fix them if possible
    if automol.geom.connected(geo):

        # Remove the imaginary mode
        geo, imag_fix_needed = remove_imag(
            geo, spc_info, mod_thy_info, thy_run_fs,
            run_fs, kickoff_size, kickoff_backward, kickoff_mode=0,
            overwrite=overwrite)

        # Recheck connectivity for imag-checked geometry
        if geo is not None:

            conf_found = True
            if automol.geom.connected(geo):

                print('\nSaving structure as the first conformer...')
                locs = [autofile.schema.generate_new_conformer_id()]
                job = elstruct.Job.OPTIMIZATION
                filesys.save_struct(
                    run_fs, cnf_save_fs, locs, job, mod_thy_info,
                    zma_locs=(0,), in_zma_fs=False,
                    cart_to_zma=imag_fix_needed)

            else:

                print('Saving disconnected species...')
                _, opt_ret = es_runner.read_job(
                    job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
                structure.instab.write_instab(
                    zma_init, zma,
                    thy_save_fs, mod_thy_info[1:4],
                    opt_ret,
                    zma_locs=(0,),
                    save_cnf=True
                )

        else:

            print('\n No geom found...')
            conf_found = False

    else:

        print('Saving disconnected species...')
        conf_found = False
        _, opt_ret = es_runner.read_job(
            job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
        structure.instab.write_instab(
            zma_init, zma,
            thy_save_fs, mod_thy_info[1:4],
            opt_ret,
            zma_locs=(0,),
            save_cnf=True
        )

    # Save geom in thy filesys if a good geom is found
    if conf_found:
        thy_save_fs[-1].create(mod_thy_info[1:4])
        thy_save_path = thy_save_fs[-1].path(mod_thy_info[1:4])
        thy_save_fs[-1].file.geometry.write(geo, mod_thy_info[1:4])

        print('Saving reference geometry')
        print(" - Save path: {}".format(thy_save_path))

    return conf_found


def _init_geom_opt(zma_init, spc_info, mod_thy_info,
                   run_fs, thy_run_fs,
                   opt_script_str, overwrite, **opt_kwargs):
    """ Generate initial geometry via optimization from either reference
    geometries or from inchi
    """

    # Set up the filesystem
    thy_run_fs[-1].create(mod_thy_info[1:4])
    thy_run_path = thy_run_fs[-1].path(mod_thy_info[1:4])

    # Call the electronic structure optimizer
    run_fs = autofile.fs.run(thy_run_path)
    success, ret = execute_job(
        job=elstruct.Job.OPTIMIZATION,
        script_str=opt_script_str,
        run_fs=run_fs,
        geom=zma_init,
        spc_info=spc_info,
        thy_info=mod_thy_info,
        overwrite=overwrite,
        **opt_kwargs
    )

    geo, zma = None, None
    if success:
        print('Succesful reference geometry optimization')
        inf_obj, _, out_str = ret
        prog = inf_obj.prog
        geo = elstruct.reader.opt_geometry(prog, out_str)
        zma = elstruct.reader.opt_zmatrix(prog, out_str)

    return geo, zma


def conformer_sampling(zma, spc_info,
                       mod_thy_info, thy_save_fs,
                       cnf_run_fs, cnf_save_fs,
                       script_str, overwrite,
                       saddle=False, nsamp_par=(False, 3, 3, 1, 50, 50),
                       tors_names='',
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
    print('tors_names', tors_names)
    print('tors_range_dct', tors_range_dct)
    if not saddle:
        gra = automol.inchi.graph(ich)
        ntaudof = len(
            automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
    else:
        ntaudof = len(tors_names)
    nsamp = util.nsamp_init(nsamp_par, ntaudof)

    # Check samples and if nsamp met and no resave

    print('\nSaving any conformers in run filesys...')
    save_conformers(
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
        thy_info=mod_thy_info,
        saddle=saddle,
        rxn_class=rxn_class,
        orig_ich=ich
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
        rxn_class=rxn_class,
        orig_ich=ich
    )

    # Save information about the minimum energy conformer in top directory
    min_cnf_locs, _ = filesys.mincnf.min_energy_conformer_locators(
        cnf_save_fs, mod_thy_info)
    if min_cnf_locs:
        print('min_cnf_locs test in save_conformer:', min_cnf_locs)
        geo = cnf_save_fs[-1].file.geometry.read(min_cnf_locs)
        if not saddle:
            thy_save_fs[-1].file.geometry.write(geo, mod_thy_info[1:4])
        else:
            thy_save_fs[0].file.geometry.write(geo)

    return bool(min_cnf_locs)


def single_conformer(zma, spc_info, mod_thy_info,
                     cnf_run_fs, cnf_save_fs,
                     script_str, overwrite,
                     retryfail=True, saddle=False, **kwargs):
    """ generate single optimized geometry to be saved into a
        filesystem
    """

    # Build the filesystem
    locs = [autofile.schema.generate_new_conformer_id()]
    cnf_run_fs[-1].create(locs)
    cnf_run_path = cnf_run_fs[-1].path(locs)
    run_fs = autofile.fs.run(cnf_run_path)

    # Run the optimization
    print('Optimizing a single conformer')
    success, ret = execute_job(
        job=elstruct.Job.OPTIMIZATION,
        script_str=script_str,
        run_fs=run_fs,
        geom=zma,
        spc_info=spc_info,
        thy_info=mod_thy_info,
        overwrite=overwrite,
        frozen_coordinates=(),
        saddle=saddle,
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

        if _geo_unique(geo, ene, saved_geos, saved_enes, saddle):
            sym_id = _sym_unique(
                geo, ene, saved_geos, saved_enes)
            if sym_id is None:
                _save_unique_conformer(
                    ret, mod_thy_info, cnf_save_fs, locs,
                    saddle=saddle, zma_locs=(0,))
                saved_geos.append(geo)
                saved_enes.append(ene)
                saved_locs.append(locs)

        # Update the conformer trajectory file
        print('')
        filesys.mincnf.traj_sort(cnf_save_fs, mod_thy_info)


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

        print('\nChecking if ZMA has high repulsion...')
        bad_geom_count = 0
        while (not automol.intmol.low_repulsion_struct(zma, samp_zma) and
               bad_geom_count < 1000):
            print('  ZMA has high repulsion.')
            print('\n  Generating new sample ZMA')
            samp_zma, = automol.zmatrix.samples(zma, 1, tors_range_dct)
            bad_geom_count += 1
        print('  ZMA is fine...')

        cid = autofile.schema.generate_new_conformer_id()
        locs = [cid]

        cnf_run_fs[-1].create(locs)
        cnf_run_path = cnf_run_fs[-1].path(locs)
        run_fs = autofile.fs.run(cnf_run_path)

        print("Run {}/{}".format(samp_idx, tot_samp))
        tors_names = tuple(tors_range_dct.keys())
        if two_stage and tors_names:
            frozen_coords_lst = ((), tors_names)
            _, _ = es_runner.multi_stage_optimization(
                script_str=script_str,
                run_fs=run_fs,
                geom=samp_zma,
                spc_info=spc_info,
                thy_info=thy_info,
                frozen_coords_lst=frozen_coords_lst,
                overwrite=overwrite,
                saddle=saddle,
                retryfail=retryfail,
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
                    rxn_class='', orig_ich=''):
    """ save the conformers that have been found so far
        # Only go through save procedure if conf not in save
        # may need to get geo, ene, etc; maybe make function
    """

    saved_locs, saved_geos, saved_enes = _saved_cnf_info(
        cnf_save_fs, thy_info)

    if not saddle:
        _check_old_inchi(orig_ich, saved_geos, saved_locs, cnf_save_fs)

    if not cnf_run_fs[0].exists():
        print(" - No conformers in run filesys to save.")
    else:
        print(" - Found conformers in run filesys to save.\n")
        for locs in cnf_run_fs[-1].existing():
            cnf_run_path = cnf_run_fs[-1].path(locs)
            run_fs = autofile.fs.run(cnf_run_path)
            print("\nReading from conformer run at {}".format(cnf_run_path))

            # Read the electronic structure optimization job
            success, ret = es_runner.read_job(
                job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)

            # Assess the geometry and save it if so
            if success:
                inf_obj, _, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)
                geo = elstruct.reader.opt_geometry(prog, out_str)
                zma = elstruct.reader.opt_zmatrix(prog, out_str)
                # zma = automol.geom.zmatrix(geo)

                # Assess if geometry is properly connected
                if _geo_connected(geo, saddle):
                    if saddle:
                        ts_viable = _ts_geo_viable(
                            zma, cnf_save_fs, rxn_class, thy_info)
                        if not ts_viable:
                            continue
                    elif not _inchi_are_same(orig_ich, geo):
                        continue
                    # Assess viability of transition state conformer

                    # Determine uniqueness of conformer, save if needed
                    if _geo_unique(geo, ene, saved_geos, saved_enes, saddle):
                        # iso check breaks because of zma location
                        # if _is_proper_isomer(cnf_save_fs, zma):
                        sym_id = _sym_unique(
                            geo, ene, saved_geos, saved_enes)
                        if sym_id is None:
                            _save_unique_conformer(
                                ret, thy_info, cnf_save_fs,
                                locs, saddle=saddle)
                            saved_geos.append(geo)
                            saved_enes.append(ene)
                            saved_locs.append(locs)
                        else:
                            sym_locs = saved_locs[sym_id]
                            _save_sym_indistinct_conformer(
                                geo, cnf_save_fs, locs, sym_locs)

        # Update the conformer trajectory file
        print('')
        filesys.mincnf.traj_sort(cnf_save_fs, thy_info)


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


def _inchi_are_same(orig_ich, geo):
    """ Assess if a geometry has the same connectivity to
     saved geos evaluated in temrs of inchi
    """
    same = False
    ich = automol.geom.inchi(geo)
    assert automol.inchi.is_complete(orig_ich), (
        'the inchi {} orig_ich is not complete'.format(orig_ich))
    if ich == orig_ich:
        same = True
    if not same:
        print(" - new inchi {} not the same as old {}".format(ich, orig_ich))

    return same


def _check_old_inchi(orig_ich, seen_geos, saved_locs, cnf_save_fs):
    """
    This assumes you already have bad geos in your save
    """
    for i, geoi in enumerate(seen_geos):
        if not orig_ich == automol.geom.inchi(geoi):
            smi = automol.geom.smiles(geoi)
            print('ERROR: inchi do not match for {}'.format(smi))
            print(cnf_save_fs[-1].path(saved_locs[i]))


def _sym_unique(geo, ene, saved_geos, saved_enes, ethresh=1.0e-5):
    """ Check if a conformer is symmetrically distinct from the
        existing conformers in the filesystem
    """

    sym_idx = None
    for idx, (geoi, enei) in enumerate(zip(saved_geos, saved_enes)):
        if abs(enei - ene) < ethresh:
            unique = geomprep.is_unique_coulomb_energy(
                geo, ene, [geoi], [enei])
            if not unique:
                sym_idx = idx

    if sym_idx is not None:
        print(' - Structure is not symmetrically unique.')

    return sym_idx


def _is_proper_isomer(cnf_save_fs, zma):
    """ Check if geom is the same isomer as those in the filesys
    """
    vma = automol.zmatrix.var_(zma)
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


def _ts_geo_viable(zma, cnf_save_fs, rxn_class, mod_thy_info, zma_locs=(0,)):
    """ Perform a series of checks to assess the viability
        of a transition state geometry prior to saving
    """

    # Initialize viable
    viable = True

    # Obtain the min-ene zma and bond keys
    min_cnf_locs, cnf_save_path = filesys.mincnf.min_energy_conformer_locators(
        cnf_save_fs, mod_thy_info)
    zma_save_fs = fs.manager(cnf_save_path, 'ZMATRIX')
    ref_zma = zma_save_fs[-1].file.zmatrix.read(zma_locs)

    # Read the form and broken keys from the min conf
    frm_bnd_keys, brk_bnd_keys = tsprep.all_rxn_bnd_keys(
        cnf_save_fs, min_cnf_locs, zma_locs=zma_locs)

    # Calculate the distance of bond being formed
    cnf_geo = automol.zmatrix.geometry(zma)
    ref_geo = automol.zmatrix.geometry(ref_zma)

    cnf_dist_lst = []
    ref_dist_lst = []
    bnd_key_lst = []
    cnf_ang_lst = []
    ref_ang_lst = []
    for frm_bnd_key in frm_bnd_keys:
        print('frm_bnd_key test:', frm_bnd_key)
        frm_idx1, frm_idx2 = list(frm_bnd_key)
        cnf_dist = automol.geom.distance(cnf_geo, frm_idx1, frm_idx2)
        ref_dist = automol.geom.distance(ref_geo, frm_idx1, frm_idx2)
        cnf_dist_lst.append(cnf_dist)
        ref_dist_lst.append(ref_dist)
        bnd_key_lst.append(frm_bnd_key)

    for brk_bnd_key in brk_bnd_keys:
        brk_idx1, brk_idx2 = list(brk_bnd_key)
        cnf_dist = automol.geom.distance(cnf_geo, brk_idx1, brk_idx2)
        ref_dist = automol.geom.distance(ref_geo, brk_idx1, brk_idx2)
        cnf_dist_lst.append(cnf_dist)
        ref_dist_lst.append(ref_dist)
        bnd_key_lst.append(brk_bnd_key)

    for frm_bnd_key in frm_bnd_keys:
        for brk_bnd_key in brk_bnd_keys:
            for frm_idx in frm_bnd_key:
                for brk_idx in brk_bnd_key:
                    if frm_idx == brk_idx:
                        idx2 = frm_idx
                        idx1 = list(frm_bnd_key - frozenset({idx2}))[0]
                        idx3 = list(brk_bnd_key - frozenset({idx2}))[0]
                        cnf_ang = automol.geom.central_angle(
                            cnf_geo, idx1, idx2, idx3)
                        ref_ang = automol.geom.central_angle(
                            ref_geo, idx1, idx2, idx3)
                        cnf_ang_lst.append(cnf_ang)
                        ref_ang_lst.append(ref_ang)

    # Set the maximum allowed displacement for a TS conformer
    max_disp = 0.6
    # better to check for bond-form length in bond scission with ring forming
    if 'addition' in rxn_class:
        max_disp = 0.8
    if 'abstraction' in rxn_class:
        # this was 1.4 - SJK reduced it to work for some OH abstractions
        max_disp = 1.0

    # Check forming bond angle similar to ini config
    if 'elimination' not in rxn_class:
        for ref_angle, cnf_angle in zip(ref_ang_lst, cnf_ang_lst):
            if abs(cnf_angle - ref_angle) > .44:
                print(" - Transition State conformer has",
                      "diverged from original structure of",
                      "angle {:.3f} with angle {:.3f}".format(
                          ref_angle, cnf_angle))
                viable = False

    symbols = automol.geom.symbols(cnf_geo)
    lst_info = zip(ref_dist_lst, cnf_dist_lst, bnd_key_lst)
    for ref_dist, cnf_dist, bnd_key in lst_info:
        if 'add' in rxn_class or 'abst' in rxn_class:
            bnd_key1, bnd_key2 = min(list(bnd_key)), max(list(bnd_key))
            symb1 = symbols[bnd_key1]
            symb2 = symbols[bnd_key2]

            if bnd_key in frm_bnd_keys:
                # Check if radical atom is closer to some atom
                # other than the bonding atom
                cls = geomprep.is_atom_closest_to_bond_atom(
                    zma, bnd_key2, cnf_dist)
                if not cls:
                    print(" - Transition State conformer has",
                          "diverged from original structure of",
                          "dist {:.3f} with dist {:.3f}".format(
                              ref_dist, cnf_dist))
                    print(' - Radical atom now has a new nearest neighbor')
                    viable = False
                # check forming bond distance
                if abs(cnf_dist - ref_dist) > max_disp:
                    print(" - Transition State conformer has",
                          "diverged from original structure of",
                          "dist {:.3f} with dist {:.3f}".format(
                              ref_dist, cnf_dist))
                    viable = False

            # Check distance relative to equi. bond
            if (symb1, symb2) in bnd.LEN_DCT:
                equi_bnd = bnd.LEN_DCT[(symb1, symb2)]
            elif (symb2, symb1) in bnd.LEN_DCT:
                equi_bnd = bnd.LEN_DCT[(symb2, symb1)]
            else:
                equi_bnd = 0.0
            displace_from_equi = cnf_dist - equi_bnd
            dchk1 = abs(cnf_dist - ref_dist) > 0.1
            dchk2 = displace_from_equi < 0.2
            if dchk1 and dchk2:
                print(" - Transition State conformer has",
                      "converged to an",
                      "equilibrium structure with dist",
                      " {:.3f} comp with equil {:.3f}".format(
                          cnf_dist, equi_bnd))
                viable = False
        else:
            # check forming/breaking bond distance
            # if abs(cnf_dist - ref_dist) > 0.4:
            # max disp of 0.4 causes problems for bond scission w/ ring forming
            # not sure if setting it to 0.3 will cause problems for other cases
            if abs(cnf_dist - ref_dist) > 0.3:
                print(" - Transition State conformer has",
                      "diverged from original structure of",
                      "dist {:.3f} with dist {:.3f}".format(
                          ref_dist, cnf_dist))
                viable = False

    return viable


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
        _, inigeo = filesys.inf.cnf_fs_zma_geo(
            ini_cnf_save_fs, ini_locs)
        ini_cnf_save_path = ini_cnf_save_fs[-1].path(ini_locs)
        print('Checking structure from path {}'.format(ini_cnf_save_path))
        for locs in cnf_save_locs_lst:
            _, geo = filesys.inf.cnf_fs_zma_geo(cnf_save_fs, locs)
            if automol.geom.almost_equal_dist_matrix(inigeo, geo, thresh=.15):
                cnf_save_path = cnf_save_fs[-1].path(locs)
                print('- Similar structure found at {}'.format(cnf_save_path))
                found = True
                break

        # If no match was found, add to unique locs lst
        if not found:
            uni_ini_cnf_save_locs.append(ini_locs)

    return uni_ini_cnf_save_locs


def _save_unique_conformer(ret, thy_info, cnf_save_fs, locs,
                           saddle=False, zma_locs=(0,)):
    """ Save the conformer in the filesystem
    """

    # Set the path to the conformer save filesystem
    cnf_save_path = cnf_save_fs[-1].path(locs)

    # Unpack the ret object and obtain the prog and method
    inf_obj, inp_str, out_str = ret
    prog = inf_obj.prog
    method = inf_obj.method

    # Read the energy and geom from the output
    ene = elstruct.reader.energy(prog, method, out_str)
    geo = elstruct.reader.opt_geometry(prog, out_str)
    zma = elstruct.reader.opt_zmatrix(prog, out_str)

    # Read the tra and graph
    if saddle:
        _, ts_min_path = filesys.mincnf.min_energy_conformer_locators(
            cnf_save_fs, thy_info)
        ts_min_zma_fs = fs.manager(ts_min_path, 'ZMATRIX')
        print('ts_min_path test:', ts_min_path)
        tra = ts_min_zma_fs[-1].file.transformation.read(zma_locs)
        print('zma_locs test:', zma_locs)
        rct_gra = ts_min_zma_fs[-1].file.reactant_graph.read(zma_locs)

    # Build the conformer filesystem and save the structural info
    print(" - Geometry is unique. Saving...")
    print(" - Save path: {}".format(cnf_save_path))
    cnf_save_fs[-1].create(locs)
    cnf_save_fs[-1].file.geometry_info.write(inf_obj, locs)
    cnf_save_fs[-1].file.geometry_input.write(inp_str, locs)
    cnf_save_fs[-1].file.energy.write(ene, locs)
    cnf_save_fs[-1].file.geometry.write(geo, locs)

    # Build the zma filesystem and save the z-matrix
    zma_save_fs = fs.manager(cnf_save_path, 'ZMATRIX')
    zma_save_fs[-1].create(zma_locs)
    zma_save_fs[-1].file.geometry_info.write(inf_obj, zma_locs)
    zma_save_fs[-1].file.geometry_input.write(inp_str, zma_locs)
    zma_save_fs[-1].file.zmatrix.write(zma, zma_locs)

    # Get the tors names
    tors_names = ()
    zma_save_fs[-1].file.torsional_names.write(tors_names, zma_locs)

    # Save the tra and gra for a saddle
    if saddle:
        zma_save_fs[-1].file.transformation.write(tra, zma_locs)
        zma_save_fs[-1].file.reactant_graph.write(rct_gra, zma_locs)

    # Saving the energy to a SP filesystem
    print(" - Saving energy of unique geometry...")
    sp_save_fs = autofile.fs.single_point(cnf_save_path)
    sp_save_fs[-1].create(thy_info[1:4])
    sp_save_fs[-1].file.input.write(inp_str, thy_info[1:4])
    sp_save_fs[-1].file.info.write(inf_obj, thy_info[1:4])
    sp_save_fs[-1].file.energy.write(ene, thy_info[1:4])


def _save_sym_indistinct_conformer(geo, cnf_save_fs,
                                   cnf_tosave_locs, cnf_saved_locs):
    """ Save a structure into the SYM directory of a conformer
    """

    # Set the path to the conf filesys with sym similar
    cnf_save_path = cnf_save_fs[-1].path(cnf_saved_locs)

    # Build the sym file sys
    sym_save_fs = fs.manager(cnf_save_path, 'SYMMETRY')
    sym_save_path = cnf_save_fs[-1].path(cnf_saved_locs)
    print(" - Saving structure in a sym directory at path {}".format(
        sym_save_path))
    sym_save_fs[-1].create(cnf_tosave_locs)
    sym_save_fs[-1].file.geometry.write(geo, cnf_tosave_locs)


def _saved_cnf_info(cnf_save_fs, mod_thy_info):
    """ get the locs, geos and enes for saved conformers
    """

    saved_locs = list(cnf_save_fs[-1].existing())
    saved_geos = [cnf_save_fs[-1].file.geometry.read(locs)
                  for locs in saved_locs]
    saved_enes = []
    for locs in saved_locs:
        path = cnf_save_fs[-1].path(locs)
        sp_save_fs = autofile.fs.single_point(path)
        saved_enes.append(sp_save_fs[-1].file.energy.read(
            mod_thy_info[1:4]))

    return saved_locs, saved_geos, saved_enes
