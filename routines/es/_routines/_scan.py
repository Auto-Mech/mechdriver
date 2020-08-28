""" es_runners for coordinate scans
"""

import automol
import autofile
import elstruct
from routines.es import runner as es_runner
from routines.es._routines import sp
from routines.es._routines._fs import save_struct
from routines.es._routines._fs import _read as read_zma_geo
from lib.structure import tors as torsprep
from lib.structure import instab
from lib import filesys


def run_scan(zma, spc_info, mod_thy_info, thy_save_fs,
             coord_names, coord_grids,
             scn_run_fs, scn_save_fs, scn_typ,
             script_str, overwrite,
             update_guess=True, reverse_sweep=True,
             saddle=False,
             constraint_dct=None, retryfail=True,
             chkstab=False,
             **kwargs):
    """ run constrained optimization scan
    """

    # Build the SCANS/CSCANS filesystems
    if constraint_dct is None:
        scn_save_fs[1].create([coord_names])
        inf_obj = autofile.schema.info_objects.scan_branch(
            dict(zip(coord_names, coord_grids)))
        scn_save_fs[1].file.info.write(inf_obj, [coord_names])
    else:
        scn_save_fs[1].create([constraint_dct])
        inf_obj = autofile.schema.info_objects.scan_branch(
            dict(zip(coord_names, coord_grids)))
        scn_save_fs[1].file.info.write(inf_obj, [constraint_dct])

    # Build the grid of values
    _, grid_vals = torsprep.set_scan_dims(coord_grids)

    _run_scan(
        guess_zma=zma,
        spc_info=spc_info,
        mod_thy_info=mod_thy_info,
        thy_save_fs=thy_save_fs,
        coord_names=coord_names,
        grid_vals=grid_vals,
        scn_run_fs=scn_run_fs,
        scn_save_fs=scn_save_fs,
        scn_typ=scn_typ,
        script_str=script_str,
        overwrite=overwrite,
        retryfail=retryfail,
        update_guess=update_guess,
        saddle=saddle,
        constraint_dct=constraint_dct,
        chkstab=chkstab,
        **kwargs
    )

    if reverse_sweep:
        print('\nDoing a reverse sweep of the HR scan to catch errors...')
        _run_scan(
            guess_zma=zma,
            spc_info=spc_info,
            mod_thy_info=mod_thy_info,
            thy_save_fs=thy_save_fs,
            coord_names=coord_names,
            grid_vals=tuple(reversed(grid_vals)),
            scn_run_fs=scn_run_fs,
            scn_save_fs=scn_save_fs,
            scn_typ=scn_typ,
            script_str=script_str,
            overwrite=overwrite,
            retryfail=retryfail,
            update_guess=update_guess,
            saddle=saddle,
            constraint_dct=constraint_dct,
            chkstab=chkstab,
            **kwargs
        )


def _run_scan(guess_zma, spc_info, mod_thy_info, thy_save_fs,
              coord_names, grid_vals,
              scn_run_fs, scn_save_fs, scn_typ,
              script_str, overwrite,
              errors=(), options_mat=(),
              retryfail=True, update_guess=True,
              saddle=False, constraint_dct=None,
              chkstab=False,
              **kwargs):
    """ new run function
    """

    # Get a connected geometry from the init guess_zma for instability checks
    # conn_geo = automol.zmatrix.geometry(guess_zma)
    conn_zma = guess_zma

    # Set the frozen coordinates (set job at this point?)
    if constraint_dct is not None:
        frozen_coordinates = tuple(coord_names) + tuple(constraint_dct)
    else:
        frozen_coordinates = coord_names

    # Read the energies and Hessians from the filesystem
    for vals in grid_vals:

        # Set the locs for the scan point
        locs = [coord_names, vals]
        if constraint_dct is not None:
            locs = [constraint_dct] + locs

        # Create the filesys
        scn_run_fs[-1].create(locs)
        run_fs = autofile.fs.run(scn_run_fs[-1].path(locs))

        # Build the zma
        zma = automol.zmatrix.set_values(
            guess_zma, dict(zip(coord_names, vals)))

        # Set the job
        job = _set_job(scn_typ)

        # Run an optimization or energy job, as needed.
        geo_exists = scn_save_fs[-1].file.geometry.exists(locs)
        if not geo_exists or overwrite:
            if job == elstruct.Job.OPTIMIZATION:
                es_runner.run_job(
                    job=job,
                    script_str=script_str,
                    run_fs=run_fs,
                    geom=zma,
                    spc_info=spc_info,
                    thy_info=mod_thy_info,
                    overwrite=overwrite,
                    frozen_coordinates=frozen_coordinates,
                    errors=errors,
                    options_mat=options_mat,
                    retryfail=retryfail,
                    saddle=saddle,
                    **kwargs
                )

                # Read the output for the zma and geo
                _, opt_zma, opt_geo = read_zma_geo(run_fs, job)

                if opt_zma is not None and opt_geo is not None:

                    # Check connectivity, save instability files if needed
                    if chkstab:
                        connected = automol.geom.connected(opt_geo)
                    else:
                        connected = True

                    # If connected and update requested: update geom
                    # If disconnected: save instab files and break loop
                    if connected:
                        if update_guess:
                            guess_zma = opt_zma
                    else:
                        _, opt_ret = es_runner.read_job(job=job, run_fs=run_fs)
                        instab.write_instab(
                            conn_zma, opt_zma,
                            thy_save_fs, mod_thy_info[1:4],
                            opt_ret,
                            zma_locs=(0,),
                            save_cnf=False
                        )
                        break
                else:
                    print('No output found in file')

            elif job == elstruct.Job.ENERGY:
                es_runner.run_job(
                    job=job,
                    script_str=script_str,
                    run_fs=run_fs,
                    geom=zma,
                    spc_info=spc_info,
                    thy_info=mod_thy_info,
                    overwrite=overwrite,
                    errors=errors,
                    options_mat=options_mat,
                    retryfail=retryfail,
                    **kwargs
                )

                # Run read_job to print status message
                _, _ = es_runner.read_job(job=job, run_fs=run_fs)

                # Write initial geos in run fs as they are needed later
                run_fs[-1].file.zmatrix.write(zma, [job])
                run_fs[-1].file.geometry.write(
                    automol.zmatrix.geometry(zma), [job])


def save_scan(scn_run_fs, scn_save_fs, scn_typ,
              coo_names, mod_thy_info, in_zma_fs=False):
    """ save the scans that have been run so far
    """

    if not scn_run_fs[1].exists([coo_names]):
        print("No scan to save. Skipping...")
    else:
        locs_lst = []
        for locs in scn_run_fs[-1].existing([coo_names]):

            # Set run filesys
            run_path = scn_run_fs[-1].path(locs)
            run_fs = autofile.fs.run(run_path)
            print("Reading from scan run at {}".format(run_path))

            # Save the structure
            saved = save_struct(
                run_fs, scn_save_fs, locs, _set_job(scn_typ),
                mod_thy_info, in_zma_fs=in_zma_fs)

            # Add to locs lst if the structure is saved
            if saved:
                locs_lst.append(locs)

        # Build the trajectory file
        if locs_lst:
            _hr_traj(coo_names, scn_save_fs, mod_thy_info, locs_lst)


def save_cscan(cscn_run_fs, cscn_save_fs, scn_typ,
               constraint_dct, mod_thy_info, in_zma_fs=True):
    """ save the scans that have been run so far
    """

    if not cscn_run_fs[1].exists([constraint_dct]):
        print("No scan to save. Skipping...")
    else:
        locs_lst = []
        for locs1 in cscn_run_fs[2].existing([constraint_dct]):
            if cscn_run_fs[2].exists(locs1):
                for locs2 in cscn_run_fs[3].existing(locs1):

                    # Set run filesys
                    run_path = cscn_run_fs[-1].path(locs2)
                    run_fs = autofile.fs.run(run_path)
                    print("Reading from scan run at {}".format(run_path))

                    # Save the structure
                    saved = save_struct(
                        run_fs, cscn_save_fs, locs2, _set_job(scn_typ),
                        mod_thy_info, in_zma_fs=in_zma_fs)

                    # Add to locs lst if the structure is saved
                    if saved:
                        locs_lst.append(locs2)
            else:
                print("No scan to save. Skipping...")

        # Build the trajectory file
        if locs_lst:
            _hr_traj(constraint_dct, cscn_save_fs, mod_thy_info, locs_lst)


def _set_job(scn_typ):
    """ Determine if scan is rigid or relaxed and set the appropriate
        electronic structure job.
    """

    assert scn_typ in ('relaxed', 'rigid'), (
        '{} is not relaxed or rigid'.format(scn_typ)
    )

    if scn_typ == 'relaxed':
        job = elstruct.Job.OPTIMIZATION
    else:
        job = elstruct.Job.ENERGY

    return job


def _hr_traj(ini_locs, scn_save_fs, mod_thy_info, locs_lst):
    """ Save a hindered rotor trajectory
    """

    idxs_lst = [locs[-1] for locs in locs_lst]
    enes = []
    for locs in locs_lst:
        path = scn_save_fs[-1].path(locs)
        sp_save_fs = autofile.fs.single_point(path)
        enes.append(sp_save_fs[-1].file.energy.read(mod_thy_info[1:4]))
    geos = [scn_save_fs[-1].file.geometry.read(locs)
            for locs in locs_lst]

    traj = []
    for idxs, ene, geo in zip(idxs_lst, enes, geos):
        comment = (
            'energy: {:>15.10f}, '.format(ene) +
            'grid idxs: {}'.format(idxs)
        )
        traj.append((comment, geo))

    traj_path = scn_save_fs[1].file.trajectory.path([ini_locs])
    print("Updating scan trajectory file at {}".format(traj_path))
    scn_save_fs[1].file.trajectory.write(traj, [ini_locs])


# DERIVED FUNCTION THAT RUNS RUN_SCAN AND SAVE IN TWO DIRECTIONS #
def run_two_way_scan(ts_zma, ts_info, mod_var_scn_thy_info,
                     grid1, grid2, coord_name,
                     thy_save_fs,
                     scn_run_fs, scn_save_fs,
                     opt_script_str, overwrite,
                     update_guess=True,
                     reverse_sweep=True,
                     saddle=False,
                     constraint_dct=None,
                     retryfail=False,
                     **opt_kwargs):
    """ Run a two-part scan that goes into two directions, as for rxn path
    """

    # Setup and run the first part of the scan to shorter distances
    for grid in (grid1, grid2):
        run_scan(
            zma=ts_zma,
            spc_info=ts_info,
            mod_thy_info=mod_var_scn_thy_info,
            thy_save_fs=thy_save_fs,
            coord_names=[coord_name],
            coord_grids=[grid],
            scn_run_fs=scn_run_fs,
            scn_save_fs=scn_save_fs,
            scn_typ='relaxed',
            script_str=opt_script_str,
            overwrite=overwrite,
            update_guess=update_guess,
            reverse_sweep=reverse_sweep,
            saddle=saddle,
            constraint_dct=constraint_dct,
            retryfail=retryfail,
            chkstab=False,
            **opt_kwargs
        )

    print('\nSaving the scans...')
    if constraint_dct is None:
        save_scan(
            scn_run_fs=scn_run_fs,
            scn_save_fs=scn_save_fs,
            scn_typ='relaxed',
            coo_names=[coord_name],
            mod_thy_info=mod_var_scn_thy_info,
            in_zma_fs=True)
    else:
        save_cscan(
            cscn_run_fs=scn_run_fs,
            cscn_save_fs=scn_save_fs,
            scn_typ='relaxed',
            constraint_dct=constraint_dct,
            mod_thy_info=mod_var_scn_thy_info,
            in_zma_fs=True)


# DERIVED FUNCTION FOR MULTIREFERENCE SCAN CALCULATIONS
def multiref_rscan(ts_zma, ts_info,
                   grid1, grid2, coord_name,
                   mod_var_scn_thy_info,
                   vscnlvl_thy_save_fs,
                   scn_run_fs, scn_save_fs,
                   overwrite, update_guess=True,
                   constraint_dct=None,
                   **cas_kwargs):
    """ run constrained optimization scan
    """

    # Set the opt script string and build the opt_kwargs
    [prog, method, _, _] = mod_var_scn_thy_info
    _, opt_script_str, _, opt_kwargs = es_runner.qchem_params(
        prog, method)
    opt_kwargs.update(cas_kwargs)

    # Run the scans
    run_two_way_scan(
        ts_zma, ts_info, mod_var_scn_thy_info,
        grid1, grid2, coord_name,
        vscnlvl_thy_save_fs,
        scn_run_fs, scn_save_fs,
        opt_script_str, overwrite,
        update_guess=update_guess,
        reverse_sweep=False,
        saddle=False,
        constraint_dct=constraint_dct,
        retryfail=False,
        **opt_kwargs
    )


# CALCULATE INFINITE SEPARATION ENERGY #
def molrad_inf_sep_ene(rct_info, rcts_cnf_fs,
                       mod_thy_info, overwrite):
    """ Calculate the inf spe ene for a mol-rad ene
    """
    sp_script_str, _, kwargs, _ = es_runner.qchem_params(
        *mod_thy_info[0:2])
    inf_sep_ene = reac_sep_ene(
        rct_info, rcts_cnf_fs,
        mod_thy_info, overwrite, sp_script_str, **kwargs)

    return inf_sep_ene


def radrad_inf_sep_ene(hs_info, ref_zma,
                       rct_info, rcts_cnf_fs,
                       var_sp2_thy_info,
                       hs_var_sp1_thy_info, hs_var_sp2_thy_info,
                       geo, geo_run_path, geo_save_path,
                       scn_save_fs, inf_locs,
                       overwrite=False,
                       **cas_kwargs):
    """ Obtain the infinite separation energy from the multireference energy
        at a given reference point, the high-spin low-spin splitting at that
        reference point, and the high level energy for the high spin state
        at the reference geometry and for the fragments
        scn = thy for optimizations
        sp1 = low-spin single points
        sp2 = high-spin single points for inf sep

        inf = spc0 + spc1 - hs_sr_e + hs_mr_ene
    """

    # Initialize infinite sep energy
    inf_sep_ene = -1.0e12

    # Calculate the energies for the two cases
    hs_inf = (hs_var_sp2_thy_info, hs_var_sp1_thy_info)
    for idx, thy_info in enumerate(hs_inf):

        if idx == 0:
            print(" - Running high-spin single reference energy ...")
        else:
            print(" - Running high-spin multi reference energy ...")

        # Calculate the single point energy
        script_str, _, kwargs, _ = es_runner.qchem_params(*thy_info[0:2])
        cas_kwargs.update(kwargs)

        sp.run_energy(ref_zma, geo, hs_info, thy_info,
                      scn_save_fs, geo_run_path, geo_save_path, inf_locs,
                      script_str, overwrite, highspin=True, **cas_kwargs)

        # Read the energty from the filesystem
        hs_save_fs, _ = filesys.build.high_spin_from_prefix(
            geo_save_path, thy_info)
        if not hs_save_fs[-1].file.energy.exists(thy_info[1:4]):
            print('ERROR: High-spin energy job failed: ',
                  'energy is needed to evaluate infinite separation energy')
            ene = None
        else:
            print(" - Reading high-spin energy from filesystem...")
            ene = hs_save_fs[-1].file.energy.read(thy_info[1:4])

        if idx == 0:
            hs_sr_ene = ene
        else:
            hs_mr_ene = ene

    # get the single reference energy for each of the reactant configurations
    sp_script_str, _, kwargs, _ = es_runner.qchem_params(
        *var_sp2_thy_info[0:2])
    reac_ene = reac_sep_ene(
        rct_info, rcts_cnf_fs,
        var_sp2_thy_info, overwrite, sp_script_str, **kwargs)

    # Calculate the infinite seperation energy
    all_enes = (reac_ene, hs_sr_ene, hs_mr_ene)
    if all(ene is not None for ene in all_enes):
        inf_sep_ene = reac_ene - hs_sr_ene + hs_mr_ene
    else:
        inf_sep_ene = None

    return inf_sep_ene


def reac_sep_ene(rct_info, rcts_cnf_fs, thy_info,
                 overwrite, sp_script_str, **kwargs):
    """ Calculate the sum of two reactants
    """

    # get the single reference energy for each of the reactant configurations
    spc_enes = []
    for (cnf_run_fs, cnf_save_fs, min_locs), inf in zip(rcts_cnf_fs, rct_info):

        # Set the modified thy info
        mod_thy_info = filesys.inf.modify_orb_restrict(inf, thy_info)

        # Build filesys
        cnf_save_paths = filesys.build.cnf_paths_from_locs(
            cnf_save_fs, min_locs)
        zma_fs, _ = filesys.build.zma_fs_from_prefix(
            cnf_save_paths[0])

        # Read the geometry and set paths
        zma = zma_fs[-1].file.zmatrix.read([0])
        geo = cnf_save_fs[-1].file.geometry.read(min_locs)
        geo_run_path = cnf_run_fs[-1].path(min_locs)
        geo_save_path = cnf_save_fs[-1].path(min_locs)

        # Build the single point filesys objects
        sp_save_fs, _ = filesys.build.sp_from_prefix(
            cnf_save_paths[0], mod_thy_info)

        # Calculate the save single point energy
        sp.run_energy(zma, geo, inf, mod_thy_info,
                      cnf_save_fs, geo_run_path, geo_save_path,
                      min_locs, sp_script_str, overwrite,
                      **kwargs)
        exists = sp_save_fs[-1].file.energy.exists(mod_thy_info[1:4])
        if not exists:
            print('No ene found')
            ene = None
        else:
            print('Reading energy')
            ene = sp_save_fs[-1].file.energy.read(mod_thy_info[1:4])

        # Append ene to list
        spc_enes.append(ene)

    # Analyze the energies in the list
    inf_ene = 0.0
    for ene, inf in zip(spc_enes, rct_info):
        if ene is not None:
            inf_ene += ene
        else:
            print('ERROR: Single reference energy job fails',
                  'for {}: '.format(inf),
                  'Energy needed to evaluate infinite separation energy')
            inf_ene = None

    return inf_ene
