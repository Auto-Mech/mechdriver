""" es_runners for coordinate scans
"""

import automol
import autofile
import elstruct
from mechlib import filesys
from mechlib.amech_io import printer as ioprinter
from mechroutines.es.runner._run import execute_job
from mechroutines.es.runner._run import read_job


def execute_scan(zma, spc_info, mod_thy_info,
                 coord_names, coord_grids,
                 scn_run_fs, scn_save_fs, scn_typ,
                 script_str, overwrite,
                 update_guess=True, reverse_sweep=False,
                 saddle=False,
                 constraint_dct=None, retryfail=True,
                 **kwargs):
    """ Run and save the scan
    """

    # Need a resave option
    _fin = _scan_finished(
        coord_names, coord_grids, scn_save_fs, constraint_dct=constraint_dct)

    if not _fin:
        run_scan(
            zma, spc_info, mod_thy_info,
            coord_names, coord_grids,
            scn_run_fs, scn_save_fs, scn_typ,
            script_str, overwrite,
            update_guess=update_guess, reverse_sweep=reverse_sweep,
            saddle=saddle,
            constraint_dct=constraint_dct, retryfail=retryfail,
            **kwargs)

        save_scan(
            scn_run_fs=scn_run_fs,
            scn_save_fs=scn_save_fs,
            scn_typ=scn_typ,
            coord_names=coord_names,
            constraint_dct=constraint_dct,
            mod_thy_info=mod_thy_info)


def run_scan(zma, spc_info, mod_thy_info,
             coord_names, coord_grids,
             scn_run_fs, scn_save_fs, scn_typ,
             script_str, overwrite,
             update_guess=True, reverse_sweep=True,
             saddle=False,
             constraint_dct=None, retryfail=True,
             **kwargs):
    """ run constrained optimization scan
    """

    # Build the SCANS/CSCANS filesystems
    if constraint_dct is None:
        coord_locs = coord_names
    else:
        coord_locs = constraint_dct

    scn_save_fs[1].create([coord_locs])
    inf_obj = autofile.schema.info_objects.scan_branch(
        dict(zip(coord_locs, coord_grids)))
    scn_save_fs[1].file.info.write(inf_obj, [coord_locs])

    # Build the grid of values
    mixed_grid_vals = automol.pot.coords(coord_grids)
    if not reverse_sweep:
        all_grid_vals = [mixed_grid_vals]
    else:
        all_grid_vals = [mixed_grid_vals, tuple(reversed(mixed_grid_vals))]

    for idx, grid_vals in enumerate(all_grid_vals):
        if idx == 1:
            print('\nDoing a reverse sweep of the scan to catch errors...')
        _run_scan(
            guess_zma=zma,
            spc_info=spc_info,
            mod_thy_info=mod_thy_info,
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
            **kwargs
        )


def _run_scan(guess_zma, spc_info, mod_thy_info,
              coord_names, grid_vals,
              scn_run_fs, scn_save_fs, scn_typ,
              script_str, overwrite,
              errors=(), options_mat=(),
              retryfail=True, update_guess=True,
              saddle=False, constraint_dct=None,
              **kwargs):
    """ new run function
    """

    # Get a connected geometry from the init guess_zma for instability checks
    # conn_geo = automol.zmatrix.geometry(guess_zma)
    # conn_zma = guess_zma

    # Set the frozen coordinates (set job at this point?)
    if constraint_dct is not None:
        frozen_coordinates = tuple(coord_names) + tuple(constraint_dct)
    else:
        frozen_coordinates = coord_names

    # Set the job
    job = _set_job(scn_typ)

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
        zma = automol.zmat.set_values_by_name(
            guess_zma, dict(zip(coord_names, vals)),
            angstrom=False, degree=False)

        # Run an optimization or energy job, as needed.
        geo_exists = scn_save_fs[-1].file.geometry.exists(locs)
        if not geo_exists or overwrite:
            if job == elstruct.Job.OPTIMIZATION:
                success, ret = execute_job(
                    job=job,
                    script_str=script_str,
                    run_fs=run_fs,
                    geo=zma,
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
                if success:
                    opt_zma = filesys.save.read_job_zma(ret, init_zma=zma)
                    if update_guess:
                        guess_zma = opt_zma

            elif job == elstruct.Job.ENERGY:
                _, _ = execute_job(
                    job=job,
                    script_str=script_str,
                    run_fs=run_fs,
                    geo=zma,
                    spc_info=spc_info,
                    thy_info=mod_thy_info,
                    overwrite=overwrite,
                    errors=errors,
                    options_mat=options_mat,
                    retryfail=retryfail,
                    **kwargs
                )

                # Write initial geos in run fs as they are needed later
                run_fs[-1].file.zmatrix.write(zma, [job])
                run_fs[-1].file.geometry.write(
                    automol.zmat.geometry(zma), [job])


def save_scan(scn_run_fs, scn_save_fs, scn_typ,
              coord_names, constraint_dct,
              mod_thy_info):
    """ save the scan
    """

    ioprinter.info_message(
        'Saving any newly run HR scans in run filesys...', newline=1)

    # Set job
    job = _set_job(scn_typ)

    # Set locs for scan
    coord_locs, save_locs = scan_locs(
        scn_run_fs, coord_names, constraint_dct=constraint_dct)

    if not scn_run_fs[1].exists([coord_locs]):
        print("No scan to save. Skipping...")
    else:
        locs_lst = []
        for locs in save_locs:

            # Set run filesys
            run_path = scn_run_fs[-1].path(locs)
            run_fs = autofile.fs.run(run_path)
            print("Reading from scan run at {}".format(run_path))

            # Save the structure
            success, ret = read_job(job, run_fs)
            if success:
                # Need to get the init zma structure in here
                # could write init zma to run filesys; wont work retro
                # get init zma readers?
                init_geo = run_fs[-1].file.geometry.read([job])
                init_zma = run_fs[-1].file.zmatrix.read([job])
                filesys.save.scan_point_structure(
                    ret, scn_save_fs, locs, mod_thy_info[1:], job,
                    init_zma=init_zma, init_geo=init_geo)
                locs_lst.append(locs)

        # Build the trajectory file
        if locs_lst:
            _write_traj(coord_locs, scn_save_fs, mod_thy_info, locs_lst)


def scan_locs(scn_fs, coord_names, constraint_dct=None):
    """ Get all of the scan locs for a given tors name for either a
        SCAN or CSCAN filesys
    """

    if constraint_dct is None:
        coord_locs = coord_names
        scn_locs = scn_fs[-1].existing([coord_locs])
    else:
        coord_locs = constraint_dct
        scn_locs = ()
        # scn_locs = scn_fs[3].existing()
        scn_locs = ()
        for locs1 in scn_fs[2].existing([coord_locs]):
            if scn_fs[2].exists(locs1):
                for locs2 in scn_fs[3].existing(locs1):
                    scn_locs += (locs2,)

    return scn_locs


def _scan_finished(coord_names, coord_grids, scn_save_fs, constraint_dct=None):
    """ See if the scan needs to be run

        maybe return the grid that is not finished?
    """

    run_finished = True

    grid_vals = automol.pot.coords(coord_grids)
    for vals in grid_vals:

        # Set the locs for the scan point
        locs = [coord_names, vals]
        if constraint_dct is not None:
            locs = [constraint_dct] + locs

        # Check if ZMA (other info?) exists
        if not scn_save_fs[-1].file.zmatrix.exists(locs):
            run_finished = False
            break

    if run_finished:
        ioprinter.message('Scan saved previously at {}'.format(
            scn_save_fs[0].path()))
    else:
        ioprinter.message('Need to run scans')

    return run_finished


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


def _write_traj(ini_locs, scn_save_fs, mod_thy_info, locs_lst):
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
        idx_lst = ['{0:.2f}'.format(idx) for idx in idxs]
        idx_str = ','.join(idx_lst)
        comment = (
            'energy: {:>15.10f}, '.format(ene) +
            'grid idxs: {}'.format(idx_str)
        )
        traj.append((geo, comment))

    traj_path = scn_save_fs[1].file.trajectory.path([ini_locs])
    print("Updating scan trajectory file at {}".format(traj_path))
    scn_save_fs[1].file.trajectory.write(traj, [ini_locs])


# DERIVED FUNCTION THAT RUNS RUN_SCAN AND SAVE IN TWO DIRECTIONS #
def run_two_way_scan(ts_zma, ts_info, mod_var_scn_thy_info,
                     grid1, grid2, coord_name,
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

    for grid in (grid1, grid2):
        execute_scan(
            zma=ts_zma,
            spc_info=ts_info,
            mod_thy_info=mod_var_scn_thy_info,
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
            **opt_kwargs
        )
