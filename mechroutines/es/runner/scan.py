""" es_runners for coordinate scans
"""

import itertools
import automol
import autofile
import elstruct
from mechroutines.es.runner import scan
from mechroutines.es.runner._run import execute_job
from mechroutines.es._routines import sp
from mechlib.structure import tors as torsprep
from mechlib.structure import instab
from mechlib import filesys
from mechlib.submission import qchem_params
from mechlib.amech_io import printer as ioprinter


def execute_scan(zma, spc_info, mod_thy_info, thy_save_fs,
                 coord_names, coord_grids,
                 scn_run_fs, scn_save_fs, scn_typ,
                 script_str, overwrite,
                 update_guess=True, reverse_sweep=True,
                 saddle=False,
                 constraint_dct=None, retryfail=True,
                 chkstab=False,
                 **kwargs):
    """ Run and save the scan
    """

    run_scan(zma, spc_info, mod_thy_info, thy_save_fs,
             coord_names, coord_grids,
             scn_run_fs, scn_save_fs, scn_typ,
             script_str, overwrite,
             update_guess=update_guess, reverse_sweep=reverse_sweep,
             saddle=saddle,
             constraint_dct=constraint_dct, retryfail=retryfail,
             chkstab=chkstab,
             **kwargs)

    ioprinter.info_message(
        'Saving any newly run HR scans in run filesys...',
        newline=1)
    if constraint_dct is None:
        save_scan(
            scn_run_fs=scn_run_fs,
            scn_save_fs=scn_save_fs,
            scn_typ=scn_typ,
            coo_names=coord_names,
            mod_thy_info=mod_thy_info,
            in_zma_fs=True)
    else:
        save_cscan(
            cscn_run_fs=scn_run_fs,
            cscn_save_fs=scn_save_fs,
            scn_typ=scn_typ,
            constraint_dct=constraint_dct,
            mod_thy_info=mod_thy_info,
            in_zma_fs=True)


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
        print('names', coord_names)
        print('grids', coord_grids)
        inf_obj = autofile.schema.info_objects.scan_branch(
            dict(zip(coord_names, coord_grids)))
        scn_save_fs[1].file.info.write(inf_obj, [coord_names])
    else:
        scn_save_fs[1].create([constraint_dct])
        inf_obj = autofile.schema.info_objects.scan_branch(
            dict(zip(coord_names, coord_grids)))
        scn_save_fs[1].file.info.write(inf_obj, [constraint_dct])

    # Build the grid of values
    mixed_grid_vals = automol.pot.coords(coord_grids)
    print('mixed grid_vals\n', mixed_grid_vals)

    _run_scan(
        guess_zma=zma,
        spc_info=spc_info,
        mod_thy_info=mod_thy_info,
        thy_save_fs=thy_save_fs,
        coord_names=coord_names,
        grid_vals=mixed_grid_vals,
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
            grid_vals=tuple(reversed(mixed_grid_vals)),
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
        print('locs', locs)
        if constraint_dct is not None:
            locs = [constraint_dct] + locs
        print('locs', locs)

        # Create the filesys
        scn_run_fs[-1].create(locs)
        run_fs = autofile.fs.run(scn_run_fs[-1].path(locs))

        # Build the zma
        zma = automol.zmat.set_values_by_name(
            guess_zma, dict(zip(coord_names, vals)),
            angstrom=False, degree=False)

        # Set the job
        job = _set_job(scn_typ)

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
                _, opt_zma, opt_geo = filesys.read_zma_geo(run_fs, job)

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
                        print('WARNING: Structure seems to be unstable...')
                        # _, opt_ret = read_job(job=job, run_fs=run_fs)
                        # instab.write_instab(
                        #     conn_zma, opt_zma,
                        #     thy_save_fs, mod_thy_info[1:4],
                        #     opt_ret,
                        #     zma_locs=(0,),
                        #     save_cnf=False
                        # )
                        # break
                else:
                    print('No output found in file')

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
            saved = filesys.save_struct(
                run_fs, scn_save_fs, locs, _set_job(scn_typ),
                mod_thy_info, in_zma_fs=in_zma_fs)

            # Add to locs lst if the structure is saved
            if saved:
                locs_lst.append(locs)

        # Build the trajectory file
        if locs_lst:
            _write_traj(coo_names, scn_save_fs, mod_thy_info, locs_lst)


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
                    saved = filesys.save_struct(
                        run_fs, cscn_save_fs, locs2, _set_job(scn_typ),
                        mod_thy_info, in_zma_fs=in_zma_fs)

                    # Add to locs lst if the structure is saved
                    if saved:
                        locs_lst.append(locs2)
            else:
                print("No scan to save. Skipping...")

        # Build the trajectory file
        if locs_lst:
            _write_traj(constraint_dct, cscn_save_fs, mod_thy_info, locs_lst)


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
    _, opt_script_str, _, opt_kwargs = qchem_params(
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
