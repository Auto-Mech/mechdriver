""" es_runners for coordinate scans
"""

import automol
import elstruct
import autofile
from routines.es import runner as es_runner
from routines.es._routines._fs import save_struct
from routines.es._routines._fs import check_isomer
from lib.structure import tors as torsprep


def run_scan(zma, spc_info, mod_thy_info,
             coord_names, coord_grids,
             scn_run_fs, scn_save_fs,
             script_str, overwrite,
             update_guess=True, reverse_sweep=True, saddle=False,
             constraint_dct=None, retryfail=True,
             **kwargs):
    """ run constrained optimization scan
    """

    # Check if ZMA matches one in filesys
    check_isomer(zma, scn_save_fs)

    # for now, running only one-dimensional hindered rotor scans
    scn_save_fs[1].create([coord_names])
    inf_obj = autofile.schema.info_objects.scan_branch(
        {coord_names: coord_grids})  # WRONG
    scn_save_fs[1].file.info.write(inf_obj, [coord_names])

    # Build run prefixses?
    _run_scan(
        guess_zma=zma,
        spc_info=spc_info,
        mod_thy_info=mod_thy_info,
        coord_names=coord_names,
        coord_grids=coord_grids,
        scn_run_fs=scn_run_fs,
        scn_save_fs=scn_save_fs,
        script_str=script_str,
        overwrite=overwrite,
        errors=(),
        options_mat=(),
        retryfail=retryfail,
        update_guess=update_guess,
        saddle=saddle,
        constraint_dct=constraint_dct
    )

    if reverse_sweep:
        print('\nDoing a reverse sweep of the HR scan to catch errors...')
        _run_scan(
            script_str=script_str,
            run_prefixes=list(reversed(run_prefixes)),
            scn_save_fs=scn_save_fs,
            guess_zma=zma,
            coo_name=coo_names[0],
            grid_idxs=list(reversed(grid_idxs)),
            grid_vals=list(reversed(grid_vals[0])),
            spc_info=spc_info,
            thy_info=thy_info,
            overwrite=overwrite,
            update_guess=update_guess,
            saddle=saddle,
            constraint_dct=constraint_dct,
            **kwargs
        )


def _run_scan(guess_zma, spc_info, mod_thy_info,
              coord_names, coord_grids,
              scn_run_fs, scn_save_fs,
              script_str, overwrite,
              errors=(), options_mat=(),
              retryfail=True, update_guess=True, saddle=False,
              constraint_dct=None,
              **kwargs):
    """ new run function
    """

    # Set the frozen coordinates (set job at this point?)
    if constraint_dct is not None:
        frozen_coordinates = coord_names + list(constraint_dct)
    else:
        frozen_coordinates = coord_names

    # Set the appropriate job based on the frozen_coordinates
    if set(frozen_coordinates) != set(automol.zmatrix.coordinates(guess_zma)):
        job = elstruct.Job.OPTIMIZATION
    else:
        job = elstruct.Job.ENERGY

    # Read the energies and Hessians from the filesystem
    _, grid_vals = torsprep.set_hr_dims(coord_grids)
    for vals in zip(grid_vals):

        # Set the locs for the scan point
        locs = [coord_names, vals]
        if constraint_dct is not None:
            locs.append(constraint_dct)

        # Create the filesys
        scn_run_fs[-1].create(locs)
        run_fs = autofile.fs.run(scn_run_fs[-1].path(locs))

        # Build the zma
        zma = automol.zmatrix.set_values(guess_zma, dict(coord_names, vals))

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

                ret = es_runner.read_job(job=job, run_fs=run_fs)
                if ret is not None:
                    inf_obj, _, out_str = ret
                    prog = inf_obj.prog
                    opt_zma = elstruct.reader.opt_zmatrix(prog, out_str)
                    if update_guess:
                        guess_zma = opt_zma

            elif job == elstruct.Job.ENERGY:
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

                ret = es_runner.read_job(job=job, run_fs=run_fs)


def save_scan(scn_run_fs, scn_save_fs, coo_names, mod_thy_info, job):
    """ save the scans that have been run so far
    """

    if not scn_run_fs[1].exists([coo_names]):
        print("No scan to save. Skipping...")
    else:
        locs_lst = []
        for locs in scn_run_fs[-1].existing([coo_names]):
            if not isinstance(locs[1][0], float):
                continue

            # Set and print the path
            run_path = scn_run_fs[-1].path(locs)
            print("Reading from scan run at {}".format(run_path))

            # Build run fs and save the structure
            run_fs = autofile.fs.run(run_path)
            saved = save_struct(run_fs, scn_save_fs, locs, job,
                                mod_thy_info, inzma, in_zma_fs=True)

            # Add to locs lst if the structure is saved
            if saved:
                locs_lst.append(locs)

        # Build the trajectory file
        if locs_lst:
            _hr_traj(coo_names, scn_save_fs, locs_lst)


def save_cscan(cscn_run_fs, cscn_save_fs, coo_names, mod_thy_info, job):
    """ save the scans that have been run so far
    """

    if cscn_run_fs[1].exists([coo_names]):

        locs_lst = []
        for locs1 in cscn_run_fs[2].existing([coo_names]):
            if cscn_run_fs[2].exists(locs1):
                for locs2 in cscn_run_fs[3].existing(locs1):

                    # Set and print the path
                    run_path = cscn_run_fs[-1].path(locs2)
                    print("Reading from scan run at {}".format(run_path))

                    # Build run fs and save the structure
                    run_fs = autofile.fs.run(run_path)
                    saved = save_struct(run_fs, cscn_save_fs, locs2, job,
                                        mod_thy_info, inzma, in_zma_fs=True)

                    # Add to locs lst if the structure is saved
                    if saved:
                        locs_lst.append(locs2)

        # Build the trajectory file
        if locs_lst:
            _hr_traj(coo_names, cscn_save_fs, locs_lst)

    else:
        print("No cscan to save. (1) Skipping...")


def _hr_traj(coord_names, scn_save_fs, locs_lst):
    """ Save a hindered rotor trajectory
    """

    idxs_lst = [locs[-1] for locs in locs_lst]
    enes = [scn_save_fs[-1].file.energy.read(locs)
            for locs in locs_lst]
    geos = [scn_save_fs[-1].file.geometry.read(locs)
            for locs in locs_lst]

    traj = []
    for idxs, ene, geo in zip(idxs_lst, enes, geos):
        comment = (
            'energy: {:>15.10f}, '.format(ene) +
            'grid idxs: {}'.format(idxs)
        )
        traj.append((comment, geo))

    traj_path = scn_save_fs[1].file.trajectory.path([coord_names])
    print("Updating scan trajectory file at {}".format(traj_path))
    scn_save_fs[1].file.trajectory.write(traj, [coord_names])
