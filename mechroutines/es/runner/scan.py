""" Library to perform sequences of electronic structure calculations
    along a molecular coordinate and save the resulting information to
    SCAN or CSAN layers of the SAVE filesystem.
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
    """ Wrapper to perform both the run_scan and save_scan functions for
        a series of electronic structure calculations along some
        internal coordinates defined in a Z-Matrix.

        Prior to launching any electronic structure jobs, function will
        assess if Z-matrices exist at all points of the input grids in
        the SCAN/CSAN filesystem layer of the SAVE filesystem.

        :param zma: Z-Matrix to start scan from
        :type zma: automol.zmtat object
        :param spc_info:
        :type spc_info:
        :param mod_thy_info:
        :type mod_thy_info:
        :param coord_names: names of the scan coordinates
        :type coord_names: tuple(tuple(str))
        :param coord_grids: values of all the scan coordinates
        :type coord_grids: tuple(tuple(float))
        :param scn_run_fs: SCAN/CSCAN object with RUN filesys prefix
        :type scn_run_fs: autofile.fs.scan or autofile.fs.cscan object
        :param scn_save_fs: SCAN/CSCAN object with SAVE filesys prefix
        :type scn_save_fs: autofile.fs.scan or autofile.fs.cscan object
        :param scn_typ: label for scan type ('relaxed' or 'rigid')
        :type scn_typ: str
        :param overwrite: overwrite existing input file with new one and rerun
        :type overwrite: bool
        :param update_guess: start optimization at point of scan using 
            optimized structure from the previous point of the scan
        :type update_guess: bool
        :param reverse_sweep: attempt to run scan in both directions
            of the coordinate defined by the grid
        :type reverse_sweep: bool
        :param saddle: perform saddle-point optimization
        :type saddle: bool
        :param constraint_dct: values of coordinates to constrain during scan
        :type constraint_dct: dict[str: float]
        :param retryfail: re-run the job if failed job found in RUN filesys
        :type retryfail: bool
        :param kwargs: addiitonal options for electronic structure job
        :type kwargs: dict
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

        :param zma: Z-Matrix to start scan from
        :type zma: automol.zmtat object
        :param spc_info:
        :type spc_info:
        :param mod_thy_info:
        :type mod_thy_info:
        :param coord_names: names of the scan coordinates
        :type coord_names: tuple(tuple(str))
        :param coord_grids: values of all the scan coordinates
        :type coord_grids: tuple(tuple(float))
        :param scn_run_fs: SCAN/CSCAN object with run filesys prefix
        :type scn_run_fs: autofile.fs.scan or autofile.fs.cscan object
        :param scn_save_fs: SCAN/CSCAN object with save filesys prefix
        :type scn_save_fs: autofile.fs.scan or autofile.fs.cscan object
        :param scn_typ: label for scan type ('relaxed' or 'rigid')
        :type scn_typ: str
        :param overwrite: overwrite existing input file with new one and rerun
        :type overwrite: bool
        :param update_guess: start optimization at point of scan using 
            optimized structure from the previous point of the scan
        :type update_guess: bool
        :param reverse_sweep: attempt to run scan in both directions
            of the coordinate defined by the grid
        :type reverse_sweep: bool
        :param saddle: perform saddle-point optimization
        :type saddle: bool
        :param constraint_dct: values of coordinates to constrain during scan
        :type constraint_dct: dict[str: float]
        :param retryfail: re-run the job if failed job found in RUN filesys
        :type retryfail: bool
        :param kwargs: addiitonal options for electronic structure job
        :type kwargs: dict
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

        :param guess_zma: Z-Matrix to start scan from
        :type guess_zma: automol.zmtat object
        :param spc_info:
        :type spc_info:
        :param mod_thy_info:
        :type mod_thy_info:
        :param coord_names: names of the scan coordinates
        :type coord_names: tuple(tuple(str))
        :param grid_vals: values of all the scan coordinates
        :type grid_vals: same as coord_grids???
        :param scn_run_fs: SCAN/CSCAN object with run filesys prefix
        :type scn_run_fs: autofile.fs.scan or autofile.fs.cscan object
        :param scn_save_fs: SCAN/CSCAN object with save filesys prefix
        :type scn_save_fs: autofile.fs.scan or autofile.fs.cscan object
        :param scn_typ: label for scan type ('relaxed' or 'rigid')
        :type scn_typ: str
        :param script_str:
        :type script_str: str
        :param overwrite: overwrite existing input file with new one and rerun
        :type overwrite: bool
        :param errors: list of error message types to search output for
        :type errors: tuple(str)
        :param options_mat: varopis options to run job with
        :type options_mat: tuple(dict[str: str])
        :param update_guess: start optimization at point of scan using 
            optimized structure from the previous point of the scan
        :type update_guess: bool
        :param reverse_sweep: attempt to run scan in both directions
            of the coordinate defined by the grid
        :type reverse_sweep: bool
        :param saddle: perform saddle-point optimization
        :type saddle: bool
        :param constraint_dct: values of coordinates to constrain during scan
        :type constraint_dct: dict[str: float]
        :param retryfail: re-run the job if failed job found in RUN filesys
        :type retryfail: bool
        :param kwargs: addiitonal options for electronic structure job
        :type kwargs: dict
    """

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


def save_scan(scn_run_fs, scn_save_fs, scn_typ,
              coord_names, constraint_dct,
              mod_thy_info):
    """ Search for output of electronic structure calculation along the scan
        coordinate that exist in the SCAN/CSCAN layer of the run filesys. Then
        parse out the required information and save it into formatted files
        in the SCAN/CSAN layer of the save filesys.

        :param scn_run_fs: SCAN/CSCAN object with run filesys prefix
        :type scn_run_fs: autofile.fs.scan or autofile.fs.cscan object
        :param scn_save_fs: SCAN/CSCAN object with save filesys prefix
        :type scn_save_fs: autofile.fs.scan or autofile.fs.cscan object
        :param coord_names: names of the scan coordinates
        :type coord_names: tuple(tuple(str))
        :param constraint_dct: values of coordinates to constrain during scan
        :type constraint_dct: dict[str: float]
    """

    ioprinter.info_message(
        'Saving any newly run HR scans in run filesys...', newline=1)

    # Set job
    job = _set_job(scn_typ)

    # Set locs for scan
    coord_locs, save_locs = _scan_locs(
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
                if run_fs[-1].file.zmatrix.exists([job]):
                    init_zma = run_fs[-1].file.zmatrix.read([job])
                filesys.save.scan_point_structure(
                    ret, scn_save_fs, locs, mod_thy_info[1:], job,
                    init_zma=init_zma, init_geo=None)
                locs_lst.append(locs)

        # Build the trajectory file
        if locs_lst:
            _write_traj(coord_locs, scn_save_fs, mod_thy_info, locs_lst)


def _scan_locs(scn_save_fs, coord_names, constraint_dct=None):
    """  Determine the locs for all of the directories that currently
         exist in the SCAN/CSAN layer of the SAVE filesystem.

        :param coord_names: names of the scan coordinates
        :type coord_names: tuple(tuple(str))
        :param scn_save_fs: SCAN/CSCAN object with save filesys prefix
        :type scn_save_fs: autofile.fs.scan or autofile.fs.cscan object
        :param constraint_dct: values of coordinates to constrain during scan
        :type constraint_dct: dict[str: float]
    """

    if constraint_dct is None:
        coord_locs = coord_names
        scn_locs = scn_save_fs[-1].existing([coord_locs])
    else:
        coord_locs = constraint_dct
        scn_locs = ()
        # scn_locs = scan_save_fs[3].existing()
        scn_locs = ()
        for locs1 in scn_save_fs[2].existing([coord_locs]):
            if scn_save_fs[2].exists(locs1):
                for locs2 in scn_save_fs[3].existing(locs1):
                    scn_locs += (locs2,)

    return coord_locs, scn_locs


def _scan_finished(coord_names, coord_grids, scn_save_fs, constraint_dct=None):
    """ Assesses if the scan calculations requested by the user have been
        completed by assessing if Z-Matrices exist in the filesystem for
        all grid values of the scan coordinates.

        :param coord_names: names of the scan coordinates
        :type coord_names: tuple(tuple(str))
        :param coord_grids: values of all the scan coordinates
        :type coord_grids: tuple(tuple(float))
        :param scn_save_fs: SCAN/CSCAN object with save filesys prefix
        :type scn_save_fs: autofile.fs.scan or autofile.fs.cscan object
        :param constraint_dct: values of coordinates to constrain during scan
        :type constraint_dct: dict[str: float]
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
    """ Set the appropriate job label defined in the elstruct
        package using the input scan type.

        :param scn_typ: label for scan type ('relaxed' or 'rigid')
        :type scn_typ: str
        :rtype: str
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
    """ Read the geometries and energies of all optimized geometries
        in the SCAN/CSCAN filesystem and collate them into am .xyz
        trajectory file, which is thereafter written into a file inthe
        filesystem.

        :param ini_locs: locs for high level SCAN/CSCAN to save traj file
        :type ini_locs: dict[]
        :param scn_save_fs: SCAN/CSCAN object with save filesys prefix
        :type scn_save_fs: autofile.fs.scan or autofile.fs.cscan object
        :param mod_thy_info: thy info object with
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
        Wrapper to the execute_scan to run in two directions
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
