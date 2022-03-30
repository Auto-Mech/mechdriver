""" irc x
"""

import automol.reac
import autofile
import elstruct
from mechlib.reaction import grid as rxngrid
from mechlib.amech_io import printer as ioprinter
from mechroutines.es import runner as es_runner
from mechroutines.es.runner import qchem_params


# Intrinsic Reaction Coordinates
def execute_irc(zma, ts_info,
                mod_ini_thy_info, ini_method_dct,
                ini_scn_run_fs, ini_scn_save_fs,
                es_keyword_dct,
                directions=(elstruct.Job.IRCF, elstruct.Job.IRCR)):
    """ Run and save the IRC
    """

    coord_name = 'IRC'

    overwrite = es_keyword_dct['overwrite']
    retryfail = es_keyword_dct['retryfail']

    # Set up run filesys
    run_fs = autofile.fs.run(ini_scn_run_fs[1].path([coord_name]))

    # Run and Read the IRC in the forward and reverse direction
    for direction in directions:
        script_str, kwargs = qchem_params(
            ini_method_dct, job=direction,
            geo=automol.zmat.geometry(zma), spc_info=ts_info)
        run_irc(
            zma,
            direction,
            coord_name,
            run_fs,
            ini_scn_save_fs,
            ts_info,
            mod_ini_thy_info,
            overwrite,
            retryfail,
            script_str,
            **kwargs
        )
        success, _ = es_runner.read_job(
            job=direction,
            run_fs=run_fs,
        )
        if success:
            save_irc(
                direction,
                coord_name,
                run_fs,
                ini_scn_save_fs,
                mod_ini_thy_info
            )

    return success


def run_irc(zma, irc_job, coord_name, run_fs, ini_scn_save_fs,
            ts_info, mod_ini_thy_info, overwrite, retryfail,
            opt_script_str, **opt_kwargs):
    """ Run the irc job
    """

    def _irc_ran(ini_scn_save_fs, coord_name, irc_job):
        """ See if coords are available
        """

        coords = ini_scn_save_fs[-1].existing([coord_name])
        if irc_job == elstruct.Job.IRCF:
            ran_coords = [coord[1][0] for coord in coords if coord[1][0] > 0.0]
        else:
            ran_coords = [coord[1][0] for coord in coords if coord[1][0] < 0.0]

        return bool(ran_coords)

    # Maybe check for positive coords
    if not _irc_ran(ini_scn_save_fs, coord_name, irc_job):
        print('No IRC calculation in save filesystem')
        opt_success, _ = es_runner.read_job(
            job=irc_job,
            run_fs=run_fs,
        )
        need_irc = not opt_success
    else:
        print('Found IRC directory at '
              f'{ini_scn_save_fs[1].path([coord_name])}')
        need_irc = False

    if need_irc:
        print('Running IRC calculation...')
        es_runner.run_job(
            job=irc_job,
            script_str=opt_script_str,
            run_fs=run_fs,
            geo=zma,
            spc_info=ts_info,
            thy_info=mod_ini_thy_info,
            overwrite=overwrite,
            retryfail=retryfail,
            **opt_kwargs
        )


def save_irc(irc_job, coord_name,
             run_fs, ini_scn_save_fs, mod_ini_thy_info):
    """ Read IRC output and store data in filesystem
    """

    opt_success, opt_ret = es_runner.read_job(
        job=irc_job,
        run_fs=run_fs,
    )
    locs_lst = []
    if opt_success is not None:

        # Read the IRC output file
        inf_obj, inp_str, out_str = opt_ret
        prog = inf_obj.prog
        geos, gras, hessians = elstruct.reader.irc_points(prog, out_str)
        coord_vals, enes = elstruct.reader.irc_path(prog, out_str)

        # Write the data for each geom along IRC to the filesystem
        save_path = ini_scn_save_fs[1].path([coord_name])
        print(" - Saving...")
        print(f" - Save path: {save_path}")
        locs_lst = []
        for idx, val in enumerate(coord_vals):

            # Set locs idx; for reverse, ignore SadPt and flip idx to negative
            locs_idx = idx
            if irc_job == elstruct.Job.IRCR:
                if locs_idx == 0:
                    continue

            # Scale the coordinates so rounding to .2f number is non-zero
            locs = [coord_name, [val*100.0]]
            locs_lst.append(locs)

            # Save files
            ini_scn_save_fs[-1].create(locs)
            ini_scn_save_fs[-1].file.energy.write(enes[idx], locs)
            ini_scn_save_fs[-1].file.geometry.write(geos[idx], locs)
            ini_scn_save_fs[-1].file.geometry_input.write(inp_str, locs)
            ini_scn_save_fs[-1].file.geometry_info.write(inf_obj, locs)
            if gras:
                ini_scn_save_fs[-1].file.gradient.write(gras[idx], locs)
                ini_scn_save_fs[-1].file.gradient_info.write(inf_obj, locs)
            if hessians:
                ini_scn_save_fs[-1].file.hessian.write(hessians[idx], locs)
                ini_scn_save_fs[-1].file.hessian_info.write(inf_obj, locs)

            scn_save_path = ini_scn_save_fs[-1].path(locs)
            sp_save_fs = autofile.fs.single_point(scn_save_path)
            sp_save_fs[-1].create(mod_ini_thy_info[1:4])
            sp_save_fs[-1].file.input.write(inp_str, mod_ini_thy_info[1:4])
            sp_save_fs[-1].file.info.write(inf_obj, mod_ini_thy_info[1:4])
            sp_save_fs[-1].file.energy.write(enes[idx], mod_ini_thy_info[1:4])

    update_traj_file(coord_name, ini_scn_save_fs, mod_ini_thy_info)

    return locs_lst


def update_traj_file(coord_name, ini_scn_save_fs, mod_ini_thy_info):
    """ Update the full IRC trajectory file based on what is in SAVE
        filesystem
    """
    saved_locs = ini_scn_save_fs[-1].existing()
    if saved_locs:
        es_runner.scan.write_traj(
            coord_name, ini_scn_save_fs, mod_ini_thy_info, sorted(saved_locs)
        )


def launch_point_zmatrices(ts_dct, mod_thy_info,
                           scn_alg, scn_fs, cnf_fs, cnf_locs):
    """ Determine the point to launch an IRC from

        Try to find saddle point at inplvl
        Then search for the max


        'auto': use sadpt, then max series
        'sadpt': sadpt
        'series': max series
    """

    if 'sadpt' in scn_alg:
        _, cnf_save_fs = cnf_fs
        zma_locs = (ts_dct['zma_idx'],)
        zma_fs = autofile.fs.zmatrix(cnf_save_fs[-1].path(cnf_locs))
        if zma_fs[-1].file.zmatrix.exists(zma_locs):
            geo_path = zma_fs[-1].file.zmatrix.path(zma_locs)
            ioprinter.info_message(
                ' - Z-Matrix found.')
            ioprinter.info_message(
                f' - Reading Z-Matrix from path {geo_path}')
            irc_zmas = (zma_fs[-1].file.zmatrix.read(zma_locs),)
    elif 'max' in scn_alg:
        _, scn_save_fs = scn_fs
        zma, zrxn = ts_dct['zma'], ts_dct['zrxn']

        scan_inf = automol.reac.build_scan_info(zrxn, zma)
        coord_names, constraint_dct, coord_grids, _ = scan_inf

        irc_zmas = rxngrid.grid_maximum_zmatrices(
            zrxn.class_, zma, coord_grids, coord_names, scn_save_fs,
            mod_thy_info, constraint_dct, series='full-n1')

    print('irc zmas', irc_zmas)

    import sys
    sys.exit()

    return irc_zmas
