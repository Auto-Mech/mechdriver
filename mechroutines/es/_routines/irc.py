"""
IRC calcs
"""

import elstruct
import autofile
from routines.es import runner as es_runner
from lib.submission import qchem_params


def scan(zma, ts_info, coord_name,
         mod_ini_thy_info, ini_method_dct,
         ini_scn_save_fs, geo_run_path,
         overwrite):
    """ Run the IRC
    """

    # Set up run filesys
    run_fs = autofile.fs.run(geo_run_path)

    # Run and Read the IRC in the forward and reverse direction
    for irc_job in (elstruct.Job.IRCF, elstruct.Job.IRCR):

        # Set up the script
        script_str, kwargs = qchem_params(
            ini_method_dct, job=irc_job)

        run_irc(
            zma,
            irc_job,
            coord_name,
            run_fs,
            ini_scn_save_fs,
            ts_info,
            mod_ini_thy_info,
            overwrite,
            script_str,
            **kwargs
        )
        save_irc(
            irc_job,
            coord_name,
            run_fs,
            ini_scn_save_fs,
            mod_ini_thy_info
        )


def run_irc(zma, irc_job, coord_name, run_fs, ini_scn_save_fs,
            ts_info, mod_ini_thy_info, overwrite,
            opt_script_str, **opt_kwargs):
    """ Run the irc job
    """

    # Maybe check for positive coords
    if not _irc_ran(ini_scn_save_fs, coord_name, irc_job):
        print('No IRC calculation in save filesystem')
        opt_success, _ = es_runner.read_job(
            job=irc_job,
            run_fs=run_fs,
        )
        need_irc = bool(not opt_success)
    else:
        print('Found IRC directory at {}'.format(
            ini_scn_save_fs[1].path([coord_name])))
        need_irc = False

    if need_irc:
        print('Running IRC calculation')
        es_runner.run_job(
            job=irc_job,
            script_str=opt_script_str,
            run_fs=run_fs,
            geom=zma,
            spc_info=ts_info,
            thy_info=mod_ini_thy_info,
            overwrite=overwrite,
            **opt_kwargs,
        )
    else:
        print('Skipping IRC calculation')


def save_irc(irc_job, coord_name,
             run_fs, ini_scn_save_fs, mod_ini_thy_info):
    """ Read IRC output and store data in filesystem
    """

    opt_success, opt_ret = es_runner.read_job(
        job=irc_job,
        run_fs=run_fs,
    )
    if opt_success is not None:

        # Read the IRC output file
        inf_obj, inp_str, out_str = opt_ret
        prog = inf_obj.prog
        geos, gras, hessians = elstruct.reader.irc_points(prog, out_str)
        coord_vals, enes = elstruct.reader.irc_path(prog, out_str)

        # Write the data for each geom along IRC to the filesystem
        save_path = ini_scn_save_fs[1].path([coord_name])
        print(" - Saving...")
        print(" - Save path: {}".format(save_path))
        locs_lst = []
        for idx, val in enumerate(coord_vals):

            # Set locs idx; for reverse, ignore SadPt and flip idx to negative
            locs_idx = idx
            if irc_job == elstruct.Job.IRCR:
                if locs_idx == 0:
                    continue
                # val *= -1  coord vals negative from elstruct fxn for g09

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
        
        # Build the trajectory file
        if locs_lst:
            es_runner.write_traj(
                coord_name, ini_scn_save_fs, mod_ini_thy_info, locs_lst
            )


def _irc_ran(ini_scn_save_fs, coord_name, irc_job):
    """ See if coords are available
    """

    coords = ini_scn_save_fs[-1].existing([coord_name])
    if irc_job == elstruct.Job.IRCF:
        ran_coords = [coord[1][0] for coord in coords if coord[1][0] > 0.0]
    else:
        ran_coords = [coord[1][0] for coord in coords if coord[1][0] < 0.0]

    return bool(ran_coords)
