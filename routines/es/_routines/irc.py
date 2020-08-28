"""
IRC calcs
"""

import elstruct
import autofile
from routines.es import runner as es_runner


def scan(zma, ts_info, mod_ini_thy_info, coord_name,
         ini_scn_save_fs, geo_run_path,
         overwrite, opt_script_str, **opt_kwargs):
    """ Run the IRC
    """

    # Set up run filesys
    run_fs = autofile.fs.run(geo_run_path)

    # Run and Read the IRC in the forward and reverse direction
    for irc_job in (elstruct.Job.IRCF, elstruct.Job.IRCR):
        run_irc(
            zma,
            irc_job,
            coord_name,
            run_fs,
            ini_scn_save_fs,
            ts_info,
            mod_ini_thy_info,
            overwrite,
            opt_script_str,
            **opt_kwargs
        )
        save_irc(
            irc_job,
            run_fs,
            ini_scn_save_fs,
            coord_name,
        )


def run_irc(zma, irc_job, coord_name, run_fs, ini_scn_save_fs,
            ts_info, mod_ini_thy_info, overwrite,
            opt_script_str, **opt_kwargs):
    """ Run the irc job
    """

    # Maybe check for positive coords
    if not _irc_ran(ini_scn_save_fs, coord_name, irc_job):
        print('No IRC calculation in save filesystem')
        print('Running IRC calculation')
        need_irc = True
    else:
        print('Found IRC directory at {}'.format(
            ini_scn_save_fs[1].path([coord_name])))
        print('Skipping IRC calculation')
        need_irc = False

    if need_irc:
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


def save_irc(irc_job, run_fs, ini_scn_save_fs, coord_name):
    """ Read IRC output and store data in filesystem
    """

    opt_success, opt_ret = es_runner.read_job(
        job=irc_job,
        run_fs=run_fs,
    )
    if opt_success is not None:
        inf_obj, _, out_str = opt_ret
        prog = inf_obj.prog
        geos, gras, hessians = elstruct.reader.irc_points(prog, out_str)
        enes = elstruct.reader.irc_energies(prog, out_str)
        coord_vals = elstruct.reader.irc_coordinates(prog, out_str)

        # Write the IRC inf file and input file string
        # scn_save_fs[1].file.irc_info.write(inf_obj, [coord_name])
        # scn_save_fs[1].file.irc_input.write(inp_str, [coord_name])

        # Write the IRC coords and enes to a yaml file
        # irc_inf_obj = autofile.schema.info_objects.irc(
        #    idxs=irc_idxs, coords=coords)
        # scn_save_fs[1].file.info.write(irc_inf_obj, [coord_name])

        # Write the data for each geom along IRC to the filesystem
        save_path = ini_scn_save_fs[1].path([coord_name])
        print(" - Saving...")
        print(" - Save path: {}".format(save_path))
        for idx, val in enumerate(coord_vals):

            # Set locs idx; for reverse, ignore SadPt and flip idx to negative
            locs_idx = idx
            if irc_job == elstruct.Job.IRCR:
                if locs_idx == 0:
                    continue
                val *= -1

            # Scale the coordinates so that there are all zeros in the .2f number
            locs = [coord_name, [val*100.0]]

            # Save files
            ini_scn_save_fs[-1].create(locs)
            ini_scn_save_fs[-1].file.energy.write(enes[idx], locs)
            ini_scn_save_fs[-1].file.geometry.write(geos[idx], locs)
            ini_scn_save_fs[-1].file.gradient.write(gras[idx], locs)
            ini_scn_save_fs[-1].file.hessian.write(hessians[idx], locs)


def _irc_ran(ini_scn_save_fs, coord_name, irc_job):
    """ See if coords are available
    """

    print('irc coord name', coord_name)
    coords = ini_scn_save_fs[-1].existing([coord_name])
    if irc_job == elstruct.Job.IRCF:
        ran_coords = [coord[1][0] for coord in coords if coord[1][0] > 0.0]
    else:
        ran_coords = [coord[1][0] for coord in coords if coord[1][0] < 0.0]

    return bool(ran_coords)
