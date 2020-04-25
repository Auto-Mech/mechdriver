"""
IRC calcs
"""

import automol
import elstruct
import autofile
from lib.filesystem import orb as fsorb
from lib.runner import par as runpar
from lib.runner import driver as rundriver


def scan(zma, ts_info, mod_thy_info, coo_name, irc_idxs,
         scn_save_fs, scn_run_fs, geo_run_path,
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
            coo_name,
            irc_idxs,
            run_fs,
            scn_save_fs,
            ts_info,
            mod_thy_info,
            overwrite,
            opt_script_str,
            **opt_kwargs
        )
        save_irc(
            irc_job,
            run_fs,
            scn_run_fs,
            scn_save_fs,
            coo_name,
            irc_idxs
        )


def run_irc(zma, irc_job, coo_name, irc_idxs, run_fs, scn_save_fs,
            ts_info, scn_thy_info, overwrite,
            opt_script_str, **opt_kwargs):
    """ Run the irc job
    """

    irc_ran = True
    for idx in irc_idxs: 
        locs = [[coo_name], [idx]]
        if not scn_save_fs[-1].file.geometry.exists(locs):
            irc_ran = False
    
    if not irc_ran:
        print("No IRC ran...")
        # set irc options here for now
        opt_kwargs['job_options'] = ['calcall', 'stepsize=3', 'maxpoints=4']

        # Run the calculations
        rundriver.run_job(
            job=irc_job,
            script_str=opt_script_str,
            run_fs=run_fs,
            geom=zma,
            spc_info=ts_info,
            thy_level=scn_thy_info,
            overwrite=overwrite,
            **opt_kwargs,
            )
    else:
        print("IRC ran with all points requested")


def save_irc(irc_job, run_fs, scn_run_fs, scn_save_fs,
             coo_name, irc_idxs):
    """ Read IRC output and store data in filesystem
    """

    opt_ret = rundriver.read_job(
        job=irc_job,
        run_fs=run_fs,
    )
    if opt_ret is not None:
        inf_obj, inp_str, out_str = opt_ret
        prog = inf_obj.prog
        geos, gras, hessians = elstruct.reader.irc_points(prog, out_str)
        enes = elstruct.reader.irc_energies(prog, out_str)
        coords = elstruct.reader.irc_coordinates(prog, out_str)

        # Write the IRC inf file and input file string
        # scn_save_fs[1].file.irc_info.write(inf_obj, [coo_name])
        # scn_save_fs[1].file.irc_input.write(inp_str, [coo_name])

        # Write the IRC coords and enes to a yaml file
        # irc_inf_obj = autofile.system.info.irc(idxs=irc_idxs, coords=coords)
        # scn_save_fs[1].file.info.write(irc_inf_obj, [coo_name])

        # Write the data for each geom along IRC to the filesystem
        save_path = scn_save_fs[1].path([[coo_name]])
        print(" - Saving...")
        print(" - Save path: {}".format(save_path))
        for idx, _ in enumerate(geos):

            # Set locs idx; for reverse, ignore SadPt and flip idx to negative
            locs_idx = idx
            if irc_job == elstruct.Job.IRCR:
                if locs_idx == 0:
                    continue
                locs_idx *= -1
            locs = [[coo_name], [locs_idx]]
            scn_save_fs[-1].create(locs)
            scn_save_fs[-1].file.energy.write(enes[idx], locs)
            scn_save_fs[-1].file.energy.write(enes[idx], locs)
            scn_save_fs[-1].file.geometry.write(geos[idx], locs)
            scn_save_fs[-1].file.gradient.write(gras[idx], locs)
            scn_save_fs[-1].file.hessian.write(hessians[idx], locs)
