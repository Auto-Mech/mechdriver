"""
IRC calcs
"""

import automol
import elstruct
import autofile
from lib.filesystem import orb as fsorb
from lib.runner import par as runpar
from lib.runner import driver as rundriver


def irc_scan(zma, ts_info, mod_thy_info, coo_name, irc_idxs,
             scn_save_fs, scn_run_fs, run_fs,
             overwrite, opt_script_str, **opt_kwargs):
    """ Run the IRC
    """

    # Run and Read the IRC in the forward and reverse direction
    for irc_job in (elstruct.Job.IRCF, elstruct.Job.IRCR):
        run_irc(
            zma,
            irc_job,
            run_fs,
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


def run_irc(zma, irc_job, run_fs, ts_info, scn_thy_info, overwrite,
            opt_script_str, **opt_kwargs):
    """ Run the irc job
    """

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


def save_irc(irc_job, run_fs, scn_run_fs, scn_save_fs,
             coo_name, irc_idxs):
    """ Read IRC output and store data in filesystem
    """

    # if not scn_run_fs[1].exists([coo_name]):
    #     print("No IRC to save. Skipping...")
    # else:
    opt_ret = rundriver.read_job(
        job=irc_job,
        run_fs=run_fs,
    )
    if opt_ret is not None:
        inf_obj, _, out_str = opt_ret
        # inf_obj, inp_str, out_str = opt_ret
        prog = inf_obj.prog
        geos, gras, hessians = elstruct.reader.irc_points(prog, out_str)
        enes = elstruct.reader.irc_energies(prog, out_str)
        # coords = elstruct.reader.irc_coordinates(prog, out_str)

        # Write the IRC inf file and input file string
        # scn_save_fs[-1].file.geometry_info.write(inf_obj, locs)
        # scn_save_fs[1].file.geometry_input.write(inp_str, locs)

        # Write the IRC coords and enes to a yaml file
        # inf_obj = autofile.system.info.scan_branch({coo_name: irc_idxs})
        # scn_save_fs[1].file.info.write(inf_obj, [coo_name])

        # Write the data for each geom along IRC to the filesystem
        # save_path = scn_save_fs[1].path([coo_name])
        # print(" - Saving...")
        # print(" - Save path: {}".format(save_path))
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


def irc_sp(ts_dct, scn_thy_level, sp_thy_level,
           irc_idxs, overwrite):
    """ Run single points on top of each geom along the IRC
    """

    # Set up coordinate name
    coo_name = 'RC'

    # Set up info objects needed for the run
    ts_info = (ts_dct['ich'], ts_dct['chg'], ts_dct['mul'])
    orb_restr = fsorb.orbital_restriction(ts_info, sp_thy_level)
    sp_level = sp_thy_level[0:3]
    sp_level.append(orb_restr)

    # Build the TS scan file system
    scn_save_fs, scn_run_fs, _, _ = ts_scn_fs(
        ts_dct, ts_info, scn_thy_level)

    # Set up script and kwargs for the irc run
    script_str, _, kwargs, _ = runpar.run_qchem_par(*sp_thy_level[0:2])

    for idx in irc_idxs:
        locs = [[coo_name], [idx]]

        # Set up single point filesys
        scn_run_fs[-1].create(locs)
        scn_run_path = scn_run_fs[-1].path(locs)
        scn_save_path = scn_save_fs[-1].path(locs)
        print('scn_run_path')
        print(scn_run_path)

        # Set up the run filesys for the job
        sp_run_fs = autofile.fs.single_point(scn_run_path)
        sp_save_fs = autofile.fs.single_point(scn_save_path)
        sp_run_fs[-1].create(sp_level[1:4])
        sp_run_path = sp_run_fs[-1].path(sp_level[1:4])
        run_fs = autofile.fs.run(sp_run_path)

        # Read the geometry from the save filesys
        exists = sp_save_fs[-1].file.energy.exists(sp_level[1:4])
        if not exists or overwrite:
            geo = scn_save_fs[-1].file.geometry.read(locs)
            rundriver.run_job(
                job='energy',
                script_str=script_str,
                run_fs=run_fs,
                geom=geo,
                spc_info=ts_info,
                thy_level=sp_thy_level,
                overwrite=overwrite,
                **kwargs,
            )

        ret = rundriver.read_job(
            job='energy',
            run_fs=run_fs,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            print(" - Reading energy from output...")
            ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)

            print(" - Saving energy...")
            sp_save_fs[-1].create(sp_level[1:4])
            sp_save_fs[-1].file.input.write(inp_str, sp_level[1:4])
            sp_save_fs[-1].file.info.write(inf_obj, sp_level[1:4])
            sp_save_fs[-1].file.energy.write(ene, sp_level[1:4])
