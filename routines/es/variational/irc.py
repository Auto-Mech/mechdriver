"""
IRC calcs
"""

<<<<<<< HEAD
import automol
=======
>>>>>>> 9bf6c45022eda7c7b6875f3290b6564dd75cf7b9
import elstruct
import autofile
from lib.filesystem import path as fpath
from lib.filesystem import minc as fsmin
from lib.filesystem import orb as fsorb
from lib.runner import par as runpar
from lib.runner import driver as rundriver


def irc_opt(ts_dct, thy_info, irc_idxs, overwrite):
    """ Run the IRC
    """

    # Set up coordinate name
    coo_name = 'RC'

    # Set up info objects needed for the run
    ts_info = (ts_dct['ich'], ts_dct['chg'], ts_dct['mul'])
    orb_restr = fsorb.orbital_restriction(ts_info, thy_info)
    scn_thy_level = thy_info[0:3]
    scn_thy_level.append(orb_restr)

    # Build the TS scan file system
    scn_save_fs, scn_run_fs, run_fs, cnf_save_fs = ts_scn_fs(
        ts_dct, ts_info, thy_info)

    # Get the geometry from the filesystem
    zma = cnf_save_fs[-1].file.zmatrix.read(
        fsmin.min_energy_conformer_locators(cnf_save_fs))

    # Run and Read the IRC in the forward and reverse direction
    for irc_direction in ('forward', 'reverse'):

        # Set elstruct job type based on the direction
        if irc_direction == 'forward':
            job = elstruct.Job.IRCF
        else:
            job = elstruct.Job.IRCR

        # Run and read the IRC jobs
        run_irc(
            irc_direction,
            job,
            run_fs,
            zma,
            ts_info,
            scn_thy_level,
            overwrite
        )

        save_irc(
            irc_direction,
            job,
            run_fs,
            scn_run_fs,
            scn_save_fs,
            coo_name,
            irc_idxs
        )


def run_irc(irc_direction, job, run_fs, zma, ts_info, scn_thy_info, overwrite):
    """ Run the irc job
    """

    # Set up script and kwargs for the irc run
    _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
        *scn_thy_info[0:2])

    # set irc options here for now
    opt_kwargs['job_options'] = ['calcall', 'stepsize=3', 'maxpoints=4']

    # Run the calculations
    rundriver.run_job(
        job=job,
        script_str=opt_script_str,
        run_fs=run_fs,
        geom=zma,
        spc_info=ts_info,
        thy_level=scn_thy_info,
        overwrite=overwrite,
        irc_direction=irc_direction,  # remove run_job if keeping dif job types
        **opt_kwargs,
        )


def save_irc(irc_direction, job,
             run_fs, scn_run_fs, scn_save_fs,
             coo_name, irc_idxs):
    """ Read IRC output and store data in filesystem
    """

    # if not scn_run_fs[1].exists([coo_name]):
    #     print("No IRC to save. Skipping...")
    # else:
    opt_ret = rundriver.read_job(
        job=job,
        run_fs=run_fs,
    )
    if opt_ret is not None:
        inf_obj, _, out_str = opt_ret
        # inf_obj, inp_str, out_str = opt_ret
        prog = inf_obj.prog
        geos, gras, hessians = elstruct.reader.irc_points(prog, out_str)
        enes = elstruct.reader.irc_energies(prog, out_str)
<<<<<<< HEAD

        print(irc_direction)
        print('enes')
        for x in enes:
            print(x)
            print('\n')
        print('geos')
        for x in geos:
            print(automol.geom.string(x))
            print('\n')
=======
>>>>>>> 9bf6c45022eda7c7b6875f3290b6564dd75cf7b9
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

<<<<<<< HEAD
            # Set locs idx; for reverse, ignore SadPt and flip idx to negative
            locs_idx = idx
            if irc_direction == 'reverse':
                if locs_idx == 0:
                    continue
                locs_idx *= -1
            locs = [[coo_name], [locs_idx]]
=======
            # For reverse direction, ignore TS and flip idx to negative
            if irc_direction == 'reverse':
                if idx == 0:
                    continue
                idx *= -1
            locs = [[coo_name], [idx]]
>>>>>>> 9bf6c45022eda7c7b6875f3290b6564dd75cf7b9
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


def ts_scn_fs(ts_dct, ts_info, thy_info):
    """ set up the scan fs for a ts
    """
    # Set up info objects needed for the run
    orb_restr = fsorb.orbital_restriction(ts_info, thy_info)
    ref_level = thy_info[0:3]
    ref_level.append(orb_restr)

    # Build the conformer filesytems Get the file system
    [_, _, rxn_run_path, rxn_save_path] = ts_dct['rxn_fs']
    _, _, ts_run_path, ts_save_path = fpath.get_ts_fs(
        rxn_run_path, rxn_save_path, ref_level)
    cnf_run_fs = autofile.fs.conformer(ts_run_path)
    cnf_save_fs = autofile.fs.conformer(ts_save_path)

    # Set up the scan filesys for the minimum-ene conformer of the TS
    min_cnf_locs = fsmin.min_energy_conformer_locators(cnf_save_fs)
    cnf_run_path = cnf_run_fs[-1].path(min_cnf_locs)
    cnf_save_path = cnf_save_fs[-1].path(min_cnf_locs)
    scn_save_fs = autofile.fs.scan(cnf_save_path)
    scn_run_fs = autofile.fs.scan(cnf_run_path)

    # Set up the IRC run filesystem
    run_fs = autofile.fs.run(cnf_run_path)

    return scn_save_fs, scn_run_fs, run_fs, cnf_save_fs
