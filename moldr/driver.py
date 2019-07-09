""" drivers
"""
import functools
import numpy
import automol
import elstruct
import autofile
import moldr.runner


# conformer sampling
def run_conformers(zma, charge, mult, method, basis, orb_restr,
                   nsamp, tors_range_dct, run_prefix, save_prefix, script_str,
                   prog, **kwargs):
    """ run sampling algorithm to find conformers
    """
    if not tors_range_dct:
        print("No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1

    afs = autofile.fs.conformer()
    afs.conf_trunk.dir.create(save_prefix)

    vma = automol.zmatrix.var_(zma)
    if afs.conf_trunk.file.vmatrix.exists(save_prefix):
        existing_vma = afs.conf_trunk.file.vmatrix.read(save_prefix)
        assert vma == existing_vma

    if afs.conf_trunk.file.info.exists(save_prefix):
        inf_obj = afs.conf_trunk.file.info.read(save_prefix)
        nsamp = max(nsamp - inf_obj.nsamp, 0)
        print("Found previous saved run. Adjusting nsamp.")
        print("    New nsamp is {:d}.".format(nsamp))
    else:
        inf_obj = autofile.system.info.conformer_trunk(
            nsamp=0, tors_ranges=tors_range_dct)
    afs.conf_trunk.file.vmatrix.write(vma, save_prefix)
    afs.conf_trunk.file.info.write(inf_obj, save_prefix)

    samp_zmas = automol.zmatrix.samples(zma, nsamp, tors_range_dct)
    alocs = [[autofile.system.generate_new_conformer_id()]
             for _ in range(nsamp)]

    for idx, (aloc, samp_zma) in enumerate(zip(alocs, samp_zmas)):
        afs.conf.dir.create(run_prefix, aloc)

        run_path = afs.conf.dir.path(run_prefix, aloc)

        print("Run {}/{}".format(idx+1, nsamp))
        run_job(
            job=elstruct.Job.OPTIMIZATION,
            script_str=script_str,
            prefix=run_path,
            geom=samp_zma,
            charge=charge,
            mult=mult,
            method=method,
            basis=basis,
            orb_restr=orb_restr,
            prog=prog,
            **kwargs
        )


def save_conformers(run_prefix, save_prefix):
    """ save the conformers that have been found so far
    """
    afs = autofile.fs.conformer()

    saved_geos = [afs.conf.file.geometry.read(save_prefix, alocs)
                  for alocs in afs.conf.dir.existing(save_prefix)]

    if not afs.conf_trunk.dir.exists(run_prefix):
        print("No conformers to save. Skipping...")
    else:
        for alocs in afs.conf.dir.existing(run_prefix):
            run_path = afs.conf.dir.path(run_prefix, alocs)
            print("Reading from conformer run at {}".format(run_path))

            ret = read_job(job=elstruct.Job.OPTIMIZATION, prefix=run_path)
            if ret:
                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)

                geo = elstruct.reader.opt_geometry(prog, out_str)

                if not _unique_coulomb_spectrum(geo, saved_geos):
                    print(" - Geometry is not unique. Skipping...")
                else:
                    save_path = afs.conf.dir.path(save_prefix, alocs)
                    print(" - Geometry is unique. Saving...")
                    print(" - Save path: {}".format(save_path))

                    afs.conf.dir.create(save_prefix, alocs)
                    afs.conf.file.geometry_info.write(
                        inf_obj, save_prefix, alocs)
                    afs.conf.file.geometry_input.write(
                        inp_str, save_prefix, alocs)
                    afs.conf.file.energy.write(ene, save_prefix, alocs)
                    afs.conf.file.geometry.write(geo, save_prefix, alocs)

                    saved_geos.append(geo)

        # update the number of samples in the trunk information file
        alocs_lst = afs.conf.dir.existing(save_prefix)
        trunk_inf_obj = afs.conf_trunk.file.info.read(save_prefix)
        trunk_inf_obj.nsamp = len(alocs_lst)
        afs.conf_trunk.file.info.write(trunk_inf_obj, save_prefix)

        # update the conformer trajectory file
        if alocs_lst:
            enes = [afs.conf.file.energy.read(save_prefix, alocs)
                    for alocs in alocs_lst]
            geos = [afs.conf.file.geometry.read(save_prefix, alocs)
                    for alocs in alocs_lst]

            traj = []
            for ene, geo in sorted(zip(enes, geos), key=lambda x: x[0]):
                comment = 'energy: {:>15.10f}'.format(ene)
                traj.append((comment, geo))

            traj_path = afs.conf_trunk.file.trajectory.path(save_prefix)
            print("Updating conformer trajectory file at {}".format(traj_path))
            afs.conf_trunk.file.trajectory.write(traj, save_prefix)


def _unique_coulomb_spectrum(geo, seen_geos, rtol=2e-5):
    return not any(
        automol.geom.almost_equal_coulomb_spectrum(geo, seen_geo, rtol)
        for seen_geo in seen_geos)


# constrained optimization scans
def run_scan(zma, charge, mult, method, basis, orb_restr,
             grid_dct, run_prefix, save_prefix, script_str,
             prog, update_guess=True, **kwargs):
    """ run constrained optimization scan
    """
    if len(grid_dct) > 1:
        raise NotImplementedError

    afs = autofile.fs.scan()
    afs.scan_trunk.dir.create(save_prefix)

    vma = automol.zmatrix.var_(zma)
    if afs.scan_trunk.file.vmatrix.exists(save_prefix):
        existing_vma = afs.scan_trunk.file.vmatrix.read(save_prefix)
        assert vma == existing_vma

    # for now, running only one-dimensional hindered rotor scans
    ((coo_name, grid_vals),) = grid_dct.items()
    afs.scan_branch.dir.create(save_prefix, [[coo_name]])

    if afs.scan_branch.file.info.exists(save_prefix, [[coo_name]]):
        inf_obj = afs.scan_branch.file.info.read(save_prefix, [[coo_name]])
        existing_grid_dct = dict(inf_obj.grids)
        assert numpy.allclose([grid_vals], existing_grid_dct[coo_name])

    inf_obj = autofile.system.info.scan_branch({coo_name: grid_vals})
    afs.scan_branch.file.info.write(inf_obj, save_prefix, [[coo_name]])

    npoint = len(grid_vals)
    grid_idxs = tuple(range(npoint))

    for grid_idx in grid_idxs:
        afs.scan.dir.create(run_prefix, [[coo_name], [grid_idx]])

    prefixes = tuple(afs.scan.dir.path(run_prefix, [[coo_name], [grid_idx]])
                     for grid_idx in grid_idxs)

    _run_1d_scan(
        script_str=script_str,
        prefixes=prefixes,
        guess_zma=zma,
        coo_name=coo_name,
        grid_idxs=grid_idxs,
        grid_vals=grid_vals,
        charge=charge,
        mult=mult,
        method=method,
        basis=basis,
        orb_restr=orb_restr,
        prog=prog,
        update_guess=update_guess,
        **kwargs
    )

    _run_1d_scan(
        script_str=script_str,
        prefixes=list(reversed(prefixes)),
        guess_zma=zma,
        coo_name=coo_name,
        grid_idxs=list(reversed(grid_idxs)),
        grid_vals=list(reversed(grid_vals)),
        charge=charge,
        mult=mult,
        method=method,
        basis=basis,
        orb_restr=orb_restr,
        prog=prog,
        update_guess=update_guess,
        **kwargs
    )


def _run_1d_scan(script_str, prefixes,
                 guess_zma, coo_name, grid_idxs, grid_vals,
                 charge, mult, method, basis, orb_restr, prog,
                 errors=(), options_mat=(), retry_failed=True,
                 update_guess=True,
                 **kwargs):
    npoints = len(grid_idxs)
    assert len(grid_vals) == len(prefixes) == npoints
    for grid_idx, grid_val, prefix in zip(grid_idxs, grid_vals, prefixes):
        print("Point {}/{}".format(grid_idx+1, npoints))
        zma = automol.zmatrix.set_values(guess_zma, {coo_name: grid_val})

        run_job(
            job=elstruct.Job.OPTIMIZATION,
            script_str=script_str,
            prefix=prefix,
            geom=zma,
            charge=charge,
            mult=mult,
            method=method,
            basis=basis,
            orb_restr=orb_restr,
            prog=prog,
            frozen_coordinates=[coo_name],
            errors=errors,
            options_mat=options_mat,
            retry_failed=retry_failed,
            **kwargs
        )

        ret = read_job(job=elstruct.Job.OPTIMIZATION, prefix=prefix)
        if update_guess and ret is not None:
            inf_obj, _, out_str = ret
            prog = inf_obj.prog
            method = inf_obj.method
            guess_zma = elstruct.reader.opt_zmatrix(prog, out_str)


def save_scan(run_prefix, save_prefix, coo_names):
    """ save the scans that have been run so far
    """
    if len(coo_names) > 1:
        raise NotImplementedError

    afs = autofile.fs.scan()
    if not afs.scan_branch.dir.exists(run_prefix, [coo_names]):
        print("No scan to save. Skipping...")
    else:
        alocs_lst = []

        for rlocs in afs.scan.dir.existing(run_prefix, [coo_names]):
            alocs = [coo_names] + rlocs
            run_path = afs.scan.dir.path(run_prefix, alocs)
            print("Reading from scan run at {}".format(run_path))

            ret = read_job(job=elstruct.Job.OPTIMIZATION, prefix=run_path)
            if ret:
                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)
                geo = elstruct.reader.opt_geometry(prog, out_str)

                save_path = afs.scan.dir.path(save_prefix, alocs)
                print(" - Saving...")
                print(" - Save path: {}".format(save_path))

                afs.scan.dir.create(save_prefix, alocs)
                afs.scan.file.geometry_info.write(inf_obj, save_prefix, alocs)
                afs.scan.file.geometry_input.write(inp_str, save_prefix, alocs)
                afs.scan.file.energy.write(ene, save_prefix, alocs)
                afs.scan.file.geometry.write(geo, save_prefix, alocs)

                alocs_lst.append(alocs)

        if alocs_lst:
            idxs_lst = [alocs[-1] for alocs in alocs_lst]
            enes = [afs.scan.file.energy.read(save_prefix, alocs)
                    for alocs in alocs_lst]
            geos = [afs.scan.file.geometry.read(save_prefix, alocs)
                    for alocs in alocs_lst]

            traj = []
            for idxs, ene, geo in zip(idxs_lst, enes, geos):
                comment = 'energy: {:>15.10f}, grid idxs: {}'.format(ene, idxs)
                traj.append((comment, geo))

            traj_path = afs.scan_branch.file.trajectory.path(
                save_prefix, [coo_names])
            print("Updating scan trajectory file at {}".format(traj_path))
            afs.scan_branch.file.trajectory.write(
                traj, save_prefix, [coo_names])


# centralized job runners and readers
JOB_ERROR_DCT = {
    elstruct.Job.ENERGY: elstruct.Error.SCF_NOCONV,
    elstruct.Job.GRADIENT: elstruct.Error.SCF_NOCONV,
    elstruct.Job.HESSIAN: elstruct.Error.SCF_NOCONV,
    elstruct.Job.OPTIMIZATION: elstruct.Error.OPT_NOCONV,
}

JOB_RUNNER_DCT = {
    elstruct.Job.ENERGY: functools.partial(
        moldr.runner.options_matrix_run, elstruct.writer.energy),
    elstruct.Job.GRADIENT: functools.partial(
        moldr.runner.options_matrix_run, elstruct.writer.gradient),
    elstruct.Job.HESSIAN: functools.partial(
        moldr.runner.options_matrix_run, elstruct.writer.hessian),
    elstruct.Job.OPTIMIZATION: moldr.runner.feedback_optimization,
}


def run_job(job, script_str, prefix,
            geom, charge, mult, method, basis, orb_restr, prog,
            errors=(), options_mat=(), retry_failed=True,
            **kwargs):
    """ run an elstruct job by name
    """
    assert job in JOB_RUNNER_DCT
    assert job in JOB_ERROR_DCT

    afs = autofile.fs.run()

    run_path = afs.run.dir.path(prefix, [job])
    if not afs.run.file.info.exists(prefix, [job]):
        do_run = True
        print(" - Running {} job at {}".format(job, run_path))
    else:
        inf_obj = afs.run.file.info.read(prefix, [job])
        if inf_obj.status == autofile.system.RunStatus.FAILURE:
            print(" - Found failed {} job at {}".format(job, run_path))
            if retry_failed:
                print(" - Retrying...")
                do_run = True
            else:
                do_run = False
        else:
            do_run = False
            if inf_obj.status == autofile.system.RunStatus.SUCCESS:
                print(" - Found completed {} job at {}".format(job, run_path))
            else:
                print(" - Found running {} job at {}".format(job, run_path))
            print(" - Skipping...")

    if do_run:
        # create the run directory
        afs.run.dir.create(prefix, [job])

        run_path = afs.run.dir.path(prefix, [job])

        status = autofile.system.RunStatus.RUNNING
        inf_obj = autofile.system.info.run(
            job=job, prog=prog, method=method, basis=basis, status=status)
        inf_obj.utc_start_time = autofile.system.info.utc_time()
        afs.run.file.info.write(inf_obj, prefix, [job])

        runner = JOB_RUNNER_DCT[job]

        print(" - Starting the run...")
        inp_str, out_str = runner(
            script_str, run_path,
            geom=geom, charge=charge, mult=mult, method=method, basis=basis,
            orb_restricted=orb_restr, prog=prog, errors=errors,
            options_mat=options_mat,
            **kwargs
        )

        inf_obj.utc_end_time = autofile.system.info.utc_time()

        if is_successful_output(out_str, job, prog):
            afs.run.file.output.write(out_str, prefix, [job])
            print(" - Run succeeded.")
            status = autofile.system.RunStatus.SUCCESS
        else:
            print(" - Run failed.")
            status = autofile.system.RunStatus.FAILURE
        inf_obj.status = status
        afs.run.file.info.write(inf_obj, prefix, [job])
        afs.run.file.input.write(inp_str, prefix, [job])


def read_job(job, prefix):
    """ read from an elstruct job by name
    """
    ret = None

    afs = autofile.fs.run()
    path = afs.run.dir.path(prefix, [job])
    print(" - Reading from {} job at {}".format(job, path))
    if not afs.run.file.output.exists(prefix, [job]):
        print(" - No output file found. Skipping...")
    else:
        assert afs.run.file.info.exists(prefix, [job])
        assert afs.run.file.input.exists(prefix, [job])
        inf_obj = afs.run.file.info.read(prefix, [job])
        inp_str = afs.run.file.input.read(prefix, [job])
        out_str = afs.run.file.output.read(prefix, [job])
        prog = inf_obj.prog

        if is_successful_output(out_str, job, prog):
            print(" - Found succesfull output. Reading...")
            ret = (inf_obj, inp_str, out_str)
        else:
            print(" - Output has an error message. Skipping...")

    return ret


def is_successful_output(out_str, job, prog):
    """ is this a successful output string?
    """
    assert job in JOB_ERROR_DCT
    error = JOB_ERROR_DCT[job]

    ret = False
    if elstruct.reader.has_normal_exit_message(prog, out_str):
        if not elstruct.reader.has_error_message(prog, error, out_str):
            ret = True

    return ret
