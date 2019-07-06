""" drivers
"""
import functools
import automol
import elstruct
import autofile
import moldr.runner


def run_conformers(zma, charge, mult, method, basis, orb_restr,
                   nsamp, tors_range_dct, run_prefix, save_prefix, script_str,
                   prog, **kwargs):
    """ run sampling algorithm to find conformers
    """
    if not tors_range_dct:
        print("No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1

    afs = autofile.fs.conformer()

    vma = automol.zmatrix.var_(zma)
    if afs.conf_trunk.dir.exists(save_prefix):
        _vma = afs.conf_trunk.file.vmatrix.read(save_prefix)
        assert vma == _vma
        inf_obj = afs.conf_trunk.file.info.read(save_prefix)
        nsamp = max(nsamp - inf_obj.nsamp, 0)
        print("Found previous saved run. Adjusting nsamp.")
        print("    New nsamp is {:d}.".format(nsamp))
    else:
        afs.conf_trunk.dir.create(save_prefix)
        inf_obj = autofile.system.info.conformer_trunk(
            nsamp=0, tors_ranges=tors_range_dct)
    afs.conf_trunk.file.vmatrix.write(vma, save_prefix)
    afs.conf_trunk.file.info.write(inf_obj, save_prefix)

    samp_zmas = automol.zmatrix.samples(zma, nsamp, tors_range_dct)
    alocs = [[autofile.system.generate_new_conformer_id()]
             for _ in range(nsamp)]

    for idx, (aloc, samp_zma) in enumerate(zip(alocs, samp_zmas)):
        afs.conf.dir.create(run_prefix, aloc)

        path = afs.conf.dir.path(run_prefix, aloc)

        print("Run {}/{}".format(idx+1, nsamp))
        run_job(
            job=elstruct.Job.OPTIMIZATION,
            script_str=script_str,
            prefix=path,
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
        print("Nothing to save. Skipping...")
    else:
        for alocs in afs.conf.dir.existing(run_prefix):
            run_path = afs.conf.dir.path(run_prefix, alocs)
            print("Reading from run at {}".format(run_path))

            ret = read_job(
                job=elstruct.Job.OPTIMIZATION,
                prefix=run_path)
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
