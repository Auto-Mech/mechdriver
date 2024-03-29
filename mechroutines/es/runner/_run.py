""" Centralized job runners and readers for electronic structure calcualtions
"""

import functools
import elstruct
import autofile
import automol
from . import _seq as optseq


BASE_ERRS = (elstruct.Error.SCF_NOCONV,
             elstruct.Error.MCSCF_NOCONV,
             elstruct.Error.CC_NOCONV,
             elstruct.Error.LIN_DEP_BASIS)
JOB_ERROR_DCT = {
    elstruct.Job.ENERGY: BASE_ERRS,
    elstruct.Job.GRADIENT: BASE_ERRS,
    elstruct.Job.HESSIAN: BASE_ERRS,
    elstruct.Job.VPT2: BASE_ERRS,
    elstruct.Job.MOLPROP: BASE_ERRS,
    elstruct.Job.OPTIMIZATION: BASE_ERRS+(elstruct.Error.OPT_NOCONV,),
    elstruct.Job.IRCF: BASE_ERRS+(elstruct.Error.IRC_NOCONV,),
    elstruct.Job.IRCR: BASE_ERRS+(elstruct.Error.IRC_NOCONV,)
}

JOB_SUCCESS_DCT = {
    elstruct.Job.ENERGY: elstruct.Success.SCF_CONV,
    elstruct.Job.GRADIENT: elstruct.Success.SCF_CONV,
    elstruct.Job.HESSIAN: elstruct.Success.SCF_CONV,
    elstruct.Job.VPT2: elstruct.Success.SCF_CONV,
    elstruct.Job.MOLPROP: elstruct.Success.SCF_CONV,
    elstruct.Job.OPTIMIZATION: elstruct.Success.OPT_CONV,
    elstruct.Job.IRCF: elstruct.Success.IRC_CONV,
    elstruct.Job.IRCR: elstruct.Success.IRC_CONV,
}

JOB_RUNNER_DCT = {
    elstruct.Job.ENERGY: functools.partial(
        optseq.options_matrix_run, elstruct.writer.energy),
    elstruct.Job.GRADIENT: functools.partial(
        optseq.options_matrix_run, elstruct.writer.gradient),
    elstruct.Job.HESSIAN: functools.partial(
        optseq.options_matrix_run, elstruct.writer.hessian),
    elstruct.Job.VPT2: functools.partial(
        optseq.options_matrix_run, elstruct.writer.vpt2),
    elstruct.Job.MOLPROP: functools.partial(
        optseq.options_matrix_run, elstruct.writer.molecular_properties),
    elstruct.Job.OPTIMIZATION: optseq.options_matrix_optimization,
    elstruct.Job.IRCF: functools.partial(
        optseq.options_matrix_run, elstruct.writer.irc),
    elstruct.Job.IRCR: functools.partial(
        optseq.options_matrix_run, elstruct.writer.irc),
}


def execute_job(job, script_str, run_fs,
                geo, spc_info, thy_info,
                zrxn=None,
                errors=(), options_mat=(),
                retryfail=True, feedback=False,
                frozen_coordinates=(), freeze_dummy_atoms=True,
                overwrite=False,
                **kwargs):
    """ Both ruBoth runs and reads electrouct jobs
    """

    run_job(job, script_str, run_fs,
            geo, spc_info, thy_info,
            zrxn=zrxn,
            errors=errors,
            options_mat=options_mat,
            retryfail=retryfail,
            feedback=feedback,
            frozen_coordinates=frozen_coordinates,
            freeze_dummy_atoms=freeze_dummy_atoms,
            overwrite=overwrite,
            **kwargs)

    success, ret = read_job(job, run_fs)

    return success, ret


def run_job(job, script_str, run_fs,
            geo, spc_info, thy_info,
            zrxn=None,
            errors=(), options_mat=(), retryfail=True, feedback=False,
            frozen_coordinates=(), freeze_dummy_atoms=True, overwrite=False,
            **kwargs):
    """ Run an electronic structure job in the specified RUN filesys layer
        by calling the elstruct package to write the input file with the
        information and executing the job with the specified script string.

        Will first look into the RUN filesys and will either rewrite the
        input and rerun the job if requested.

        :param geo: input molecular geometry or Z-Matrix
        :type geo:
        :param errors: list of error message types to search output for
        :type errors: tuple(str)
        :param options_mat: varopis options to run job with
        :type options_mat: tuple(dict[str: str])
        :param retryfail: re-run the job if failed job found in RUN filesys
        :type retryfail: bool
        :param feedback: update geom with job from previous sequence
        :type feedback: bool
        :param frozen_coordinates: Z-matrix coordinate names to freeze in opts
        :type frozen_coordinates: tuple(str)
        :param freeze_dummy_atoms: freeze any coords defined by dummy atoms
        :type freeze_dummy_atoms: bool
        :param overwrite: overwrite existing input file with new one and rerun
        :type overwrite: bool
        :param kwargs: additional options for electronic structure job
        :type kwargs: dict[str]
    """

    assert job in JOB_RUNNER_DCT
    assert job in JOB_ERROR_DCT
    assert job in JOB_SUCCESS_DCT

    run_path = run_fs[-1].path([job])
    run_fs[-1].create([job])
    run_path = run_fs[-1].path([job])
    if overwrite:
        do_run = True
        print(f" - Running {job} job at {run_path}")
    else:
        if not run_fs[-1].file.info.exists([job]):
            do_run = True
            print(f" - Running {job} job at {run_path}")
        else:
            inf_obj = run_fs[-1].file.info.read([job])
            if inf_obj.status == autofile.schema.RunStatus.FAILURE:
                print(f" - Found failed {job} job at {run_path}")
                if retryfail:
                    print(" - Retrying...")
                    do_run = True
                else:
                    print(" - Skipping failed job, per user request...")
                    do_run = False
            else:
                do_run = False
                if inf_obj.status == autofile.schema.RunStatus.SUCCESS:
                    print(f" - Found completed {job} job at {run_path}")
                else:
                    print(f" - Found running {job} job at {run_path}")
                    print(" - Skipping...")

    if do_run:
        # Create the run directory
        status = autofile.schema.RunStatus.RUNNING
        prog = thy_info[0]
        method = thy_info[1]
        basis = thy_info[2]
        inf_obj = autofile.schema.info_objects.run(
            job=job, prog=prog, version='',
            method=method, basis=basis, status=status)
        inf_obj.utc_start_time = autofile.schema.utc_time()
        run_fs[-1].file.info.write(inf_obj, [job])

        # Write the initial geo/zma
        _write_input_geo(geo, job, run_fs)

        # Set job runner based on user request; set special options as needed
        runner = JOB_RUNNER_DCT[job]

        if job == elstruct.Job.OPTIMIZATION:
            runner = functools.partial(
                runner, feedback=feedback,
                frozen_coordinates=frozen_coordinates,
                freeze_dummy_atoms=freeze_dummy_atoms)
        inp_str, out_str = runner(
            script_str, run_path, geo=geo, chg=spc_info[1],
            mul=spc_info[2], method=thy_info[1], basis=thy_info[2],
            orb_type=thy_info[3], prog=thy_info[0], zrxn=zrxn,
            errors=errors, options_mat=options_mat, **kwargs
        )

        inf_obj.utc_end_time = autofile.schema.utc_time()
        prog = inf_obj.prog
        if is_successful_output(out_str, job, prog):
            run_fs[-1].file.output.write(out_str, [job])
            print(" - Run succeeded.")
            status = autofile.schema.RunStatus.SUCCESS
        else:
            # Added writing output at point even for fail
            # Need to check if this is bad. But read_job changes
            # should address this hopefully
            run_fs[-1].file.output.write(out_str, [job])
            print(" - Run failed.")
            status = autofile.schema.RunStatus.FAILURE
        version = elstruct.reader.program_version(prog, out_str)
        inf_obj.version = version
        inf_obj.status = status
        run_fs[-1].file.info.write(inf_obj, [job])
        run_fs[-1].file.input.write(inp_str, [job])


def read_job(job, run_fs):
    """ Searches for an output file for the specified electronic
        structure job in the specified RUN filesytem. If an output file
        is found, it is parsed for job success messages. If successful,
        function returns job input, output and autofile job info object.

        :param job: label for job formatted to elstruct package definitions
        :type job: str
        :param run_fs: filesystem object for the run filesys where job is run
        :type run_fs: autofile.fs.run object
        :rtype: (bool, (autofile.info_object object???, str, str))
    """

    inf_exists = run_fs[-1].file.info.exists([job])
    inp_exists = run_fs[-1].file.input.exists([job])
    out_exists = run_fs[-1].file.output.exists([job])

    if inf_exists and inp_exists and out_exists:
        inf_obj = run_fs[-1].file.info.read([job])
        inp_str = run_fs[-1].file.input.read([job])
        out_str = run_fs[-1].file.output.read([job])
        prog = inf_obj.prog
        ret = (inf_obj, inp_str, out_str)

        success = bool(is_successful_output(out_str, job, prog))
        if success:
            print(" - Reading successful output...")
    else:
        if not out_exists:
            print(" - No output file found.")
        if not inp_exists:
            print(" - No input file found.")
        if not inf_exists:
            print(" - No info file found.")
        print("Skipping...")
        success = False
        ret = None

    return success, ret


def is_successful_output(out_str, job, prog):
    """ Parses the output string of the electronic structure job
        and calls the appropraite elstruct status readers to assess
        if program has exited normally, contains approprate success messages
        for the job, and precludes error messages signifying job failure.

        :param out_str: string for job output file
        :type out_str: str
        :param job: label for job formatted to elstruct package definitions
        :type job: str
        :param prog: name of the electronic structure program to run
        :type prog: str
        :rtype: bool
    """

    assert job in JOB_ERROR_DCT
    assert job in JOB_SUCCESS_DCT
    errors = JOB_ERROR_DCT[job]
    success = JOB_SUCCESS_DCT[job]

    ret = False
    if elstruct.reader.has_normal_exit_message(prog, out_str):
        for error in errors:
            conv = elstruct.reader.check_convergence_messages(
                prog, error, success, out_str)
        if conv:
            ret = True
        else:
            print(" - Output has an error message. Skipping...")
    else:
        print(' - Output does not contain normal exit message. Skipping...')

    return ret


# Helpers
def _write_input_geo(geo, job, run_fs):
    """ Writes input molecular structures into the specified RUN filesystem.

        Can only reliably a Z-Matrix in the filesystem if the input `geo`
        object is a Z-Matrix and not a geometry (due to conversion issues).

        :param geo: input molecular geometry or Z-Matrix
        :type geo:
        :param job: label for job formatted to elstruct package definitions
        :type job: str
        :param run_fs: filesystem object for the run filesys where job is run
        :type run_fs: autofile.fs.run object
    """

    if automol.zmat.is_valid(geo):
        run_fs[-1].file.zmatrix.write(geo, [job])
        run_fs[-1].file.geometry.write(
            automol.zmat.geometry(geo), [job])
    else:
        run_fs[-1].file.geometry.write(geo, [job])
