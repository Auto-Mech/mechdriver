""" Centralized job runners and readers for electronic structure calcualtions
"""

import functools
import elstruct
import autofile
from . import _optseq as optseq


JOB_ERROR_DCT = {
    elstruct.Job.ENERGY: elstruct.Error.SCF_NOCONV,
    elstruct.Job.GRADIENT: elstruct.Error.SCF_NOCONV,
    elstruct.Job.HESSIAN: elstruct.Error.SCF_NOCONV,
    elstruct.Job.VPT2: elstruct.Error.SCF_NOCONV,
    elstruct.Job.OPTIMIZATION: elstruct.Error.OPT_NOCONV,
    elstruct.Job.IRCF: elstruct.Error.IRC_NOCONV,
    elstruct.Job.IRCR: elstruct.Error.IRC_NOCONV,
}

JOB_SUCCESS_DCT = {
    elstruct.Job.ENERGY: elstruct.Success.SCF_CONV,
    elstruct.Job.GRADIENT: elstruct.Success.SCF_CONV,
    elstruct.Job.HESSIAN: elstruct.Success.SCF_CONV,
    elstruct.Job.VPT2: elstruct.Success.SCF_CONV,
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
    elstruct.Job.OPTIMIZATION: optseq.options_matrix_optimization,
    elstruct.Job.IRCF: functools.partial(
        optseq.options_matrix_run, elstruct.writer.irc),
    elstruct.Job.IRCR: functools.partial(
        optseq.options_matrix_run, elstruct.writer.irc),
}


def run_job(
        job, script_str, run_fs,
        geom, spc_info, thy_info,
        errors=(), options_mat=(), retryfail=True, feedback=False,
        frozen_coordinates=(), freeze_dummy_atoms=True, overwrite=False,
        irc_direction=None,
        **kwargs):
    """ run an elstruct job by name
    """
    assert job in JOB_RUNNER_DCT
    assert job in JOB_ERROR_DCT
    assert job in JOB_SUCCESS_DCT

    run_path = run_fs[-1].path([job])
    run_fs[-1].create([job])
    run_path = run_fs[-1].path([job])
    if overwrite:
        do_run = True
        print(" - Running {} job at {}".format(job, run_path))
    else:
        if not run_fs[-1].file.info.exists([job]):
            do_run = True
            print(" - Running {} job at {}".format(job, run_path))
        else:
            inf_obj = run_fs[-1].file.info.read([job])
            if inf_obj.status == autofile.schema.RunStatus.FAILURE:
                print(" - Found failed {} job at {}".format(job, run_path))
                if retryfail:
                    print(" - Retrying...")
                    do_run = True
                else:
                    print(" - Skipping failed job, per user request...")
                    do_run = False
            else:
                do_run = False
                if inf_obj.status == autofile.schema.RunStatus.SUCCESS:
                    print(" - Found completed {} job at {}"
                          .format(job, run_path))
                else:
                    print(" - Found running {} job at {}"
                          .format(job, run_path))
                    print(" - Skipping...")

    if do_run:
        # create the run directory
        status = autofile.schema.RunStatus.RUNNING
        prog = thy_info[0]
        method = thy_info[1]
        basis = thy_info[2]
        inf_obj = autofile.schema.info.run(
            job=job, prog=prog, version='',
            method=method, basis=basis, status=status)
        inf_obj.utc_start_time = autofile.schema.utc_time()
        run_fs[-1].file.info.write(inf_obj, [job])

        # Set job runner based on user request; set special options as needed
        runner = JOB_RUNNER_DCT[job]

        if job == elstruct.Job.OPTIMIZATION:
            runner = functools.partial(
                runner, feedback=feedback,
                frozen_coordinates=frozen_coordinates,
                freeze_dummy_atoms=freeze_dummy_atoms)
        elif job in (elstruct.Job.IRCF, elstruct.Job.IRCR):
            if job == elstruct.Job.IRCF:
                irc_direction = 'forward'
            else:
                irc_direction = 'reverse'
            runner = functools.partial(
                runner, irc_direction=irc_direction)

        inp_str, out_str = runner(
            script_str, run_path, geom=geom, chg=spc_info[1],
            mul=spc_info[2], method=thy_info[1], basis=thy_info[2],
            orb_type=thy_info[3], prog=thy_info[0],
            errors=errors, options_mat=options_mat, **kwargs
        )

        inf_obj.utc_end_time = autofile.schema.utc_time()
        prog = inf_obj.prog
        if is_successful_output(out_str, job, prog):
            run_fs[-1].file.output.write(out_str, [job])
            print(" - Run succeeded.")
            status = autofile.schema.RunStatus.SUCCESS
        else:
            print(" - Run failed.")
            status = autofile.schema.RunStatus.FAILURE
        version = elstruct.reader.program_version(prog, out_str)
        inf_obj.version = version
        inf_obj.status = status
        run_fs[-1].file.info.write(inf_obj, [job])
        run_fs[-1].file.input.write(inp_str, [job])


def read_job(job, run_fs):
    """ read from an elstruct job by name
    """
    ret = None

    run_path = run_fs[-1].path([job])

    print(" - Reading from {} job at {}".format(job, run_path))
    if not run_fs[-1].file.output.exists([job]):
        print(" - No output file found. Skipping...")
    else:
        assert run_fs[-1].file.info.exists([job])
        assert run_fs[-1].file.input.exists([job])
        inf_obj = run_fs[-1].file.info.read([job])
        inp_str = run_fs[-1].file.input.read([job])
        out_str = run_fs[-1].file.output.read([job])
        prog = inf_obj.prog

        if is_successful_output(out_str, job, prog):
            print(" - Found successful output. Reading...")
            ret = (inf_obj, inp_str, out_str)
        else:
            print(" - Output has an error message. Skipping...")

    return ret


def is_successful_output(out_str, job, prog):
    """ is this a successful output string?
    """
    assert job in JOB_ERROR_DCT
    assert job in JOB_SUCCESS_DCT
    error = JOB_ERROR_DCT[job]
    success = JOB_SUCCESS_DCT[job]

    ret = False
    if elstruct.reader.has_normal_exit_message(prog, out_str):
        conv = elstruct.reader.check_convergence_messages(
            prog, error, success, out_str)
        if conv:
            ret = True

    return ret
