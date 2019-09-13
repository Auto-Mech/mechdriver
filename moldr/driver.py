""" drivers
"""
import functools
from qcelemental import constants as qcc
import elstruct
import autofile
from moldr import runner

WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

# centralized job runners and readers
JOB_ERROR_DCT = {
    elstruct.Job.ENERGY: elstruct.Error.SCF_NOCONV,
    elstruct.Job.GRADIENT: elstruct.Error.SCF_NOCONV,
    elstruct.Job.HESSIAN: elstruct.Error.SCF_NOCONV,
    elstruct.Job.OPTIMIZATION: elstruct.Error.OPT_NOCONV,
}

JOB_RUNNER_DCT = {
    elstruct.Job.ENERGY: functools.partial(
        runner.options_matrix_run, elstruct.writer.energy),
    elstruct.Job.GRADIENT: functools.partial(
        runner.options_matrix_run, elstruct.writer.gradient),
    elstruct.Job.HESSIAN: functools.partial(
        runner.options_matrix_run, elstruct.writer.hessian),
    elstruct.Job.OPTIMIZATION: runner.options_matrix_optimization,
}


def run_job(
        job, script_str, run_fs,
        geom, spc_info, thy_level,
        errors=(), options_mat=(), retry_failed=True, feedback=False,
        frozen_coordinates=(), freeze_dummy_atoms=True, overwrite=False,
        **kwargs):
    """ run an elstruct job by name
    """
    assert job in JOB_RUNNER_DCT
    assert job in JOB_ERROR_DCT

    run_fs.leaf.create([job])
    run_path = run_fs.leaf.path([job])

    if overwrite:
        do_run = True
        print(" - Running {} job at {}".format(job, run_path))
    else:
        if not run_fs.leaf.file.info.exists([job]):
            do_run = True
            print(" - Running {} job at {}".format(job, run_path))
        else:
            inf_obj = run_fs.leaf.file.info.read([job])
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
                    print(" - Found completed {} job at {}"
                          .format(job, run_path))
                else:
                    print(" - Found running {} job at {}"
                          .format(job, run_path))
                    print(" - Skipping...")

    # I think the len(geom) test was supposed to stop the call for gradients and hessians for atoms
    # It didn't work properly for two reasons:
    # (i) if the geom is passed as a z-matrix, len(geom) is 2 not 1
    # (ii) it was checking for UPPERCASE job words and its or logic made it always true
    # (ii) if I fix those problems then if you turn off the optimization, then the conformer directories are not created
    # which creates evern more problems. For now I just turned these tests off, since gaussians seems to run fine for atom
    # gradients and hessians, etc.
#    print('do_run_0 test:', do_run)
#    print('geom test:',geom)
#    print(len(geom))
#    if len(geom) == 1:
#        print(len(geom))
#        print(job)
#        if job in ('hessian', 'optimization', 'gradient', 'vpt2'):
#            print(job)
#            do_run = False
#            print('do_run_1 test:', do_run)

#    print ('do_run test:', do_run)
    if do_run:
        # create the run directory
        status = autofile.system.RunStatus.RUNNING
        prog = thy_level[0]
        method = thy_level[1]
        basis = thy_level[2]
        inf_obj = autofile.system.info.run(
            job=job, prog=prog, method=method, basis=basis, status=status)
        inf_obj.utc_start_time = autofile.system.info.utc_time()
        run_fs.leaf.file.info.write(inf_obj, [job])

        runner = JOB_RUNNER_DCT[job]
        if job == elstruct.Job.OPTIMIZATION:
            runner = functools.partial(
                runner, feedback=feedback,
                frozen_coordinates=frozen_coordinates,
                freeze_dummy_atoms=freeze_dummy_atoms)

        inp_str, out_str = runner(
            script_str, run_path, geom=geom, chg=spc_info[1],
            mul=spc_info[2], method=thy_level[1], basis=thy_level[2],
            orb_restricted=thy_level[3], prog=thy_level[0],
            errors=errors, options_mat=options_mat, **kwargs
        )

        inf_obj.utc_end_time = autofile.system.info.utc_time()
        prog = inf_obj.prog

        if is_successful_output(out_str, job, prog):
            run_fs.leaf.file.output.write(out_str, [job])
            print(" - Run succeeded.")
            status = autofile.system.RunStatus.SUCCESS
        else:
            print(" - Run failed.")
            status = autofile.system.RunStatus.FAILURE
        inf_obj.status = status
        run_fs.leaf.file.info.write(inf_obj, [job])
        run_fs.leaf.file.input.write(inp_str, [job])
        print('finished run_job')


def read_job(job, run_fs):
    """ read from an elstruct job by name
    """
    ret = None

    run_path = run_fs.leaf.path([job])

    print(" - Reading from {} job at {}".format(job, run_path))
    if not run_fs.leaf.file.output.exists([job]):
        print(" - No output file found. Skipping...")
    else:
        assert run_fs.leaf.file.info.exists([job])
        assert run_fs.leaf.file.input.exists([job])
        inf_obj = run_fs.leaf.file.info.read([job])
        inp_str = run_fs.leaf.file.input.read([job])
        out_str = run_fs.leaf.file.output.read([job])
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
    error = JOB_ERROR_DCT[job]

    ret = False
    if elstruct.reader.has_normal_exit_message(prog, out_str):
        if not elstruct.reader.has_error_message(prog, error, out_str):
            ret = True

    return ret


