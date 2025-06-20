""" status checkers
"""
import autoparse.pattern as app
import autoparse.find as apf
import elstruct.par


def has_normal_exit_message(output_str):
    """ does this output string have a normal exit message?
    """
    pattern = app.escape('Normal termination of mrcc.')
    return apf.has_match(pattern, output_str, case=False)


def _has_scf_nonconvergence_error_message(output_str):
    """ does this output string have an SCF non-convergence message?
    """
    pattern = app.padded(app.NEWLINE).join([
        app.escape('THE SCF ITERATION HAS NOT CONVERGED,'),
        app.escape('IN MAXIMAL NUMBER OF STEPS SET BY USER!')
    ])
    return apf.has_match(pattern, output_str, case=False)


ERROR_READER_DCT = {
    elstruct.par.Error.SCF_NOCONV: _has_scf_nonconvergence_error_message,
    elstruct.par.Error.MCSCF_NOCONV: False,
    elstruct.par.Error.CC_NOCONV: False,  # not checked
    elstruct.par.Error.OPT_NOCONV: False,
    elstruct.par.Error.IRC_NOCONV: False,
    elstruct.par.Error.LIN_DEP_BASIS: False
}


def error_list():
    """ list of errors that be identified from the output file
    """
    return tuple(sorted(ERROR_READER_DCT.keys()))


def has_error_message(error, output_str):
    """ does this output string have an error message?
    """

    assert error in error_list()

    error_reader = ERROR_READER_DCT[error]
    if isinstance(error_reader, bool):
        err_val = False
    else:
        err_val = error_reader(output_str)

    return err_val


def check_convergence_messages(error, success, output_str):
    """ check if error messages should trigger job success or failure
    """
    assert error in error_list()
    # assert success in sucess_list()
    _ = success

    job_success = True
    has_error = ERROR_READER_DCT[error](output_str)
    if has_error:
        job_success = False

    return job_success
