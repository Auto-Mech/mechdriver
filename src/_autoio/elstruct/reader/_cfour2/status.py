""" status checkers
"""

import autoparse.pattern as app
import autoparse.find as apf
import elstruct.par


# Exit message for the program
def has_normal_exit_message(output_str):
    """ Assess whether the output file string contains the
        normal program exit message.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    pattern = ('The final electronic energy is' +
               app.one_or_more(app.WILDCARD) +
               'This computation required')

    return apf.has_match(pattern, output_str, case=False)


# Parsers for convergence success messages
def _has_scf_convergence_message(output_str):
    """ Assess whether the output file string contains the
        message signaling successful convergence of the SCF procedure.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    pattern = 'SCF has converged.'

    return apf.has_match(pattern, output_str, case=False)


def _has_cc_convergence_message(output_str):
    """ Assess whether the output file string contains the
        message signaling the failure of the coupled-cluster iterations.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    pattern = app.escape('timing for (T)')
    match = apf.has_match(pattern, output_str, case=False)
    if match is None:
        pattern = (
            'A miracle come to pass. The CC iterations have converged.')
        match = apf.has_match(pattern, output_str, case=False)

    return match


def _has_opt_convergence_message(output_str):
    """ Assess whether the output file string contains the
        message signaling successful convergence of the geometry optimization.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    pattern = app.escape(
        'Convergence criterion satisfied.  Optimization completed.')

    return apf.has_match(pattern, output_str, case=False)


# Parsers for various error messages
def _has_scf_nonconvergence_error_message(output_str):
    """ Assess whether the output file string contains the
        message signaling the failure of the SCF procedure.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    pattern = 'SCF failed to converge'

    return apf.has_match(pattern, output_str, case=False)


def _has_cc_nonconvergence_error_message(output_str):
    """ Assess whether the output file string contains the
        message signaling the failure of the coupled-cluster iterations.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    pattern = app.escape('CC did not converge !!!')

    return apf.has_match(pattern, output_str, case=False)


def _has_opt_nonconvergence_error_message(output_str):
    """ Assess whether the output file string contains the
        message signaling the failure of the geometry optimization.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    pattern = app.escape('*Maximum number of optimization steps exceeded.')

    return apf.has_match(pattern, output_str, case=False)


ERROR_READER_DCT = {
    elstruct.par.Error.SCF_NOCONV: _has_scf_nonconvergence_error_message,
    elstruct.par.Error.MCSCF_NOCONV: False,
    elstruct.par.Error.CC_NOCONV: _has_cc_nonconvergence_error_message,
    elstruct.par.Error.OPT_NOCONV: _has_opt_nonconvergence_error_message,
    elstruct.par.Error.IRC_NOCONV: False,
    elstruct.par.Error.LIN_DEP_BASIS: False  # not checked
}
SUCCESS_READER_DCT = {
    elstruct.par.Success.SCF_CONV: _has_scf_convergence_message,
    elstruct.par.Success.CC_CONV: _has_cc_convergence_message,
    elstruct.par.Success.OPT_CONV: _has_opt_convergence_message,
}


def error_list():
    """ Constructs a list of errors that be identified from the output file.
    """
    return tuple(sorted(ERROR_READER_DCT.keys()))


def success_list():
    """ Constructs a list of successes that be identified from the output file.
    """
    return tuple(sorted(SUCCESS_READER_DCT.keys()))


def has_error_message(error, output_str):
    """ Assess whether the output file string contains error messages
        for any of the procedures in the job.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    assert error in error_list()

    error_reader = ERROR_READER_DCT[error]
    if isinstance(error_reader, bool):
        err_val = False
    else:
        err_val = error_reader(output_str)

    return err_val


def check_convergence_messages(error, success, output_str):
    """ Assess whether the output file string contains messages
        denoting all of the requested procedures in the job have converged.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    assert error in error_list()
    assert success in success_list()

    if has_error_message(error, output_str):
        job_success = SUCCESS_READER_DCT[success](output_str)
    else:
        job_success = True

    return job_success
