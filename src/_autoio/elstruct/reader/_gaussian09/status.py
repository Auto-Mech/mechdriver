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

    pattern = app.escape('Normal termination of Gaussian 09')

    return apf.has_match(pattern, output_str, case=False)


# Parsers for convergence success messages
def _has_scf_convergence_message(output_str):
    """ Assess whether the output file string contains the
        message signaling successful convergence of the SCF procedure.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    scf_str1 = (
        'Initial convergence to {} achieved.  Increase integral accuracy.' +
        app.LINE_FILL + app.NEWLINE + app.LINE_FILL + app.escape('SCF Done:')
    ).format(app.EXPONENTIAL_FLOAT_D)
    scf_str2 = app.escape('Rotation gradient small -- convergence achieved.')

    pattern = app.one_of_these([scf_str1, scf_str2])

    return apf.has_match(pattern, output_str, case=False)


def _has_opt_convergence_message(output_str):
    """ Assess whether the output file string contains the
        message signaling successful convergence of the geometry optimization.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    pattern = (
        app.escape('Optimization completed.') +
        app.LINE_FILL + app.NEWLINE + app.LINE_FILL +
        app.escape('-- Stationary point found.')
    )

    return apf.has_match(pattern, output_str, case=False)


def _has_irc_convergence_message(output_str):
    """ Assess whether the output file string contains the
        message signaling successful convergence of the
        Intrinsic Reaction Coordinate search.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    pattern = app.escape('Reaction path calculation complete.')

    return apf.has_match(pattern, output_str, case=False)


# Parsers for various error messages
def _has_scf_nonconvergence_error_message(output_str):
    """ Assess whether the output file string contains the
        message signaling the failure of the SCF procedure.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    pattern = app.padded(app.NEWLINE).join([
        app.escape('Convergence criterion not met.'),
        app.escape('SCF Done:')
    ])

    return apf.has_match(pattern, output_str, case=False)


def _has_opt_nonconvergence_error_message(output_str):
    """ Assess whether the output file string contains the
        message signaling the failure of the geometry optimization.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    pattern = app.padded(app.NEWLINE).join([
        app.escape('Optimization stopped.'),
        app.escape('-- Number of steps exceeded,')
    ])

    return apf.has_match(pattern, output_str, case=False)


def _has_irc_nonconvergence_error_message(output_str):
    """ Assess whether the output file string contains the
        message signaling the failure of the
        Intrinsic Reaction Coordinate search.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    pattern = app.escape('Maximum number of corrector steps exceeded')

    return apf.has_match(pattern, output_str, case=False)


ERROR_READER_DCT = {
    elstruct.par.Error.SCF_NOCONV: _has_scf_nonconvergence_error_message,
    elstruct.par.Error.MCSCF_NOCONV: False,
    elstruct.par.Error.CC_NOCONV: False,  # not checked
    elstruct.par.Error.OPT_NOCONV: _has_opt_nonconvergence_error_message,
    elstruct.par.Error.IRC_NOCONV: _has_irc_nonconvergence_error_message,
    elstruct.par.Error.LIN_DEP_BASIS: False  # not checked
}
SUCCESS_READER_DCT = {
    elstruct.par.Success.SCF_CONV: _has_scf_convergence_message,
    elstruct.par.Success.OPT_CONV: _has_opt_convergence_message,
    elstruct.par.Success.IRC_CONV: _has_irc_convergence_message
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
