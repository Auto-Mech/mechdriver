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

    pattern = 'Molpro calculation terminated'
    # pattern = app.padded(app.NEWLINE).join([
    #     app.escape('Molpro calculation terminated'),
    #    app.escape('Variable memory released')])

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
    scf_str3 = app.escape('The problem occurs in Multi')
    scf_str4 = app.escape('The problem occurs in cipro')

    pattern = app.one_of_these([scf_str1, scf_str2, scf_str3, scf_str4])

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


# Parsers for various error messages
def _has_scf_nonconvergence_error_message(output_str):
    """ Assess whether the output file string contains the
        message signaling the failure of the SCF procedure.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    pattern = app.escape('No convergence') + app.not_followed_by(
        app.padded('in max. number of iterations'))

    return apf.has_match(pattern, output_str, case=False)


def _has_mcscf_nonconvergence_error_message(output_str):
    """ Assess whether the output file string contains the
        message signaling the failure of the MCSCF procedure.

        When an MCSCF WF calculation fails to converge at
        some step in opt it is not clear how to save the
        optimization, so we give up on opt

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    pattern = app.one_of_these([
        app.escape('The problem occurs in Multi'),
        app.escape('The problem occurs in cipro')
    ])
    return apf.has_match(pattern, output_str, case=True)


def _has_opt_nonconvergence_error_message(output_str):
    """ Assess whether the output file string contains the
        message signaling the failure of the geometry optimization.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    # First check for scf convergence issues
    if _has_scf_nonconvergence_error_message(output_str):
        error = True
    else:
        pattern = app.escape('No convergence in max. number of iterations')
        error = apf.has_match(pattern, output_str, case=False)
    return error


def _has_linear_dependence_error_message(output_str):
    """ Assess whether the output file string contains the
        message signaling Hartee-Fock failure as a result of
        linear dependencies in the basis set.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """
    pattern = 'ERROR: BASIS LINEARLY DEPENDENT OR WRONG S'
    return apf.has_match(pattern, output_str, case=False)


ERROR_READER_DCT = {
    elstruct.par.Error.SCF_NOCONV: _has_scf_nonconvergence_error_message,
    elstruct.par.Error.MCSCF_NOCONV: _has_mcscf_nonconvergence_error_message,
    elstruct.par.Error.CC_NOCONV: False,  # not checked
    elstruct.par.Error.OPT_NOCONV: _has_opt_nonconvergence_error_message,
    elstruct.par.Error.IRC_NOCONV: False,
    elstruct.par.Error.LIN_DEP_BASIS: _has_linear_dependence_error_message
}
SUCCESS_READER_DCT = {
    elstruct.par.Success.SCF_CONV: _has_scf_convergence_message,
    elstruct.par.Success.OPT_CONV: _has_opt_convergence_message
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

    if not job_success:
        print(f'  - Error found: {error}')

    return job_success
