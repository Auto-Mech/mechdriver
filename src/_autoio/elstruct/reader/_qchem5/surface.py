""" potential energy surface information readers
"""

import autoread as ar
import autoparse.pattern as app
import autoparse.find as apf


def gradient(output_str):
    """ Reads the molecular gradient (in Cartesian coordinates) from
        the output file string. Returns the gradient in atomic units.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(tuple(float))
    """

    # Try to read the gradient from several differnt methods
    ptt1 = (
        'Gradient of SCF Energy' + app.NEWLINE +
        app.SPACES + app.series(app.INTEGER, app.SPACES) + app.NEWLINE)
    ptt2 = (
        'Full Analytical Gradient' + app.LINE_FILL + app.NEWLINE +
        app.SPACES + app.series(app.INTEGER, app.SPACES) + app.NEWLINE)
    start_ptt = app.one_of_these([ptt1, ptt2])

    grad = ar.matrix.read(
        output_str,
        val_ptt=app.FLOAT,
        start_ptt=start_ptt,
        line_start_ptt=app.UNSIGNED_INTEGER)

    # Try and read a general tensor from a numerical gradient
    if grad is None:
        grad = _general_xyz_tensor(output_str)

    return grad


def hessian(output_str):
    """ Reads the molecular Hessian (in Cartesian coordinates) from
        the output file string. Returns the Hessian in atomic units.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(tuple(float))
    """

    comp_ptt = app.UNSIGNED_INTEGER
    mat = ar.matrix.read(
        output_str,
        val_ptt=app.FLOAT,
        start_ptt=(
            app.escape('Final Hessian.') +
            app.lpadded(app.NEWLINE)),
        block_start_ptt=(app.series(comp_ptt, app.LINESPACES) +
                         app.padded(app.NEWLINE)),
        line_start_ptt=comp_ptt,
        tril=False)

    if mat is not None:
        mat = tuple(map(tuple, mat))

    return mat


def harmonic_frequencies(output_str):
    """ Reads the harmonic vibrational frequencies from
        the output file string. Returns the frequencies in cm-1.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(float)
    """

    pattern = 'Frequency:' + app.capturing(app.LINE_FILL)
    captures = apf.all_captures(pattern, output_str)
    if captures is not None:
        freqs = ()
        for capture in captures:
            vals = capture.split()
            for val in vals:
                freqs += (float(val),)
    else:
        freqs = None
    return freqs


def _general_xyz_tensor(output_str):
    """ Reads the general tensor that could be used for lots of things

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(float)
    """

    init_ptt = (
        'FINAL TENSOR RESULT:' + app.NEWLINE +
        'Order ' + app.INTEGER + ', Length ' + app.INTEGER)
    atom_ptt = app.padded(
        'Atom' + app.SPACES + 'X' + app.SPACES + 'Y' + app.SPACES + 'Z')
    start_ptt = init_ptt + app.NEWLINE + atom_ptt + app.NEWLINE

    _tensor = ar.matrix.read(
        output_str,
        val_ptt=app.EXPONENTIAL_FLOAT,
        start_ptt=start_ptt,
        line_start_ptt=app.UNSIGNED_INTEGER)

    return _tensor
