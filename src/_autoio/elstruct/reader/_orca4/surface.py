""" gradient and hessian readers
"""

import autoread as ar
import autoparse.pattern as app


def gradient(output_string):
    """ read gradient from the output string
    """

    grad = ar.matrix.read(
        output_string,
        start_ptt=app.padded(app.NEWLINE).join([
            app.padded(app.escape('CARTESIAN GRADIENT'), app.NONNEWLINE),
            app.LINE, app.LINE, '']),
        line_start_ptt=app.LINESPACES.join([
            app.UNSIGNED_INTEGER,
            app.one_or_more(app.LETTER),
            ':']))

    return grad


def hessian(output_string):
    """ read hessian from the output string
    """

    comp_ptt = app.UNSIGNED_INTEGER
    mat = ar.matrix.read(
        output_string,
        val_ptt=app.EXPONENTIAL_FLOAT,
        start_ptt=app.padded(app.NEWLINE).join([
            app.escape('$hessian'),
            app.LINE, '']),
        block_start_ptt=(app.series(comp_ptt, app.LINESPACES) +
                         app.padded(app.NEWLINE)),
        line_start_ptt=comp_ptt,
        tril=False)

    if mat is not None:
        mat = tuple(map(tuple, mat))

    return mat
