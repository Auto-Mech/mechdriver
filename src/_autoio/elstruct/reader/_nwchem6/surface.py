""" gradient and hessian readers
"""

# import numpy
import autoread as ar
import autoparse.pattern as app


def gradient(output_str):
    """ Reads the molecular gradient (in Cartesian coordinates) from
        the output file string. Returns the gradient in atomic units.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(tuple(float))
    """

    grad = ar.matrix.read(
        output_str,
        start_ptt=app.padded(app.NEWLINE).join([
            app.padded(app.escape('ENERGY GRADIENTS'), app.NONNEWLINE),
            app.LINE, app.LINE, app.LINE, '']),
        line_start_ptt=app.LINESPACES.join([
            app.UNSIGNED_INTEGER,
            app.one_or_more(app.LETTER),
            app.FLOAT,
            app.FLOAT,
            app.FLOAT])
        )

    return grad

# def hessian(output_str):
#    """ Reads the molecular Hessian (in Cartesian coordinates) from
#        the output file string. Returns the Hessian in atomic units.
#
#        :param output_str: string of the program's output file
#        :type output_str: str
#        :rtype: tuple(tuple(float))
#    """
#     comp_ptt = app.UNSIGNED_INTEGER
#     mat = ar.matrix.read(
#         output_str,
#         start_ptt=app.padded(app.NEWLINE).join([
#             app.padded(app.escape(
#                 'MASS-WEIGHTED PROJECTED HESSIAN (Hartree/Bohr/Bohr/Kamu)'),
#                 app.NONNEWLINE),
#             app.LINE, app.LINE, app.LINE, '']),
#         block_start_ptt=(
#            app.series(comp_ptt, app.LINESPACES) +
#                       app.padded(app.NEWLINE) +
#            app.padded('----- ----- ----- ----- -----', app.LINESPACES) +
#                       app.padded(app.NEWLINE)),
#         line_start_ptt=comp_ptt,
#         val_ptt=app.EXPONENTIAL_FLOAT_D,
#         tril=True)
#
#     mat = tuple(map(tuple, mat))
#     return mat
