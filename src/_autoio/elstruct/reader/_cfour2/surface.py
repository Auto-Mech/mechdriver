""" gradient and hessian readers
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

    # Grab a block of text containing the gradient
    block_ptt = ('Molecular gradient' +
                 app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
                 'Molecular gradient norm')
    block = apf.last_capture(block_ptt, output_str)

    if block is not None:
        # Trim the block to start it at the gradient lines
        blank_count = 0
        for i, line in enumerate(block.splitlines()):
            if line.strip() == '':
                blank_count += 1
                if blank_count == 3:
                    grad_start = i
                    break
        trim_block = '\n'.join(block.splitlines()[grad_start:])

        # Grab the gradient from the trimmed block string
        grad = ar.matrix.read(
            trim_block,
            line_start_ptt=app.LINESPACES.join([
                app.LETTER,
                app.escape('#') + app.UNSIGNED_INTEGER,
                app.maybe(app.UNSIGNED_INTEGER)]))
    else:
        grad = None

    return grad

# def hessian(output_str):
#     """ Reads the molecular Hessian (in Cartesian coordinates) from
#         the output file string. Returns the Hessian in atomic units.
#
#         :param output_str: string of the program's output file
#         :type output_str: str
#         :rtype: tuple(tuple(float))
#     """
#     try:
#         comp_ptt = app.one_of_these(['X', 'Y', 'Z']) + app.UNSIGNED_INTEGER
#         mat = ar.matrix.read(
#             output_str,
#             start_ptt=(app.escape('The second derivative matrix:') +
#                        app.lpadded(app.NEWLINE)),
#             block_start_ptt=(app.series(comp_ptt, app.LINESPACES) +
#                              app.padded(app.NEWLINE)),
#             line_start_ptt=comp_ptt,
#             tril=True)
#     except TypeError:
#         comp_ptt = app.UNSIGNED_INTEGER
#         mat = ar.matrix.read(
#             output_str,
#             val_ptt=app.EXPONENTIAL_FLOAT_D,
#             start_ptt=(
#                 app.escape('Force constants in Cartesian coordinates:') +
#                 app.lpadded(app.NEWLINE)),
#             block_start_ptt=(app.series(comp_ptt, app.LINESPACES) +
#                              app.padded(app.NEWLINE)),
#             line_start_ptt=comp_ptt,
#             tril=True)
#
#     mat = tuple(map(tuple, mat))
#     return mat
