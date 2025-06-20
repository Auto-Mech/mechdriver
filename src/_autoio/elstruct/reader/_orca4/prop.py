""" property readers
"""

from autoparse import pattern as app
from autoparse import find as apf


def dipole_moment(output_str):
    """ Reads the xyz-coordinates of an SCF-computed static dipole moment
        from the output file string. Returns the dipole moment in Debye.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    pattern = app.LINESPACES.join([
        'Total Dipole Moment',
        ':',
        app.capturing(app.FLOAT),
        app.capturing(app.FLOAT),
        app.capturing(app.FLOAT)])
    captures = apf.last_capture(pattern, output_str)

    vals = captures if captures is not None else []
    if vals:
        vals = [float(val) for val in vals]
    else:
        vals = None

    return vals
#
#
# def polarizability(output_str):
#     """
#     Reads the static polarizability
#     """
#     tensor = ar.matrix.read(
#         output_str,
#         start_ptt=app.padded(
#                 'The raw cartesian tensor (atomic units)') +
#                 app.NEWLINE)
#     return tensor
