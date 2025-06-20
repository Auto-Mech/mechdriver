""" property readers
"""

from autoparse import pattern as app
from autoparse import find as apf
import autoread as ar


def dipole_moment(output_str):
    """ Reads the xyz-coordinates of an SCF-computed static dipole moment
        from the output file string. Returns the dipole moment in Debye.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    ptt = (app.escape('Dipole Moment: [D]') + app.NEWLINE +
           app.padded('X:') + app.capturing(app.FLOAT) +
           app.padded('Y:') + app.capturing(app.FLOAT) +
           app.padded('Z:') + app.capturing(app.FLOAT) +
           app.padded('Total:') + app.FLOAT)
    captures = apf.last_capture(ptt, output_str)

    vals = captures if captures is not None else ()
    if vals:
        vals = tuple(float(val) for val in vals)
    else:
        vals = None

    return vals


def polarizability(output_str):
    """ Reads the xyz-components of a dipole polarizability tensor
        from the output file string. Returns the polarizability in _.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    start_ptt = app.padded(app.NEWLINE).join([
        (app.escape('Dipole Polarizability (Length Gauge)') + app.LINE_FILL),
        app.LINE, app.LINE, app.LINE, app.LINE, app.LINE, app.LINE, ''])

    polar = ar.matrix.read(
        output_str,
        start_ptt=start_ptt,
        line_start_ptt=app.INTEGER)

    return polar
