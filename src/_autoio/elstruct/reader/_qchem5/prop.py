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

    pattern = (app.escape('Dipole Moment (Debye)') +
               app.LINE_FILL + app.NEWLINE +
               app.padded('X') + app.capturing(app.FLOAT) +
               app.padded('Y') + app.capturing(app.FLOAT) +
               app.padded('Z') + app.capturing(app.FLOAT))
    captures = apf.last_capture(pattern, output_str)

    vals = captures if captures is not None else []
    vals = tuple(float(val) for val in vals) if vals else None

    return vals
