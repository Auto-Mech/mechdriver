""" molecular geometry and structure readers
"""

import autoread as ar
import autoparse.pattern as app
import automol


def opt_geometry(output_str):
    """ Reads the optimized molecular geometry (in Cartesian coordinates) from
        the output file string. Returns the geometry in Bohr.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: automol molecular geometry data structure
    """

    opt_str = ('Output coordinates in a.u. ' +
               'scale by  1.000000000 to convert to a.u.)')

    symbs, xyzs = ar.geom.read(
        output_str,
        start_ptt=app.padded(app.NEWLINE).join([
            app.escape(opt_str),
            app.LINE, app.LINE, app.LINE, '']),
        line_start_ptt=app.UNSIGNED_INTEGER,
        line_sep_ptt=app.FLOAT,)
    geo = automol.geom.from_data(symbs, xyzs, angstrom=False)

    return geo
