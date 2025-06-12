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

    symbs, xyzs = ar.geom.read(
        output_str,
        start_ptt=app.padded(app.NEWLINE).join([
            app.padded(
                app.escape('CARTESIAN COORDINATES (ANGSTROEM)'),
                app.NONNEWLINE),
            app.LINE, '']))

    geo = automol.geom.from_data(symbs, xyzs, angstrom=True)

    return geo
