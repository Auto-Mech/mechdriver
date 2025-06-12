""" molecular geometry and structure readers
"""

import autoread as ar
import autoparse.pattern as app
import autoparse.find as apf
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
            app.escape('Final (previous) structure:'), app.LINE, '']))

    geo = automol.geom.from_data(symbs, xyzs, angstrom=True)

    return geo


def opt_zmatrix(output_str):
    """ Reads the optimized Z-Matrix from the output file string.
        Returns the Z-Matrix in Bohr and Radians.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: automol molecular geometry data structure
    """

    # Read the matrix and the values from the output
    symbs, key_mat, name_mat, val_mat = ar.zmat.read(
        output_str,
        start_ptt=(
            app.padded(app.escape('Geometry (in Angstrom),'), app.NONNEWLINE) +
            2 * app.padded(app.NEWLINE)))

    # Call the automol constructor
    if all(x is not None for x in (symbs, key_mat, name_mat, val_mat)):
        zma = automol.zmat.from_data(
            symbs, key_mat, val_mat, name_mat,
            one_indexed=True, angstrom=True, degree=True)
    else:
        zma = None

    return zma


def inp_zmatrix(output_str):
    """ Reads the Z-Matrix specified in the input from the output file string.
        Returns the Z-Matrix in Bohr and Radians.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: automol molecular geometry data structure
    """

    # Block with init geom
    block_ptt = (
        'molecule' + app.SPACES + app.escape('{') +
        app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
        app.escape('}'))
    block = apf.first_capture(block_ptt, output_str)

    # Read the matrix and the values from the output
    start_ptt = app.LINE_FILL + app.NEWLINE + app.LINE_FILL + app.NEWLINE
    symbs, key_mat, name_mat, val_mat = ar.zmat.read(
        block, start_ptt=start_ptt)

    # Call the automol constructor
    if all(x is not None for x in (symbs, key_mat, name_mat, val_mat)):
        zma = automol.zmat.from_data(
            symbs, key_mat, val_mat, name_mat,
            one_indexed=True, angstrom=True, degree=True)
    else:
        zma = None

    return zma
