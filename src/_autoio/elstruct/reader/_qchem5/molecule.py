""" molecular geometry and structure readers
"""

from phydat import ptab
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

    geo = None

    # Read the block of text with the geometry
    opt_block_ptt = (
        app.escape('**  OPTIMIZATION CONVERGED  **') +
        app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
        'Z-matrix Print:'
    )
    opt_block = apf.last_capture(opt_block_ptt, output_str)

    # Read the geometry info
    if opt_block is not None:
        nums, xyzs = ar.geom.read(
            output_str,
            start_ptt=app.padded(app.NEWLINE).join([
                app.escape('Coordinates (Angstroms)'), app.LINE, '']),
            line_start_ptt=app.UNSIGNED_INTEGER)

    if all(x is not None for x in (nums, xyzs)):
        symbs = tuple(map(ptab.to_symbol, nums))
        geo = automol.geom.from_data(symbs, xyzs, angstrom=True)

    return geo


def opt_zmatrix(output_str):
    """ Reads the optimized Z-Matrix from the output file string.
        Returns the Z-Matrix in Bohr and Radians.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: automol Z-Matrix data structure
    """

    # Read the inital zmatrix that has the coordinates defined
    start_ptt = (
        app.escape('$molecule') + app.NEWLINE +
        # app.INTEGER + app.SPACES + app.INTEGER + app.NEWLINE +
        app.LINE_FILL + app.NEWLINE
    )
    symbs, key_mat, name_mat = ar.vmat.read(
        output_str,
        start_ptt=start_ptt,
        last=False)

    # Read the optimized values from the final z-matrix
    start_ptt = (
        'Z-matrix Print:' + app.NEWLINE +
        app.escape('$molecule') + app.NEWLINE +
        # app.INTEGER + app.SPACES + app.INTEGER + app.NEWLINE +
        app.LINE_FILL + app.NEWLINE
    )
    _, _, val_mat = ar.vmat.read(
        output_str,
        start_ptt=start_ptt,
        name_ptt=app.FLOAT)

    # Call the automol constructor
    if all(x is not None for x in (symbs, key_mat, name_mat, val_mat)):
        zma = automol.zmat.from_data(
            symbs, key_mat, val_mat, name_mat,
            one_indexed=True, angstrom=True, degree=True)
    else:
        zma = None

    return zma
