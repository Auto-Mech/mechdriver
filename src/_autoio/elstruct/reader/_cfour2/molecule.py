""" molecular geometry and structure readers
"""

import numbers
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
            app.escape('Coordinates (in bohr)'),
            app.LINE, app.LINE, '']),
        line_sep_ptt=app.UNSIGNED_INTEGER,)
    geo = automol.geom.from_data(symbs, xyzs, angstrom=False)

    return geo


def opt_zmatrix(output_str):
    """ Reads the optimized Z-Matrix from the output file string.
        Returns the Z-Matrix in Bohr and Radians.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: automol molecular geometry data structure
    """

    # complicated string patterns for the initial matrix read
    mat_ptt = app.padded(app.NEWLINE).join([
        app.LINESPACES.join([app.escape('Input from ZMAT file'),
                             app.escape('*')]),
        app.LINE, app.LINE, ''])

    nam_ptt = (app.LETTER +
               app.one_or_more(
                   app.one_of_these([app.LETTER, app.UNDERSCORE, app.DIGIT])) +
               app.maybe(app.escape('*')))

    # read the matrix from the beginning of the output
    symbs, key_mat, name_mat = ar.vmat.read(
        output_str,
        start_ptt=mat_ptt,
        symb_ptt=ar.par.Pattern.ATOM_SYMBOL + app.maybe(app.UNSIGNED_INTEGER),
        key_ptt=app.one_of_these([app.UNSIGNED_INTEGER, app.VARIABLE_NAME]),
        name_ptt=nam_ptt,
        last=False)

    # Remove any asterisks(*) from the entries in the name matrix
    if all(x is not None for x in (symbs, key_mat, name_mat)):
        name_mat = tuple(tuple(
            name.replace('*', '')
            if name is not None else None for name in name_row)
                         for name_row in name_mat)

        # complicated string patterns for the value dictionary read
        start_ptt = app.padded(app.NEWLINE).join(
            [app.padded('Final ZMATnew file', app.NONNEWLINE)] +
            [app.LINE for i in range(len(symbs)+3)] + [''])

        # read the values from the end of the output
        if len(symbs) == 1:
            # val_dct = {}
            val_mat = ((None, None, None),)
        else:
            val_dct = ar.setval.read(
                output_str,
                start_ptt=start_ptt,
                entry_sep_ptt='=',
                last=True)
            val_mat = ar.setval.convert_dct_to_matrix(val_dct, name_mat)

        # for the case when variable names are used instead of integer keys:
        # (otherwise, does nothing)
        key_dct = dict(map(reversed, enumerate(symbs)))
        key_dct[None] = 0
        key_mat = [
            [key_dct[val]+1 if not isinstance(val, numbers.Real) else val
             for val in row] for row in key_mat]
        symb_ptt = app.STRING_START + app.capturing(ar.par.Pattern.ATOM_SYMBOL)
        symbs = [apf.first_capture(symb_ptt, symb) for symb in symbs]

        # call the automol constructor
        zma = automol.zmat.from_data(
            symbs, key_mat, val_mat, name_mat,
            one_indexed=True, angstrom=True, degree=True)
    else:
        zma = None

    return zma
