""" molecular geometry and structure readers
"""
import numbers
import numpy
import autoread as ar
import autoparse.pattern as app
import automol
import autoparse.find as apf


MOLPRO_ENTRY_START_PATTERN = (
        'SETTING' + app.SPACES +
        app.one_of_these(['R', 'A', 'D']) +
        app.one_or_more(app.INTEGER)
)
# 'SETTING' + app.not_followed_by(app.padded('MOLPRO_ENERGY')),
# 'SETTING' + app.not_followed_by('MOLPRO_ENERGY'),
# 'SETTING' + app.not_followed_by(app.padded('SPIN')),
# 'SETTING' + app.not_followed_by(app.padded('CHARGE'))


def opt_geometry(output_str):
    """ Reads the optimized molecular geometry (in Cartesian coordinates) from
        the output file string. Returns the geometry in Bohr.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: automol molecular geometry data structure
    """

    ptt = app.padded(app.NEWLINE).join([
        app.escape('Current geometry (xyz format, in Angstrom)'),
        '',
        app.UNSIGNED_INTEGER,
        (app.one_or_more(app.NONNEWLINE) + app.SPACES +
         'ENERGY=' + app.FLOAT),
        ''
    ])

    symbs, xyzs = ar.geom.read(
        output_str,
        start_ptt=ptt)
    geo = automol.geom.from_data(symbs, xyzs, angstrom=True)

    return geo


def hess_geometry(output_str):
    """ Reads the optimized molecular geometry (in Cartesian coordinates) from
        the output file string that is associated with a Hessian calculation
        so that the two are in the same coordinate system.
        Returns the geometry in Bohr.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: automol molecular geometry data structure
    """

    symbs, xyzs = ar.geom.read(
        output_str,
        start_ptt=app.padded(app.NEWLINE).join([
            app.escape('ATOMIC COORDINATES'),
            app.LINE, app.LINE, app.LINE, '']),
        line_start_ptt=app.UNSIGNED_INTEGER,
        line_sep_ptt=app.FLOAT,)
    geo = automol.geom.from_data(symbs, xyzs, angstrom=False)

    return geo


def opt_zmatrix(output_str):
    """ Reads the optimized Z-Matrix from the output file string.
        Returns the Z-Matrix in Bohr and Radians.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: automol molecular geometry data structure
    """

    # Reads the matrix from the beginning of the output
    # If geom found at beginning, no opt zmat should exist
    start_ptt = app.maybe(app.SPACES).join([
        'geometry', app.escape('='), app.escape('{'), ''])

    symbs, xyzs = ar.geom.read(output_str, start_ptt=start_ptt)
    if all(x is None for x in (symbs, xyzs)):
        symbs, key_mat, name_mat = ar.vmat.read(
            output_str,
            start_ptt=start_ptt,
            entry_start_ptt=app.maybe(','),
            last=False,
            case=False)
    else:
        symbs, key_mat, name_mat = None, None, None

    # Read the initial z-matrix values from the beginning out the output
    if all(x is not None for x in (symbs, key_mat, name_mat)):
        if len(symbs) == 1:
            val_dct = {}
        else:
            val_dct = ar.setval.read(
                output_str,
                # entry_start_ptt=MOLPRO_ENTRY_START_PATTERN,
                entry_start_ptt='SETTING',
                name_ptt=(
                    app.one_of_these(['R', 'A', 'D']) +
                    app.one_or_more(app.INTEGER)),
                val_ptt=app.NUMBER,
                last=True,
                case=False)

        names = sorted(set(numpy.ravel(name_mat)) - {None})
        caps_names = list(map(str.upper, names))
        name_dct = dict(zip(caps_names, names))
        assert set(caps_names) <= set(val_dct)
        val_dct = {name_dct[caps_name]: val_dct[caps_name]
                   for caps_name in caps_names}

        # Read optimized z-matrix values from the end of the output
        var_string = app.one_of_these([
            app.padded('Optimized variables'),
            app.padded('Current variables')
        ])
        opt_val_dct = ar.setval.read(
            output_str,
            start_ptt=var_string + app.NEWLINE,
            entry_end_ptt=app.one_of_these(['ANGSTROM', 'DEGREE']),
            last=True,
            case=False)
        opt_val_dct = {name_dct[name]: val_dct
                       for name, val_dct in opt_val_dct.items()}
        assert set(opt_val_dct) <= set(val_dct)
        val_dct.update(opt_val_dct)

        val_mat = ar.setval.convert_dct_to_matrix(val_dct, name_mat)

        # Call the automol constructor
        zma = automol.zmat.from_data(
            symbs, key_mat, val_mat, name_mat,
            one_indexed=True, angstrom=True, degree=True)
    else:
        zma = None

    return zma


def inp_zmatrix(inp_str):
    """ Reads the input z-matrix from the input file string
        Returns the Z-Matrix in Bohr and Radians.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: automol molecular geometry data structure
    """

     # Remove the constants line (not sure if needed, but retaining anyways)
    inp_str = inp_str.replace('  Constants:\n', '') 
    
    # Reads the matrix from the beginning of the input
    symbs, key_mat, name_mat = ar.vmat.read(
        inp_str,
        start_ptt='geometry' + app.padded(app.escape('=')) + app.escape('{') \
            + app.NEWLINE,
        symb_ptt=(ar.par.Pattern.ATOM_SYMBOL +
                  app.not_followed_by(app.SPACES + app.FLOAT) +
                  app.maybe(app.UNSIGNED_INTEGER)),
        key_ptt=app.one_of_these([app.UNSIGNED_INTEGER, app.VARIABLE_NAME]),
        line_end_ptt=app.maybe(app.UNSIGNED_INTEGER),
        last=False)

    # Reads the values from the input
    if all(x is not None for x in (symbs, key_mat, name_mat)):
        if len(symbs) == 1:
            # val_dct = {}
            val_mat = ((None, None, None),)
        else:
            val_dct = ar.setval.read(
                inp_str,
                start_ptt='geometry' + app.padded(app.escape('=')) + \
                    app.escape('{') + app.NEWLINE + \
                    app.one_or_more(app.one_or_more(app.WILDCARD2) + \
                    app.NEWLINE) + app.rpadded(app.escape('}')) + \
                    app.NEWLINE,
                entry_sep_ptt=app.maybe(app.escape('=')),  # either '=' or ''
                entry_start_ptt='',
                sep_ptt=app.maybe(app.LINESPACES).join([
                    app.NEWLINE]),
                last=True)
            val_mat = ar.setval.convert_dct_to_matrix(val_dct, name_mat)

        # Check for the pattern
        # For the case when variable names are used instead of integer keys:
        # (otherwise, does nothing)
        key_dct = dict(map(reversed, enumerate(symbs)))
        key_dct[None] = 0 
        key_mat = [ 
            [key_dct[val]+1 if not isinstance(val, numbers.Real) else val 
             for val in row] for row in key_mat]
        symb_ptt = app.STRING_START + app.capturing(ar.par.Pattern.ATOM_SYMBOL)
        symbs = [apf.first_capture(symb_ptt, symb) for symb in symbs]

        # Call the automol constructor
        zma = automol.zmat.from_data(
            symbs, key_mat, val_mat, name_mat,
            one_indexed=True, angstrom=True, degree=True)
    else:
        zma = None

    return zma

