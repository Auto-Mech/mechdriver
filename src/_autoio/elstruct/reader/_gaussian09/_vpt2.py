""" Obtain information from VPT2 calculations
"""

from phydat import phycon
import automol
import autoread as ar
import autoparse.pattern as app
import autoparse.find as apf


def anharmonic_frequencies(output_str):
    """ Reads the anharmonic vibrational frequencies from
        the output file string. Returns the frequencies in cm-1.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(float)
    """

    block = apf.last_capture(
        (app.escape('Fundamental Bands (DE w.r.t. Ground State)') +
         app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
         app.escape('Overtones (DE w.r.t. Ground State)')), output_str)

    pattern = (
        app.INTEGER +
        app.escape('(1') + app.maybe(app.escape(',')) +
        app.maybe(app.SIGN) +
        app.maybe(app.INTEGER) + app.escape(')') +
        app.SPACE +
        app.maybe(app.one_or_more(app.LOWERCASE_LETTER)) +
        app.one_or_more(app.SPACE) +
        app.FLOAT +
        app.one_or_more(app.SPACE) +
        app.capturing(app.FLOAT)
    )

    # Get list of anharmonic frequency values

    anharm_freq = [float(val)
                   for val in apf.all_captures(pattern, block)]

    return sorted(anharm_freq)


def anharmonic_zpve(output_str):
    """ Reads the VPT2-computed anharmonic zero-point vibrational energy from
        the output file string. Returns the energy in Hartrees.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    anharm_zpve_pattern = (
        'Total Anharm' +
        app.one_or_more(app.SPACE) +
        ':' +
        app.SPACE +
        'cm-1' +
        app.SPACE +
        '=' +
        app.one_or_more(app.SPACE) +
        app.FLOAT +
        app.SPACE +
        ';' +
        app.SPACE +
        'Kcal/mol' +
        app.SPACE +
        '=' +
        app.one_or_more(app.SPACE) +
        app.FLOAT +
        app.SPACE +
        ';' +
        app.SPACE +
        'KJ/mol' +
        app.SPACE +
        '=' +
        app.one_or_more(app.SPACE) +
        app.capturing(app.FLOAT)
    )

    # Retrieve the anharm ZPVE
    anh_zpve = apf.last_capture(anharm_zpve_pattern, output_str)

    # Convert the ZPVE units
    anh_zpve = float(anh_zpve.replace('D', 'E'))
    anh_zpve *= phycon.KJ2EH

    return anh_zpve


def anharmonicity_matrix(output_str):
    """ Reads the VPT2-computed anharmonicity matrix from the output file string.
        Returns the matrix in _.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(tuple(float))
    """

    start_string = 'Total Anharmonic X Matrix (in cm^-1)'
    comp_ptt = app.UNSIGNED_INTEGER
    mat = ar.matrix.read(
        output_str,
        val_ptt=app.EXPONENTIAL_FLOAT_D,
        start_ptt=app.padded(app.NEWLINE).join([
            app.padded(app.escape(start_string), app.NONNEWLINE),
            app.LINE, '']),
        block_start_ptt=(app.series(comp_ptt, app.LINESPACES) +
                         app.padded(app.NEWLINE)),
        line_start_ptt=comp_ptt,
        tril=True)

    mat = tuple(tuple(float(val) for val in row)
                for row in mat)

    return mat


def vibrorot_alpha_matrix(output_str):
    """ Reads the VPT2-computed Vibrational-Rotational Alpha matrix
        from the output file string. Returns the matrix in _.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(tuple(float))
    """

    begin_string = 'Vibro-Rot alpha Matrix (in cm^-1)'
    end_string = app.escape('Q(') + app.maybe(app.SPACE) + app.UNSIGNED_INTEGER + app.escape(')')

    vib_rot_mat = ar.matrix.read(
        output_str,
        start_ptt=app.padded(app.NEWLINE).join([
            app.padded(app.escape(begin_string), app.NONNEWLINE),
            app.LINE, app.LINE, '']),
        line_start_ptt=end_string)

    return vib_rot_mat


def centrifugal_distortion_constants(output_str):
    """ Reads the VPT2-computed quartic centrifugal distortion constants
        from the output file string. Returns the constants in _.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(tuple(float))
    """

    # Set patterns for all molecule types and symmetries
    block = apf.last_capture(
        ('Quartic Centrifugal Distortion Constants Tau Prime' +
         app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
         'Asymmetric Top Reduction'),
        output_str)
    if not block:
        block = apf.last_capture(
            ('Quartic Centrifugal Distortion Constants Tau Prime' +
             app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
             'Constants in the Symmetrically Reduced Hamiltonian'),
            output_str)
    if not block:
        block = apf.last_capture(
            ('Quartic Centrifugal Distortion Constants Tau Prime' +
             app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
             'Rotational l-type doubling constants'),
            output_str)

    # Read values
    pattern = (
        'TauP' +
        app.SPACE +
        app.capturing(app.one_or_more(app.LOWERCASE_LETTER)) +
        app.SPACES +
        app.capturing(app.EXPONENTIAL_FLOAT_D) +
        app.SPACES +
        app.EXPONENTIAL_FLOAT_D
    )

    cent_dist_const = [[lbl, float(val.replace('D', 'E'))]
                       for (lbl, val) in apf.all_captures(pattern, block)]

    return cent_dist_const


def cubic_force_constants(output_str):
    """ Reads the cubic force constants
        from the output file string. Returns the constants in _.
        Hartree*amu(-3/2)*Bohr(-3)

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(tuple(float))
    """

    block = apf.last_capture(
        ('CUBIC FORCE CONSTANTS IN NORMAL MODES' +
         app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
         'QUARTIC FORCE CONSTANTS IN NORMAL MODES'),
        output_str)

    pattern = (
        app.capturing(app.INTEGER) +
        app.SPACES +
        app.capturing(app.INTEGER) +
        app.SPACES +
        app.capturing(app.INTEGER) +
        app.SPACES +
        app.FLOAT +
        app.SPACES +
        app.FLOAT +
        app.SPACES +
        app.capturing(app.FLOAT)
    )

    caps = apf.all_captures(pattern, block)
    if caps:
        cfc_mat = _fc_mat(caps)
    else:
        cfc_mat = None

    return cfc_mat


def quartic_force_constants(output_str):
    """ Reads the quartic force constants
        from the output file string. Returns the constants in _.
        Hartree*amu(2)*Bohr(-4)

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(tuple(float))
    """

    block = apf.last_capture(
        ('QUARTIC FORCE CONSTANTS IN NORMAL MODES' +
         app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
         'Input to Restart Anharmonic Calculations'),
        output_str)
    if block is None:
        block = apf.last_capture(
            ('QUARTIC FORCE CONSTANTS IN NORMAL MODES' +
             app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
             'Input for POLYMODE'),
            output_str)

    pattern = (
        app.capturing(app.INTEGER) +
        app.SPACES +
        app.capturing(app.INTEGER) +
        app.SPACES +
        app.capturing(app.INTEGER) +
        app.SPACES +
        app.capturing(app.INTEGER) +
        app.SPACES +
        app.FLOAT +
        app.SPACES +
        app.FLOAT +
        app.SPACES +
        app.capturing(app.FLOAT)
    )

    caps = apf.all_captures(pattern, block)
    if caps:
        qfc_mat = _fc_mat(caps)
    else:
        qfc_mat = None

    return qfc_mat


def _fc_mat(fc_caps):
    """ caps: (
            ((idx1, idx2, ..., idxn), val1),
            ((idx1, idx2, ..., idxn), val2),
            ...,
            ((idx1, idx2, ..., idxn), valn),
    """

    # Convert the types of the force constant data
    fc_idxs, fc_vals = [], []
    for caps in fc_caps:
        fc_idxs.append(tuple(int(val)-1 for val in caps[:-1]))
        fc_vals.append(float(caps[-1]))

    fc_mat = automol.util.tensor.build_full_array(
        fc_idxs, fc_vals, fill_perms=True)

    return fc_mat


def vpt2(output_str):
    """ Reads out various pieces data from the output string of a
        VPT2 calculations.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: dict[str: obj]
    """

    anharm_dct = {
        'freqs': anharmonic_frequencies(output_str),
        'zpve': anharmonic_zpve(output_str),
        'x_mat': anharmonicity_matrix(output_str),
        'vibrot_mat': vibrorot_alpha_matrix(output_str),
        'cent_dist_const': centrifugal_distortion_constants(output_str),
        'cubic_fc': cubic_force_constants(output_str),
        'quartic_fc': quartic_force_constants(output_str)
    }

    return anharm_dct
