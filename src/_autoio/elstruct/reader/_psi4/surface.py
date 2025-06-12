""" gradient and hessian readers
"""

import automol
import autoread as ar
import autoparse.pattern as app
import autoparse.find as apf


def gradient(output_str):
    """ Reads the molecular gradient (in Cartesian coordinates) from
        the output file string. Returns the gradient in atomic units.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(tuple(float))
    """

    comp_ptt = app.UNSIGNED_INTEGER

    start_ptts = (
        (app.padded(app.NEWLINE).join([
            app.escape('## Gradient (Symmetry 0) ##'),
            app.LINE, '', app.LINE, '', ''])),
        (app.padded(app.NEWLINE).join([
            app.escape('-Total gradient:'),
            app.LINE, app.LINE, ''])),
        (app.padded(app.NEWLINE).join([
            app.escape('-Total Gradient:'),
            app.LINE, app.LINE, '']))
    )

    for ptt in start_ptts:
        grad = ar.matrix.read(
            output_str,
            start_ptt=ptt,
            line_start_ptt=comp_ptt)
        if grad is not None:
            break

    return grad


def hessian(output_str):
    """ Reads the molecular Hessian (in Cartesian coordinates) from
        the output file string. Returns the Hessian in atomic units.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(tuple(float))
    """

    comp_ptt = app.UNSIGNED_INTEGER
    hess = ar.matrix.read(
        output_str,
        start_ptt=app.padded(app.NEWLINE).join([
            app.escape('## Hessian (Symmetry 0) ##'), app.LINE, '']),
        block_start_ptt=app.padded(app.NEWLINE).join([
            '', app.series(comp_ptt, app.LINESPACES), '', '']),
        line_start_ptt=comp_ptt)

    return hess


def harmonic_frequencies(output_str):
    """ Reads the harmonic vibrational frequencies from
        the output file string. Returns the frequencies in cm-1.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(float)
    """

    pattern = app.escape('Freq [cm^-1]') + app.capturing(app.LINE_FILL)

    captures = apf.all_captures(pattern, output_str)
    if captures is not None:
        freqs = ()
        for capture in captures:
            vals = capture.split()
            for val in vals:
                if 'i' not in val:
                    freqs += (float(val),)
                else:
                    freqs += (-1.0*float(val.replace('i', '')),)
    else:
        freqs = None

    return freqs


def irc_points(output_str):
    """ Reads the geometries, gradients, and Hessians at each point along the
        Intrinsic Reaction Coordinate from the output string.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: (geom data structure, tuple(tuple(float)), tuple(tuple(float)))
    """

    # Set pattern to find the end of each IRC optimization
    pattern = app.escape(
        '@IRC  **** Point ' +
        app.INTEGER +
        ' on IRC path is optimized ****'
    )
    block = apf.last_capture(
        (pattern +
         app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
         app.escape('    Back-transformation to cartesian coordinates...')),
        output_str)

    # Set pattern for grabbing the geometry from the block
    geo_head_ptt = (
        app.escape('@IRC    Cartesian Geometry (in Angstrom)') +
        app.LINE_FILL +
        '\n'
    )

    # Grab all of the optimized geometries
    captures = apf.all_captures(pattern, block)
    if captures is not None:
        geoms = []
        for string in captures:
            syms, xyzs = ar.geom.read(
                string,
                start_ptt=geo_head_ptt,
                line_start_ptt=app.escape('@IRC')
            )
            geoms.append(automol.geom.from_data(syms, xyzs, angstrom=True))
    else:
        geoms = []

    # Set the gradients and hessians to empty lists since they MAY not be run
    grads, hessians = [], []

    return geoms, grads, hessians


def irc_path(output_str):
    """ Reads the coordinates and electronic energies (relative to saddple point)
        of the Intrinsic Reaction Coordinate.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(automol geom data structure)
    """

    # coordinates
    block = apf.last_capture(
        (app.escape('@IRC              ****     IRC Steps     ****') +
         app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
         app.escape('---Fragment 1 Intrafragment Coordinates---')),
        output_str)

    pattern = (
        app.escape('@IRC') + app.SPACES + app.INTEGER + app.SPACES +
        app.FLOAT + app.SPACES +
        app.capturing(app.FLOAT) + app.SPACES +
        app.FLOAT + app.SPACES +
        app.LINE_FILL
    )

    captures = apf.all_captures(pattern, block)
    if captures is not None:
        # Remove duplicates that may appear because of Psi4 output printing
        unique_coords = []
        for coord in captures:
            if coord not in unique_coords:
                unique_coords.append(coord)
        coords = [float(coord) for coord in unique_coords]
    else:
        coords = None

    # energies
    block = apf.last_capture(
        (app.escape('@IRC            ****      IRC Report      ****') +
         app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
         app.escape('@IRC              ****     IRC Steps     ****')),
        output_str)

    pattern = (
        app.escape('@IRC') + app.SPACES + app.INTEGER + app.SPACES +
        app.capturing(app.FLOAT) + app.SPACES +
        app.FLOAT
    )

    captures = apf.all_captures(pattern, block)
    if captures is not None:
        energies = [float(capture) for capture in captures]
    else:
        energies = None

    return (coords, energies)
