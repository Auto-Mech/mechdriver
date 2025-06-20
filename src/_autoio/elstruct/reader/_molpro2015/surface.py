""" gradient and hessian readers
"""

import numpy
from phydat import phycon
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

    head_ptt = ('Atom' + app.SPACES +
                app.escape('dE/dx') + app.SPACES +
                app.escape('dE/dy') + app.SPACES +
                app.escape('dE/dz'))
    grad = ar.matrix.read(
        output_str,
        start_ptt=app.padded(app.NEWLINE).join([
            app.padded(head_ptt, app.NONNEWLINE),
            app.LINE, '']),
        line_start_ptt=app.UNSIGNED_INTEGER)

    return grad


def hessian(output_str):
    """ Reads the molecular Hessian (in Cartesian coordinates) from
        the output file string. Returns the Hessian in atomic units.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(tuple(float))
    """

    comp_ptt = (
        app.one_or_more(app.LETTER) +
        app.one_of_these(['X', 'Y', 'Z']) +
        app.UNSIGNED_INTEGER
    )
    mat = ar.matrix.read(
        output_str,
        start_ptt=(
            app.escape('Force Constants (Second Derivatives of the Energy) ') +
            app.escape('in [a.u.]') +
            app.lpadded(app.NEWLINE)),
        block_start_ptt=(app.series(comp_ptt, app.LINESPACES) +
                         app.padded(app.NEWLINE)),
        line_start_ptt=comp_ptt,
        tril=True)

    if mat is not None:
        mat = tuple(map(tuple, mat))

    return mat


def harmonic_frequencies(output_str):
    """ Reads the harmonic vibrational frequencies from
        the output file string. Returns the frequencies in cm-1.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(float)
    """

    pattern = (app.escape('Wavenumbers [cm-1]') +
               app.capturing(app.LINE_FILL))
    captures = apf.all_captures(pattern, output_str)
    if captures is not None:
        freqs = []
        for capture in captures:
            vals = capture.split()
            for val in vals:
                freqs.append(float(val))
        freqs = tuple(freqs)
    else:
        freqs = None

    imag_block = apf.last_capture(
        (app.escape('Imaginary Vibration  Wavenumber') +
         app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
         app.escape('Vibration        Wavenumber')),
        output_str)
    if imag_block is not None:
        ptt = app.INTEGER + app.SPACES + app.capturing(app.FLOAT)
        imag_freqs = apf.all_captures(ptt, imag_block)
        if imag_freqs is not None:
            imag_freqs = tuple(-1*float(val) for val in imag_freqs)
            freqs = imag_freqs + freqs

    return freqs


def normal_coordinates(output_str):
    """ Reads the displacement along the normal modes (in Cartesian coordinates)
        from the output string. Returns the coordinates in Bohr.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(tuple(float))
    """

    block = apf.last_capture(
        ('Normal Modes' +
         app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
         'Frequencies dumped to record'),
        output_str)

    if block is not None:

        start_ptt = app.escape('Intensities [relative]') + app.LINE_FILL + '\n'
        line_start_ptt = app.one_or_more(app.NONSPACE)
        mats = ar.matrix.read_all(
            block,
            start_ptt=start_ptt,
            line_start_ptt=line_start_ptt
        )
        mats = tuple(numpy.array(mat) * phycon.ANG2BOHR for mat in mats)

        nmodes = ()
        for mat in mats:
            for column in mat.transpose():
                nrow = len(column) // 3
                ncol = 3

                nmode = numpy.zeros(shape=(nrow, ncol))
                for i in range(nrow):
                    for j in range(ncol):
                        nmode[i, j] = column[ncol*i+j]

                nmodes += (nmode,)

    else:
        nmodes = None

    return nmodes


if __name__ == '__main__':
    with open('run.out', encoding='utf-8') as fobj:
        OUTSTR = fobj.read()
    print(harmonic_frequencies(OUTSTR))
    NCS = normal_coordinates(OUTSTR)
    for x in NCS:
        print(x)
