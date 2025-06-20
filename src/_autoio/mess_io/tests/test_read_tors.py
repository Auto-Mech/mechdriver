"""
 tests torsional mode reader
"""

import os
import numpy
from ioformat import pathtools
import mess_io.reader


PATH = os.path.dirname(os.path.realpath(__file__))
OUT_PATH = os.path.join(PATH, 'data', 'out')

OUT_STR = pathtools.read_file(OUT_PATH, 'freq.out')


def test__freqs():
    """ test mess_io.reader.tors.analytic_frequencies
        test mess_io.reader.tors.grid_minimum_frequencies
    """

    ref_analyt_freqs = (97.5241, 111.919, 97.5351, 255.494, 255.076)
    ref_gridmin_freqs = (100.299, 108.73, 100.312, 251.535, 251.531)

    analyt_freqs = mess_io.reader.tors.analytic_frequencies(OUT_STR)
    gridmin_freqs = mess_io.reader.tors.grid_minimum_frequencies(OUT_STR)

    assert numpy.allclose(analyt_freqs, ref_analyt_freqs)
    assert numpy.allclose(gridmin_freqs, ref_gridmin_freqs)


def test__zpves():
    """ test mess_io.reader.tors.zpves
    """

    ref_zpves = (0.0002253782706237488,
                 0.00025139700130794593,
                 0.00022540376824676267,
                 0.0005721284139963076,
                 0.0005719403690265805)

    zpves = mess_io.reader.tors.zero_point_vibrational_energies(OUT_STR)
    assert numpy.allclose(zpves, ref_zpves)

if __name__ == '__main__':
    test__freqs()
    test__zpves()