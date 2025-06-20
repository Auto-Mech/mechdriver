""" test intder_io.reader
"""

import os
import numpy
from ioformat import pathtools
import intder_io.reader


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')

OUT_STR = pathtools.read_file(DAT_PATH, 'intder.out')


def test__coords():
    """ test intder_io.reader.internal_coordinates
    """

    ref_intl_coords = (
        ('STRE', (1, 0)),
        ('STRE', (2, 0)),
        ('STRE', (3, 0)),
        ('STRE', (4, 0)),
        ('STRE', (5, 1)),
        ('STRE', (6, 1)),
        ('STRE', (7, 1)),
        ('STRE', (8, 5)),
        ('BEND', (2, 0, 1)),
        ('BEND', (3, 0, 1)),
        ('BEND', (4, 0, 1)),
        ('BEND', (5, 1, 0)),
        ('BEND', (6, 1, 0)),
        ('BEND', (7, 1, 0)),
        ('BEND', (8, 5, 1)),
        ('TORS', (3, 0, 1, 2)),
        ('TORS', (4, 0, 1, 2)),
        ('TORS', (5, 1, 0, 2)),
        ('TORS', (6, 1, 0, 5)),
        ('TORS', (7, 1, 0, 5)),
        ('TORS', (8, 5, 1, 0))
    )

    intl_coords = intder_io.reader.internal_coordinates(OUT_STR)

    assert ref_intl_coords == intl_coords


def test__ted_assignments():
    """ test intder_io.reader.ted_assignments
    """

    ref_ted = (
        (3.7, {4: 68.9, 2: 5.5, -12: 5.3, -6: 5.0}),
        (162.1, {0: 59.4, 11: 11.3, 10: 9.9, -20: 5.4}),
        (279.0, {2: 28.6, 1: 19.9, -10: 19.8, -4: 12.9}),
        (459.4, {20: 36.9, -7: 33.5, -17: 13.1, -5: 4.6}),
        (715.2, {17: 45.6, 11: 19.9, 15: 8.7, -19: 7.5}),
        (858.6, {3: 30.4, -12: 17.0, -6: 14.2, 11: 7.8}),
        (1086.6, {11: 42.2, 0: 20.6, -3: 6.9, -17: 6.6}),
        (1279.2, {3: 17.6, 6: 14.1, -19: 12.7, 0: 10.2}),
        (1319.7, {19: 26.3, 5: 21.2, -13: 11.3, 11: 9.5}),
        (1410.1, {15: 40.8, -5: 15.6, -17: -10.1, 10: 9.9}),
        (1502.6, {2: 34.5, -3: 23.6, -7: 9.4, -1: 9.0}),
        (1520.8, {6: 38.8, -12: 20.8, -5: 9.6, 19: 9.0}),
        (1541.5, {10: 34.2, 1: 19.2, -12: 14.7, -8: 6.7}),
        (1631.5, {7: 25.0, -1: 19.7, -15: 14.1, 20: 13.0}),
        (1830.7, {15: 15.6, 18: 15.5, 5: 14.4, -17: 12.3}),
        (2673.9, {9: 48.6, -18: 18.0, 15: 13.8, 5: 10.1}),
        (2909.4, {8: 29.4, -18: 24.2, 5: 12.6, -15: 12.3}),
        (3084.4, {9: 27.3, 18: 22.2, 8: 16.9, 17: 10.7}),
        (3145.5, {13: 38.2, -16: 23.8, 19: 16.7, 8: 10.8}),
        (3188.0, {16: 61.1, 13: 22.5, 19: 9.3, -17: -4.7}),
        (3758.4, {14: 83.5, 20: 12.0})
    )

    ted = intder_io.reader.ted_assignments(OUT_STR)

    for rmode, mode in zip(ref_ted, ted):
        assert numpy.isclose(rmode[0], mode[0])
        rdct, dct = rmode[1], mode[1]
        assert set(rdct.keys()) == set(dct.keys())
        assert numpy.allclose(tuple(rdct.values()), tuple(dct.values()))
