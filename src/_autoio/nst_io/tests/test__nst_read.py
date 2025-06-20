""" test nst_io.reader
"""

import os
import numpy
import automol.geom
from ioformat import pathtools
import nst_io.reader


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')

OPT_OUT_STR = pathtools.read_file(DAT_PATH, 'opt.out')
ROT_OUT_STR = pathtools.read_file(DAT_PATH, 'rot.out')
FREQ_OUT_STR = pathtools.read_file(DAT_PATH, 'freq.out')


def test__msx_geometry():
    """ test nst_io.reader.optimized_msx_geometry
    """

    ref_opt_geo = (
        ('N', (0.48496662750727754, -0.5028027113553117, 0.0)),
        ('N', (2.472105383571757, 0.28846982341581795, 0.0)),
        ('O', (-2.617832759628675, 0.172984360844439, 0.0)))

    opt_geo = nst_io.reader.optimized_msx_geometry(OPT_OUT_STR)

    assert automol.geom.almost_equal_dist_matrix(ref_opt_geo, opt_geo)


def test__rotated_geometry():
    """ test nst_io.reader.rotated_geometry
    """

    ref_rot_geo = (
        ('N', (-0.5017302360246954, -0.28916123340201094, 0.0)),
        ('N', (-2.5910876099662348, 0.1751418972083358, 0.0)),
        ('O', (2.707672893605562, 0.09982064296437704, 0.0)))

    rot_geo = nst_io.reader.rotated_geometry(ROT_OUT_STR)

    assert automol.geom.almost_equal_dist_matrix(ref_rot_geo, rot_geo)


def test__frequencies():
    """ test nst_io.reader.msx_vibrational_frequencies
    """

    ref_freqs = (173.50368745416804, 2020.7934568716717)

    freqs = nst_io.reader.msx_vibrational_frequencies(FREQ_OUT_STR)

    assert numpy.allclose(ref_freqs, freqs)


if __name__ == '__main__':
    test__msx_geometry()
    test__rotated_geometry()
    test__frequencies()
