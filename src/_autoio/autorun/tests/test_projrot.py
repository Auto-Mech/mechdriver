""" tests projrot
"""

import os
import tempfile
import numpy
import automol
from ioformat import pathtools
import autorun


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')

GEO = automol.geom.from_string(
    pathtools.read_file(DAT_PATH, 'c2h5oh.xyz'))
GRAD = pathtools.read_numpy_file(DAT_PATH, 'c2h5oh.grad')
HESS = pathtools.read_numpy_file(DAT_PATH, 'c2h5oh.hess')

GEOS = (GEO,)
GRADS = (GRAD,)
HESSIANS = (HESS,)
PROJROT_ROT_STR = pathtools.read_file(DAT_PATH, 'c2h5oh.prot')


def test__frequencies():
    """ test autorun.projrot.frequencies
    """

    ref_rtproj_freqs = (
        269.39, 325.04, 429.05, 818.35, 913.01, 1084.71, 1116.71,
        1169.23, 1307.13, 1408.92, 1434.56, 1459.29, 1527.27, 1531.68,
        1555.98, 3020.12, 3062.32, 3121.51, 3142.58, 3160.32, 3841.83)
    ref_hrproj_freqs = (
        428.01, 815.69, 912.95, 1084.7, 1116.7, 1169.2, 1306.78,
        1408.87, 1434.54, 1459.25, 1527.25, 1531.57, 1555.97, 3020.11,
        3062.31, 3121.51, 3142.57, 3160.32, 3841.82)

    script_str = autorun.SCRIPT_DCT['projrot']

    with tempfile.TemporaryDirectory(dir=PATH) as run_dir:
        freq_inf = autorun.projrot.frequencies(
            script_str, run_dir, GEOS, GRADS, HESSIANS,
            rotors_str=PROJROT_ROT_STR)
        rtproj_freqs, hrproj_freqs, rt_imag_freq, hr_imag_freq = freq_inf

        assert numpy.allclose(rtproj_freqs, ref_rtproj_freqs)
        assert numpy.allclose(hrproj_freqs, ref_hrproj_freqs)
        assert not rt_imag_freq
        assert not hr_imag_freq
