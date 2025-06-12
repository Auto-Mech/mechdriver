""" tests onedmine
"""

import os
import tempfile
import numpy
import autorun.onedmin


PATH = os.path.dirname(os.path.realpath(__file__))
TMP_DIR = tempfile.mkdtemp()
print('dir', TMP_DIR)


TGT_GEO = (
    ('C', (0.000, 0.00, 0.00)),
    ('H', (-1.2769423144257581, 1.613505423002926, -0.1610503348317155)),
    ('H', (-0.7457064171885416, -1.579904976024016, -1.0989739020896072)),
    ('H', (0.16053595786521754, -0.5552367534488531, 1.9813794625419021)),
    ('H', (1.8621125191304513, 0.521636700036555, -0.7213555059339084)))
BATH_GEO = (('He', (0.0, 0.0, 0.0)),)

EXP6_THY_INFO = ('exp6', None, None, None)
G09_THY_INFO = ('gaussian09', 'mp2', 'cc-pvdz', 'R')
CHARGE = 0
MULT = 1

NSAMP = 3
NJOBS = 3
SMIN = 3.779
SMAX = 9.448

SPIN_METHOD = 1
RANSEEDS = (153214316, 129539275, 174930586)


def test__exp6():
    """ test autorun.onedmin.lennard_jones_params
    """

    ref_lj_sig = (6.462050766831826, 6.503530255285625, 5.7414604006722465,
                  6.462050766831826, 6.503530255285625, 5.7414604006722465)
    ref_lj_eps = (19.05864, 18.71206, 27.96192, 19.05864, 18.71206, 27.96192)

    script_str = None
    lj_sig, lj_eps = autorun.onedmin.lennard_jones_params(
        script_str, TMP_DIR, NSAMP, NJOBS,
        TGT_GEO, BATH_GEO, EXP6_THY_INFO, CHARGE, MULT,
        smin=SMIN, smax=SMAX, spin_method=1, ranseeds=RANSEEDS)

    assert numpy.allclose(ref_lj_sig, lj_sig)
    assert numpy.allclose(ref_lj_eps, lj_eps)
