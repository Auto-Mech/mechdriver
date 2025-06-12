""" test autorun._multiprog
"""

import os
import tempfile
import numpy
import automol
from ioformat import pathtools
import autorun


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')

NAME = 'Methane'
FORMULA = {'C': 1, 'H': 4}
HFORM0 = -0.02535174495
ENTHALPYT = 0.
BREAKT = 1000.
PF_STR = pathtools.read_file(DAT_PATH, 'pf.dat')

GEO = automol.geom.from_string(
    pathtools.read_file(DAT_PATH, 'c2h5oh.xyz'))
GRAD = pathtools.read_numpy_file(DAT_PATH, 'c2h5oh.grad')
HESS = pathtools.read_numpy_file(DAT_PATH, 'c2h5oh.hess')
PROJROT_ROT_STR = pathtools.read_file(DAT_PATH, 'c2h5oh.prot')
MESS_ROT_STR = pathtools.read_file(DAT_PATH, 'c2h5oh.mrot')


def test__thermo():
    """ test
    """

    ref_hform298 = -0.028389069688367432

    formula_str = automol.form.string(FORMULA)

    thermp_script_str = autorun.SCRIPT_DCT['thermp']
    pac99_script_str = autorun.SCRIPT_DCT['pac99'].format(formula_str)

    with tempfile.TemporaryDirectory(dir=PATH) as run_dir:
        hform298, nasa_poly = autorun.thermo(
            thermp_script_str, pac99_script_str, run_dir,
            PF_STR, NAME, FORMULA, HFORM0,
            enthalpyt=ENTHALPYT, breakt=BREAKT, convert=True)

    assert numpy.isclose(hform298, ref_hform298)
    assert nasa_poly == pathtools.read_file(DAT_PATH, 'ch4nasa.ckin')


def test__projected_frequencies():
    """ test for sadpt and species
        good to have a test with needing version 2
    """

    ref_proj_freqs = (
        428.01, 815.69, 912.95, 1084.7, 1116.7, 1169.2, 1306.78,
        1408.87, 1434.54, 1459.25, 1527.25, 1531.57, 1555.97,
        3020.11, 3062.31, 3121.51, 3142.57, 3160.32, 3841.82)
    ref_proj_zpe = 0.07996395706617601
    ref_harm_freqs = (
        269.39, 325.04, 429.05, 818.35, 913.01, 1084.71, 1116.71,
        1169.23, 1307.13, 1408.92, 1434.56, 1459.29, 1527.27,
        1531.68, 1555.98, 3020.12, 3062.32, 3121.51, 3142.58,
        3160.32, 3841.83)
    ref_tors_freqs = (264.0, 322.0)

    mess_script_str = autorun.SCRIPT_DCT['messpf']
    projrot_script_str = autorun.SCRIPT_DCT['projrot']

    with tempfile.TemporaryDirectory(dir=PATH) as run_dir:
        proj_inf = autorun.projected_frequencies(
            mess_script_str, projrot_script_str, run_dir,
            MESS_ROT_STR, PROJROT_ROT_STR,
            GEO, GEO, HESS,
            saddle=False)
        proj_freqs, proj_imag, proj_zpe, harm_freqs, tors_freqs = proj_inf

        assert numpy.allclose(proj_freqs, ref_proj_freqs)
        assert proj_imag is None
        assert numpy.isclose(proj_zpe, ref_proj_zpe)
        assert numpy.allclose(harm_freqs, ref_harm_freqs)
        assert numpy.allclose(tors_freqs, ref_tors_freqs)


if __name__ == '__main__':
    test__thermo()
    test__projected_frequencies()
