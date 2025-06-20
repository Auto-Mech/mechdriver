""" Test the calculator/thermo.py functions
"""

import numpy as np
from mechanalyzer.calculator import thermo


TEMPS = np.array([1000, 1500, 2000])
BAD_TEMPS = np.array([1000, 1500, 7000])
SPC_NASA7_DCT = {
    'N2O': [
        'L 7/88', 'N   1O   1          ',
        'G',
        [200.0, 6000.0, 1000.0],
        ([0.48230729E+01, 0.26270251E-02, -0.95850872E-06,
          0.16000712E-09, -0.97752302E-14, 0.80734047E+04, -0.22017208E+01],
         [0.22571502E+01, 0.11304728E-01, -0.13671319E-04, 0.96819803E-08,
          -0.29307182E-11, 0.87417746E+04, 0.10757992E+02])
    ]
}

# Correct values
CORR_H = np.array([27678.841886015183, 34523.21839952253, 41721.44174812351])
CORR_CP = np.array([13.198655572491381, 14.104187903369416, 14.63922007474708])
CORR_S = np.array([66.20082629652916, 71.73872123759558, 75.87663200436498])
CORR_G = np.array([-38521.984410513964, -73084.86345687084, -110031.82226060645])


def test__valid_temps():
    """ Test the thermo calculator for temps within the valid range
    """
    spc_therm_dct = thermo.create_spc_therm_dct(SPC_NASA7_DCT, TEMPS)
    calc_h = spc_therm_dct['N2O'][1]
    calc_cp = spc_therm_dct['N2O'][2]
    calc_s = spc_therm_dct['N2O'][3]
    calc_g = spc_therm_dct['N2O'][4]
    assert np.allclose(calc_h, CORR_H, rtol=1e-3)
    assert np.allclose(calc_cp, CORR_CP, rtol=1e-3)
    assert np.allclose(calc_s, CORR_S, rtol=1e-3)
    assert np.allclose(calc_g, CORR_G, rtol=1e-3)


def test__invalid_temps():
    """ Test the thermo calculator for temps outside the valid range
    """
    spc_therm_dct = thermo.create_spc_therm_dct(SPC_NASA7_DCT, BAD_TEMPS)
    calc_h = spc_therm_dct['N2O'][1]
    calc_cp = spc_therm_dct['N2O'][2]
    calc_s = spc_therm_dct['N2O'][3]
    calc_g = spc_therm_dct['N2O'][4]
    assert np.isnan(calc_h[2])
    assert np.isnan(calc_cp[2])
    assert np.isnan(calc_s[2])
    assert np.isnan(calc_g[2])


if __name__ == '__main__':
    test__valid_temps()
    test__invalid_temps()
