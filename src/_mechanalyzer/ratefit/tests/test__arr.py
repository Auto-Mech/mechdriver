""" Test arrhenius fitter
"""

import numpy
from ratefit.fit import arr
from ratefit.fit import err

# Define stuff for testing single fitting
# KTS_1 and TEMPS_1 are from experimental data on N2O + Ar = N2 + O + Ar
KTS_1 = numpy.array([
    3.78e6, 5.22e6, 1.12e7, 1.77e7, 1.83e7,
    2.76e7, 3.50e7, 6.58e7, 1.33e8, 2.21e8,
    3.08e8, 6.73e8, 1.37e9, 2.44e9, 5.40e9])
TEMPS_1 = numpy.array([
    1546, 1576, 1640, 1681, 1688, 1729, 1754, 1821, 1894, 1954,
    1997, 2112, 2230, 2326, 2476])
KTP_DCT_1 = {'high': (TEMPS_1, KTS_1)}

# Define stuff for testing double fitting
# KTS_2 and TEMPS_2 are from the expressions given by Joshi and Wang for
# CO + OH = CO2 + H
TEMPS_2 = numpy.array([
    80, 85, 90, 95, 100, 105, 110, 115, 120, 125,
    130, 135, 140, 145, 150, 155, 160, 165, 170, 175,
    180, 185, 190, 195, 200, 215, 225, 240, 255, 270,
    285, 300, 400, 500, 600, 700, 800, 900, 1000, 1100,
    1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100,
    2200, 2300, 2400, 2500])
KTS_2 = numpy.array([
    4.21E+10, 4.55E+10, 4.87E+10, 5.17E+10, 5.45E+10,
    5.71E+10, 5.95E+10, 6.17E+10, 6.37E+10, 6.56E+10,
    6.74E+10, 6.90E+10, 7.04E+10, 7.18E+10, 7.31E+10,
    7.42E+10, 7.53E+10, 7.63E+10, 7.72E+10, 7.81E+10,
    7.89E+10, 7.96E+10, 8.03E+10, 8.09E+10, 8.15E+10,
    8.31E+10, 8.39E+10, 8.51E+10, 8.60E+10, 8.69E+10,
    8.77E+10, 8.84E+10, 9.29E+10, 9.88E+10, 1.07E+11,
    1.18E+11, 1.32E+11, 1.48E+11, 1.66E+11, 1.87E+11,
    2.11E+11, 2.37E+11, 2.65E+11, 2.96E+11, 3.29E+11,
    3.65E+11, 4.03E+11, 4.43E+11, 4.86E+11, 5.31E+11,
    5.78E+11, 6.28E+11, 6.80E+11, 7.34E+11])
KTP_DCT_2 = {'high': (TEMPS_2, KTS_2)}


def test_single():
    """ Test the double Arrhenius fitter
    """
    _, err_dct = arr.get_params(KTP_DCT_1)
    max_err = err.get_max_err(err_dct)
    assert max_err < 5  # %


def test_double():
    """ Test the double Arrhenius fitter
    """
    _, err_dct = arr.get_params(KTP_DCT_2)
    max_err = err.get_max_err(err_dct)
    assert max_err < 15  # %


if __name__ == '__main__':
    test_single()
    test_double()
