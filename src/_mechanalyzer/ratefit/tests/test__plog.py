""" Run tests
"""

import numpy
from ratefit.fit import plog
from ratefit.fit import err


TEMPS = numpy.linspace(400, 900, 11)
# used with KTS1 since no 400 K val
TEMPS_SHORT = numpy.linspace(450, 900, 10)

# Define rates at 0.01, 0.1, 1, 10, and 100 atm
KTS1 = numpy.array([5.45E-05, 0.000126326, 0.000302054, 0.00107086,
                    0.00553726, 0.0416265, 0.369654, 3.08305,
                    20.4691, 130.603])  # leaving off 400 K value
KTS2 = numpy.array([0.00373048, 0.0124119, 0.0307059, 0.066723,
                    0.162508, 0.489855, 2.14679, 13.2429,
                    87.8471, 462.577, 2226.93])
KTS3 = numpy.array([0.508358, 2.4336, 7.23671, 16.2493,
                    33.3004, 68.2562, 161.142, 479.797,
                    1768.32, 6495.67, 23780])
KTS4 = numpy.array([12.9012, 101.173, 430.225, 1199.49,
                    2600.21, 4880.51, 8806.79, 16446.2,
                    33401.4, 73633.2, 186461])
KTS5 = numpy.array([51.9327, 655.99, 4402.47, 18140.5,
                    51929.1, 114477, 213629, 362666,
                    589511, 946634, 1.84E+06])
KTP_DCT = {
    0.01:  (TEMPS_SHORT, KTS1),
    #0.1:   (TEMPS, KTS2), #gives 60% error
    1.0:   (TEMPS, KTS3),
    10.0:  (TEMPS, KTS4),
    100.0: (TEMPS, KTS5)}

# MORE COMPLEX FITS; CHECK FIT ACCURACY WITH SCIPY.OPTIMIZE
KTP_DCT_003 = {0.03: (numpy.array([ 700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700,
       1800, 1900, 2000]), numpy.array([1.70839926e+01, 3.59052809e+02, 3.29199792e+03, 1.65649708e+04,
       5.25326997e+04, 1.16776837e+05, 1.97924251e+05, 2.73040784e+05,
       3.22394912e+05, 3.38361666e+05, 3.23964108e+05, 2.87644090e+05,
       2.39691713e+05, 1.89607419e+05])), }
KTP_DCT_01 = {0.1: (numpy.array([ 700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700,
       1800, 1900, 2000]), numpy.array([1.17042854e+01, 3.62340701e+02, 3.96695222e+03, 2.17418827e+04,
       7.28404669e+04, 1.71444252e+05, 3.12533539e+05, 4.71391229e+05,
       6.15833072e+05, 7.19922497e+05, 7.70306454e+05, 7.65803355e+05,
       7.15134504e+05, 6.33489714e+05])), }
KTP_DCT_1 = {1.0: (numpy.array([ 700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700,
       1800, 1900, 2000]), numpy.array([6.09370479e-01, 4.82578824e+01, 1.05773202e+03, 9.79865455e+03,
       5.00060751e+04, 1.67098531e+05, 4.10455552e+05, 8.03791141e+05,
       1.33361425e+06, 1.96485747e+06, 2.65818345e+06, 3.36101673e+06,
       3.99605414e+06, 4.48473231e+06])),}
KTP_DCT_10 = {10.0: (numpy.array([ 700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700,
       1800, 1900, 2000]), numpy.array([4.26692513e-01, 1.93815729e+01, 3.37509413e+02, 3.05107943e+03,
       1.75693053e+04, 7.73320057e+04, 3.03780241e+05, 1.02326408e+06,
       2.56408368e+06, 4.67402312e+06, 6.80145399e+06, 8.93009298e+06,
       1.15236982e+07, 1.50188781e+07]))}


def test_plog():

    params, err_dct = plog.get_params(KTP_DCT, dbltol=15, dbl_iter=1)
    max_err = err.get_max_err(err_dct)
    assert max_err < 25  # error in %

def test_plog_dbltols():

    params, err_dct = plog.get_params(KTP_DCT_003, dbltol=1)
    max_err = err.get_max_err(err_dct)
    assert max_err < 0.6  # error in %
    params, err_dct = plog.get_params(KTP_DCT_01, dbltol=0.1)
    max_err = err.get_max_err(err_dct)
    assert max_err < 0.3  # error in %
    params, err_dct = plog.get_params(KTP_DCT_1, dbltol=0.1)
    max_err = err.get_max_err(err_dct)
    assert max_err < 0.5  # error in %
    params, err_dct = plog.get_params(KTP_DCT_10, dbltol=1)
    max_err = err.get_max_err(err_dct)
    assert max_err < 11.5  # error in %


if __name__ == '__main__':
    test_plog()
    test_plog_dbltols()

