""" tests fitter
"""

import numpy
from ratefit.fit import _fit as fit


# Define things for the assess_fit_method test
PDEP_KTP_DCT = {
    0.1: ((500.0, 1000.0), (1.020, 1.050)),
    1.0: ((300.0, 500.0, 1000.0), (1.0, 1.030, 1.060)),
    10.0: ((500.0, 1000.0), (1.040, 2.080)),
    100.0: ((500.0, 1300.0, 1000.0), (1.050, 1.50, 2.090)),
    'high': ((500.0, 1000.0), (1.060, 3.040))}
PDEP_KTP_DCT2 = {
    0.1: ((500.0, 1000.0), (1.020, 1.050)),
    1.0: ((500.0, 1000.0), (1.0, 1.030)),
    10.0: ((500.0, 1000.0), (1.040, 2.080)),
    100.0: ((500.0, 1000.0), (1.050, 2.090)),
    'high': ((500.0, 1000.0), (1.060, 3.040))}
PDEP_KTP_DCT3 = {  # designed to fail Chebyshev validity test
    0.1: ((600.0, 1000.0), (1.020, 1.050)),
    1.0: ((500.0, 1000.0), (1.0, 1.030)),
    10.0: ((500.0, 1000.0), (1.040, 2.080)),
    100.0: ((500.0, 1000.0), (1.050, 2.090)),
    'high': ((500.0, 1000.0), (1.060, 3.040))}
ASSESS_TEMPS = (500.0,)

# Define things for the Arrhenius test
ARR_KTS = numpy.array([
    3.78e6, 5.22e6, 1.12e7, 1.77e7, 1.83e7,
    2.76e7, 3.50e7, 6.58e7, 1.33e8, 2.21e8,
    3.08e8, 6.73e8, 1.37e9, 2.44e9, 5.40e9])
ARR_TEMPS = numpy.array([
    1546, 1576, 1640, 1681, 1688, 1729, 1754, 1821, 1894, 1954,
    1997, 2112, 2230, 2326, 2476])
ARR_KTP_DCT = {'high': (ARR_TEMPS, ARR_KTS)}

# Define things for the PLOG test
PLOG_TEMPS = numpy.linspace(400, 900, 11)
# used with KTS1 since no 400 K val
PLOG_TEMPS_SHORT = numpy.linspace(450, 900, 10)
# Define rates at 0.01, 0.1, 1, 10, and 100 atm
PLOG_KTS1 = numpy.array([5.45E-05, 0.000126326, 0.000302054, 0.00107086,
                         0.00553726, 0.0416265, 0.369654, 3.08305,
                         20.4691, 130.603])  # leaving off 400 K value
PLOG_KTS2 = numpy.array([0.00373048, 0.0124119, 0.0307059, 0.066723,
                         0.162508, 0.489855, 2.14679, 13.2429,
                         87.8471, 462.577, 2226.93])
PLOG_KTS3 = numpy.array([0.508358, 2.4336, 7.23671, 16.2493,
                         33.3004, 68.2562, 161.142, 479.797,
                         1768.32, 6495.67, 23780])
PLOG_KTP_DCT = {
    0.01: (PLOG_TEMPS_SHORT, PLOG_KTS1),
    0.1: (PLOG_TEMPS, PLOG_KTS2),
    1.0: (PLOG_TEMPS, PLOG_KTS3)}

CHEB_TEMPS = numpy.arange(300.0, 2500, 100.0)  # only goes up to 2400
CHEB_KTP_DCT = {
    'high': (CHEB_TEMPS, [0, 0]),  # should be ignored
    0.1: (CHEB_TEMPS, numpy.array(
        [1.08450103e-06, 4.02056587e-02, 2.67879037e+01, 2.18532255e+03,
         4.53398238e+04, 3.60139663e+05, 1.45018006e+06, 3.60988324e+06,
         6.38901400e+06, 8.87215267e+06, 1.03634259e+07, 1.06987289e+07,
         1.01130485e+07, 8.97902728e+06, 7.62804534e+06, 6.28500509e+06,
         5.07244048e+06, 4.03946844e+06, 3.19135182e+06, 2.51133072e+06,
         1.97421057e+06, 1.55376775e+06])),
    1.0: (CHEB_TEMPS, numpy.array(
        [1.05204297e-06, 4.37287193e-02, 2.68683407e+01, 2.33553171e+03,
         6.11167966e+04, 6.64960593e+05, 3.74422528e+06, 1.28835884e+07,
         3.07244837e+07, 5.57883099e+07, 8.26967252e+07, 1.05371941e+08,
         1.19913201e+08, 1.25396286e+08, 1.23093791e+08, 1.15262582e+08,
         1.04210676e+08, 9.18153940e+07, 7.93887558e+07, 6.77312296e+07,
         5.72544387e+07, 4.81070825e+07])),
    5.0: (CHEB_TEMPS, numpy.array(
        [1.03195059e-06, 4.51888681e-02, 2.58346404e+01, 2.19539776e+03,
         6.18976976e+04, 7.74446059e+05, 5.18193172e+06, 2.14400511e+07,
         6.14663642e+07, 1.33290864e+08, 2.33700127e+08, 3.48414202e+08,
         4.58824325e+08, 5.49334987e+08, 6.11221003e+08, 6.42712213e+08,
         6.46985110e+08, 6.29740256e+08, 5.97285483e+08, 5.55368956e+08,
         5.08649782e+08, 4.60589012e+08])),
    10.0: (CHEB_TEMPS, numpy.array(
        [1.02854867e-06, 4.54462673e-02, 2.55882776e+01, 2.14313125e+03,
         6.09581958e+04, 7.87229443e+05, 5.52751500e+06, 2.42486903e+07,
         7.41295252e+07, 1.71833810e+08, 3.22150011e+08, 5.13016160e+08,
         7.20276630e+08, 9.17230565e+08, 1.08265014e+09, 1.20439945e+09,
         1.27916208e+09, 1.31010944e+09, 1.30413070e+09, 1.26952728e+09,
         1.21446106e+09, 1.14609891e+09])),
    20.0: (CHEB_TEMPS, numpy.array(
        [1.03006701e-06, 4.54387780e-02, 2.56341148e+01, 2.12010613e+03,
         5.99626704e+04, 7.83171633e+05, 5.65725653e+06, 2.58820537e+07,
         8.33172377e+07, 2.04648252e+08, 4.08031125e+08, 6.92217829e+08,
         1.03569539e+09, 1.40463446e+09, 1.76350702e+09, 2.08321497e+09,
         2.34488271e+09, 2.53997245e+09, 2.66829042e+09, 2.73529370e+09,
         2.74957776e+09, 2.72093196e+09])),
    50.0: (CHEB_TEMPS, numpy.array(
        [1.04183916e-06, 4.49835129e-02, 2.63629986e+01, 2.16291378e+03,
         5.93540953e+04, 7.59794305e+05, 5.50077721e+06, 2.58129437e+07,
         8.69088124e+07, 2.26620973e+08, 4.84815568e+08, 8.88840854e+08,
         1.44354631e+09, 2.13010602e+09, 2.91213575e+09, 3.74488944e+09,
         4.58379302e+09, 5.39036895e+09, 6.13522323e+09, 6.79867258e+09,
         7.36986370e+09, 7.84514937e+09])),
    100.0: (CHEB_TEMPS, numpy.array(
        [1.06001201e-06, 4.42811270e-02, 2.76274082e+01, 2.27730457e+03,
         6.00191004e+04, 7.34814182e+05, 5.16987678e+06, 2.40869598e+07,
         8.21693105e+07, 2.20782995e+08, 4.93093215e+08, 9.52892469e+08,
         1.64235752e+09, 2.58362271e+09, 3.77615383e+09, 5.19915635e+09,
         6.81694453e+09, 8.58518980e+09, 1.04565997e+10, 1.23853058e+10,
         1.43297921e+10, 1.62545195e+10]))}

# Define things for the fit_rxn_ktp_dct test
RXN1 = (('H', 'O2'), ('OH', 'O'), (None,))
RXN2 = (('N', 'O2'), ('NO', 'O'), (None,))
RXN_KTP_DCT = {RXN1: ARR_KTP_DCT, RXN2: ARR_KTP_DCT}


def test_assess_fit_method():
    """ Tests the assess_fit_method function
    """

    # This dct is pressure dependent
    pdep_ktp_dct = fit.get_pdep_ktp_dct(PDEP_KTP_DCT)
    # Do some checks on the pressure-dependent dct
    fit_method1 = fit.assess_fit_method(pdep_ktp_dct, 'arr')
    assert fit_method1 == 'plog'
    fit_method2 = fit.assess_fit_method(pdep_ktp_dct, 'plog')
    assert fit_method2 == 'plog'
    fit_method3 = fit.assess_fit_method(pdep_ktp_dct, 'cheb')
    assert fit_method3 == 'plog'  # should switch to plog
    fit_method4 = fit.assess_fit_method(pdep_ktp_dct, 'troe')
    assert fit_method4 == 'troe'

    # This dct is pressure independent (due to the changed assess_temps)
    non_pdep_ktp_dct = fit.get_pdep_ktp_dct(
        PDEP_KTP_DCT, assess_temps=ASSESS_TEMPS)
    # Do some checks on the pressure-independent dct
    fit_method5 = fit.assess_fit_method(non_pdep_ktp_dct, 'arr')
    assert fit_method5 == 'arr'
    fit_method6 = fit.assess_fit_method(non_pdep_ktp_dct, 'plog')
    assert fit_method6 == 'arr'
    fit_method7 = fit.assess_fit_method(non_pdep_ktp_dct, 'cheb')
    assert fit_method7 == 'arr'
    fit_method8 = fit.assess_fit_method(non_pdep_ktp_dct, 'troe')
    assert fit_method8 == 'arr'

    # This dct is pressure dependent and also is valid for Chebyshev
    pdep_ktp_dct2 = fit.get_pdep_ktp_dct(PDEP_KTP_DCT2)

    # Do some checks on the Chebyshev-valid dct
    fit_method9 = fit.assess_fit_method(pdep_ktp_dct2, 'cheb')
    assert fit_method9 == 'cheb'

    # This dct is pressure dependent and is invalid for Chebyshev
    pdep_ktp_dct3 = fit.get_pdep_ktp_dct(PDEP_KTP_DCT3)

    # Do some checks on the Chebyshev-invalid dct
    fit_method10 = fit.assess_fit_method(pdep_ktp_dct3, 'cheb')
    assert fit_method10 == 'plog'


def test_fit_arr():
    """ Tests the fitting of a ktp_dct with the Arrhenius form
    """

    ref_arr_params = numpy.asarray((1651834009420615.0, -0.0582113, 59922.119))
    params, _ = fit.fit_ktp_dct(ARR_KTP_DCT, 'plog')  # for fun, specify plog
    assert numpy.allclose(ref_arr_params, params.arr)


def test_fit_plog():
    """ Tests the fitting of a ktp_dct with the PLOG form
    """

    ref_plog_dct = {
        0.01: ((2.6464e-233, 74.0887, -65453),),
        0.1: ((7.6153e-159, 51.2153, -41468),),
        1.0: ((2.5969e-81, 27.2407, -17616),)}

    arrfit_dct = {'dbltol': 80, 'dbl_iter': 30}  # 80 for only single fits
    # Just for fun, specify arr; it should switch to plog
    params, _ = fit.fit_ktp_dct(PLOG_KTP_DCT, 'arr', arrfit_dct=arrfit_dct)
    for pressure, ref_params in ref_plog_dct.items():
        fit_params = params.plog[pressure]
        assert numpy.allclose(ref_params[0], fit_params[0], rtol=1e-3)


def test_fit_cheb():
    """ Tests the fitting of a ktp_dct with the Chebyshev form
    """

    ref_alpha = numpy.array(
        [[1.86421309e+00, 4.20602838e-01, -5.74358452e-02, -5.45222311e-05],
         [7.61423648e+00, 7.51012552e-01, -9.23375204e-02, -8.25427040e-03],
         [-4.89211391e-01, 5.10360005e-01, -2.71105409e-02, -1.04075446e-02],
         [-3.93397030e-01, 2.67821927e-01, 1.58876205e-02, -6.32223880e-03],
         [-2.15290577e-01, 7.79192168e-02, 4.05605101e-02, 3.99721924e-03],
         [-8.40067735e-02, -2.24969601e-03, 2.50720909e-02, 5.36853083e-03]])
    ref_tmin = 300
    ref_tmax = 2400
    ref_pmin = 0.1
    ref_pmax = 100

    params, _ = fit.fit_ktp_dct(CHEB_KTP_DCT, 'cheb')

    assert numpy.allclose(params.cheb['alpha'], ref_alpha)
    assert numpy.allclose(params.cheb['tlim'], (ref_tmin, ref_tmax))
    assert numpy.allclose(params.cheb['plim'], (ref_pmin, ref_pmax))


def test_fit_rxn_ktp_dct():
    """ Tests the fitting of a rxn_ktp_dct
    """

    ref_arr_params = numpy.asarray((1651834009420615.0, -0.0582113, 59922.119))
    rxn_param_dct, _ = fit.fit_rxn_ktp_dct(RXN_KTP_DCT, 'arr')
    for params in rxn_param_dct.values():
        assert numpy.allclose(ref_arr_params, params.arr)


if __name__ == '__main__':
    test_assess_fit_method()
    test_fit_arr()
    test_fit_plog()
    test_fit_cheb()
    test_fit_rxn_ktp_dct()
