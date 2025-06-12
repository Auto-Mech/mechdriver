""" Tests the pdep functions in ratefit.fit
"""

from ratefit.fit import _fit as fit

KTP_DCT = {
    0.1:    ((500.0, 1000.0), (1.020, 1.050)),
    1.0:    ((300.0, 500.0, 1000.0), (1.0, 1.030, 1.060)),
    10.0:   ((500.0, 1000.0), (1.040, 2.080)),
    100.0:  ((500.0, 1300.0, 1000.0), (1.050, 1.50, 2.090)),
    'high': ((500.0, 1000.0), (1.060, 3.040))
}
TEMPS = (500.0,)
PVAL = 2.0


def test_assess_pdep():
    """ Tests the assess_pdep function
    """

    # Check if pdep using the default temps, 500 and 1000; should be True
    is_pdep1 = fit.assess_pdep(KTP_DCT)
    assert is_pdep1

    # Check if pdep using only 500 K; should be False
    is_pdep2 = fit.assess_pdep(KTP_DCT, assess_temps=TEMPS)
    assert not is_pdep2


def test_get_pdep_ktp_dct():
    """ Tests the get_pdep_ktp_dct function
    """

    ref1 = {
        0.1: ((500.0, 1000.0), (1.02, 1.05)),
        1.0: ((300.0, 500.0, 1000.0), (1.0, 1.03, 1.06)),
        10.0: ((500.0, 1000.0), (1.04, 2.08)),
        100.0: ((500.0, 1300.0, 1000.0), (1.05, 1.5, 2.09)),
        'high': ((500.0, 1000.0), (1.06, 3.04))}
    ref2 = {1.0: ((300.0, 500.0, 1000.0), (1.0, 1.03, 1.06))}
    ref3 = {100.0: ((500.0, 1300.0, 1000.0), (1.05, 1.5, 2.09))}

    print(fit.get_pdep_ktp_dct(KTP_DCT, assess_temps=TEMPS))
    assert fit.get_pdep_ktp_dct(KTP_DCT) == ref1
    assert fit.get_pdep_ktp_dct(KTP_DCT, assess_temps=TEMPS) == ref2
    assert fit.get_pdep_ktp_dct(KTP_DCT, assess_temps=TEMPS, pval=PVAL) == ref3


if __name__ == '__main__':
    test_assess_pdep()
    test_get_pdep_ktp_dct()
