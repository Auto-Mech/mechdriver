""" Tests the comments module
"""

import numpy
from chemkin_io.writer import comments

TEMPS = numpy.array([1000, 1500, 2000])
ERRS = numpy.array([10, 20, 30])
ERR_DCT = {'high': (TEMPS, ERRS),
           1.0:    (TEMPS, ERRS),
           10.0:   (TEMPS, ERRS)}
RXN = (('H', 'O2'), ('OH', 'O'), (None,))
RXN_ERR_DCT = {RXN: ERR_DCT}
REF_ERR_STR = (
    '! Fitting errors and ranges:\n'
    '! HPL: fit betw. 1000 and 2000 K, '
    'MeanAbsErr of 20.0%, MaxAbsErr of 30.0%\n'
    '! 1.00e+00 atm: fit betw. 1000 and 2000 K, '
    'MeanAbsErr of 20.0%, MaxAbsErr of 30.0%\n'
    '! 1.00e+01 atm: fit betw. 1000 and 2000 K, '
    'MeanAbsErr of 20.0%, MaxAbsErr of 30.0%\n')
SORT_DCT = {'cmts_top': 'This is a header comment',
            'cmts_inline': 'This is an inline comment'}
RXN_SORT_DCT = {RXN: SORT_DCT}


def test_get_err_str():
    """ Tests the get_err_str function
    """
    err_str = comments.get_err_str(ERR_DCT)
    assert err_str == REF_ERR_STR


def test_get_rxn_cmts_dct():
    """ Tests the get_rxn_cmts_dct function
    """
    # Test with only an err_dct
    rxn_cmts_dct = comments.get_rxn_cmts_dct(rxn_err_dct=RXN_ERR_DCT)
    assert rxn_cmts_dct[RXN]['header'] == ''
    assert rxn_cmts_dct[RXN]['inline'] == ''
    assert rxn_cmts_dct[RXN]['footer'] == REF_ERR_STR

    # Test with an err_dct and a sort_dct
    rxn_cmts_dct = comments.get_rxn_cmts_dct(rxn_err_dct=RXN_ERR_DCT,
                                             rxn_sort_dct=RXN_SORT_DCT)
    assert rxn_cmts_dct[RXN]['header'] == 'This is a header comment'
    assert rxn_cmts_dct[RXN]['inline'] == 'This is an inline comment'
    assert rxn_cmts_dct[RXN]['footer'] == REF_ERR_STR


if __name__ == '__main__':
    test_get_err_str()
    test_get_rxn_cmts_dct()
