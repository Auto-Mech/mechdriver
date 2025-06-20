""" reaction read write testing
"""

import numpy as np
from autoreact.params import RxnParams
from chemkin_io.writer.mechanism import reactions_block as writer


# Define miscellaneous stuff
RXN1 = (('H', 'O2'), ('OH', 'O'), (None,))
RXN2 = (('H', 'O2'), ('OH', 'O'), ('(+N2)',))
CMTS_DCT = {'header': '! Comment\n',
            'inline': '! Comment',
            'footer': '! Comment\n'}
RXN_CMTS_DCT1 = {RXN1: CMTS_DCT, 'block': '! Block comment\n\n'}
RXN_CMTS_DCT2 = {RXN2: CMTS_DCT}

# Define stuff for testing Arrhenius
ARR_DCT = {'arr_tuples': [[1E+15, 0.00, 25000]]}
ARR_PARAMS = RxnParams(arr_dct=ARR_DCT)
ARR_RXN_PARAM_DCT = {RXN1: ARR_PARAMS}
DUP_ARR_DCT = {'arr_tuples': [[1E+15, 0.00, 25000], [1E+15, 0.00, 25000]]}
DUP_ARR_PARAMS = RxnParams(arr_dct=DUP_ARR_DCT)
DUP_ARR_RXN_PARAM_DCT = {RXN1: DUP_ARR_PARAMS}
DUP_ARR_DCT2 = {'arr_tuples': [[1E+15, 0.00, 27000], [1E+15, 0.00, 25000]]}
DUP_ARR_PARAMS2 = RxnParams(arr_dct=DUP_ARR_DCT2)
DUP_ARR_RXN_PARAM_DCT2 = {RXN1: DUP_ARR_PARAMS2}
COLLID = {'N2': 1.4, 'AR': 1.0}
COLLID_ARR_DCT = {'arr_tuples': [[1E+15, 0.00, 25000]], 'arr_collid': COLLID}
COLLID_ARR_PARAMS = RxnParams(arr_dct=COLLID_ARR_DCT)
COLLID_ARR_RXN_PARAM_DCT = {RXN1: COLLID_ARR_PARAMS}

# rxns to test "=>" sign if 1) both fw and bw rxns present, 2) >2 products
RXNA = (('H', 'O2'), ('OH', 'O'), (None,))
RXNB = (('OH', 'O'), ('H', 'O2'), (None,))
RXNC = (('H', 'O2'), ('O', 'O', 'H'), (None,))
RXN_ALL_DCT = {RXNA: ARR_PARAMS, RXNB: ARR_PARAMS, RXNC: ARR_PARAMS}

# Define stuff for testing PLOG
PLOG_DCT = {
    0.1: [[1E+15, 0.00, 25000]],
    1.0: [[1E+15, 0.00, 25000], [1E+15, 0.00, 25000]],
    100.0: [[1E+15, 0.00, 25000]]}
PLOG_PARAMS = RxnParams(plog_dct=PLOG_DCT)
PLOG_RXN_PARAM_DCT = {RXN1: PLOG_PARAMS}

# Define stuff for testing Chebyshev
CHEB_DCT = {
    'tlim': [500.0, 2000.0],
    'plim': [0.03, 100.0],
    'alpha': np.array([
        [-1.620e+01, -1.183e-01, -5.423e-02, -1.476e-02],
        [2.578e+00, 1.614e-01, 7.359e-02, 1.880e-02],
        [1.068e-01, -7.235e-02, -2.733e-02, -3.778e-03],
        [3.955e-02, 1.207e-02, 3.402e-04, -2.695e-03],
        [8.557e-03, 4.345e-03, 3.670e-03, 1.608e-03],
        [8.599e-04, -1.758e-03, -7.502e-04, 7.396e-07]])}
CHEB_PARAMS = RxnParams(cheb_dct=CHEB_DCT)
CHEB_RXN_PARAM_DCT = {RXN2: CHEB_PARAMS}

# Define stuff for testing Troe
TROE_DCT = {'highp_arr': [[1e12, 1.5, 50000]],
            'lowp_arr': [[1e12, 1.5, 50000]],
            'troe_params': [1.5, 8000, 100, 1000],
            'collid': {'AR': 1.4, 'N2': 1.7}}
TROE_PARAMS = RxnParams(troe_dct=TROE_DCT)
TROE_RXN_PARAM_DCT = {RXN2: TROE_PARAMS}

# Define stuff for testing Lindemann
LIND_DCT = {'highp_arr': [[1e12, 1.5, 50000]],
            'lowp_arr': [[1e12, 1.5, 50000]],
            'collid': {'AR': 1.4, 'N2': 1.7}}
LIND_PARAMS = RxnParams(lind_dct=LIND_DCT)
LIND_RXN_PARAM_DCT = {RXN2: LIND_PARAMS}

# Define stuff for testing duplicates of two PLOGs
DUP_PLOG_PARAMS = RxnParams(plog_dct=PLOG_DCT)
DUP_PLOG_PARAMS.combine_objects(DUP_PLOG_PARAMS)  # combine with itself
DUP_PLOG_RXN_PARAM_DCT = {RXN2: DUP_PLOG_PARAMS}

# Define stuff for testing duplicates of PLOG and Chebyshev
DUP_PLOG_CHEB_PARAMS = RxnParams(plog_dct=PLOG_DCT)
DUP_PLOG_CHEB_PARAMS.combine_objects(CHEB_PARAMS)  # combine with Cheb
DUP_PLOG_CHEB_RXN_PARAM_DCT = {RXN2: DUP_PLOG_CHEB_PARAMS}

# Define stuff for testing duplicates of two Chebyshevs
DUP_CHEB_PARAMS = RxnParams(cheb_dct=CHEB_DCT)
DUP_CHEB_PARAMS.combine_objects(DUP_CHEB_PARAMS)  # combine with itself
DUP_CHEB_RXN_PARAM_DCT = {RXN2: DUP_CHEB_PARAMS}

# Define stuff for testing duplicates of two Troes
DUP_TROE_PARAMS = RxnParams(troe_dct=TROE_DCT)
DUP_TROE_PARAMS.combine_objects(DUP_TROE_PARAMS)  # combine with itself
DUP_TROE_RXN_PARAM_DCT = {RXN2: DUP_TROE_PARAMS}

# Define stuff for testing duplicates of two Lindemanns
DUP_LIND_PARAMS = RxnParams(lind_dct=LIND_DCT)
DUP_LIND_PARAMS.combine_objects(DUP_LIND_PARAMS)  # combine with itself
DUP_LIND_RXN_PARAM_DCT = {RXN2: DUP_LIND_PARAMS}


def test_orderalphabetically():
    """ reorder reactions alphabetically when writing
    """
    ref_ckin_str = (
          'REACTIONS     CAL/MOLE     MOLES\n\n'
          'H + O2 => O + O + H          1.000E+15     0.000    25000\n\n'
          'H + O2 => OH + O             1.000E+15     0.000    25000\n\n'
          'OH + O => H + O2             1.000E+15     0.000    25000\n\n\n\n'
          'END\n\n')
    ckin_str = writer(RXN_ALL_DCT, sortrxns=True)
    assert ckin_str == ref_ckin_str

def test_headersymbol():
    """ Test if => is written correctly
    """
    ref_ckin_str = (
         'REACTIONS     CAL/MOLE     MOLES\n\n'
         'H + O2 => OH + O             1.000E+15     0.000    25000\n\n'
         'OH + O => H + O2             1.000E+15     0.000    25000\n\n'
         'H + O2 => O + O + H          1.000E+15     0.000    25000\n\n\n\n'
         'END\n\n')
    ckin_str = writer(RXN_ALL_DCT)
    assert ckin_str == ref_ckin_str

def test_arr():
    """ Tests the Arrhenius writer
    """
    ref_ckin_str1 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        '! Block comment\n\n'
        '! Comment\n'
        'H + O2 = OH + O          1.000E+15     0.000    25000  ! Comment\n'
        '! Comment\n\n\n\n'
        'END\n\n')
    ref_ckin_str21 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        '! Block comment\n\n'
        '! Comment\n'
        'H + O2 = OH + O          2.000E+15     0.000    25000  ! Comment\n'
        '! Comment\n\n\n\n'
        'END\n\n')
    ref_ckin_str22 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        '! Block comment\n\n'
        '! Comment\n'
        'H + O2 = OH + O          1.000E+15     0.000    27000  ! Comment\n'
        '  DUP\n'
        'H + O2 = OH + O          1.000E+15     0.000    25000\n'
        '  DUP\n'
        '! Comment\n\n\n\n'
        'END\n\n')
    ref_ckin_str3 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H + O2 = OH + O          1.000E+15     0.000    25000\n'
        '  N2/1.400/   AR/1.000/   \n\n\n\n'  # note three spaces before \n
        'END\n\n')
    ckin_str1 = writer(ARR_RXN_PARAM_DCT, rxn_cmts_dct=RXN_CMTS_DCT1)
    ckin_str21 = writer(DUP_ARR_RXN_PARAM_DCT, rxn_cmts_dct=RXN_CMTS_DCT1)
    ckin_str22 = writer(DUP_ARR_RXN_PARAM_DCT2, rxn_cmts_dct=RXN_CMTS_DCT1)
    ckin_str3 = writer(COLLID_ARR_RXN_PARAM_DCT)
    assert ckin_str1 == ref_ckin_str1
    assert ckin_str21 == ref_ckin_str21
    assert ckin_str22 == ref_ckin_str22
    assert ckin_str3 == ref_ckin_str3


def test_plog():
    """ Tests the PLOG writer
    """
    ref_ckin_str = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H + O2 = OH + O          1.000E+15     0.000    25000\n'
        '  PLOG /     1.000E-01   1.000E+15     0.000    25000 /\n'
        '  PLOG /     1.000E+00   1.000E+15     0.000    25000 /\n'
        '  PLOG /     1.000E+00   1.000E+15     0.000    25000 /\n'
        '  PLOG /     1.000E+02   1.000E+15     0.000    25000 /\n\n\n\nEND\n\n')
    ckin_str = writer(PLOG_RXN_PARAM_DCT)
    assert ckin_str == ref_ckin_str


def test_cheb():
    """ Tests the Chebyshev writer
    """
    ref_ckin_str = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H + O2 (+N2) = OH + O (+N2)          1.000E+00     0.000        0\n'
        '  TCHEB/      500.00     2000.00 /\n'
        '  PCHEB/        0.03      100.00 /\n'
        '  CHEB /           6           4 /\n'
        '  CHEB /  -1.620E+01  -1.183E-01  -5.423E-02  -1.476E-02 /\n'
        '  CHEB /   2.578E+00   1.614E-01   7.359E-02   1.880E-02 /\n'
        '  CHEB /   1.068E-01  -7.235E-02  -2.733E-02  -3.778E-03 /\n'
        '  CHEB /   3.955E-02   1.207E-02   3.402E-04  -2.695E-03 /\n'
        '  CHEB /   8.557E-03   4.345E-03   3.670E-03   1.608E-03 /\n'
        '  CHEB /   8.599E-04  -1.758E-03  -7.502E-04   7.396E-07 /\n\n\n\n'
        'END\n\n')
    ckin_str = writer(CHEB_RXN_PARAM_DCT)
    assert ckin_str == ref_ckin_str


def test_troe():
    """ Tests the Troe writer
    """
    ref_ckin_str = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H + O2 (+N2) = OH + O (+N2)          1.000E+12     1.500    50000\n'
        '  LOW  /                             1.000E+12     1.500    50000  /\n'
        '  TROE /   1.500E+00   8.000E+03   1.000E+02   1.000E+03 /\n'
        '  AR/1.400/   N2/1.700/   \n\n\n\nEND\n\n')
    ckin_str = writer(TROE_RXN_PARAM_DCT)
    assert ckin_str == ref_ckin_str


def test_lind():
    """ Tests the Lindemann writer
    """
    ref_ckin_str = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H + O2 (+N2) = OH + O (+N2)          1.000E+12     1.500    50000\n'
        '  LOW  /                             1.000E+12     1.500    50000  /\n'
        '  AR/1.400/   N2/1.700/   \n\n\n\nEND\n\n')
    ckin_str = writer(LIND_RXN_PARAM_DCT)
    assert ckin_str == ref_ckin_str


def test_dups():
    """ Tests the writing of duplicates
    """
    ref_ckin_str1 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H + O2 (+N2) = OH + O (+N2)          1.000E+15     0.000    25000\n'
        '  PLOG /                 1.000E-01   1.000E+15     0.000    25000 /\n'
        '  PLOG /                 1.000E+00   1.000E+15     0.000    25000 /\n'
        '  PLOG /                 1.000E+00   1.000E+15     0.000    25000 /\n'
        '  PLOG /                 1.000E+02   1.000E+15     0.000    25000 /\n'
        '  DUP\n'
        'H + O2 (+N2) = OH + O (+N2)          1.000E+15     0.000    25000\n'
        '  PLOG /                 1.000E-01   1.000E+15     0.000    25000 /\n'
        '  PLOG /                 1.000E+00   1.000E+15     0.000    25000 /\n'
        '  PLOG /                 1.000E+00   1.000E+15     0.000    25000 /\n'
        '  PLOG /                 1.000E+02   1.000E+15     0.000    25000 /\n'
        '  DUP\n\n\n\nEND\n\n')
    ref_ckin_str2 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H + O2 (+N2) = OH + O (+N2)          1.000E+15     0.000    25000\n'
        '  PLOG /                 1.000E-01   1.000E+15     0.000    25000 /\n'
        '  PLOG /                 1.000E+00   1.000E+15     0.000    25000 /\n'
        '  PLOG /                 1.000E+00   1.000E+15     0.000    25000 /\n'
        '  PLOG /                 1.000E+02   1.000E+15     0.000    25000 /\n'
        '  DUP\n'
        'H + O2 (+N2) = OH + O (+N2)          1.000E+00     0.000        0\n'
        '  TCHEB/      500.00     2000.00 /\n'
        '  PCHEB/        0.03      100.00 /\n'
        '  CHEB /           6           4 /\n'
        '  CHEB /  -1.620E+01  -1.183E-01  -5.423E-02  -1.476E-02 /\n'
        '  CHEB /   2.578E+00   1.614E-01   7.359E-02   1.880E-02 /\n'
        '  CHEB /   1.068E-01  -7.235E-02  -2.733E-02  -3.778E-03 /\n'
        '  CHEB /   3.955E-02   1.207E-02   3.402E-04  -2.695E-03 /\n'
        '  CHEB /   8.557E-03   4.345E-03   3.670E-03   1.608E-03 /\n'
        '  CHEB /   8.599E-04  -1.758E-03  -7.502E-04   7.396E-07 /\n'
        '  DUP\n\n\n\nEND\n\n')
    ref_ckin_str3 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H + O2 (+N2) = OH + O (+N2)          1.000E+00     0.000        0\n'
        '  TCHEB/      500.00     2000.00 /\n'
        '  PCHEB/        0.03      100.00 /\n'
        '  CHEB /           6           4 /\n'
        '  CHEB /  -1.620E+01  -1.183E-01  -5.423E-02  -1.476E-02 /\n'
        '  CHEB /   2.578E+00   1.614E-01   7.359E-02   1.880E-02 /\n'
        '  CHEB /   1.068E-01  -7.235E-02  -2.733E-02  -3.778E-03 /\n'
        '  CHEB /   3.955E-02   1.207E-02   3.402E-04  -2.695E-03 /\n'
        '  CHEB /   8.557E-03   4.345E-03   3.670E-03   1.608E-03 /\n'
        '  CHEB /   8.599E-04  -1.758E-03  -7.502E-04   7.396E-07 /\n'
        '  DUP\n'
        'H + O2 (+N2) = OH + O (+N2)          1.000E+00     0.000        0\n'
        '  TCHEB/      500.00     2000.00 /\n'
        '  PCHEB/        0.03      100.00 /\n'
        '  CHEB /           6           4 /\n'
        '  CHEB /  -1.620E+01  -1.183E-01  -5.423E-02  -1.476E-02 /\n'
        '  CHEB /   2.578E+00   1.614E-01   7.359E-02   1.880E-02 /\n'
        '  CHEB /   1.068E-01  -7.235E-02  -2.733E-02  -3.778E-03 /\n'
        '  CHEB /   3.955E-02   1.207E-02   3.402E-04  -2.695E-03 /\n'
        '  CHEB /   8.557E-03   4.345E-03   3.670E-03   1.608E-03 /\n'
        '  CHEB /   8.599E-04  -1.758E-03  -7.502E-04   7.396E-07 /\n'
        '  DUP\n\n\n\nEND\n\n')

    ref_ckin_str4 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H + O2 (+N2) = OH + O (+N2)          1.000E+12     1.500    50000\n'
        '  LOW  /                             1.000E+12     1.500    50000  /\n'
        '  TROE /   1.500E+00   8.000E+03   1.000E+02   1.000E+03 /\n'
        '  AR/1.400/   N2/1.700/   \n'
        '  DUP\n'
        'H + O2 (+N2) = OH + O (+N2)          1.000E+12     1.500    50000\n'
        '  LOW  /                             1.000E+12     1.500    50000  /\n'
        '  TROE /   1.500E+00   8.000E+03   1.000E+02   1.000E+03 /\n'
        '  AR/1.400/   N2/1.700/   \n'
        '  DUP\n\n\n\nEND\n\n')

    ref_ckin_str5 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H + O2 (+N2) = OH + O (+N2)          1.000E+12     1.500    50000\n'
        '  LOW  /                             1.000E+12     1.500    50000  /\n'
        '  AR/1.400/   N2/1.700/   \n'
        '  DUP\n'
        'H + O2 (+N2) = OH + O (+N2)          1.000E+12     1.500    50000\n'
        '  LOW  /                             1.000E+12     1.500    50000  /\n'
        '  AR/1.400/   N2/1.700/   \n'
        '  DUP\n\n\n\nEND\n\n')

    ckin_str1 = writer(DUP_PLOG_RXN_PARAM_DCT)
    ckin_str2 = writer(DUP_PLOG_CHEB_RXN_PARAM_DCT)
    ckin_str3 = writer(DUP_CHEB_RXN_PARAM_DCT)
    ckin_str4 = writer(DUP_TROE_RXN_PARAM_DCT)
    ckin_str5 = writer(DUP_LIND_RXN_PARAM_DCT)
    assert ckin_str1 == ref_ckin_str1
    assert ckin_str2 == ref_ckin_str2
    assert ckin_str3 == ref_ckin_str3
    assert ckin_str4 == ref_ckin_str4
    assert ckin_str5 == ref_ckin_str5


if __name__ == '__main__':
    test_headersymbol()
    test_orderalphabetically()
    test_arr()
    test_plog()
    test_cheb()
    test_troe()
    test_lind()
    test_dups()
