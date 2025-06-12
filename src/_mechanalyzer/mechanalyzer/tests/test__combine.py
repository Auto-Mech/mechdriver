""" Tests the procedures for combining mechanisms
"""

import numpy
from mechanalyzer.calculator import combine
from autoreact.params import RxnParams


ARR_TUPLES1 = [[1e13, 0, 10000]]
ARR_TUPLES2 = [[2e13, 0, 20000]]

RXN_PARAM_DCT1 = {
    (('H2', 'O'), ('OH', 'H'), (None,)): RxnParams({'arr_tuples': ARR_TUPLES1}),
    (('H', 'O2'), ('OH', 'O'), (None,)): RxnParams({'arr_tuples': ARR_TUPLES1})}
RXN_PARAM_DCT2 = {
    (('HV', 'O2V'), ('OHV', 'OV'), (None,)): RxnParams({'arr_tuples': ARR_TUPLES2}),
    (('H2V', 'OV'), ('OHV', 'OHV'), (None,)): RxnParams({'arr_tuples': ARR_TUPLES2})}
RXN_PARAM_DCT3 = {
    (('NH3',), ('NH2', 'H'), (None,)): RxnParams({'arr_tuples': ARR_TUPLES1}),
    (('cN2H2',), ('NNH', 'H'), (None,)): RxnParams({'arr_tuples': ARR_TUPLES1}),
    (('H', 'NNH'), ('tN2H2',), (None,)): RxnParams({'arr_tuples': ARR_TUPLES1})}
RXN_PARAM_DCT4 = {
    (('NH3',), ('NH2', 'H'), (None,)): RxnParams({'arr_tuples': ARR_TUPLES2}),
    (('NNH',), ('N2', 'H'), (None,)): RxnParams({'arr_tuples': ARR_TUPLES2}),
    (('H', 'NNH'), ('N2H2',), (None,)): RxnParams({'arr_tuples': ARR_TUPLES2})}


NASA7_PARAMS1 = [
    '1 7/88', 'N   1O   1          ',  'G', [200.0, 6000.0, 1000.0],
    ([0.48230729E+01, 0.26270251E-02, -0.95850872E-06,
      0.16000712E-09, -0.97752302E-14, 0.80734047E+04, -0.22017208E+01],
     [0.22571502E+01, 0.11304728E-01, -0.13671319E-04, 0.96819803E-08,
      -0.29307182E-11, 0.87417746E+04, 0.10757992E+02])]
NASA7_PARAMS2 = [
    '2 7/88', 'N   1O   1          ',  'G', [200.0, 6000.0, 1000.0],
    ([0.48230729E+01, 0.26270251E-02, -0.95850872E-06,
      0.16000712E-09, -0.97752302E-14, 0.80734047E+04, -0.22017208E+01],
     [0.22571502E+01, 0.11304728E-01, -0.13671319E-04, 0.96819803E-08,
      -0.29307182E-11, 0.87417746E+04, 0.10757992E+02])]

SPC_NASA7_DCT1 = {
    'H': NASA7_PARAMS1,
    'OH': NASA7_PARAMS1,
    'O': NASA7_PARAMS1,
    'H2': NASA7_PARAMS1,
    'O2': NASA7_PARAMS1,
    'O(S)': NASA7_PARAMS1}
SPC_NASA7_DCT2 = {
    'HV': NASA7_PARAMS2,
    'OHV': NASA7_PARAMS2,
    'OV': NASA7_PARAMS2,
    'O2V': NASA7_PARAMS2,
    'H2V': NASA7_PARAMS2,
    'HO2V': NASA7_PARAMS2,
    'O(S)V': NASA7_PARAMS2}
SPC_NASA7_DCT3 = {
    'cN2H2': NASA7_PARAMS1,
    'tN2H2': NASA7_PARAMS1,
    'NH3': NASA7_PARAMS1,
    'NH2': NASA7_PARAMS1,
    'H': NASA7_PARAMS1,
    'NNH': NASA7_PARAMS1}
SPC_NASA7_DCT4 = {
    'N2H2': NASA7_PARAMS2,
    'NH3': NASA7_PARAMS2,
    'NH2': NASA7_PARAMS2,
    'H': NASA7_PARAMS2,
    'NNH': NASA7_PARAMS2,
    'N2': NASA7_PARAMS2}

MECH_SPC_DCT1 = {
    'H': {'inchi': 'InChI=1S/H', 'mult': 2, 'charge': 0, 'exc_flag': 0,
          'fml': {'H':1}},
    'OH': {'inchi': 'InChI=1S/HO/h1H', 'mult': 2, 'charge': 0, 'exc_flag': 0,
           'fml': {'H':1, 'O':1}},
    'O': {'inchi': 'InChI=1S/O', 'mult': 3, 'charge': 0, 'exc_flag': 0,
          'fml': {'O':1}},
    'H2': {'inchi': 'InChI=1S/H2/h1H', 'mult': 1, 'charge': 0, 'exc_flag': 0,
           'fml': {'H':2}},
    'O2': {'inchi': 'InChI=1S/O2/c1-2', 'mult': 1, 'charge': 0, 'exc_flag': 0,
           'fml': {'O':2}},
    'O(S)': {'inchi': 'InChI=1S/O', 'mult': 1, 'charge': 0, 'exc_flag': 0,
             'fml': {'O':1}}}
MECH_SPC_DCT2 = {
    'HV': {'inchi': 'InChI=1S/H', 'mult': 2, 'charge': 0, 'exc_flag': 0,
           'fml': {'H':1}},
    'OHV': {'inchi': 'InChI=1S/HO/h1H', 'mult': 2, 'charge': 0, 'exc_flag': 0,
            'fml': {'H':1, 'O':1}},
    'OV': {'inchi': 'InChI=1S/O', 'mult': 3, 'charge': 0, 'exc_flag': 0,
           'fml': {'O':1}},
    'O2V': {'inchi': 'InChI=1S/O2/c1-2', 'mult': 1, 'charge': 0, 'exc_flag': 0, 
            'fml': {'O':2}},
    'H2V': {'inchi': 'InChI=1S/H2/h1H', 'mult': 1, 'charge': 0, 'exc_flag': 0,
            'fml': {'H':2}},
    'HO2V': {'inchi': 'InChI=1S/HO2/c1-2/h1H', 'mult': 2, 'charge': 0, 'exc_flag': 0, 
             'fml': {'H': 1,'O':2}},
    'O(S)V': {'inchi': 'InChI=1S/O', 'mult': 1, 'charge': 0, 'exc_flag': 0,
              'fml': {'O':1}}}
MECH_SPC_DCT3 = {
    'cN2H2': {'inchi': 'InChI=1S/H2N2/c1-2/h1-2H/b2-1+', 'mult': 1, 'charge': 0, 'exc_flag': 0,
              'fml': {'H':2, 'N':2}},
    'tN2H2': {'inchi': 'InChI=1S/H2N2/c1-2/h1-2H/b2-1-', 'mult': 1, 'charge': 0, 'exc_flag': 0,
              'fml': {'H':2, 'N':2}},
    'NH3': {'inchi': 'InChI=1S/H3N/h1H3', 'mult': 1, 'charge': 0, 'exc_flag': 0,
            'fml': {'H':3, 'N':1}},
    'NH2': {'inchi': 'InChI=1S/H2N/h1H2', 'mult': 2, 'charge': 0, 'exc_flag': 0,
            'fml': {'H':2, 'N':1}},
    'H': {'inchi': 'InChI=1S/H', 'mult': 2, 'charge': 0, 'exc_flag': 0,
          'fml': {'H':1}},
    'NNH': {'inchi': 'InChI=1S/HN2/c1-2/h1H', 'mult': 2, 'charge': 0, 'exc_flag': 0,
            'fml': {'H':1, 'N':2}}}
MECH_SPC_DCT4 = {
    'N2H2': {'inchi': 'InChI=1S/H2N2/c1-2/h1-2H', 'mult': 1, 'charge': 0, 'exc_flag': 0,
             'fml': {'H':2, 'N':2}},
    'NH3': {'inchi': 'InChI=1S/H3N/h1H3', 'mult': 1, 'charge': 0, 'exc_flag': 0,
            'fml': {'H':3, 'N':1}},
    'NH2': {'inchi': 'InChI=1S/H2N/h1H2', 'mult': 2, 'charge': 0, 'exc_flag': 0,
            'fml': {'H':2, 'N':1}},
    'H': {'inchi': 'InChI=1S/H', 'mult': 2, 'charge': 0, 'exc_flag': 0,
          'fml': {'H':1}},
    'NNH': {'inchi': 'InChI=1S/HN2/c1-2/h1H', 'mult': 2, 'charge': 0, 'exc_flag': 0,
            'fml': {'H':1, 'N':2}},
    'N2': {'inchi': 'InChI=1S/N2/c1-2', 'mult': 1, 'charge': 0, 'exc_flag': 0,
           'fml': {'H':1, 'N':2}}}



def test_comb_mechs():
    """ Tests the compare.comb_mechs function
    """

    comb_rxn_param_dct, comb_spc_nasa7_dct, _ = combine.comb_mechs(
        RXN_PARAM_DCT1, RXN_PARAM_DCT2, SPC_NASA7_DCT1, SPC_NASA7_DCT2,
        MECH_SPC_DCT1, MECH_SPC_DCT2)

    for rxn, params in comb_rxn_param_dct.items():
        _, prds, _ = rxn
        if prds in (('OH', 'H'), ('OH', 'O')):
            assert numpy.allclose(params.arr[0], ARR_TUPLES1[0])
        elif prds == ('OH', 'OH'):
            assert numpy.allclose(params.arr[0], ARR_TUPLES2[0])
        else:
            raise NameError(f'An improper reaction, {rxn}, was found')

    for spc, nasa7_params in comb_spc_nasa7_dct.items():
        if spc != 'HO2V':
            assert nasa7_params[0][0] == '1'
        else:  # if on HO2V
            assert nasa7_params[0][0] == '2'


def test_comb_stereo():
    """ Tests the compare.comb_mechs function with stereo involved (i.e., in
        mech1 but not in mech2)
    """

    comb_rxn_param_dct, comb_spc_nasa7_dct, _ = combine.comb_mechs(
        RXN_PARAM_DCT3, RXN_PARAM_DCT4, SPC_NASA7_DCT3, SPC_NASA7_DCT4,
        MECH_SPC_DCT3, MECH_SPC_DCT4, ste_mech1_only=True)

    for rxn, params in comb_rxn_param_dct.items():
        rcts, _, _ = rxn
        if rcts == ('NNH',):  # one reaction comes from the second dct
            assert numpy.allclose(params.arr[0], ARR_TUPLES2[0])
        else:
            assert numpy.allclose(params.arr[0], ARR_TUPLES1[0]), f'{rxn}'

    for spc in comb_spc_nasa7_dct.keys():
        assert spc != 'N2H2'


if __name__ == '__main__':
    #test_comb_mechs()
    test_comb_stereo()
