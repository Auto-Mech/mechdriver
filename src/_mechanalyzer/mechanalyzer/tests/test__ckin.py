""" Tests the ckin.py file
"""

import os
import numpy
from mechanalyzer.parser import ckin_ as ckin


# Set paths
PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')

# File locations
MECH_FILENAMES = ['mech3.txt', 'mech3.txt']
THERMO_FILENAMES = ['thermo3.txt', 'thermo3.txt']

# Values, etc.
TEMPS_LST = [numpy.array([500, 1000, 1500])]
TEMPS = TEMPS_LST[0]  # used for thermo calls, which take a single Numpy array
PRESSURES = [1, 10]
CORRECT_RATES = numpy.array([2.166e+07, 4.746e+10, 6.164e+11])  # Hong
CORRECT_RXN = ((('H', 'O2'), ('OH', 'O'), (None,)),)
CORRECT_PARAMS = [1.04e14, 0, 15286]
CORRECT_SPCS = ('O',)
CORRECT_THERM = numpy.array([60591.29919473, 63105.13563443, 65600.33145539])
MECH_STR = 'REACTIONS\nH+O2=OH+O 1.04e14 0 15286\nEND'


def test_load_rxn_ktp_dcts():
    """ Tests the load_rxn_ktp_dcts function
    """

    rxn_ktp_dcts = ckin.load_rxn_ktp_dcts(MECH_FILENAMES, DAT_PATH, TEMPS_LST,
                                          PRESSURES)
    assert len(rxn_ktp_dcts) == 2
    for rxn_ktp_dct in rxn_ktp_dcts:
        assert tuple(rxn_ktp_dct.keys()) == CORRECT_RXN
        for ktp_dct in rxn_ktp_dct.values():
            assert numpy.allclose(ktp_dct[10][1], CORRECT_RATES, rtol=1e-2)


def test_load_rxn_param_dcts():
    """ Tests the load_rxn_param_dcts function
    """

    rxn_param_dcts = ckin.load_rxn_param_dcts(MECH_FILENAMES, DAT_PATH)
    assert len(rxn_param_dcts) == 2
    for rxn_param_dct in rxn_param_dcts:
        assert tuple(rxn_param_dct.keys()) == CORRECT_RXN
        for params in rxn_param_dct.values():
            arr_tuples = params.arr
            for arr_tuple in arr_tuples:
                assert numpy.allclose(arr_tuple, CORRECT_PARAMS)


def test_load_spc_therm_dcts():
    """ Tests the load_spc_therm_dcts function
    """
    spc_therm_dcts = ckin.load_spc_therm_dcts(THERMO_FILENAMES, DAT_PATH, TEMPS)
    assert len(spc_therm_dcts) == 2
    for spc_therm_dct in spc_therm_dcts:
        assert tuple(spc_therm_dct.keys()) == CORRECT_SPCS
        assert numpy.allclose(spc_therm_dct['O'][1], CORRECT_THERM, rtol=1e-2)


def test_load_spc_nasa7_dcts():
    """ Tests the load_spc_nasa7_dcts function
    """
    spc_nasa7_dcts = ckin.load_spc_nasa7_dcts(THERMO_FILENAMES, DAT_PATH)
    assert len(spc_nasa7_dcts) == 2
    for spc_nasa7_dct in spc_nasa7_dcts:
        assert tuple(spc_nasa7_dct.keys()) == CORRECT_SPCS


def test_parse_rxn_ktp_dct():
    """ Tests the parse_rxn_ktp_dct function
    """

    rxn_ktp_dct = ckin.parse_rxn_ktp_dct(MECH_STR, TEMPS_LST, PRESSURES)
    assert tuple(rxn_ktp_dct.keys()) == CORRECT_RXN
    for ktp_dct in rxn_ktp_dct.values():
        assert numpy.allclose(ktp_dct[10][1], CORRECT_RATES, rtol=1e-2)


if __name__ == '__main__':
    test_load_rxn_ktp_dcts()
    test_load_rxn_param_dcts()
    test_load_spc_therm_dcts()
    test_load_spc_nasa7_dcts()
    test_parse_rxn_ktp_dct()
