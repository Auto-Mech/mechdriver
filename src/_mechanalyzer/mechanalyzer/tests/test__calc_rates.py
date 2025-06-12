"""
Test the mechanalyzer.calculator.rates functions
"""

import numpy as np
from autoreact.params import RxnParams
from mechanalyzer.calculator import rates


# Misc definitions used for various tests
PRESSURES = [0.316, 1.0, 10.0, 100.0, 'high']
PRESSURES_NO_HIGH = [0.316, 1.0, 10.0, 100.0]
TEMPS = [np.array([1000.0, 1500.0, 2000.0])]
TEMPS2 = [np.array([300.0, 500.0, 1000.0])]  # only for Chebyshev test
TEMPS3 = [np.array([1000.0, 1500.0, 2000.0]),
          np.array([1000.0, 1500.0, 2000.0]),
          np.array([1000.0, 1500.0, 2000.0]),
          np.array([1000.0, 1500.0, 2000.0, 2500.0]),
          np.array([1000.0, 1500.0, 2000.0, 2500.0])]  # for PLOG test
LOW_P_RXN = (('N2O',), ('N2', 'O'), ('+M',))
HIGH_P_RXN = (('N2O',), ('N2', 'O'), ('(+M)',))
LOW_P_PARAMS = [[1.04E+15, 0, 59810]]  # for Troe and Lindemann
HIGH_P_PARAMS = [[1.26E+12, 0, 62620]]  # for Troe and Lindemann

# Define stuff for testing Arrhenius
ARR_DCT = {'arr_tuples': [[1.04E+15, 0, 59810]]}
ARR_PARAMS = RxnParams(arr_dct=ARR_DCT)
ARR_RXN_PARAM_DCT = {LOW_P_RXN: ARR_PARAMS}
ARRHENIUS_KTS = np.array([8.8273E+1, 2.0086E+6, 3.0299E+8])

# Define stuff for testing PLOG
PLOG_DCT = {
    0.1: [[1.04E+15, 0, 59810]],
    1: [[1.04E+16, 0, 59810]],
    10: [[1.04E+17, 0, 59810]],
    100: [[1.04E+18, 0, 59810]]}
PLOG_PARAMS = RxnParams(plog_dct=PLOG_DCT)
PLOG_RXN_PARAM_DCT = {
    HIGH_P_RXN: PLOG_PARAMS}
# Rates at 0.316 atm
PLOG_0_3ATM_KTS = np.array([2.78943384e+02, 6.34722925e+06, 9.57454718e+08])

# Define stuff for testing Chebyshev
CHEB_DCT = {
    'tlim': [300.0, 2500.0],
    'plim': [1.0, 100.0],
    'alpha': np.array([
        [1.0216E+01, -1.1083E+00, -1.9807E-01],  # from Bozzelli (2002)
        [7.8325E-01, 1.1609E+00, 1.1762E-01],    # C2H5 + O2 = C2H4 + HO2
        [-9.5707E-02, 1.0928E-01, 1.1551E-01],
        [-8.0290E-02, -1.0978E-01, 3.7074E-04],
        [-1.4830E-02, -6.0589E-02, -2.8056E-02],
        [6.9203E-03, -9.7259E-03, -1.3556E-02],
        [7.6648E-03, 6.6865E-03, -8.8244E-04]]),
    'one_atm_arr': [[1, 0, 0]]}  # 1 atm info not used to calculate anything
CHEB_PARAMS = RxnParams(cheb_dct=CHEB_DCT)
CHEB_RXN_PARAM_DCT = {HIGH_P_RXN: CHEB_PARAMS}
CHEB_100ATM_KTS = np.array(
    [1.23879008e+07, 2.77379527e+08, 2.31421081e+10])

# Define stuff for testing Troe
TROE_COEFFS = [1.18, 1E-30, 7900]
TROE_DCT = {'highp_arr': HIGH_P_PARAMS, 'lowp_arr': LOW_P_PARAMS,
            'troe_params': TROE_COEFFS}
TROE_PARAMS = RxnParams(troe_dct=TROE_DCT)
TROE_RXN_PARAM_DCT = {HIGH_P_RXN: TROE_PARAMS}
TROE_10ATM_KTS = np.array([7.76752254e-03, 1.37906401e+02, 1.62555727e+04])

# Define stuff for testing Lindemann
LIND_DCT = {'highp_arr': HIGH_P_PARAMS, 'lowp_arr': LOW_P_PARAMS}
LIND_PARAMS = RxnParams(lind_dct=LIND_DCT)
LIND_RXN_PARAM_DCT = {HIGH_P_RXN: LIND_PARAMS}
LIND_10ATM_KTS = np.array([7.6096e-03, 1.3922e+02, 1.6753e+04])

# Define stuff for testing duplicate Arrhenius
DUP_ARR_DCT = {'arr_tuples': [[1.04E+15, 0, 59810], [1.04E+15, 0, 59810]]}
DUP_ARR_PARAMS = RxnParams(arr_dct=DUP_ARR_DCT)
DUP_ARR_RXN_PARAM_DCT = {LOW_P_RXN: DUP_ARR_PARAMS}

# Define stuff for testing duplicate PLOG
DUP_PLOG_PARAMS = RxnParams(plog_dct=PLOG_DCT)
DUP_PLOG_PARAMS.combine_objects(DUP_PLOG_PARAMS)  # combine with itself
DUP_PLOG_RXN_PARAM_DCT = {HIGH_P_RXN: DUP_PLOG_PARAMS}


def test_arr():
    """ Test the Arrhenius calculator
    """
    rxn_ktp_dct = rates.eval_rxn_param_dct(
        ARR_RXN_PARAM_DCT, TEMPS, PRESSURES)
    calc_rates = rates.read_rxn_ktp_dct(rxn_ktp_dct, LOW_P_RXN, 'high', 'rates')
    assert np.allclose(calc_rates, ARRHENIUS_KTS, rtol=1e-3)


def test_plog():
    """ Test the PLOG calculator
    """
    rxn_ktp_dct = rates.eval_rxn_param_dct(
        PLOG_RXN_PARAM_DCT, TEMPS3, PRESSURES)
    calc_rates = rates.read_rxn_ktp_dct(rxn_ktp_dct, HIGH_P_RXN, 0.316, 'rates')
    assert np.allclose(calc_rates, PLOG_0_3ATM_KTS, rtol=1e-3)


def test_cheb():
    """ Test the Chebyshev calculator
    """
    rxn_ktp_dct = rates.eval_rxn_param_dct(
        CHEB_RXN_PARAM_DCT, TEMPS2, PRESSURES)
    calc_rates = rates.read_rxn_ktp_dct(rxn_ktp_dct, HIGH_P_RXN, 100.0, 'rates')
    assert np.allclose(calc_rates, CHEB_100ATM_KTS, rtol=1e-3)


def test_troe():
    """ Test the Troe calculator
    """
    rxn_ktp_dct = rates.eval_rxn_param_dct(
        TROE_RXN_PARAM_DCT, TEMPS, PRESSURES)
    calc_rates = rates.read_rxn_ktp_dct(rxn_ktp_dct, HIGH_P_RXN, 10.0, 'rates')
    assert np.allclose(calc_rates, TROE_10ATM_KTS, rtol=1e-3)


def test_lind():
    """ Test the Lindemann calculator
    """
    rxn_ktp_dct = rates.eval_rxn_param_dct(
        LIND_RXN_PARAM_DCT, TEMPS, PRESSURES)
    calc_rates = rates.read_rxn_ktp_dct(rxn_ktp_dct, HIGH_P_RXN, 10.0, 'rates')
    assert np.allclose(calc_rates, LIND_10ATM_KTS, rtol=1e-3)


def test_dup_arrhenius():
    """ Test the Arrhenius calculator for a duplicate reaction
    """
    rxn_ktp_dct = rates.eval_rxn_param_dct(
        DUP_ARR_RXN_PARAM_DCT, TEMPS, PRESSURES)
    calc_rates = rates.read_rxn_ktp_dct(rxn_ktp_dct, LOW_P_RXN, 'high', 'rates')
    assert np.allclose(calc_rates, 2*ARRHENIUS_KTS, rtol=1e-3)


def test_dup_plog():
    """ Test the PLOG calculator for a duplicate reaction
    """
    rxn_ktp_dct = rates.eval_rxn_param_dct(
        DUP_PLOG_RXN_PARAM_DCT, TEMPS, PRESSURES)
    calc_rates = rates.read_rxn_ktp_dct(rxn_ktp_dct, HIGH_P_RXN, 0.316, 'rates')
    assert np.allclose(calc_rates, 2*PLOG_0_3ATM_KTS, rtol=1e-3)


def test_check_p_t():
    """ Test the enforcement of the P and T array rules
    """

    temps1 = [np.array([1000, 1500]), np.array([1000, 1500])]
    pressures1 = [1, 10]
    _ = rates.check_p_t(temps1, pressures1)

    temps2 = [np.array([1000, 1500])]
    pressures2 = [1, 10]
    temps2 = rates.check_p_t(temps2, pressures2)
    for temp_arr in temps2:
        assert np.allclose(temp_arr, [1000, 1500])

    temps3 = (np.array([1000, 1500]), np.array([2000, 2500]))
    pressures3 = [1, 10]
    temps3 = rates.check_p_t(temps3, pressures3)
    for temp_arr in temps3:
        assert np.allclose(temp_arr, [1000, 1500]) or \
            np.allclose(temp_arr, [2000, 2500])


def test_read_rxn_ktp_dct():
    """ Test the rxn_ktp_dct reader
    """
    rxn_ktp_dct = {
        (('N2O',), ('N2', 'O'), ('+M',)): {
            'high': (np.array([1000., 1500., 2000.]),
                     np.array([10.0, 100.0, 1000.0]))
            }
    }
    calc_rates = rates.read_rxn_ktp_dct(rxn_ktp_dct, LOW_P_RXN, 'high', 'rates')
    temps = rates.read_rxn_ktp_dct(rxn_ktp_dct, LOW_P_RXN, 'high', 'temps')
    assert np.allclose(calc_rates, np.array([10.0, 100.0, 1000.0]))
    assert np.allclose(temps, np.array([1000., 1500., 2000.]))


def test_remove_high():
    """ Test the remove_high function
    """

    # Test the case where 'high' is included
    temps_lst = rates.check_p_t(TEMPS, PRESSURES)  # adds duplicate temps
    temps_lst, pressures = rates.remove_high(temps_lst, PRESSURES)
    assert np.allclose(pressures, PRESSURES_NO_HIGH)
    assert len(temps_lst) == len(PRESSURES_NO_HIGH)

    # Test the case where 'high' is not included
    temps_lst = rates.check_p_t(TEMPS, PRESSURES_NO_HIGH)  # adds duplicate temps
    temps_lst, pressures = rates.remove_high(temps_lst, PRESSURES_NO_HIGH)
    assert np.allclose(pressures, PRESSURES_NO_HIGH)
    assert len(temps_lst) == len(PRESSURES_NO_HIGH)


if __name__ == '__main__':
    test_arr()
    test_plog()
    test_cheb()
    test_troe()
    test_lind()
    test_dup_arrhenius()
    test_dup_plog()
    test_check_p_t()
    test_read_rxn_ktp_dct()
    test_remove_high()
