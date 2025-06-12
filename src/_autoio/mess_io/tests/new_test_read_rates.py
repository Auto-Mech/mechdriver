""" test mess_io.reader.rates
"""

import os
import numpy
import automol.util.dict_
from ioformat import pathtools
import mess_io.reader


PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH = os.path.join(PATH, 'data', 'inp')
OUT_PATH = os.path.join(PATH, 'data', 'out')

KTP_INP_STR = pathtools.read_file(INP_PATH, 'example.inp')
KTP_OUT_STR = pathtools.read_file(OUT_PATH, 'rate.out')
KTP_BAD_OUT_STR = pathtools.read_file(OUT_PATH, 'new_rate.out')
KTP_OUT_BAR_STR = pathtools.read_file(OUT_PATH, 'rate.out_bar')
KTP_OUT_TORR_STR = pathtools.read_file(OUT_PATH, 'rate.out_torr')
KE_OUT_STR = pathtools.read_file(OUT_PATH, 'ke.out')
KE_PED_OUT_DBL = pathtools.read_file(OUT_PATH, 'ke_ped_c3h8_h.out')
KTP_OUT_HABS_STR = pathtools.read_file(OUT_PATH, 'HABS.out')

KE_OUT_SGL_NEW = pathtools.read_file(OUT_PATH, 'dmm_habs_ke.out')
KTP_OUT_STR_WFAKE = pathtools.read_file(OUT_PATH, 'dmm_wext.out')
# Set the REACTANT and PRODUCT
REACTANT = 'W1'
PRODUCT = 'P1'

# Define things for testing the filtering function
KTP_DCT1 = {
    # Weird negative, non-physical k(T)s and undefined k(T)s
    0.01: ((200., 400., 600., 800., 1000., 1200.),
           (-1.0e3, 2.0e3, -3.0e3, None, None, None)),
    # Just negative k(T) at first position and undefined k(T)s
    0.1: ((200., 400., 600., 800., 1000., 1200.),
          (-1.0e4, 2.0e4, 3.0e4, None, None, None)),
    # Positive k(T) and undefined k(T)s
    1.0: ((200., 400., 600., 800., 1000., 1200.),
          (1.0e5, 2.0e5, 3.0e5, None, None, None)),
    # Positive k(T)s
    10.0: ((200., 400., 600., 800., 1000., 1200.),
           (1.0e6, 2.0e6, 3.0e6, 4.0e6, 5.0e6, 6.0e6)),
    # Positive k(T)s
    100.0: ((200., 400., 600., 800., 1000., 1200.),
            (1.0e7, 2.0e7, 3.0e7, 4.0e7, 5.0e7, 6.0e7))}
# Handles very small rate constants that may not wish to be modelled
KTP_DCT2 = {
    0.01: ((200., 400., 600., 800., 1000., 1200.),
           (1.0e-25, 5.0e-25, 9.0e-24, 1.5e-23, 6.5e-22, 1.0e-20)),
    1.0: ((200., 400., 600., 800., 1000., 1200.),
          (8.0e-25, 3.0e-24, 2.0e-22, 6.0e-20, 2.5e-18, 5.0e-16)), }
# Dict that should come back empty after filter
KTP_DCT3 = {
    0.01: ((200., 400., 600., 800., 1000., 1200.),
           (-1.0e-25, -5.0e-25, -9.0e-24, -1.5e-23, -6.5e-22, -1.0e-20)),
    1.0: ((200., 400., 600., 800., 1000., 1200.),
          (-8.0e-25, -3.0e-24, -2.0e-22, -6.0e-20, -2.5e-18, -5.0e-16)), }


def test__ktp():
    """ test mess_io.reader.new_rates.ktp_dct
    """

    ref_ktp_vals = [
       [3.62e-14, 4.14e-12, 2.13e-10, 5.90e-09, 1.00e-07, 1.15e-06,
        9.65e-06, 6.21e-05, 3.21e-04, 1.37e-03, 5.04e-03, 4.59e-02,
        2.78e-01, 1.23e+00, 4.25e+00, 1.21e+01, 2.93e+01, 6.26e+01,
        1.21e+02, 2.13e+02, 3.53e+02, 8.30e+02, 2.28e+03],
       [3.63e-13, 3.94e-11, 1.93e-09, 5.12e-08, 8.35e-07, 9.24e-06,
        7.46e-05, 4.64e-04, 2.32e-03, 9.70e-03, 3.46e-02, 3.01e-01,
        1.75e+00, 7.47e+00, 2.50e+01, 6.90e+01, 1.63e+02, 3.40e+02,
        6.41e+02, 1.11e+03, 1.80e+03, 4.09e+03, 1.07e+04],
       [2.49e-12, 2.67e-10, 1.29e-08, 3.36e-07, 5.38e-06, 5.85e-05,
        4.64e-04, 2.84e-03, 1.40e-02, 5.74e-02, 2.02e-01, 1.71e+00,
        9.65e+00, 4.02e+01, 1.31e+02, 3.55e+02, 8.24e+02, 1.69e+03,
        3.12e+03, 5.31e+03, 8.49e+03, 1.86e+04, 4.66e+04],
       [1.34e-11, 1.48e-09, 7.27e-08, 1.92e-06, 3.11e-05, 3.40e-04,
        2.71e-03, 1.66e-02, 8.19e-02, 3.36e-01, 1.18e+00, 9.93e+00,
        5.59e+01, 2.31e+02, 7.49e+02, 2.01e+03, 4.62e+03, 9.36e+03,
        1.72e+04, 2.89e+04, 4.57e+04, 9.80e+04, 2.37e+05],
       [3.96e-11, 4.59e-09, 2.35e-07, 6.43e-06, 1.07e-04, 1.20e-03,
        9.81e-03, 6.14e-02, 3.08e-01, 1.28e+00, 4.56e+00, 3.93e+01,
        2.25e+02, 9.43e+02, 3.10e+03, 8.37e+03, 1.94e+04, 3.95e+04,
        7.26e+04, 1.23e+05, 1.94e+05, 4.14e+05, 9.86e+05],
       [8.39e-11, 1.04e-08, 5.62e-07, 1.62e-05, 2.83e-04, 3.32e-03,
        2.81e-02, 1.82e-01, 9.43e-01, 4.05e+00, 1.48e+01, 1.33e+02,
        7.95e+02, 3.45e+03, 1.17e+04, 3.24e+04, 7.66e+04, 1.59e+05,
        2.97e+05, 5.09e+05, 8.12e+05, 1.76e+06, 4.22e+06],
       [1.23e-10, 1.61e-08, 9.18e-07, 2.78e-05, 5.10e-04, 6.25e-03,
        5.51e-02, 3.71e-01, 1.99e+00, 8.85e+00, 3.34e+01, 3.19e+02,
        2.00e+03, 9.08e+03, 3.19e+04, 9.18e+04, 2.24e+05, 4.78e+05,
        9.14e+05, 1.60e+06, 2.59e+06, 5.77e+06, 1.42e+07],
       [1.51e-10, 2.06e-08, 1.23e-06, 3.92e-05, 7.52e-04, 9.65e-03,
        8.90e-02, 6.26e-01, 3.51e+00, 1.62e+01, 6.35e+01, 6.52e+02,
        4.36e+03, 2.10e+04, 7.78e+04, 2.35e+05, 5.98e+05, 1.33e+06,
        2.63e+06, 4.74e+06, 7.92e+06, 1.84e+07, 4.72e+07],
       [1.72e-10, 2.44e-08, 1.53e-06, 5.14e-05, 1.05e-03, 1.45e-02,
        1.45e-01, 1.11e+00, 6.77e+00, 3.44e+01, 1.49e+02, 1.88e+03,
        1.57e+04, 9.54e+04, 4.49e+05, 1.73e+06, 5.64e+06, 1.60e+07,
        4.07e+07, 9.39e+07, 1.99e+08, 7.36e+08, 3.54e+09]
       ]


    ktp = mess_io.reader.new_rates.get_ktp(
        KTP_OUT_STR, REACTANT, PRODUCT)
    assert numpy.allclose(ktp, ref_ktp_vals)

    # Read files that have units that are in bar, torr instead of atm
    # vals should be same as above; just changed units in output string
    # for testing
    ktp_bar = mess_io.reader.new_rates.get_ktp(
        KTP_OUT_BAR_STR, REACTANT, PRODUCT)
    ktp_torr = mess_io.reader.new_rates.get_ktp(
        KTP_OUT_TORR_STR, REACTANT, PRODUCT)

    assert numpy.allclose(ref_ktp_vals, ktp_bar)
    assert numpy.allclose(ref_ktp_vals, ktp_torr)

def test__filter():
    """ test mess_io.reader.new_rates.ktp_dct with values that need to be 
        filtered
    """

    ref_ktp_vals = [
    [numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan,
     numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan,
     numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan,
     numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan],
    [numpy.nan, numpy.nan, numpy.nan, numpy.nan, 8.35e-07, 9.24e-06,
     7.46e-05, 4.64e-04, 2.32e-03, 9.70e-03, 3.46e-02, 3.01e-01,
     1.75e+00, 7.47e+00, 2.50e+01, 6.90e+01, 1.63e+02, 3.40e+02,
     6.41e+02, 1.11e+03, 1.80e+03, numpy.nan, numpy.nan],
    [numpy.nan, numpy.nan, numpy.nan, numpy.nan, 5.38e-06, 5.85e-05,
     4.64e-04, 2.84e-03, 1.40e-02, 5.74e-02, 2.02e-01, 1.71e+00,
     9.65e+00, 4.02e+01, 1.31e+02, 3.55e+02, 8.24e+02, 1.69e+03,
     3.12e+03, 5.31e+03, 8.49e+03, numpy.nan, numpy.nan],
    [numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, 3.40e-04,
     2.71e-03, 1.66e-02, 8.19e-02, 3.36e-01, 1.18e+00, 9.93e+00,
     5.59e+01, 2.31e+02, 7.49e+02, 2.01e+03, 4.62e+03, 9.36e+03,
     1.72e+04, 2.89e+04, 4.57e+04, numpy.nan, numpy.nan],
    [numpy.nan, numpy.nan, numpy.nan, numpy.nan, 1.07e-04, 1.20e-03,
     9.81e-03, 6.14e-02, 3.08e-01, 1.28e+00, 4.56e+00, 3.93e+01,
     2.25e+02, 9.43e+02, 3.10e+03, 8.37e+03, 1.94e+04, 3.95e+04,
     7.26e+04, 1.23e+05, numpy.nan, numpy.nan, numpy.nan],
    [numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan,
     numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan,
     numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan,
     numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan],
    [numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan,
     numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan,
     numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan,
     numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan],
    [numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan,
     numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan,
     numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan,
     numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan],
    [numpy.nan, numpy.nan, numpy.nan, numpy.nan, 1.05e-03, 1.45e-02,
     1.45e-01, 1.11e+00, 6.77e+00, 3.44e+01, 1.49e+02, 1.88e+03,
     1.57e+04, 9.54e+04, 4.49e+05, 1.73e+06, 5.64e+06, 1.60e+07,
     4.07e+07, 9.39e+07, 1.99e+08, numpy.nan, numpy.nan]]

    ktp = mess_io.reader.new_rates.get_ktp(
        KTP_BAD_OUT_STR, REACTANT, PRODUCT, True, 700, 2e+3, 0.1, 10)
    assert numpy.allclose(
            ktp, ref_ktp_vals, rtol=1e-05, atol=1e-08, equal_nan=True)

def test__ke_dct():
    """ test mess_io.reader.new_rates.ke_dct
    """

    ref_ke_dct = {
        0.0: 0.0, 0.2: 3.5e-17, 0.4: 8.5e-17, 
        0.6: 1.35e-16, 0.8: 1.85e-16, 
        1.0: 2.35e-16, 1.2: 2.85e-16, 1.4: 3.35e-16, 
        1.6: 3.85e-16, 1.8: 4.35e-16, 2.0: 4.85e-16, 
        2.2: 5.35e-16, 2.4: 5.85e-16, 2.6: 6.35e-16, 
        2.8: 6.85e-16, 3.0: 7.35e-16, 3.2: 7.85e-16, 
        3.4: 8.35e-16, 3.6: 8.85e-16, 3.8: 9.35e-16, 
        4.0: 9.85e-16, 4.2: 1.04e-15, 4.4: 1.09e-15, 
        4.6: 1.14e-15, 4.8: 1.19e-15, 5.0: 1.24e-15, 
        5.2: 1.29e-15, 5.4: 1.34e-15, 5.6: 1.39e-15, 
        5.8: 1.44e-15, 6.0: 1.49e-15, 6.2: 1.54e-15, 
        6.4: 1.59e-15, 6.6: 1.64e-15, 6.8: 1.69e-15, 
        7.0: 1.74e-15, 7.2: 1.79e-15, 7.4: 1.84e-15, 
        7.6: 1.89e-15, 7.8: 1.94e-15, 8.0: 1.99e-15, 
        8.2: 2.04e-15, 8.4: 2.09e-15, 8.6: 2.14e-15, 
        8.8: 2.19e-15, 9.0: 2.24e-15, 9.2: 2.29e-15, 
        9.4: 2.34e-15, 9.6: 2.39e-15, 9.8: 2.44e-15, 
        10.0: 2.49e-15, 10.2: 2.54e-15, 10.4: 2.59e-15, 
        10.6: 2.64e-15, 10.8: 2.69e-15, 11.0: 2.74e-15, 
        11.2: 2.79e-15, 11.4: 2.84e-15, 11.6: 2.89e-15, 
        11.8: 2.94e-15, 12.0: 2.99e-15
	}

    ke_dct = mess_io.reader.new_rates.ke_dct(
        KE_OUT_STR, REACTANT, PRODUCT)

    assert set(ke_dct.keys()) == set(ref_ke_dct.keys())
    for ene, ratek in ke_dct.items():
        assert numpy.isclose(ratek, ref_ke_dct[ene])


def test__tp():
    """ test mess_io.reader.rates.pressures
        test mess_io.reader.rates.temperatures
    """

    ref_inp_temps = (600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0,
                     2000.0, 2200.0, 2400.0, 2600.0, 2800.0, 3000.0)
    ref_out_temps = (500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 
                     850.0, 900.0, 950.0, 1000.0, 1100.0, 1200.0, 1300.0, 
                     1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 
                     2000.0, 2200.0, 2500.0)
    ref_inp_press = (0.01, 0.1, 1.0, 10.0, 100.0, numpy.inf)
    ref_out_press = (0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, numpy.inf)

    inp_temps, inp_tunit = mess_io.reader.new_rates.get_temps(
        KTP_INP_STR, mess_file='inp')
    out_temps, out_tunit = mess_io.reader.new_rates.get_temps(
        KTP_OUT_STR, mess_file='out')

    inp_press, inp_punit = mess_io.reader.new_rates.get_press(
        KTP_INP_STR, mess_file='inp')
    out_press, out_punit = mess_io.reader.new_rates.get_press(
        KTP_OUT_STR, mess_file='out')

    assert numpy.allclose(ref_inp_temps, inp_temps)
    assert numpy.allclose(ref_out_temps, out_temps)
    assert inp_tunit == out_tunit == 'K'

    assert ref_inp_press[-1] == inp_press[-1]
    assert ref_out_press[-1] == out_press[-1]
    assert numpy.allclose(ref_inp_press[:-1], inp_press[:-1])
    assert numpy.allclose(ref_out_press[:-1], out_press[:-1])
    assert inp_punit == out_punit == 'atm'


def test__rxns_labels():
    """ test mess_io.reader.new_rates.reactions
    """

    ref_rxns1 = (
        (('W1',), ('P1',), (None,)), 
    	(('W1',), ('P2',), (None,)), 
	    (('W1',), ('P3',), (None,)), 
    	(('W1',), ('Loss',), (None,)), 
	    (('W1',), ('Capture',), (None,)), 
	    (('P1',), ('W1',), (None,)), 
	    (('P1',), ('P2',), (None,)), 
	    (('P1',), ('P3',), (None,)), 
	    (('P1',), ('Loss',), (None,)), 
        (('P1',), ('Capture',), (None,)), 
        (('P2',), ('W1',), (None,)), 
        (('P2',), ('P1',), (None,)), 
        (('P2',), ('P3',), (None,)), 
        (('P2',), ('Loss',), (None,)), 
        (('P2',), ('Capture',), (None,)), 
        (('P3',), ('W1',), (None,)), 
        (('P3',), ('P1',), (None,)), 
        (('P3',), ('P2',), (None,)), 
        (('P3',), ('Loss',), (None,)), 
        (('P3',), ('Capture',), (None,))
	    )

    assert ref_rxns1 == mess_io.reader.new_rates.reactions(
        KTP_OUT_STR)

def test__dos_rovib():
    """ test mess_io.reader.rates.dos_rovib
    """
    dos_df = mess_io.reader.new_rates.dos_rovib(KE_PED_OUT_DBL)
    # check dos info
    assert list(dos_df.columns) == ['CH3CH2CH2', 'H2', 'CH3CHCH3']
    assert numpy.allclose(dos_df.loc[0.4].values, numpy.array(
        [1.099400e+05, 2.86923, 7.682720e+04]))


test__ktp()
test__filter()
test__ke_dct()
test__tp()
test__rxns_labels()
test__dos_rovib()
