""" Tests the writing of the energy transfer section
"""

import os

from ioformat import pathtools
import mess_io.writer


PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH = os.path.join(PATH, 'data', 'inp')


# Core Data
GEO1 = (
    ('C', (-1.4283035320563338, 0.013425343735546437, -0.030302158896694683)),
    ('C', (1.4283027358735494, -0.013425597530894248, 0.0303022919384165)),
    ('H', (-2.1972722614281355, -0.19229727219177065, 1.8778380427620682)),
    ('H', (-2.121310184939721, 1.792702413487708, -0.8231106338374065)),
    ('H', (-2.1448124562913287, -1.5396513482615042, -1.191852168914227)),
    ('H', (2.1448121742707795, 1.539654946791746, 1.1918517388178247)),
    ('H', (2.1972712765396953, 0.1922944277301287, -1.8778395029874426)),
    ('H', (2.121312248031497, -1.7927029137609576, 0.8231123911174519)))
GEO2 = (
    ('O', (-2.4257421043, 6.0797867203, 0.0000000000)),
    ('H', (-0.8184039622, 5.1215269111, 0.0000000000)))
SYM_FACTOR1 = 3.0
SYM_FACTOR2 = 1.0
INTERP_EMAX = 1500
QUANT_LVL_EMAX = 11
POT_SURF_FILE = 'surf.dat'
FLUX_FILE = 'flux.dat'
STOICH = 'C2O1H7'
POT_PREFACTOR = 12.0
POT_EXP = 4.0
TSTLVL = 'ej'

# Rotor Data
HR_GROUP = (2, 3, 4)
HR_AXIS = (0, 1)
HR_SYMMETRY = 1.0
HR_ID = 'D5'
HR_GRID_SIZE = 100
HR_MASS_EXP_SIZE = 5
THERM_POW_MAX = 50.0
LVL_ENE_MAX = 1000.0
POT_EXP_SIZE = 7
HMIN = 15
HMAX = 110
UMBR_GROUP = (3, 4, 5)
UMBR_PLANE = (3, 4, 5)
UMBR_REF = 3
ONEDPOT = {(0.,): 0.00, (30.,): 2.91, (60.,): 9.06, (90.,): 12.63,
           (120.,): 9.97, (150.,): 3.51, (180.,): 0.03, (210.,): 3.49,
           (240.,): 9.96, (270.,): 12.63, (300.,): 9.08, (330.,): 2.93}

# MDHR Data
FLST = (10., 50., 100., 200., 300., 400., 500., 600., 700., 800., 900.)
TWODPOT = {(0, 0): 1.0, (0, 1): 2.0, (0, 2): 3.0, (0, 3): 4.0,
           (1, 1): 6.0, (1, 2): 7.0, (1, 3): 8.0,
           (2, 2): 2.1, (2, 3): 3.1,
           (3, 3): 7.1}
TWODFREQ = {(0, 0): FLST, (0, 1): FLST, (0, 2): FLST, (0, 3): FLST,
            (1, 1): FLST, (1, 2): FLST, (1, 3): FLST,
            (2, 2): FLST, (2, 3): FLST,
            (3, 3): FLST}
THREEDPOT = {(0, 0, 0): 1.0, (0, 0, 1): 2.0,
             (0, 0, 2): 3.0, (0, 0, 3): 4.0,
             (0, 1, 0): 1.0, (0, 1, 1): 2.0,
             (0, 1, 2): 3.0, (0, 1, 3): 4.0,
             (0, 2, 0): 1.0, (0, 2, 1): 2.0,
             (0, 2, 2): 3.0, (0, 2, 3): 4.0,
             (0, 3, 0): 1.0, (0, 3, 1): 2.0,
             (0, 3, 2): 3.0, (0, 3, 3): 4.0}
THREEDFREQ = {(0, 0, 0): FLST, (0, 0, 1): FLST,
              (0, 0, 2): FLST, (0, 0, 3): FLST,
              (0, 1, 0): FLST, (0, 1, 1): FLST,
              (0, 1, 2): FLST, (0, 1, 3): FLST,
              (0, 2, 0): FLST, (0, 2, 1): FLST,
              (0, 2, 2): FLST, (0, 2, 3): FLST,
              (0, 3, 0): FLST, (0, 3, 1): FLST,
              (0, 3, 2): FLST, (0, 3, 3): FLST}
FOURDPOT = {(0, 0, 0, 0): 1.0, (0, 0, 0, 1): 2.0,
            (0, 0, 0, 2): 3.0, (0, 0, 0, 3): 4.0,
            (0, 0, 1, 0): 1.0, (0, 0, 1, 1): 2.0,
            (0, 0, 1, 2): 3.0, (0, 0, 1, 3): 4.0,
            (0, 0, 2, 0): 1.0, (0, 0, 2, 1): 2.0,
            (0, 0, 2, 2): 3.0, (0, 0, 2, 3): 4.0,
            (0, 0, 3, 0): 1.0, (0, 0, 3, 1): 2.0,
            (0, 0, 3, 2): 3.0, (0, 0, 3, 3): 4.0,
            (0, 1, 0, 0): 1.0, (0, 1, 0, 1): 2.0,
            (0, 1, 0, 2): 3.0, (0, 1, 0, 3): 4.0,
            (0, 1, 1, 0): 1.0, (0, 1, 1, 1): 2.0,
            (0, 1, 1, 2): 3.0, (0, 1, 1, 3): 4.0,
            (0, 1, 2, 0): 1.0, (0, 1, 2, 1): 2.0,
            (0, 1, 2, 2): 3.0, (0, 1, 2, 3): 4.0,
            (0, 1, 3, 0): 1.0, (0, 1, 3, 1): 2.0,
            (0, 1, 3, 2): 3.0, (0, 1, 3, 3): 4.0,
            (0, 2, 0, 0): 1.0, (0, 2, 0, 1): 2.0,
            (0, 2, 0, 2): 3.0, (0, 2, 0, 3): 4.0,
            (0, 2, 1, 0): 1.0, (0, 2, 1, 1): 2.0,
            (0, 2, 1, 2): 3.0, (0, 2, 1, 3): 4.0,
            (0, 2, 2, 0): 1.0, (0, 2, 2, 1): 2.0,
            (0, 2, 2, 2): 3.0, (0, 2, 2, 3): 4.0,
            (0, 2, 3, 0): 1.0, (0, 2, 3, 1): 2.0,
            (0, 2, 3, 2): 3.0, (0, 2, 3, 3): 4.0,
            (0, 3, 0, 0): 1.0, (0, 3, 0, 1): 2.0,
            (0, 3, 0, 2): 3.0, (0, 3, 0, 3): 4.0,
            (0, 3, 1, 0): 1.0, (0, 3, 1, 1): 2.0,
            (0, 3, 1, 2): 3.0, (0, 3, 1, 3): 4.0,
            (0, 3, 2, 0): 1.0, (0, 3, 2, 1): 2.0,
            (0, 3, 2, 2): 3.0, (0, 3, 2, 3): 4.0,
            (0, 3, 3, 0): 1.0, (0, 3, 3, 1): 2.0,
            (0, 3, 3, 2): 3.0, (0, 3, 3, 3): 4.0}
FOURDFREQ = {(0, 0, 0, 0): FLST, (0, 0, 0, 1): FLST,
             (0, 0, 0, 2): FLST, (0, 0, 0, 3): FLST,
             (0, 0, 1, 0): FLST, (0, 0, 1, 1): FLST,
             (0, 0, 1, 2): FLST, (0, 0, 1, 3): FLST,
             (0, 0, 2, 0): FLST, (0, 0, 2, 1): FLST,
             (0, 0, 2, 2): FLST, (0, 0, 2, 3): FLST,
             (0, 0, 3, 0): FLST, (0, 0, 3, 1): FLST,
             (0, 0, 3, 2): FLST, (0, 0, 3, 3): FLST,
             (0, 1, 0, 0): FLST, (0, 1, 0, 1): FLST,
             (0, 1, 0, 2): FLST, (0, 1, 0, 3): FLST,
             (0, 1, 1, 0): FLST, (0, 1, 1, 1): FLST,
             (0, 1, 1, 2): FLST, (0, 1, 1, 3): FLST,
             (0, 1, 2, 0): FLST, (0, 1, 2, 1): FLST,
             (0, 1, 2, 2): FLST, (0, 1, 2, 3): FLST,
             (0, 1, 3, 0): FLST, (0, 1, 3, 1): FLST,
             (0, 1, 3, 2): FLST, (0, 1, 3, 3): FLST,
             (0, 2, 0, 0): FLST, (0, 2, 0, 1): FLST,
             (0, 2, 0, 2): FLST, (0, 2, 0, 3): FLST,
             (0, 2, 1, 0): FLST, (0, 2, 1, 1): FLST,
             (0, 2, 1, 2): FLST, (0, 2, 1, 3): FLST,
             (0, 2, 2, 0): FLST, (0, 2, 2, 1): FLST,
             (0, 2, 2, 2): FLST, (0, 2, 2, 3): FLST,
             (0, 2, 3, 0): FLST, (0, 2, 3, 1): FLST,
             (0, 2, 3, 2): FLST, (0, 2, 3, 3): FLST,
             (0, 3, 0, 0): FLST, (0, 3, 0, 1): FLST,
             (0, 3, 0, 2): FLST, (0, 3, 0, 3): FLST,
             (0, 3, 1, 0): FLST, (0, 3, 1, 1): FLST,
             (0, 3, 1, 2): FLST, (0, 3, 1, 3): FLST,
             (0, 3, 2, 0): FLST, (0, 3, 2, 1): FLST,
             (0, 3, 2, 2): FLST, (0, 3, 2, 3): FLST,
             (0, 3, 3, 0): FLST, (0, 3, 3, 1): FLST,
             (0, 3, 3, 2): FLST, (0, 3, 3, 3): FLST}

# Tunnel Data
IMAG_FREQ = 2000.0
WELL_DEPTH1 = 10.0
WELL_DEPTH2 = 20.0
CUTOFF_ENE = 3000.0
TUNNEL_FILE = 'tunnel.dat'


def test__core_rigidrotor_writer():
    """ test mess_io.writer.core_rigidrotor
    """

    core_rigrot1_str = mess_io.writer.core_rigidrotor(
        GEO1, SYM_FACTOR1)

    core_rigrot2_str = mess_io.writer.core_rigidrotor(
        GEO1, SYM_FACTOR1,
        interp_emax=INTERP_EMAX)

    assert core_rigrot1_str == pathtools.read_file(
        INP_PATH, 'core_rigrot1.inp')
    assert core_rigrot2_str == pathtools.read_file(
        INP_PATH, 'core_rigrot2.inp')


def test__core_multirotor_writer():
    """ test mess_io.writer.rotor_internal
        test mess_io.writer.mdhr_data
        test mess_io.writer.core_multirotor
    """

    # Internal Rotor Sections
    rotor_int1_str = mess_io.writer.rotor_internal(
        HR_GROUP, HR_AXIS, HR_SYMMETRY, HR_GRID_SIZE, HR_MASS_EXP_SIZE)
    rotor_int2_str = mess_io.writer.rotor_internal(
        HR_GROUP, HR_AXIS, HR_SYMMETRY, HR_GRID_SIZE, HR_MASS_EXP_SIZE,
        pot_exp_size=POT_EXP_SIZE,
        hmin=HMIN,
        hmax=HMAX,
        geo=GEO1,
        rotor_id=HR_ID)

    assert rotor_int1_str == pathtools.read_file(
        INP_PATH, 'rotor_int1.inp')
    assert rotor_int2_str == pathtools.read_file(
        INP_PATH, 'rotor_int2.inp')

    mdhr_dat_1d_str = mess_io.writer.mdhr_data(
        ONEDPOT)
    mdhr_dat_2dfr_str = mess_io.writer.mdhr_data(
        TWODPOT, freqs=TWODFREQ, nrot=3)
    mdhr_dat_3dfr_str = mess_io.writer.mdhr_data(
        THREEDPOT, freqs=THREEDFREQ, nrot=3)
    mdhr_dat_4dfr_str = mess_io.writer.mdhr_data(
        FOURDPOT, freqs=FOURDFREQ, nrot=3)

    assert mdhr_dat_1d_str == pathtools.read_file(
        INP_PATH, 'mdhr_dat_1d.inp')
    assert mdhr_dat_2dfr_str == pathtools.read_file(
        INP_PATH, 'mdhr_dat_2dfr.inp')
    assert mdhr_dat_3dfr_str == pathtools.read_file(
        INP_PATH, 'mdhr_dat_3dfr.inp')
    assert mdhr_dat_4dfr_str == pathtools.read_file(
        INP_PATH, 'mdhr_dat_4dfr.inp')

    # MultiRotor Core Sections
    core_multirot1_str = mess_io.writer.core_multirotor(
        GEO1, SYM_FACTOR1, POT_SURF_FILE, rotor_int1_str)
    core_multirot2_str = mess_io.writer.core_multirotor(
        GEO1, SYM_FACTOR1, POT_SURF_FILE, rotor_int1_str,
        interp_emax=INTERP_EMAX,
        quant_lvl_emax=QUANT_LVL_EMAX)

    assert core_multirot1_str == pathtools.read_file(
        INP_PATH, 'core_multirot1.inp')
    assert core_multirot2_str == pathtools.read_file(
        INP_PATH, 'core_multirot2.inp')


def test__core_phasespace_writer():
    """ test mess_io.writer.core_phasespace
    """

    # Use the writer to create a string for each of the core sections
    core_pst1_str = mess_io.writer.core_phasespace(
        GEO1, GEO2, SYM_FACTOR2, STOICH)

    core_pst2_str = mess_io.writer.core_phasespace(
        GEO1, GEO2, SYM_FACTOR2, STOICH,
        pot_prefactor=POT_PREFACTOR,
        pot_exp=POT_EXP,
        tstlvl=TSTLVL)

    assert core_pst1_str == pathtools.read_file(
        INP_PATH, 'core_pst1.inp')
    assert core_pst2_str == pathtools.read_file(
        INP_PATH, 'core_pst2.inp')


def test__core_rotd_writer():
    """ test mess_io.writer.core_rotd
    """

    core_rotd_str = mess_io.writer.core_rotd(
        SYM_FACTOR1, FLUX_FILE, STOICH)

    assert core_rotd_str == pathtools.read_file(
        INP_PATH, 'core_rotd.inp')


def test__rotor_hindered_writer():
    """ test mess_io.writer.rotor_hindered
    """

    rotor_hind1_str = mess_io.writer.rotor_hindered(
        HR_GROUP, HR_AXIS, HR_SYMMETRY, ONEDPOT)

    rotor_hind2_str = mess_io.writer.rotor_hindered(
        HR_GROUP, HR_AXIS, HR_SYMMETRY, ONEDPOT,
        hmin=HMIN,
        hmax=HMAX,
        lvl_ene_max=LVL_ENE_MAX,
        therm_pow_max=THERM_POW_MAX,
        geo=GEO1,
        rotor_id=HR_ID)

    rotor_hind3_str = mess_io.writer.rotor_hindered(
        HR_GROUP, HR_AXIS, HR_SYMMETRY, ONEDPOT,
        potential_form='fourier')

    rotor_hind4_str = mess_io.writer.rotor_hindered(
        HR_GROUP, HR_AXIS, HR_SYMMETRY, ONEDPOT,
        hmin=HMIN,
        hmax=HMAX,
        lvl_ene_max=LVL_ENE_MAX,
        therm_pow_max=THERM_POW_MAX,
        geo=GEO1,
        rotor_id=HR_ID,
        potential_form='fourier')

    assert rotor_hind1_str == pathtools.read_file(
        INP_PATH, 'rotor_hind1.inp')
    assert rotor_hind2_str == pathtools.read_file(
        INP_PATH, 'rotor_hind2.inp')
    assert rotor_hind3_str == pathtools.read_file(
        INP_PATH, 'rotor_hind3.inp')
    assert rotor_hind4_str == pathtools.read_file(
        INP_PATH, 'rotor_hind4.inp')


def test__umbrella_writer():
    """ test mess_io.writer.umbrella
    """

    umbrella1_str = mess_io.writer.umbrella_mode(
        UMBR_GROUP, UMBR_PLANE, UMBR_REF, ONEDPOT)
    umbrella2_str = mess_io.writer.umbrella_mode(
        UMBR_GROUP, UMBR_PLANE, UMBR_REF, ONEDPOT,
        geo=GEO1)

    assert umbrella1_str == pathtools.read_file(
        INP_PATH, 'umbrella1.inp')
    assert umbrella2_str == pathtools.read_file(
        INP_PATH, 'umbrella2.inp')


def test__tunnel_eckart_writer():
    """ test mess_io.writer.tunnel_eckart
    """

    tunnel_eckart_str = mess_io.writer.tunnel_eckart(
        IMAG_FREQ, WELL_DEPTH1, WELL_DEPTH2)

    assert tunnel_eckart_str == pathtools.read_file(
        INP_PATH, 'tunnel_eckart.inp')


def test__tunnel_read_writer():
    """ test mess_io.writer.tunnel_read
    """

    tunnel_read1_str = mess_io.writer.tunnel_read(
        IMAG_FREQ, TUNNEL_FILE)
    tunnel_read2_str = mess_io.writer.tunnel_read(
        IMAG_FREQ, TUNNEL_FILE, cutoff_energy=CUTOFF_ENE)

    assert tunnel_read1_str == pathtools.read_file(
        INP_PATH, 'tunnel_read1.inp')
    assert tunnel_read2_str == pathtools.read_file(
        INP_PATH, 'tunnel_read2.inp')

if __name__ == '__main__':
    test__core_multirotor_writer()
    test__core_phasespace_writer()
    test__core_rigidrotor_writer()
    test__core_rotd_writer()
    test__rotor_hindered_writer()
    test__tunnel_eckart_writer()
    test__tunnel_read_writer()
    test__umbrella_writer()