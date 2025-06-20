"""
Tests writing the input for a Monte Carlo sampling routine
"""

import os
from ioformat import pathtools
import mess_io.writer


PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH = os.path.join(PATH, 'data', 'inp')


GEO = (('C', (-4.0048955763, -0.3439866053, -0.0021431734)),
       ('O', (-1.3627056155, -0.3412713280, 0.0239463418)),
       ('H', (-4.7435343957, 1.4733340928, 0.7491098889)),
       ('H', (-4.7435373042, -1.9674678465, 1.1075144307)),
       ('H', (-4.6638955748, -0.5501793084, -1.9816675556)),
       ('H', (-0.8648060003, -0.1539639444, 1.8221471090)))
GRAD = ((-4.0048955763, -0.3439866053, -0.0021431734),
        (-1.3627056155, -0.3412713280, 0.0239463418),
        (-4.7435343957, 1.4733340928, 0.7491098889),
        (-4.7435373042, -1.9674678465, 1.1075144307),
        (-4.6638955748, -0.5501793084, -1.9816675556),
        (-0.8648060003, -0.1539639444, 1.8221471090))
HESS = ((0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.959, 0.0, 0.0, -0.452, 0.519, 0.0, -0.477, -0.23),
        (0.0, 0.0, 0.371, 0.0, 0.222, -0.555, 0.0, -0.279, -0.128),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, -0.479, 0.279, 0.0, 0.455, -0.256, 0.0, -0.017, 0.051),
        (0.0, 0.251, -0.185, 0.0, -0.247, 0.607, 0.0, -0.012, 0.09),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, -0.479, -0.279, 0.0, -0.003, -0.263, 0.0, 0.494, 0.279),
        (0.0, -0.251, -0.185, 0.0, 0.025, 0.947, 0.0, 0.292, 0.137))
GEOS = [GEO for i in range(21)]
GRADS = [GRAD for i in range(21)]
HESSES = [HESS for i in range(21)]
ENES = tuple(float(val) for val in range(21))
SYM_FACTOR = 3.0
ELEC_LEVELS = ((2, 0.00),)
FORMULA = 'CH3OH'
DATA_FILE_NAME = 'tau.dat'
REF_CONFIG_FILE_NAME = 'ref.dat'
GROUND_ENE = 32.724
REF_ENE = 0.00
FLUX_IDX = (5, 1, 0, 2)
FLUX_SPAN = 2.0944
EXCLUDED_VOLUME_FACTOR = 0.85679
USE_CM_SHIFT = True


# For output
TEMPS = (100.0, 200.0, 300.0)
LOGQ = (0.100, 0.200, 0.300)
DQ_DT = (0.111, 0.222, 0.333)
D2Q_DT2 = (0.777, 0.888, 0.999)


def test__flux_mode():
    """ test mess_io.writer.fluxional_mode
    """

    # Write the fluxional mode string
    flux_mode1_str = mess_io.writer.fluxional_mode(
        FLUX_IDX)
    flux_mode2_str = mess_io.writer.fluxional_mode(
        FLUX_IDX, span=FLUX_SPAN)

    assert flux_mode1_str == pathtools.read_file(
        INP_PATH, 'flux_mode1.inp')
    assert flux_mode2_str == pathtools.read_file(
        INP_PATH, 'flux_mode2.inp')


def test__species():
    """ test mess_io.writer.mc_species
    """

    flux_mode_str = pathtools.read_file(INP_PATH, 'flux_mode1.inp')

    mc_spc1_str = mess_io.writer.monte_carlo_species(
        GEO, SYM_FACTOR, ELEC_LEVELS,
        flux_mode_str, DATA_FILE_NAME)
    mc_spc2_str = mess_io.writer.monte_carlo_species(
        GEO, SYM_FACTOR, ELEC_LEVELS,
        flux_mode_str, DATA_FILE_NAME,
        ref_config_file_name=REF_CONFIG_FILE_NAME,
        ground_ene=GROUND_ENE,
        reference_ene=REF_ENE,
        excluded_volume_factor=EXCLUDED_VOLUME_FACTOR,
        use_cm_shift=USE_CM_SHIFT)

    assert mc_spc1_str == pathtools.read_file(INP_PATH, 'mc_spc1.inp')
    assert mc_spc2_str == pathtools.read_file(INP_PATH, 'mc_spc2.inp')


def test__dat():
    """ test mess_io.writer.mc_data
    """

    mc_dat1_str = mess_io.writer.monte_carlo_data(
        GEOS, ENES)
    mc_dat2_str = mess_io.writer.monte_carlo_data(
        GEOS, ENES,
        grads=GRADS,
        hessians=HESSES)

    assert mc_dat1_str == pathtools.read_file(INP_PATH, 'mc_dat1.inp')
    assert mc_dat2_str == pathtools.read_file(INP_PATH, 'mc_dat2.inp')

if __name__ == '__main__':
    test__dat()
    test__flux_mode()
    test__dat()