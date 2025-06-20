""" test mess_io.reader.ped
"""

import os
import numpy as np
from ioformat import pathtools
import mess_io


PATH = os.path.dirname(os.path.realpath(__file__))
IPATH = os.path.join(PATH, 'data', 'inp')
OPATH = os.path.join(PATH, 'data', 'out')
PED_INP_SGL = pathtools.read_file(IPATH, 'me_ktp_ped_ch2o_h.inp')
PED_OUT_SGL = pathtools.read_file(OPATH, 'thf_pyro_ped_ch2o_h.out')
PED_INP_DBL = pathtools.read_file(IPATH, 'me_ktp_ped_c3h8_h.inp')
PED_OUT_DBL = pathtools.read_file(OPATH, 'thf_pyro_ped_c3h8_h.out')
PED_INP_HOT = pathtools.read_file(IPATH, 'me_ktp_pedhot_c6h7.inp')
PED_OUT_HOT = pathtools.read_file(OPATH, 'ped_c6h7.out')
PED_INP_SGL_NEW = pathtools.read_file(IPATH, 'dmm_habs_ped.inp')
PED_OUT_SGL_NEW = pathtools.read_file(OPATH, 'dmm_habs_ped.out')

def test_ped_names():
    """ test mess_io.reader.ped.ped_names
    """

    pedspecies, pedoutput = mess_io.reader.ped.ped_names(PED_INP_SGL)
    assert pedspecies == (('CH2O+H', 'HCO+H2'),)
    assert pedoutput == 'thf_pyro_ped_ch2o_h.out'

    pedspecies, pedoutput = mess_io.reader.ped.ped_names(PED_INP_DBL)
    assert pedspecies == (('C3H8+H', 'CH3CH2CH2+H2'),
                          ('C3H8+H', 'CH3CHCH3+H2'))
    assert pedoutput == 'thf_pyro_ped_c3h8_h.out'

    pedspecies, pedoutput = mess_io.reader.ped.ped_names(PED_INP_HOT)
    assert pedspecies == (('C6H6+H', 'FULVENE+H'),)
    assert pedoutput == 'ped_c6h7.out'


def test_ped_get_ped():
    """ test mess_io.reader.ped.get_ped
    """
    # PESHOT
    energy_dct3 = {'C5H4CH3': 0.0, 'C5H5CH2-1': 1.36, 'C5H5CH2': 23.38,
                   'C5H5CH2-2': 6.78, 'W5': 15.85, 'W6': -2.95, 'FULVENE+H': 50.0, 'C6H6+H': 18.46}
    energy_dct3, _, _, _ = mess_io.reader.pes(PED_INP_HOT)
    ped_dct3 = mess_io.reader.ped.get_ped(
        PED_OUT_HOT, energy_dct3)

    assert np.isclose((ped_dct3[(('C5H4CH3',), ('C6H6', 'H'), (None,))]
                       [1.][500.].iloc[0].iloc[-1]), 0.04, atol=1e-2)
    
    # PES1
    energy_dct1 = {'W0': -0.8, 'CH2O+H': 0.0, 'HCO+H2': -16.7,
                   'B0': 0.0, 'B1': 5.2}
    ped_dct1 = mess_io.reader.ped.get_ped(
        PED_OUT_SGL, energy_dct1)

    assert np.isclose(
        ped_dct1[(('CH2O', 'H'), ('HCO', 'H2'), (None,))
                 ][1][500].index[1], 16.8339,
        atol=1e-3, rtol=1e-3)

    # PES2
    energy_dct2 = {'W0': 0.0, 'C3H8+H': 0.0, 'CH3CH2CH2+H2': -3.53, 'CH3CHCH3+H2': -6.58,
                   'B0': 0.0, 'B1': 10.19, 'B2': 7.67}
    ped_dct2 = mess_io.reader.ped.get_ped(
        PED_OUT_DBL, energy_dct2)

    assert np.isclose(
        ped_dct2[(('C3H8','H'), ('CH3CH2CH2','H2'), (None,))
                 ][1][500].index[1], 3.685142,
        atol=1e-3, rtol=1e-3)
    assert np.isclose(
        ped_dct2[(('C3H8','H'), ('CH3CHCH3','H2'), (None,))
                 ][1][500].index[1], 6.735142,
        atol=1e-3, rtol=1e-3)

    # PES3 NEWFORMAT mess 2024, skip Fake wells
    energy_dct4 = {'FakeW-OH+DMM': -3.0, 'FakeW-H2O+DMM-R1': -23.39,
                   'OH+DMM': 0.0, 'H2O+DMM-R1': -20.39,}
    ped_dct4 = mess_io.reader.ped.get_ped(PED_OUT_SGL_NEW, energy_dct4)
    assert np.isclose(ped_dct4[(('OH', 'DMM'), ('H2O', 'DMM-R1'), (None,))
                               ][1][1000].iloc[50], 2.2e-6, atol=1e-7)

if __name__ == '__main__':
    test_ped_names()
    test_ped_get_ped()
