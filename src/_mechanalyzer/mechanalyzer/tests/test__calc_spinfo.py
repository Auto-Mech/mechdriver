""" test mechanalyzer.calculator.statmodels
    mess_io required for input generation - too large input strings otherwise
    statmodels.ped_models is a class, tested in test__prompt.py
"""

import os
import numpy as np
from ioformat import pathtools
from ioformat import remove_comment_lines
import autoparse.pattern as app
import mess_io
from mechanalyzer.calculator.spinfo_frommess import get_info


PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH_SGL = os.path.join(PATH, 'data', 'prompt', 'CH2O_H')
INP_PATH_DBL = os.path.join(PATH, 'data', 'prompt', 'C3H8_OH')
PED_INP_SGL = pathtools.read_file(INP_PATH_SGL, 'me_ktp_ped.inp')
PED_INP_DBL = pathtools.read_file(INP_PATH_DBL, 'me_ktp_ped.inp')
GET_SPECIES =  pathtools.read_file(PATH, 'data/getspecies.inp')
PED_INP_SGL = remove_comment_lines(
    PED_INP_SGL, delim_pattern=app.escape('#'))

geom0A = (('C', (3.7204537, -0.8199244, 0.17718337)),
            ('C', (1.1214556, -0.2664461, -0.809017)), 
            ('H', (5.1261625, 0.32049152, -0.80527407)), 
            ('H', (3.832703, -0.39406428, 2.1863515)), 
            ('H', (4.1944613, -2.805293, -0.108317584)), 
            ('C', (-1.1213611, -0.2657656, 0.8088847)), 
            ('O', (-0.00011342155, 2.075293, -0.00018903591)), 
            ('H', (0.82217395, -0.65714556, -2.8023632)), 
            ('C', (-3.7204158, -0.82005674, -0.17689982)), 
            ('H', (-0.82215506, -0.6565595, 2.8023062)), 
            ('H', (-4.1935163, -2.8057468, 0.108090736)), 
            ('H', (-3.8332896, -0.3936862, -2.1859357)), 
            ('H', (-5.1263895, 0.31931946, 0.8063327)), 
            ('C', (0.011115313, -1.2618715, 9.782401)), 
            ('O', (0.011115313, 1.3222874, 9.782401)), 
            ('H', (-2.000775, -1.8087524, 9.782401)), 
            ('H', (0.8460492, -2.076087, 11.4915695)), 
            ('H', (0.8460492, -2.076087, 8.073252)))

freq0A =np.array([10.,20.,30.,40.,50.,60.,244.,280.,455.,469.,739.,813.,
894.,959.,1020.,1028.,1113.,1124.,1158.,1170.,1262.,1340.,1383.,1387.,
1443.,1454.,1464.,1468.,1488.,2918.,2919.,2976.,2979.,2980.,2980.,3000.,
3003.,730.,948.,1094.,1352.,1353.,1489.,2814.,2884.,2925.], dtype=np.float32)

geom_C4H8ORvEsWvAA0 = (('C', (3.7204537, -0.73232514, 0.17716447)), 
                       ('C', (1.1214556, -0.17884688, -0.809017)), 
                       ('H', (5.1261625, 0.40807185, -0.80527407)), 
                       ('H', (3.8326843, -0.30646503, 2.1863515)), 
                       ('H', (4.1944423, -2.7177126, -0.108336486)), 
                       ('C', (-1.1213611, -0.17818525, 0.8088658)), 
                       ('O', (-0.00011342155, 2.1628923, -0.00018903591)), 
                       ('H', (0.82217395, -0.56956524, -2.802382)), 
                       ('C', (-3.7204158, -0.73247635, -0.17691872)), 
                       ('H', (-0.82215506, -0.5689792, 2.8022873)), 
                       ('H', (-4.193535, -2.7181664, 0.108071834)), 
                       ('H', (-3.8332896, -0.30610585, -2.1859548)), 
                       ('H', (-5.1263895, 0.4068998, 0.8063327)),)

freq_C4H8ORvEsWvAA0 = np.array([ 244.,280.,455.,469.,739.,813.,894.,959.,1020.,1028.,1113.,1124.,
1158.,1170.,1262.,1340.,1383.,1387.,1443.,1454.,1464.,1468.,1488.,2918.,
2919.,2976.,2979.,2980.,2980.,3000.,3003.], dtype=np.float32)

geom_CH3OS58cwB = (('C', (0.020113422, -1.0881096, -0.0)), 
                    ('O', (0.020113422, 1.4960493, 0.0)), 
                    ('H', (-1.9917771, -1.6349906, -0.0)), 
                    ('H', (0.8550662, -1.9023441, 1.7091494)), 
                    ('H', (0.8550662, -1.9023441, -1.7091494)))

freq_CH3OS58cwB = np.array([730.,948.,1094.,1352.,1353.,1489.,2814.,2884.,2925.], dtype=np.float32)


def test_get_dof_info():
    """ test mechanalyzer.calculator.spinfo_frommess.get_info
    """
    # PES0 
    species_blocks_pes0 =  mess_io.reader.get_species(GET_SPECIES)
    #unimol species
    dof0A = get_info(species_blocks_pes0['FakeW-C4H8ORvEsWvAA0+CH3O-S58cwB'])
    assert np.allclose(list(dof0A.loc['FakeW-C4H8ORvEsWvAA0+CH3O-S58cwB',['n_atoms', 'vib dof','rot dof', 'mw', 'symmetry']].values), [
            18., 48., 3., 0.103, 1], atol=1e-4)
    for i, el in enumerate(dof0A.loc['FakeW-C4H8ORvEsWvAA0+CH3O-S58cwB', 'geometry']):
        assert np.allclose(el[1], geom0A[i][1])
    assert np.allclose(dof0A.loc['FakeW-C4H8ORvEsWvAA0+CH3O-S58cwB', 'freqs'], freq0A)
    # bimol species
    dof0B = get_info(species_blocks_pes0['C4H8ORvEsWvAA0+CH3O-S58cwB'])
    assert np.allclose(list(dof0B.loc['C4H8ORvEsWvAA0',['n_atoms', 'vib dof','rot dof', 'mw', 'symmetry']].values), [
            13., 33., 3., 0.072, 1], atol=1e-4)   
    assert np.allclose(list(dof0B.loc['C4H8ORvEsWvAA0',['n_atoms', 'vib dof','rot dof', 'mw', 'symmetry']].values), [
            13., 33., 3., 0.072, 1], atol=1e-4)
    assert np.allclose(list(dof0B.loc['CH3O-S58cwB',['n_atoms', 'vib dof','rot dof', 'mw', 'symmetry']].values), [
            5., 9., 3., 0.031, 1], atol=1e-4)
    assert np.allclose(list(dof0B.loc['TS',['n_atoms', 'vib dof','rot dof', 'mw',]].values), [
            18., 47., 3., 0.103], atol=1e-4) 

    for i, el in enumerate(dof0B.loc['C4H8ORvEsWvAA0', 'geometry']):
        assert np.allclose(el[1], geom_C4H8ORvEsWvAA0[i][1])
    assert np.allclose(dof0B.loc['C4H8ORvEsWvAA0', 'freqs'], freq_C4H8ORvEsWvAA0)
    for i, el in enumerate(dof0B.loc['CH3O-S58cwB', 'geometry']):
        assert np.allclose(el[1], geom_CH3OS58cwB[i][1])
    assert np.allclose(dof0B.loc['CH3O-S58cwB', 'freqs'], freq_CH3OS58cwB)   
    
    # PES1
    species_blocks_ped1 = mess_io.reader.get_species(PED_INP_SGL)
    dof1 = get_info(
        species_blocks_ped1['R'])

    assert np.allclose(list(dof1.loc['HCO',['n_atoms', 'vib dof','rot dof', 'mw']].values), [
               3., 3., 3., 0.029], atol=1e-4)
    assert np.allclose(list(dof1.loc['H2',['n_atoms', 'vib dof','rot dof', 'mw']].values), [
               2.0, 1.0, 2.0, 0.002], atol=1e-4)
    assert np.allclose(list(dof1.loc['TS',['n_atoms', 'vib dof','rot dof', 'mw']].values), [
               5.0, 8.0, 3.0, 0.031], atol=1e-4)

    # PES2
    species_blocks_ped2 = mess_io.reader.get_species(PED_INP_DBL)
    dof2 = get_info(
        species_blocks_ped2['CH3CH2CH2+H2O'])

    assert np.allclose(list(dof2.loc['CH3CH2CH2',['vib dof','rot dof', 'mw']].values), [
               24., 3., 0.043], atol=1e-4, rtol=1e-3)
    assert np.allclose(list(dof2.loc['H2O',['vib dof','rot dof', 'mw']].values), [
               3., 3., 0.018], atol=1e-4)


if __name__ == '__main__':
    test_get_dof_info()
