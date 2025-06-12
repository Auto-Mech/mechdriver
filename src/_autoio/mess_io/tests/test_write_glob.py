""" test the writing of the global keyword section for reactions and messpf
"""

import os
import numpy
import pandas
from ioformat import pathtools
import mess_io.writer


PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH = os.path.join(PATH, 'data', 'inp')
OUT_PATH = os.path.join(PATH, 'data', 'out')


TEMPS = (100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0)
PRESSURES = (0.01, 0.1, 1.0, 10.0, 100.0)
TEMP_STEP = 50
NTEMPS = 40
REL_TEMP_INC = 0.002
ATOM_DIST_MIN = 1.05
EXCESS_ENE_TEMP = 10000.0
WELL_EXTEND = 20000.0

# For messhr input
GEO = (
    ('C', (0.0, 0.0, 0.0)),
    ('C', (0.0, 0.0, 2.886481067591814)),
    ('H', (0.0, 1.9289757370835872, -0.7491252306539922)),
    ('H', (1.6705367921659746, -0.9644784199111662, -0.7491252306539922)),
    ('H', (-1.6705367921659746, -0.9644784199111662, -0.7491252306539922)),
    ('H', (1.6705367921659746, 0.9644784199111662, 3.6356062982458064)),
    ('H', (-0.0, -1.9289757370835872, 3.6356062982458064)),
    ('H', (-1.6705367921659746, 0.9644784199111662, 3.6356062982458064)))
HR_STR = """Rotor  Hindered   # D5
  Group                  3   4   5
  Axis                   1   2
  Symmetry               3
  Potential[kcal/mol]    8
    0.00    0.45    1.55    2.66    3.13    2.66
    1.55    0.45
End  ! HindRot\n"""


# For messpf output
FORMULA = 'CH4O'
THM_TEMPS = numpy.array([
    100., 200., 300., 400., 500., 600., 700., 800., 900., 1000.,
    1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000.,
    2100., 2200., 2300., 2400., 2500., 2600., 2700., 2800., 2900., 3000.,
    298.2
])
DATA = numpy.array([
    [60.1983, 0.0302761, -0.000290427, 125.642, 6.26152],
    [62.3581, 0.0163207, -6.51731e-05, 130.404, 7.79248],
    [63.7557, 0.0123154, -2.28447e-05, 134.037, 10.5982],
    [64.9053, 0.0109421, -6.74453e-06, 137.677, 15.2508],
    [65.9811, 0.010705, 1.21039e-06, 141.754, 21.8742],
    [67.0664, 0.0110813, 5.9895e-06, 146.486, 30.7096],
    [68.2106, 0.011859, 9.42064e-06, 152.044, 42.1657],
    [69.4485, 0.0129453, 1.22352e-05, 158.587, 56.7204],
    [70.8085, 0.014296, 1.47438e-05, 166.278, 74.8681],
    [72.3157, 0.0158886, 1.70879e-05, 175.279, 97.1045],
    [73.9938, 0.0177105, 1.93376e-05, 185.754, 123.924],
    [75.8652, 0.0197543, 2.15298e-05, 197.865, 155.822],
    [77.9519, 0.0220153, 2.36856e-05, 211.779, 193.291],
    [80.2754, 0.0244906, 2.58173e-05, 227.657, 236.825],
    [82.8571, 0.0271782, 2.79324e-05, 245.665, 286.916],
    [85.7181, 0.0300767, 3.00358e-05, 265.967, 344.056],
    [88.8794, 0.033185, 3.21303e-05, 288.727, 408.737],
    [92.3621, 0.0365025, 3.4218e-05, 314.109, 481.447],
    [96.1869, 0.0400285, 3.63002e-05, 342.276, 562.677],
    [100.375, 0.0437624, 3.83778e-05, 373.392, 652.913],
    [104.946, 0.0477039, 4.04514e-05, 407.622, 752.643],
    [109.922, 0.0518526, 4.25216e-05, 445.127, 862.354],
    [115.324, 0.0562081, 4.45887e-05, 486.072, 982.53],
    [121.171, 0.0607702, 4.6653e-05, 530.619, 1113.66],
    [127.485, 0.0655386, 4.87146e-05, 578.931, 1256.22],
    [134.285, 0.0705131, 5.07739e-05, 631.171, 1410.71],
    [141.594, 0.0756933, 5.2831e-05, 687.501, 1577.6],
    [149.431, 0.0810792, 5.4886e-05, 748.084, 1757.37],
    [157.817, 0.0866705, 5.6939e-05, 813.082, 1950.52],
    [166.772, 0.0924669, 5.89903e-05, 882.657, 2157.52],
    [63.7335, 0.0123569, -2.32761e-05, 133.973, 10.5319]
])
THM_VALS = pandas.DataFrame(
    data=DATA,
    index=THM_TEMPS,
    columns=('logq', 'dq_dt', 'd2q_dt2', 'svals', 'cpvals')
)


def test__global_rates_input():
    """ test mess_io.writer.global_rates_input
    """

    glob1_str_v1 = mess_io.writer.global_rates_input_v1(
        TEMPS, PRESSURES)
    glob1_str_v2 = mess_io.writer.global_rates_input_v2(
        TEMPS, PRESSURES)
    glob2_str = mess_io.writer.global_rates_input_v1(
        TEMPS, PRESSURES, well_extension=None)
    glob3_str = mess_io.writer.global_rates_input_v1(
        TEMPS, PRESSURES,
        excess_ene_temp=EXCESS_ENE_TEMP,
        chem_eig_max=0.2,
        well_extension=WELL_EXTEND)
    glob4_str = mess_io.writer.global_rates_input_v1(
        TEMPS, PRESSURES,
        calculation_method='direct')
    glob5_str = mess_io.writer.global_rates_input_v1(
        TEMPS, PRESSURES,
        float_type='quadruple')
    glob6_str = mess_io.writer.global_rates_input_v1(
        TEMPS, PRESSURES,
        ped_spc_lst=('RH_R1', 'RH_R2'))
    glob7_str = mess_io.writer.global_rates_input_v1(
        TEMPS, PRESSURES,
        hot_enes_dct={'W1': (10.0, 1.0, 30.0), 'W2': (40.0, 1.0, 60.0)})
    glob8_str = mess_io.writer.global_rates_input_v1(
        TEMPS, PRESSURES,
        micro_out_params=(0.0, 140.0, 0.1))

    assert glob1_str_v1 == pathtools.read_file(INP_PATH, 'glob_rxn1.inp')
    assert glob1_str_v2 == pathtools.read_file(INP_PATH, 'glob_rxn1_v2.inp')
    assert glob2_str == pathtools.read_file(INP_PATH, 'glob_rxn2.inp')
    assert glob3_str == pathtools.read_file(INP_PATH, 'glob_rxn3.inp')
    assert glob4_str == pathtools.read_file(INP_PATH, 'glob_rxn4.inp')
    assert glob5_str == pathtools.read_file(INP_PATH, 'glob_rxn5.inp')
    assert glob6_str == pathtools.read_file(INP_PATH, 'glob_rxn6.inp')
    assert glob7_str == pathtools.read_file(INP_PATH, 'glob_rxn7.inp')
    assert glob8_str == pathtools.read_file(INP_PATH, 'glob_rxn8.inp')


def test__global_pf_input():
    """ tests writing the section to a file for messpf
    """

    glob1_str = mess_io.writer.global_pf_input(
        temperatures=TEMPS)
    glob2_str = mess_io.writer.global_pf_input(
        temperatures=())
    glob3_str = mess_io.writer.global_pf_input(
        temperatures=(),
        temp_step=TEMP_STEP,
        ntemps=NTEMPS,
        rel_temp_inc=REL_TEMP_INC,
        atom_dist_min=ATOM_DIST_MIN)
    glob4_str = mess_io.writer.global_pf_input(
        temperatures=TEMPS,
        float_type='quadruple')

    assert glob1_str == pathtools.read_file(
        INP_PATH, 'glob_pf1.inp').rstrip()
    assert glob2_str == pathtools.read_file(
        INP_PATH, 'glob_pf2.inp').rstrip()
    assert glob3_str == pathtools.read_file(
        INP_PATH, 'glob_pf3.inp').rstrip()
    assert glob4_str == pathtools.read_file(
        INP_PATH, 'glob_pf4.inp').rstrip()


def test__messhr_inp_str():
    """ test mess_io.writer.messhr_inp_str_str
    """

    messhr_str = mess_io.writer.messhr_inp_str(GEO, HR_STR)
    assert messhr_str == pathtools.read_file(INP_PATH, 'full_hr.inp')


def test__full_str():
    """ test mess_io.writer.messrates_inp_str
        test mess_io.writer.messpf_inp_str
    """

    glob_keys_str = '<FAKE GLOB KEYS STR>'
    glob_etrans_str = '<FAKE ETRANS STR>'
    rxn_chan_str = '<FAKE REACTION CHANNEL STR>'
    pf_str = '<FAKE PF CHANNEL STR>'

    rates_str = mess_io.writer.messrates_inp_str(
        glob_keys_str, glob_etrans_str, rxn_chan_str)
    pf_str = mess_io.writer.messpf_inp_str(
        glob_keys_str, pf_str)

    assert rates_str == pathtools.read_file(INP_PATH, 'full_rates.inp')
    assert pf_str == pathtools.read_file(INP_PATH, 'full_pf.inp')


def test__pf_output():
    """ test mess_io.writer.pf_output
    """

    pf1_str = mess_io.writer.pf_output(
        FORMULA, THM_TEMPS,
        THM_VALS['logq'].values,
        THM_VALS['dq_dt'].values,
        THM_VALS['d2q_dt2'].values)

    pf2_str = mess_io.writer.pf_output(
        FORMULA, THM_TEMPS,
        THM_VALS['logq'].values,
        THM_VALS['dq_dt'].values,
        THM_VALS['d2q_dt2'].values,
        svals=THM_VALS['svals'].values,
        cpvals=THM_VALS['cpvals'].values)

    assert pf1_str == pathtools.read_file(OUT_PATH, 'pf.dat2')
    assert pf2_str == pathtools.read_file(OUT_PATH, 'pf.dat3')

if __name__ == '__main__':
    test__full_str()
    test__global_pf_input()
    test__global_rates_input()
    test__messhr_inp_str()
    test__pf_output()