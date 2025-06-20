""" test mechanalyzer.calculator prompt functionalities
"""

import os
import numpy as np
from ioformat import pathtools, remove_comment_lines
import autoparse.pattern as app
import mess_io
from mechanalyzer.calculator import spinfo_frommess
from mechanalyzer.calculator import ene_partition
from mechanalyzer.calculator import bf

PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH = os.path.join(PATH, 'data', 'prompt', 'C3H8_OH')
PED_INP = pathtools.read_file(INP_PATH, 'me_ktp_ped.inp')
PED_INP = remove_comment_lines(PED_INP, delim_pattern=app.escape('!'))
PED_INP = remove_comment_lines(PED_INP, delim_pattern=app.escape('#'))
PED_OUT = pathtools.read_file(INP_PATH, 'ped.out')
KE_PED_OUT = pathtools.read_file(INP_PATH, 'ke_ped.out')
HOT_OUT = pathtools.read_file(INP_PATH, 'me_ktp_hoten.log')

T = [400.0, 600.0, 800.0, 1200.0, 1800.0, 2000.0]
P = [0.1, 1.0, 100.0]

ENERGY_DCT = {'W0': -6.0, 'C3H8+OH': -2.2, 'CH3CH2CH2+H2O': -20.14, 'CH3CHCH3+H2O': -23.18,
              'B0': -2.2, 'B1': 0.0, 'B2': -0.89}

KTP_DCT = {
    (('C3H8', 'OH',), ('CH3CH2CH2', 'H2O',), (None,)): {
        1.0: ((400.0, 600.0, 800.0, 1200.0, 1800.0, 2000.0),
              (5.14854e-13, 1.49542e-12, 2.95517e-12,
               7.33206e-12, 1.72203e-11, 2.12587e-11))
    },
    (('C3H8', 'OH',), ('CH3CHCH3', 'H2O',), (None,)): {
        1.0: ((400.0, 600.0, 800.0, 1200.0, 1800.0, 2000.0),
              (1.8988e-13, 4.35041e-13, 7.7321e-13,
               1.73144e-12, 3.78948e-12, 4.61016e-12))
    }
}

LABELS = list(KTP_DCT.keys())


FRAG_REACS = ('C3H8', 'OH')

HOT_FRAG_DCT = {
    'CH3CH2CH2': ('CH3CH2CH2',),
    'CH3CHCH3': ('CH3CHCH3',),
    'C2H4+CH3': ('C2H4', 'CH3'),
    'CH3CHCH2+H': ('CH3CHCH2', 'H')
}

HOTSPECIES = {'CH3CH2CH2': 3.19, 'CH3CHCH3': 0.0}


# test different models
def test_equip_simple():
    """ test statmodels.pedmodels.equip_simple
    """

    dof_dct, ped_dct, _, _, _ = _read_data()

    ped_df_frag1_dct = ene_partition.ped_frag1(
        ped_dct[(('C3H8', 'OH',), ('CH3CH2CH2', 'H2O',), (None,))
                ], 'CH3CH2CH2', 'H2O', 'equip_simple',
        dof_info=dof_dct[(('C3H8', 'OH',), ('CH3CH2CH2', 'H2O',), (None,))])

    ped_600 = ped_df_frag1_dct[1.0][600]
    ped_1200 = ped_df_frag1_dct[1.0][1200]

    assert np.isclose((ped_600.iloc[100]), 0.1043, atol=1e-4, rtol=1e-4)
    assert np.isclose((ped_1200.iloc[100]), 0.02105, atol=1e-4, rtol=1e-4)
    assert np.isclose(np.trapz(ped_600.values, x=ped_600.index), 1)
    assert np.isclose(np.trapz(ped_1200.values, x=ped_1200.index), 1)


def test_equip_phi():
    """ test statmodels.pedmodels.equip_phi
    """

    dof_dct, ped_dct, _, _, _ = _read_data()

    ped_df_frag1_dct = ene_partition.ped_frag1(
        ped_dct[(('C3H8', 'OH',), ('CH3CH2CH2', 'H2O',), (None,))
                ], 'CH3CH2CH2', 'H2O', 'equip_phi',
        dof_info=dof_dct[(('C3H8', 'OH',), ('CH3CH2CH2', 'H2O',), (None,))])
    ped_600 = ped_df_frag1_dct[1.0][600]
    ped_1200 = ped_df_frag1_dct[1.0][1200]

    assert np.isclose((ped_600.iloc[166]), 0.045, atol=1e-3, rtol=1e-2)
    assert np.isclose((ped_1200.iloc[166]), 0.0376, atol=1e-3, rtol=1e-2)
    assert np.isclose(np.trapz(ped_600.values, x=ped_600.index), 1)
    assert np.isclose(np.trapz(ped_1200.values, x=ped_1200.index), 1)


def test_beta_phi1a():
    """ test statmodels.pedmodels.beta_phi1a
    """

    dof_dct, ped_dct, _, _, _ = _read_data()

    ped_df_frag1_dct = ene_partition.ped_frag1(
        ped_dct[(('C3H8', 'OH',), ('CH3CHCH3', 'H2O',), (None,))
                ], 'CH3CHCH3', 'H2O', 'beta_phi1a',
        dof_info=dof_dct[(('C3H8', 'OH',), ('CH3CHCH3', 'H2O',), (None,))])
    ped_400 = ped_df_frag1_dct[1.0][400]
    ped_800 = ped_df_frag1_dct[1.0][800]

    assert np.isclose((ped_400.iloc[166]), 0.001, atol=1e-3, rtol=1e-2)
    assert np.isclose((ped_800.iloc[166]), 0.0589, atol=1e-3, rtol=1e-2)
    assert np.isclose(np.trapz(ped_400.values, x=ped_400.index), 1)
    assert np.isclose(np.trapz(ped_800.values, x=ped_800.index), 1)


def test_beta_phi2a():
    """ test statmodels.pedmodels.beta_phi2a
    """

    dof_dct, ped_dct, _, _, _ = _read_data()

    ped_df_frag1_dct = ene_partition.ped_frag1(
        ped_dct[(('C3H8', 'OH',), ('CH3CHCH3', 'H2O',), (None,))
                ], 'CH3CHCH3', 'H2O', 'beta_phi2a',
        dof_info=dof_dct[(('C3H8', 'OH',), ('CH3CHCH3', 'H2O',), (None,))])
    ped_400 = ped_df_frag1_dct[1.0][400]
    ped_800 = ped_df_frag1_dct[1.0][800]

    assert np.isclose((ped_400.iloc[166]), 0.00045, atol=1e-4, rtol=1e-2)
    assert np.isclose((ped_800.iloc[166]), 0.0522, atol=1e-3, rtol=1e-2)
    assert np.isclose(np.trapz(ped_400.values, x=ped_400.index), 1)
    assert np.isclose(np.trapz(ped_800.values, x=ped_800.index), 1)


def test_beta_phi3a():
    """ test statmodels.pedmodels.beta_phi3a
    """

    dof_dct, ped_dct, _, _, _ = _read_data()

    ped_df_frag1_dct = ene_partition.ped_frag1(
        ped_dct[(('C3H8', 'OH',), ('CH3CHCH3', 'H2O',), (None,))
                ], 'CH3CHCH3', 'H2O', 'beta_phi3a',
        dof_info=dof_dct[(('C3H8', 'OH',), ('CH3CHCH3', 'H2O',), (None,))])
    ped_400 = ped_df_frag1_dct[1.0][400]
    ped_800 = ped_df_frag1_dct[1.0][800]

    assert np.isclose((ped_400.iloc[166]), 0.001, atol=1e-3, rtol=1e-2)
    assert np.isclose((ped_800.iloc[166]), 0.0589, atol=1e-3, rtol=1e-2)
    assert np.isclose(np.trapz(ped_400.values, x=ped_400.index), 1)
    assert np.isclose(np.trapz(ped_800.values, x=ped_800.index), 1)


def test_rovib_dos():
    """ test statmodels.pedmodels.beta_rovib_dos
    """

    dof_dct, ped_dct, dos_rovib, _, _ = _read_data()

    ped_df_frag1_dct = ene_partition.ped_frag1(
        ped_dct[(('C3H8', 'OH',), ('CH3CHCH3', 'H2O',), (None,))
                ], 'CH3CHCH3', 'H2O', 'rovib_dos',
        dos_df=dos_rovib, dof_info=dof_dct[(
            ('C3H8', 'OH',), ('CH3CHCH3', 'H2O',), (None,))],
    )
    ped_1800 = ped_df_frag1_dct[1.0][1800]
    ped_2000 = ped_df_frag1_dct[1.0][2000]

    assert np.isclose((ped_1800.iloc[166]), 0.01985, atol=1e-3, rtol=1e-2)
    assert np.isclose((ped_2000.iloc[166]), 0.01724, atol=1e-3, rtol=1e-2)
    assert np.isclose(np.trapz(ped_1800.values, x=ped_1800.index), 1)
    assert np.isclose(np.trapz(ped_2000.values, x=ped_2000.index), 1)


def test_thermal():
    """ test statmodels.pedmodels.thermal
    """

    dof_dct, ped_dct, dos_rovib, _, _ = _read_data()

    ped_df_frag1_dct = ene_partition.ped_frag1(
        ped_dct[(('C3H8', 'OH',), ('CH3CHCH3', 'H2O',), (None,))
                ], 'CH3CHCH3', 'H2O', 'thermal',
        dos_df=dos_rovib, dof_info=dof_dct[(
            ('C3H8', 'OH',), ('CH3CHCH3', 'H2O',), (None,))])
    ped_1800 = ped_df_frag1_dct[1.0][1800]
    ped_2000 = ped_df_frag1_dct[1.0][2000]

    assert np.isclose((ped_1800.iloc[166]), 0.01924792, atol=1e-5, rtol=1e-5)
    assert np.isclose((ped_2000.iloc[166]), 0.01132143, atol=1e-5, rtol=1e-5)
    assert np.isclose(np.trapz(ped_1800.values, x=ped_1800.index), 1)
    assert np.isclose(np.trapz(ped_2000.values, x=ped_2000.index), 1)


def test_bf_from_phi1a():
    """ test calculator.bf.bf_tp_dct
        calls calculator.bf.bf_tp_df_full, bf_tp_df_todct
    """

    dof_dct, ped_dct, _, hoten_dct, _ = _read_data()

    ped_df_frag1_dct = ene_partition.ped_frag1(
        ped_dct[(('C3H8', 'OH',), ('CH3CHCH3', 'H2O',), (None,))
                ], 'CH3CHCH3', 'H2O', 'beta_phi1a',
        dof_info=dof_dct[(('C3H8', 'OH',), ('CH3CHCH3', 'H2O',), (None,))])

    bf_tp_dct = bf.bf_tp_dct(
        'beta_phi1a', ped_df_frag1_dct, hoten_dct['CH3CHCH3'], 0.1)

    assert np.allclose(
        bf_tp_dct['CH3CHCH3'][1.0][1],
        np.array([1.,  0.999, 0.999, 0.888]), atol=1e-3, rtol=1e-2)
    assert np.allclose(
        bf_tp_dct['CH3CHCH3'][100.0][1],
        np.array([1., 0.9999, 0.9999, 0.9956, 0.8287, 0.6795]), atol=1e-3, rtol=1e-2)
    assert np.allclose(
        bf_tp_dct['CH3CHCH2+H'][1.0][1][-3:],
        np.array([1.115e-01,
                  9.859e-01, 9.822e-01]), atol=1e-3, rtol=1e-2)
    assert np.allclose(
        bf_tp_dct['CH3CHCH2+H'][100.0][1][-3:],
        np.array([4.325e-03, 1.67e-01,
                  3.1196e-01]), atol=1e-3, rtol=1e-2)


def test_bf_from_fne():
    """ test calculator.bf.bf_tp_dct
        calls calculator.bf.bf_tp_df_full, bf_tp_df_todct
    """
    _, _, _, _, fne_bf = _read_data()
    bf_tp_dct = bf.bf_tp_dct(
        'fne', None, None, 0.1, fne=fne_bf['CH3CHCH3'])

    assert np.allclose(
        bf_tp_dct['CH3CHCH3'][1.0][1][[0, 6, 14]],
        np.array([1., 0.99999987, 0.99798199]), atol=1e-3, rtol=1e-2)
    assert np.allclose(
        bf_tp_dct['CH3CHCH3'][100.0][1][[0, 6, 12, -5, -1]],
        np.array([1., 1., 0.99999238, 0.98542176, 0.94370273]), atol=1e-3, rtol=1e-2)
    assert np.allclose(
        bf_tp_dct['CH3CHCH2+H'][1.0][1][[16, -1, -1]],
        np.array([1.99996390e-03, 9.89494747e-01, 9.87394958e-01]), atol=1e-3, rtol=1e-2)
    assert np.allclose(
        bf_tp_dct['CH3CHCH2+H'][100.0][1][[-5, -1]],
        np.array([1.42060802e-02, 5.47827434e-02]), atol=1e-3, rtol=1e-2)


def test_new_ktp_dct():
    """ test calculator.bf.merge_bf_ktp
        calls calculator.bf.merge_bf_rates
        calls calculator.bf.rename_ktp_dct
    """

    dof_dct, ped_dct, _, hoten_dct, _ = _read_data()

    ped_df_frag1_dct = ene_partition.ped_frag1(
        ped_dct[(('C3H8', 'OH',), ('CH3CHCH3', 'H2O',), (None,))
                ], 'CH3CHCH3', 'H2O', 'equip_simple',
        dof_info=dof_dct[(('C3H8', 'OH',), ('CH3CHCH3', 'H2O',), (None,))])

    bf_tp_dct = bf.bf_tp_dct(
        'equip_simple', ped_df_frag1_dct, hoten_dct['CH3CHCH3'], 0.01)

    rxn_ktp_dct = bf.merge_bf_ktp(
        bf_tp_dct, KTP_DCT[(('C3H8', 'OH',), ('CH3CHCH3', 'H2O',), (None,))],
        (('C3H8', 'OH',), ('CH3CHCH3', 'H2O',), (None,)), HOT_FRAG_DCT)

    rxn1 = (('C3H8', 'OH'), ('CH3CHCH3', 'H2O'), (None,))
    rxn2 = (('C3H8', 'OH'), ('H2O', 'CH3CHCH2', 'H'), (None,))
    rxn3 = (('C3H8', 'OH'), ('H2O', 'C2H4', 'CH3'), (None,))

    assert np.allclose(
        rxn_ktp_dct[rxn1][1.0][1],
        np.array([1.89880000e-13, 4.35037052e-13, 7.71806509e-13, 1.43693675e-12]))
    
    assert np.allclose(
        rxn_ktp_dct[rxn2][1.0][1],
        np.array([2.51727819e-23, 3.94444858e-18, 1.39868122e-15, 2.91984766e-13,
                  3.72843123e-12, 4.51618523e-12]))

    assert np.allclose(
        rxn_ktp_dct[rxn3][1.0][1],
        np.array([5.41734589e-28, 1.39726840e-21, 3.37227513e-18, 2.35907806e-15,
                  6.10487527e-14, 9.39747611e-14]))


def _read_data():
    """ Obtain all the data needed to perform prompt tests
    """

    # get dof info
    spc_blocks_ped = mess_io.reader.get_species(PED_INP)

    dof_dct = {}
    for label in LABELS:
        prods = '+'.join(label[1])
        # NB FCT TESTED IN TEST__CALC_SPINFO
        dof_dct[label] = spinfo_frommess.get_info(
            spc_blocks_ped[prods])
    # GET PED
    ped_dct = mess_io.reader.ped.get_ped(
        PED_OUT, ENERGY_DCT)
    # GET DOS
    dos_rovib = mess_io.reader.rates.dos_rovib(KE_PED_OUT)
    # HOTEN
    hoten_dct = mess_io.reader.hoten.extract_hot_branching(
        HOT_OUT, HOTSPECIES, list(HOT_FRAG_DCT.keys()))
    # FNE
    fne_bf = mess_io.reader.hoten.extract_fne(HOT_OUT)

    return dof_dct, ped_dct, dos_rovib, hoten_dct, fne_bf


if __name__ == '__main__':
    test_equip_simple()
    test_equip_phi()
    test_beta_phi1a()
    test_beta_phi2a()
    test_beta_phi3a()
    test_thermal()
    test_bf_from_phi1a()
    test_bf_from_fne()
    test_new_ktp_dct()
    test_rovib_dos()
