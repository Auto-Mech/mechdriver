""" Test the compare.py function
"""

from mechanalyzer.calculator import compare
import numpy as np

TEMPS = np.array([500, 1000, 1500])
KTS = np.array([1e10, 1e11, 1e12])
h = np.array([-5000, -4000, -3000])
cp = np.array([-5000, -4000, -3000])
s = np.array([-5000, -4000, -3000])

SPC_IDENT_DCT1 = {
    'H': {'smiles': '', 'inchi': 'InChI=1S/H', 'inchikey': '', 'mult': 2, 'charge': 0, 'exc_flag': 0,
          'fml': {'H': 1}},
    'OH': {'smiles': '', 'inchi': 'InChI=1S/HO/h1H', 'inchikey': '', 'mult': 2, 'charge': 0,
           'exc_flag': 0, 'fml': {'H': 1, 'O': 1}},
    'O': {'smiles': '', 'inchi': 'InChI=1S/O', 'inchikey': '', 'mult': 3, 'charge': 0, 'exc_flag': 0,
          'fml': {'O': 1}},
    'H2': {'smiles': '', 'inchi': 'InChI=1S/H2/h1H', 'inchikey': '', 'mult': 1, 'charge': 0,
           'exc_flag': 0, 'fml': {'H': 2}},
    'O2': {'smiles': '', 'inchi': 'InChI=1S/O2/c1-2', 'inchikey': '', 'mult': 1, 'charge': 0,
           'exc_flag': 0, 'fml': {'O': 2}},
    'O(S)': {'smiles': '', 'inchi': 'InChI=1S/O', 'inchikey': '', 'mult': 1, 'charge': 0,
             'exc_flag': 0, 'fml': {'O': 1}},
}

SPC_IDENT_DCT2 = {
    'HV': {'smiles': '', 'inchi': 'InChI=1S/H', 'inchikey': '', 'mult': 2, 'charge': 0, 'exc_flag': 0,
           'fml': {'H': 1}},
    'OHV': {'smiles': '', 'inchi': 'InChI=1S/HO/h1H', 'inchikey': '', 'mult': 2, 'charge': 0,
            'exc_flag': 0, 'fml': {'H': 1, 'O': 1}},
    'OV': {'smiles': '', 'inchi': 'InChI=1S/O', 'inchikey': '', 'mult': 3, 'charge': 0, 'exc_flag': 0,
           'fml': {'O': 1}},
    'O2V': {'smiles': '', 'inchi': 'InChI=1S/O2/c1-2', 'inchikey': '', 'mult': 1, 'charge': 0,
            'exc_flag': 0, 'fml': {'O': 2}},
    'H2V': {'smiles': '', 'inchi': 'InChI=1S/H2/h1H', 'inchikey': '', 'mult': 1, 'charge': 0,
            'exc_flag': 0, 'fml': {'H': 2}},
    'HO2V': {'smiles': '', 'inchi': 'InChI=1S/HO2/c1-2/h1H', 'inchikey': '', 'mult': 2, 'charge': 0,
             'exc_flag': 0, 'fml': {'H': 1, 'O': 2}},
    'O(S)V': {'smiles': '', 'inchi': 'InChI=1S/O', 'inchikey': '', 'mult': 1, 'charge': 0,
              'exc_flag': 0, 'fml': {'O': 1}},
}

SPC_IDENT_DCT3 = {
    'HO2X': {'smiles': '', 'inchi': 'InChI=1S/HO2/c1-2/h1H', 'inchikey': '', 'mult': 2, 'charge': 0,
             'exc_flag': 0, 'fml': {'H': 1, 'O': 2}},
    'O2X': {'smiles': '', 'inchi': 'InChI=1S/O2/c1-2', 'inchikey': '', 'mult': 1, 'charge': 0,
            'exc_flag': 0, 'fml': {'O': 2}},
    'H2O': {'smiles': '', 'inchi': 'InChI=1S/H2O/h1H2', 'inchikey': '', 'mult': 1, 'charge': 0,
            'exc_flag': 0, 'fml': {'H': 2, 'O': 1}},
    'H2': {'smiles': '', 'inchi': 'InChI=1S/H2/h1H', 'inchikey': '', 'mult': 2, 'charge': 0,
           'exc_flag': 0, 'fml': {'H': 2}},
}

SPC_THERM_DCT1 = {
    'H': [TEMPS, h, cp, s, np.array([38112.076194796915, 22159.224380360094, 4906.62142234519])],
    'OH': [TEMPS, h, cp, s, np.array([-13437.71616429529, -38594.34206874827, -65664.834890985])],
    'O': [TEMPS, h, cp, s, np.array([40013.92863513222, 18462.12185893287, -4399.717127561729])],
    'H2': [TEMPS, h, cp, s, np.array([-16011.013763631481, -34787.09040495749, -55449.572706248])],
    'O2': [TEMPS, h, cp, s, np.array([-24919.542351752494, -52791.48721222688, -82818.01126124])],
    'O(S)': [TEMPS, h, cp, s, np.array([41013.92863513222, 19462.12185893287, -4899.7171261729])],
}

SPC_THERM_DCT2 = {
    'HV': [TEMPS, h, cp, s, np.array([38112.076194796915, 22159.224380360094, 4906.621422734519])],
    'OHV': [TEMPS, h, cp, s, np.array([-13437.71616429529, -38594.34206874827, -65664.80890985])],
    'OV': [TEMPS, h, cp, s, np.array([40013.92863513222, 18462.12185893287, -4399.717127561729])],
    'O2V': [TEMPS, h, cp, s, np.array([-24919.542351752494, -52791.48721222688, -82818.01121114])],
    'H2V': [TEMPS, h, cp, s, np.array([-16011.013763631481, -34787.09040495749, -55449.57270738])],
    'HO2V': [TEMPS, h, cp, s, np.array([-24930.759856779678, -56555.863471982615, -91105.70559])],
    'O(S)V': [TEMPS, h, cp, s, np.array([41013.92863513222, 19462.12185893287, -4899.712756729])],
}

RXN_KTP_DCT1 = {
    (('H2', 'O'), ('OH', 'H'), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('H', 'O2'), ('OH', 'O'), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('H2', 'O'), ('OH', 'OH'), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('H', 'O'), ('OH',), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('H', 'O'), ('OH',), ('(+M)',)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('H', 'O'), ('OH',), ('+O(S)',)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('H2', 'O(S)'), ('OH', 'O'), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS),
                                             10: (TEMPS, KTS)},
}

RXN_KTP_DCT2 = {
    (('H2V', 'OV'), ('OHV', 'HV'), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS),
                                              10: (TEMPS, KTS)},
    (('OV', 'OHV'), ('O2V', 'HV'), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS),
                                              10: (TEMPS, KTS)},
    (('H2V', 'O2V'), ('HO2V', 'HV'), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS),
                                                10: (TEMPS, KTS)},
    (('OHV',), ('HV', 'OV'), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('OHV',), ('HV', 'OV'), ('(+M)',)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('OHV',), ('HV', 'OV'), ('+O(S)V',)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS),
                                            10: (TEMPS, KTS)},
    (('H2V', 'O(S)V'), ('OV', 'OHV'), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS),
                                                 10: (TEMPS, KTS)},
}

RXN_KTP_DCT3 = {
    (('H', 'O2'), ('OH', 'O'), (None,)): {'high': (TEMPS, np.array([2.112e7, 4.682e10, 6.111e11]))},
}
RXN_KTP_DCT4 = {
    (('OH', 'O'), ('H', 'O2'), (None,)): {'high': (TEMPS, np.array([1.495e13, 9.23e12, 8.5e12]))},
}

CORRECT_SPC_KEYS = ('H', 'OH', 'O', 'O2', 'H2', 'HO2V', 'O(S)')
CORRECT_COMB_SPC_KEYS = ('H', 'OH', 'O', 'H2', 'O2', 'O(S)', 'HO2V')
CORRECT_COMB_SPC_KEYS2 = ('H', 'OH', 'O', 'H2', 'O2', 'O(S)', 'HO2V', 'H2O', 'H2-zz')
CORRECT_RENAMED_RXN_KEYS = (
    (('H2', 'O'), ('OH', 'H'), (None,)),
    (('O', 'OH'), ('O2', 'H'), (None,)),
    (('H2', 'O2'), ('HO2V', 'H'), (None,)),
    (('OH',), ('H', 'O'), (None,)),
#    (('OH',), ('H', 'O'), ('(+M)',)),
    (('OH',), ('H', 'O'), ('+O(S)',)),
    (('H2', 'O(S)'), ('O', 'OH'), (None,))
)

# For the test case of when rev_rates=True
# (These are seemingly out of order due to the fact that they are popped out and then readded)
CORRECT_REVERSED_RXN_KEYS = (
    (('H2', 'O'), ('OH', 'H'), (None,)),
    (('H2', 'O2'), ('HO2V', 'H'), (None,)),
    (('H', 'O2'), ('OH', 'O'), (None,)),
    (('H', 'O'), ('OH',), (None,)),
#    (('H', 'O'), ('OH',), ('(+M)',)),
    (('H', 'O'), ('OH',), ('+O(S)',)),
    (('H2', 'O(S)'), ('OH', 'O'), (None,))
)

# For the test case of when rev_rates=False
# (Only the last reaction has had the products reordered)
CORRECT_PARTIALLY_REVERSED_RXN_KEYS = (
    (('H2', 'O'), ('OH', 'H'), (None,)),
    (('O', 'OH'), ('O2', 'H'), (None,)),
    (('H2', 'O2'), ('HO2V', 'H'), (None,)),
    (('OH',), ('H', 'O'), (None,)),
#    (('OH',), ('H', 'O'), ('(+M)',)),
    (('OH',), ('H', 'O'), ('+O(S)',)),
    (('H2', 'O(S)'), ('OH', 'O'), (None,))
)

# Without removing loners
CORRECT_ALGN_RXN_KTP_KEYS = (
    (('H2', 'O'), ('OH', 'H'), (None,)),
    (('H', 'O2'), ('OH', 'O'), (None,)),
    (('H2', 'O'), ('OH', 'OH'), (None,)),
    (('H', 'O'), ('OH',), (None,)),
    (('H', 'O'), ('OH',), ('(+M)',)),
    (('H', 'O'), ('OH',), ('+O(S)',)),
    (('H2', 'O(S)'), ('OH', 'O'), (None,)),
    (('H2', 'O2'), ('HO2V', 'H'), (None,))
)

# With removing loners
CORRECT_ALGN_RXN_KTP_KEYS_NO_LONERS = (
    (('H2', 'O'), ('OH', 'H'), (None,)),
    (('H', 'O2'), ('OH', 'O'), (None,)),
    (('H', 'O'), ('OH',), (None,)),
#    (('H', 'O'), ('OH',), ('(+M)',)),
    (('H', 'O'), ('OH',), ('+O(S)',)),
    (('H2', 'O(S)'), ('OH', 'O'), (None,))
)


def test_rename_spc_dct():
    """ Test the renaming functions
    """
    rename_instr = compare.get_rename_instr(SPC_IDENT_DCT1, SPC_IDENT_DCT2)
    renamed_dct, _ = compare.rename_species(SPC_IDENT_DCT2, rename_instr, target_type='spc')
    assert tuple(renamed_dct.keys()) == CORRECT_SPC_KEYS
    assert tuple(renamed_dct.values()) == tuple(SPC_IDENT_DCT2.values())


def test_get_comb_spc_dct():
    """ Test the get_comb_spc_dct function
    """
    comb_spc_dct = compare.get_comb_mech_spc_dct(SPC_IDENT_DCT1, SPC_IDENT_DCT2)
    comb_spc_dct2 = compare.get_mult_comb_mech_spc_dct(
        [SPC_IDENT_DCT1, SPC_IDENT_DCT2, SPC_IDENT_DCT3])
    assert tuple(comb_spc_dct) == CORRECT_COMB_SPC_KEYS
    assert tuple(comb_spc_dct2) == CORRECT_COMB_SPC_KEYS2


def test_rename_spc_therm_dct():
    """ Test the rename_spc_therm_dct function
    """
    rename_instr = compare.get_rename_instr(SPC_IDENT_DCT1, SPC_IDENT_DCT2)
    renamed_dct, _ = compare.rename_species(SPC_THERM_DCT2, rename_instr, target_type='spc')
    assert tuple(renamed_dct.keys()) == CORRECT_SPC_KEYS


def test_rename_rxn_ktp_dct():
    """ Tese the rename rxn_ktp_dct function
    """
    rename_instr = compare.get_rename_instr(SPC_IDENT_DCT1, SPC_IDENT_DCT2)
    renamed_dct, _ = compare.rename_species(RXN_KTP_DCT2, rename_instr, target_type='rxn')
    assert tuple(renamed_dct.keys()) == CORRECT_RENAMED_RXN_KEYS


def test_rename_dcts():
    """ Test the rename_dcts function
    """
    rxn_ktp_dcts = [RXN_KTP_DCT1, RXN_KTP_DCT2]
    spc_therm_dcts = [SPC_THERM_DCT1, SPC_THERM_DCT2]
    spc_dcts = [SPC_IDENT_DCT1, SPC_IDENT_DCT2]

    # Rename the rxn_ktp and spc_thermo dcts
    renamed_rxn_ktp_dcts, _ = compare.rename_dcts(rxn_ktp_dcts, spc_dcts, target_type='rxn')
    renamed_spc_therm_dcts, _ = compare.rename_dcts(spc_therm_dcts, spc_dcts,
                                                    target_type='spc')
    assert tuple(renamed_rxn_ktp_dcts[1].keys()) == CORRECT_RENAMED_RXN_KEYS
    assert tuple(renamed_spc_therm_dcts[1].keys()) == CORRECT_SPC_KEYS


def test_reverse_rxn_ktp_dcts():
    """ Test the reverse_rxn_ktp_dcts function
    """
    rxn_ktp_dcts = [RXN_KTP_DCT1, RXN_KTP_DCT2]
    spc_therm_dcts = [SPC_THERM_DCT1, SPC_THERM_DCT2]
    spc_dcts = [SPC_IDENT_DCT1, SPC_IDENT_DCT2]

    # Rename the rxn_ktp and spc_thermo dcts
    renamed_rxn_ktp_dcts, _ = compare.rename_dcts(rxn_ktp_dcts, spc_dcts,
                                                  target_type='rxn')
    renamed_spc_therm_dcts, _ = compare.rename_dcts(spc_therm_dcts, spc_dcts,
                                                    target_type='spc')
    # Check functionality of the full reversal
    reversed_rxn_ktp_dcts = compare.reverse_rxn_ktp_dcts(
        renamed_rxn_ktp_dcts, renamed_spc_therm_dcts, TEMPS, rev_rates=True
    )

    # Check functionality of flipping the ordering of reactants and/or products but NOT reversing
    partially_reversed_rxn_ktp_dcts = compare.reverse_rxn_ktp_dcts(
        renamed_rxn_ktp_dcts, renamed_spc_therm_dcts, TEMPS, rev_rates=False
    )

    # Check that keys are correct after being reversed
    assert tuple(reversed_rxn_ktp_dcts[0].keys()) == tuple(RXN_KTP_DCT1.keys())
    assert tuple(reversed_rxn_ktp_dcts[1].keys()) == CORRECT_REVERSED_RXN_KEYS
    assert tuple(partially_reversed_rxn_ktp_dcts[1].keys()) == CORRECT_PARTIALLY_REVERSED_RXN_KEYS

    # Check a numerical case to make sure rates are getting reversed correctly
    # (using the rates from Hong for H+O2=OH+O)
    hong_case = compare.reverse_rxn_ktp_dcts(
        [RXN_KTP_DCT3, RXN_KTP_DCT4], renamed_spc_therm_dcts, TEMPS, rev_rates=True
    )
    assert np.allclose(
        hong_case[0][(('H', 'O2'), ('OH', 'O'), (None,))]['high'][1],
        hong_case[1][(('H', 'O2'), ('OH', 'O'), (None,))]['high'][1], rtol=1e-3
        )


def test_align_rxn_ktp_dcts():
    """ Test the align_rxn_ktp_dcts function
    """
    rxn_ktp_dcts = [RXN_KTP_DCT1, RXN_KTP_DCT2]
    spc_therm_dcts = [SPC_THERM_DCT1, SPC_THERM_DCT2]
    spc_dcts = [SPC_IDENT_DCT1, SPC_IDENT_DCT2]

    # Try without removing loners
    algn_rxn_ktp_dct = compare.get_algn_rxn_ktp_dct(
        rxn_ktp_dcts, spc_therm_dcts, spc_dcts, TEMPS, rev_rates=True, remove_loners=False,
        write_file=False
    )

    # Try with removing loners
    no_loners_algn_rxn_ktp_dct = compare.get_algn_rxn_ktp_dct(
        rxn_ktp_dcts, spc_therm_dcts, spc_dcts, TEMPS, rev_rates=True, remove_loners=True,
        write_file=False
    )

    # Check that the keys are written correctly
    assert tuple(algn_rxn_ktp_dct.keys()) == CORRECT_ALGN_RXN_KTP_KEYS
    assert tuple(no_loners_algn_rxn_ktp_dct.keys()) == CORRECT_ALGN_RXN_KTP_KEYS_NO_LONERS

    # Check that the rates were correctly reversed/not reversed for one reaction
    # Not reversed, so should be the same as KTS
    assert np.allclose(
        algn_rxn_ktp_dct[(('H', 'O2'), ('OH', 'O'), (None,))][0]['high'][1], KTS
        )
    # Reversed, so should *not* be the same as KTS
    assert not np.allclose(
        algn_rxn_ktp_dct[(('H', 'O2'), ('OH', 'O'), (None,))][1]['high'][1], KTS
        )


if __name__ == '__main__':
    test_rename_spc_dct()
    test_get_comb_spc_dct()
    test_rename_spc_therm_dct()
    test_rename_rxn_ktp_dct()
    test_rename_dcts()
    test_reverse_rxn_ktp_dcts()
    test_align_rxn_ktp_dcts()
