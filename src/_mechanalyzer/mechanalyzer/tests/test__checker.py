""" test mechanalyzer.builder.checker
"""

import numpy as np
from mechanalyzer.builder import checker


TEMPS = np.array([500, 1000, 1500])
GOOD_KTS = ([1e9, 1e10, 1e10])
NEGATIVE_KTS = np.array([-1.5, -1.8, -1.9])
LARGE_UNIMOLEC_KTS = np.array([1e12, 1e12, 1e12])
LARGE_BIMOLEC_KTS = np.array([1e16, 1e16, 1e16])
LARGE_TERMOLEC_KTS = np.array([1e23, 1e23, 1e23])
ARR_PARAMS = [0, None, None, None, None, None]
LIND_PARAMS = [0, 0, None, None, None, None]
TROE_PARAMS = [0, 0, 0, None, None, None]
CHEB_PARAMS = [0, None, None, 0, None, None]
PLOG_PARAMS = [None, None, None, None, 0, None]

# For testing sources and sinks and duplicates
RXN_PARAM_DCT1 = {
    (('H2', 'O'), ('OH', 'H'), (None,)): (ARR_PARAMS, ARR_PARAMS, ARR_PARAMS),
    (('H', 'O2'), ('OH', 'O'), (None,)): (ARR_PARAMS, PLOG_PARAMS),
    (('H2', 'O'), ('OH', 'OH'), (None,)): (LIND_PARAMS, LIND_PARAMS),
    (('H', 'O'), ('OH',), (None,)): (TROE_PARAMS, TROE_PARAMS),
    (('H', 'O'), ('OH',), ('(+M)',)): (CHEB_PARAMS, CHEB_PARAMS),
    (('H', 'O'), ('OH',), ('+O(S)',)): (PLOG_PARAMS, PLOG_PARAMS),
    (('H2', 'O(S)'), ('OH', 'O'), (None,)): (PLOG_PARAMS, PLOG_PARAMS),
    (('H2', 'O2'), ('HO2', 'H'), (None,)): (PLOG_PARAMS, PLOG_PARAMS)
}

# For testing sources and sinks and duplicates
RXN_PARAM_DCT2 = {
    (('H', 'O'), ('OH',), (None,)): (PLOG_PARAMS, PLOG_PARAMS),
    (('OH',), ('H', 'O'), (None,)): (PLOG_PARAMS, PLOG_PARAMS)
}

# For testing negative rates
RXN_KTP_DCT1 = {
    (('H2', 'O'), ('OH', 'H'), (None,)): {
        1: (TEMPS, GOOD_KTS), 10: (TEMPS, GOOD_KTS)},
    (('H', 'O2'), ('OH', 'O'), (None,)): {
        1: (TEMPS, NEGATIVE_KTS), 10: (TEMPS, GOOD_KTS)}
}

# For testing large rates
RXN_KTP_DCT2 = {
    (('OH',), ('H', 'O'), (None,)): {
        1: (TEMPS, GOOD_KTS), 10: (TEMPS, LARGE_UNIMOLEC_KTS)},
    (('HO2',), ('HO', 'O'), (None,)): {
        1: (TEMPS, GOOD_KTS), 10: (TEMPS, GOOD_KTS)},
    (('H2', 'O'), ('OH', 'H'), (None,)): {
        1: (TEMPS, GOOD_KTS), 10: (TEMPS, LARGE_BIMOLEC_KTS)},
    (('H2', 'O'), ('OH', 'OH'), (None,)): {
        1: (TEMPS, GOOD_KTS), 10: (TEMPS, GOOD_KTS)},
    (('H', 'O'), ('OH',), ('+O(S)',)): {
        1: (TEMPS, GOOD_KTS), 10: (TEMPS, LARGE_TERMOLEC_KTS)},
    (('H', 'O'), ('OH',), ('+M',)): {
        1: (TEMPS, GOOD_KTS), 10: (TEMPS, GOOD_KTS)}
}

# For testing large rates
RXN_KTP_DCT3 = {
    (('OH',), ('H', 'O'), (None,)): {
        1: (TEMPS, GOOD_KTS), 10: (TEMPS, LARGE_UNIMOLEC_KTS)},
    (('HO2',), ('HO', 'O'), (None,)): {
        1: (TEMPS, GOOD_KTS), 10: (TEMPS, GOOD_KTS)},
}

SPC_IDENT_DCT1 = {
    'H': {'smiles': '', 'inchi': 'InChI=1S/H', 'inchikey': '', 'mult': 2, 'charge': 0, 'sens': 0,
          'fml': {'H': 1}},
    'OH': {'smiles': '', 'inchi': 'InChI=1S/HO/h1H', 'inchikey': '', 'mult': 2, 'charge': 0,
           'sens': 0, 'fml': {'H': 1, 'O': 1}},
    'O': {'smiles': '', 'inchi': 'InChI=1S/O', 'inchikey': '', 'mult': 3, 'charge': 0, 'sens': 0,
          'fml': {'O': 1}},
    'H2': {'smiles': '', 'inchi': 'InChI=1S/H2/h1H', 'inchikey': '', 'mult': 1, 'charge': 0,
           'sens': 0, 'fml': {'H': 2}},
    'O2': {'smiles': '', 'inchi': 'InChI=1S/O2/c1-2', 'inchikey': '', 'mult': 1, 'charge': 0,
           'sens': 0, 'fml': {'O': 2}},
    'OHV': {'smiles': '', 'inchi': 'InChI=1S/O', 'inchikey': '', 'mult': 1, 'charge': 0,
            'sens': 0, 'fml': {'O': 1}},
}

CORRECT_NEGATIVE_KTS_STR = (
    '\nNEGATIVE RATE CONSTANTS\n\nH + O2 = OH + O\nPressure: 1 atm\n' +
    '    Temperature (K)\n    500.0       1000.0      1500.0      \n' +
    '    Rate constant\n    -1.500E+00  -1.800E+00  -1.900E+00  \n\n\n'
)

CORRECT_LARGE_KTS_STR = (
    '\nLARGE RATE CONSTANTS\n\nUnimolecular threshold: 1.0E+11 s^-1\n' +
    'Bimolecular threshold: 1.0E+15 cm^3 mol^-1 s^-1\n' +
    'Termolecular threshold: 1.0E+22 cm^6 mol^-2 s^-1\n\n' +
    'Unimolecular rate constants that exceed 1.0E+11 s^-1\n\n' +
    'OH = H + O\nPressure: 10 atm\n' +
    '    Temperature (K)\n    500.0       1000.0      1500.0      \n' +
    '    Rate constant\n    1.000E+12   1.000E+12   1.000E+12   \n\n\n' +
    'No bimolecular reactions exceed 1.0E+15 cm^3 mol^-1 s^-1\n\n' +
    'No termolecular reactions exceed 1.0E+22 cm^6 mol^-2 s^-1\n\n'
)

CORRECT_LONE_SPCS_STR = (
    '\nLONE SPECIES\n\n' +
    'These species appear in 2 or less reactions\n\n' +
    'Species  Reactions\nO2       H + O2 = OH + O, H2 + O2 = HO2 + H\n' +
    'O(S)     H2 + O(S) = OH + O\nHO2      ' +
    'H2 + O2 = HO2 + H\n\n\n'
)
CORRECT_SOURCE_SINK_STR1 = (
    '\nSOURCE AND SINK SPECIES\n\nThese species only appear as ' +
    'reactants:\nSpecies    Reactions' +
    '\nH2         H2 + O = OH + H, H2 + O = OH + OH, H2 + O(S) = OH + O, ' +
    'H2 + O2 = HO2 + H\nO(S)       H2 + O(S) = OH + O\nO2         ' +
    'H + O2 = OH + O, H2 + O2 = HO2 + H\n\n' +
    'These species only appear as products:' +
    '\nSpecies   Reactions\nHO2       H2 + O2 = HO2 + H' +
    '\nOH        H2 + O = OH + H, H + O2 = OH + O, H2 + O = OH + OH, H + O = OH, H + O(+M) = OH(+M), ' +
    'H + O+O(S) = OH+O(S), H2 + O(S) = OH + O\n\n\n'
)
CORRECT_SOURCE_SINK_STR2 = (
    '\nSOURCE AND SINK SPECIES\n\nThese species only appear as ' +
    'reactants:\nNo source species found\n\n' +
    'These species only appear as products:\n' +
    'No sink species found\n\n\n'
)
CORRECT_DUPLICATES_STR1 = (
    '\nDUPLICATE REACTIONS\n\nThese reactions have more than 2 ' +
    'rate expressions:\n' +
    '(Number of rate expressions given in parentheses)\n\n' +
    'H2 + O = OH + H     (3)\n\n\n'
)
CORRECT_DUPLICATES_STR2 = (
    '\nDUPLICATE REACTIONS\n\nThese reactions have more than 2 ' +
    'rate expressions:' +
    '\n(Number of rate expressions given in parentheses)\n\n' +
    'No reactions with more than 2 expressions found\n\n\n'
)
CORRECT_MISMATCHES_STR1 = (
    '\nMISMATCHED REACTIONS\n\nThe following reactions have ' +
    'mismatched rate expressions\nH + O2 = OH + O: Arrhenius, PLOG\n\n\n'
)

CORRECT_MISMATCHES_STR2 = (
    '\nMISMATCHED REACTIONS\n\nNo reactions with mismatching ' +
    'rate expressions found\n\n\n'
)

CORRECT_MISSING_SPC_STR = '\nSPECIES MISSING FROM CSV OR MECHANISM\n\n' \
    'These species are missing from the spc.csv file:\n' \
    'HO2\nO(S)\n\n' \
    'These species are missing from the mechanism file:\nOHV\n\n\n'


#def test__all_checks():
#    """ Test the run_all_checks function
#    """
#    k_thresholds = [1e11, 1e15, 1e22]
#    rxn_num_threshold = 2
#    _ = checker.run_all_checks(RXN_PARAM_DCT1, RXN_KTP_DCT1, k_thresholds,
#                               rxn_num_threshold)
#
#
#def test__sources_and_sinks():
#    """ Test the get_sources_and_sinks and write_sources_and_sinks functions
#    """
#    # Test the get_sources_and_sinks function for two different cases
#    sources1, sinks1 = checker.get_sources_and_sinks(RXN_PARAM_DCT1)
#    sources2, sinks2 = checker.get_sources_and_sinks(RXN_PARAM_DCT2)
#    assert list(set(list(sources1.keys())) - set(['O2', 'O(S)', 'H2'])) == []
#    assert list(set(list(sinks1.keys())) - set(['OH', 'HO2'])) == []
#    assert sources2 == sinks2 == {}
#
#    # Test the write_sources_and_sinks function
#    source_sink_str1 = checker.write_sources_and_sinks(sources1, sinks1)
#    source_sink_str2 = checker.write_sources_and_sinks(sources2, sinks2)
#    assert source_sink_str1.replace(" ", "") == CORRECT_SOURCE_SINK_STR1.replace(" ", "")
#    assert source_sink_str2.replace(" ", "") == CORRECT_SOURCE_SINK_STR2.replace(" ", "")


def test__negative_rates():
    """ Test the get_negative_kts and write_negative_kts functions
    """
    # Test the get_negative_kts function
    reac = (('H', 'O2'), ('OH', 'O'), (None,))

    negative_rxn_ktp_dct = checker.get_negative_kts(RXN_KTP_DCT1)
    assert tuple(negative_rxn_ktp_dct.keys()) == (
        (('H', 'O2'), ('OH', 'O'), (None,)),)
    assert tuple(negative_rxn_ktp_dct[reac].keys()) == (1,)

    # Test the write_negative_kts function
    negative_kts_str = checker.write_negative_kts(negative_rxn_ktp_dct)
    assert negative_kts_str.replace(" ", "") == CORRECT_NEGATIVE_KTS_STR.replace(" ", "")


def test__large_rates():
    """ Test the get_large_kts and write_large_kts functions
    """
    # Test the get_large_kts function
    thresholds = [1e11, 1e15, 1e22]
    large_rxn_ktp_dcts2 = checker.get_large_kts(RXN_KTP_DCT2, thresholds)
    assert tuple(large_rxn_ktp_dcts2[0].keys()) == (
        (('OH',), ('H', 'O'), (None,)),)
    assert tuple(large_rxn_ktp_dcts2[1].keys()) == (
        (('H2', 'O'), ('OH', 'H'), (None,)),)
    assert tuple(large_rxn_ktp_dcts2[2].keys()) == (
        (('H', 'O'), ('OH',), ('+O(S)',)),)

    # Test the write_large_kts function
    large_rxn_ktp_dcts3 = checker.get_large_kts(RXN_KTP_DCT3, thresholds)
    large_kts_str3 = checker.write_large_kts(large_rxn_ktp_dcts3, thresholds)
    assert large_kts_str3.replace(" ", "") == CORRECT_LARGE_KTS_STR.replace(" ", "")


def test__lone_species():
    """ Test the get_lone_spcs and write_lone_spcs functions
    """
    # Test the get_lone_spcs function
    threshold = 2
    lone_spcs = checker.get_lone_spcs(RXN_PARAM_DCT1, threshold)
    assert tuple(lone_spcs.keys()) == tuple(['O2', 'O(S)', 'HO2'])

    # Test the write_lone_spcs function
    lone_spcs_str = checker.write_lone_spcs(lone_spcs, threshold)
    assert lone_spcs_str.replace(" ", "") == CORRECT_LONE_SPCS_STR.replace(" ", "")


def test__duplicates():
    """ Test the get_duplicates and write_duplicates functions
    """
    # Test the get_duplicates function
    duplicate_rxns1 = checker.get_duplicates(RXN_PARAM_DCT1)
    duplicate_rxns2 = checker.get_duplicates(RXN_PARAM_DCT2)
    assert tuple(duplicate_rxns1.keys()) == (
        (('H2', 'O'), ('OH', 'H'), (None,)),)
    assert tuple(duplicate_rxns1.values()) == (3,)
    assert duplicate_rxns2 == {}

    # Test the write_duplicates function
    dup_str1 = checker.write_duplicates(duplicate_rxns1)
    dup_str2 = checker.write_duplicates(duplicate_rxns2)
    assert dup_str1.replace(" ", "") == CORRECT_DUPLICATES_STR1.replace(" ", "")
    assert dup_str2.replace(" ", "") == CORRECT_DUPLICATES_STR2.replace(" ", "")


#def test__mismatches():
#    """ Test the get_mismatches and write_mismatches functions
#    """
#    # Test the get_mismatches function
#    mismatched_rxns1 = checker.get_mismatches(RXN_PARAM_DCT1)
#    mismatched_rxns2 = checker.get_mismatches(RXN_PARAM_DCT2)
#    assert tuple(mismatched_rxns1.keys()) == (
#        (('H', 'O2'), ('OH', 'O'), (None,)),)
#    assert tuple(mismatched_rxns1.values())[0][1] == ['Arrhenius', 'PLOG']
#    assert mismatched_rxns2 == {}
#
#    # Test the write_mismatches function
#    mismatch_str1 = checker.write_mismatches(mismatched_rxns1)
#    mismatch_str2 = checker.write_mismatches(mismatched_rxns2)
#    assert mismatch_str1.replace(" ", "") == CORRECT_MISMATCHES_STR1.replace(" ", "")
#    assert mismatch_str2.replace(" ", "") == CORRECT_MISMATCHES_STR2.replace(" ", "")


def test__missing_spcs():
    """ Test the get_mismatches and write_mismatches functions
    """

    missing_from_csv, missing_from_mech = checker.get_missing_spcs(
        RXN_PARAM_DCT1, SPC_IDENT_DCT1)
    # Sort so the test always passes
    missing_from_csv = sorted(missing_from_csv)
    assert set(missing_from_csv) == set(['HO2', 'O(S)'])
    assert set(missing_from_mech) == set(['OHV'])

    missing_spcs_str = checker.write_missing_spcs(missing_from_csv,
                                                  missing_from_mech)
    assert missing_spcs_str.replace(" ", "") == CORRECT_MISSING_SPC_STR.replace(" ", "")


if __name__ == '__main__':
    test__all_checks()
    test__sources_and_sinks()
    test__negative_rates()
    test__large_rates()
    test__lone_species()
    test__duplicates()
    test__mismatches()
    test__missing_spcs()
