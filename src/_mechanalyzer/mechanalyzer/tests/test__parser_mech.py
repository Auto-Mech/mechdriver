""" test mechanalyzer.parser.mech
"""

import os
import pandas as pd
from ioformat import pathtools
from mechanalyzer.parser.mech import parse_classtype
from mechanalyzer.parser.mech import parse_sort
from mechanalyzer.parser.mech import parse_scalefactors_byclass

CWD = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(CWD, 'data')

OSCLASS_DCT_PARTIAL = {'INCEPTION_PAH-R': 'INCEPTION', 'INCEPTION_PAH-M': 'INCEPTION', 'INCEPTION_PAH-RSR': 'INCEPTION',
               'BONDFISSION_CO': 'BONDFISSION', 'MOLECULAR_DECO': 'DECOMPOSITION', 'RADICAL_DECO': 'DECOMPOSITION',
               'EL_H2': 'DECOMPOSITION', 'EL_H2O': 'DECOMPOSITION', 'ROPEN': 'DECOMPOSITION', 'BSCISS_CO': 'BETASCISSION',
               'ADD_O2_EL_HCO': 'ADD_O2', 'ADD_O2_DECO': 'ADD_O2', 'ADD_O2_ROPEN': 'ADD_O2', 'REC_O2_WELL': 'ADD_O2', 'REC_O2_EL_O': 'ADD_O2',
               'REC_O2_INSERTION': 'ADD_O2', 'REC_O2_ROPEN': 'ADD_O2', 'REC_O2_EL_OH': 'ADD_O2', 'REC_O2_EL_HO2': 'ADD_O2',
               'REC_O2_ISOM_O': 'ADD_O2', 'RO2_INSERTION': 'ADD_O2', 'RO2_EL_OH': 'ADD_O2', 'RO2_ISOM_O': 'ADD_O2', 'RO2_REC_H': 'RO2_SECONDARY',
               'RO2_DISP': 'RO2_SECONDARY', 'RO2_ABS': 'RO2_SECONDARY', 'ADD_OH_WELL': 'ADD_OH', 'ADD_OH_DECO': 'ADD_OH', 'ADD_OH_ROPEN': 'ADD_OH',
               'ADD_C2.T-M': 'GROWTH-C2', 'ADD_C2.T-R': 'GROWTH-C2', 'REC_C2.D-R': 'GROWTH-C2', 'REC_C2.T-R': 'GROWTH-C2', 'ENLARGE_C2.D-M': 'GROWTH-C2', 'ENLARGE_C2.T-M': 'GROWTH-C2',
               'ADD_C4.DT-R': 'GROWTH-C4', 'ADD_C4.D-R': 'GROWTH-C4', 'ADD_C4.D-RSR': 'GROWTH-C4', 'REC_C4.DD-R': 'GROWTH-C4',
               'REC_C4.DT-R': 'GROWTH-C4', 'ADD_C5-M': 'GROWTH-C5', 'ADD_C5O-M': 'GROWTH-C5', 'ADD_C5-RSR': 'GROWTH-C5',
               'ENLARGE_C5-RSR': 'GROWTH-C5', 'ADD_A1-M': 'GROWTH-C6', 'ADD_A1-R': 'GROWTH-C6', 'REC_A1-R': 'GROWTH-C6',
               'REC_A1O-RSR': 'GROWTH-C6', 'ENLARGE': 'GROWTH-UNIMOL', 'ENLARGE_H': 'GROWTH-H', 'UNSORTED': 'UNSORTED'}

def test__sort_readinput():
    """ test reader criteria for sorter with submech
    """
    # Read the mechanism files into strings

    try:
        sort_str = pathtools.read_file(DAT_PATH, 'sort_filterstoich_wrong.dat')
        isolate_spc, sort_lst, _ = parse_sort(sort_str)
    except ValueError as e:
        assert str(e) == 'Cannot have both keepbelow and deleteabove criteria - incompatible!'

    sort_str = pathtools.read_file(DAT_PATH, 'sort_filterstoich.dat')
    isolate_spc, sort_lst, _ = parse_sort(sort_str)
    assert isolate_spc == ['keepbelow C2H6O2']
    assert sort_lst == ['subpes', 0]

    sort_str = pathtools.read_file(DAT_PATH, 'sort_submech_deletelarge.dat')
    isolate_spc, sort_lst, _ = parse_sort(sort_str)
    assert isolate_spc == ['C2H4', 'deleteabove C3H4O2']
    assert sort_lst == ['submech_deletelarge', 0]

    sort_str = pathtools.read_file(DAT_PATH, 'sort_submech_keepsubfuel.dat')
    isolate_spc, sort_lst, _ = parse_sort(sort_str)
    assert isolate_spc == ['C2H4', 'keepbelow C2H6O2']
    assert sort_lst == ['submech_keepsubfuel', 0]

    sort_str = pathtools.read_file(DAT_PATH, 'sort_singlespecies.dat')
    isolate_spc, sort_lst, _ = parse_sort(sort_str)
    assert isolate_spc == ['C2H4', 'singlespecies']
    assert sort_lst == ['subpes', 'molecularity', 'rxn_class_broad', 0]

    sort_str = pathtools.read_file(DAT_PATH, 'sort.dat')
    isolate_spc, sort_lst, _ = parse_sort(sort_str)
    assert isolate_spc == ['C2H4']
    assert sort_lst == ['subpes', 'molecularity', 'rxn_class_broad', 0]

    sort_str = pathtools.read_file(DAT_PATH, 'sort_prompt.dat')
    isolate_spc, sort_lst, prompt_filter_dct = parse_sort(sort_str)
    dct_check = {
        'Tref': 1500.0, 'H5H3ratio': 0.0, 'DH': 0.0, 'kratio': 1e+50, 'kabs': 1e+50, 'keepfiltered': 0.0, 'lookforpromptchains': 0.0
    }
    for key, val in dct_check.items():
        assert prompt_filter_dct[key] == val

def test__parse_classtype():
    class_str = pathtools.read_file(os.path.join(DAT_PATH, 'osclass'), 'rxn_class_type.txt', remove_comments = '#')
    classtype_dct = parse_classtype(class_str)

    for key, val in OSCLASS_DCT_PARTIAL.items():
        assert val == classtype_dct[key]
"""
def test__parse_scalefactors():
    scalefactors_df = parse_scalefactors_byclass(os.path.join(DAT_PATH, 'AAA_SCALE_FACTORS.xlsx'))
    expected_output = pd.DataFrame([
        {"SPECIES": "C10H6CH3", "SPECIESTYPE": "A1-M", "REACTIONTYPE": "HABS_H",    "A": 0.70, "n": 0.0,
            "EA": 0.0,     "REFSPECIES": '', "CORRECTIONSOURCE": "CALC, LUMP", "COMMENT": ''},
        {"SPECIES": "C10H6CH3", "SPECIESTYPE": "A1-M", "REACTIONTYPE": "HABS_OH",   "A": 0.70, "n": 0.0,
            "EA": 0.0,     "REFSPECIES": '', "CORRECTIONSOURCE": "CALC, LUMP", "COMMENT": ''},
        {"SPECIES": "C10H6CH3", "SPECIESTYPE": "A1-M", "REACTIONTYPE": "HABS_O2",   "A": 0.70, "n": 0.0,
            "EA": 0.0,     "REFSPECIES": '', "CORRECTIONSOURCE": "CALC, LUMP", "COMMENT": ''},
        {"SPECIES": "C10H6CH3", "SPECIESTYPE": "A1-M", "REACTIONTYPE": "HABS_HO2",  "A": 0.70, "n": 0.0,
            "EA": 0.0,     "REFSPECIES": '', "CORRECTIONSOURCE": "CALC, LUMP", "COMMENT": ''},
        {"SPECIES": "C10H6CH3", "SPECIESTYPE": "A1-M", "REACTIONTYPE": "HABS_CH3",  "A": 0.70, "n": 0.0,
            "EA": 0.0,     "REFSPECIES": '', "CORRECTIONSOURCE": "CALC, LUMP", "COMMENT": ''},
        {"SPECIES": "C10H7CH3", "SPECIESTYPE": "A1,CH3-M", "REACTIONTYPE": "HABS_OH",   "A": 0.40,
            "n": 0.0, "EA": 0.0, "REFSPECIES": '', "CORRECTIONSOURCE": "CALC", "COMMENT": ''},
        {"SPECIES": "C10H7CH3", "SPECIESTYPE": "A1,CH3-M", "REACTIONTYPE": "HABS_CH3",  "A": 0.67,
            "n": 0.0, "EA": 0.0, "REFSPECIES": '', "CORRECTIONSOURCE": "CALC", "COMMENT": ''},
        {"SPECIES": "C10H7CH3", "SPECIESTYPE": "A1-M",     "REACTIONTYPE": "IPSO_O",    "A": 0.50,
            "n": 0.0, "EA": 1000.0, "REFSPECIES": '', "CORRECTIONSOURCE": "EST", "COMMENT": ''},
        {"SPECIES": "C10H7CH3", "SPECIESTYPE": "A1,CH3-M", "REACTIONTYPE": "IPSO_O",    "A": 0.50,
            "n": 0.0, "EA": 0.0, "REFSPECIES": '', "CORRECTIONSOURCE": "GEOM", "COMMENT": ''},
        {"SPECIES": "C10H7CH3", "SPECIESTYPE": "A1-M",     "REACTIONTYPE": "ADD_O_ISC_DECO", "A": 0.50,
            "n": 0.0, "EA": 0.0, "REFSPECIES": "C7H8", "CORRECTIONSOURCE": "GEOM", "COMMENT": ''},
    ])
    expected = expected_output.reset_index(drop=True)

    test_slice = scalefactors_df.iloc[8:18].reset_index(drop=True)
    test_slice["COMMENT"] = test_slice["COMMENT"].astype(object)
    pd.testing.assert_frame_equal(test_slice, expected)

    # Reset index of expected to match
"""

if __name__ == '__main__':
    test__sort_readinput()
    test__parse_classtype()
    #test__parse_scalefactors()


