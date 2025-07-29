""" test mechanalyzer.parser.sort for different mechanisms in 'data/'
    using different sorting options
"""
import os
import tempfile
import numpy as np
from ioformat import pathtools
from mechanalyzer.builder import sorter
from mechanalyzer.parser import mech as mparser
from mechanalyzer.parser import ckin_ as ckin_parser
from mechanalyzer.parser import new_spc as sparser

# Set Paths to test/data directory and output directory
CWD = os.path.dirname(os.path.realpath(__file__))
TMP_OUT = tempfile.mkdtemp()

# Set types for parsing mechanisms
SPC_TYPE = 'csv'
MECH_TYPE = 'chemkin'

def test__sort_submech():
    """ submech with and without "singlespecies" option
    """
    results = [[(('C2H4',), ('H2', 'H2CC'), ('(+M)',)),'! pes.subpes.NR.rxntype  1.1.1.Decomposition'],
        [(('C2H3', 'H'), ('C2H4',), ('(+M)',)),'! pes.subpes.NR.rxntype  1.1.2.Recombination H'],
        [(('C2H4', 'H'), ('C2H5',), ('(+M)',)),'! pes.subpes.NR.rxntype  2.1.2.Addition H'],
        [(('C2H4', 'H'), ('C2H3', 'H2'), (None,)),'! pes.subpes.NR.rxntype  2.2.2.H abstraction'],
        [(('CH2(S)', 'CH3'), ('C2H4', 'H'), (None,)),'! pes.subpes.NR.rxntype  2.3.2.Addition-decomposition - propagation'],
        [(('C2H5', 'H'), ('C2H6',), ('(+M)',)),'! pes.subpes.NR.rxntype  3.1.2.Recombination H'],
        [(('C2H5', 'H'), ('C2H4', 'H2'), (None,)),'! pes.subpes.NR.rxntype  3.2.2.Recombination-decomposition - termination'],
        [(('CH3', 'CH3'), ('H', 'C2H5'), (None,)),'! pes.subpes.NR.rxntype  3.3.2.Recombination-decomposition - propagation'],
        [(('C2H6', 'H'), ('C2H5', 'H2'), (None,)),'! pes.subpes.NR.rxntype  4.1.2.H abstraction'],
        [(('C2H4', 'O'), ('CH3', 'HCO'), (None,)),'! pes.subpes.NR.rxntype  5.1.2.Addition-decomposition - branching'],
        [(('C2H4', 'OH'), ('PC2H4OH',), (None,)),'! pes.subpes.NR.rxntype  6.1.2.Addition OH'],
        [(('C2H5', 'OH'), ('C2H4', 'H2O'), (None,)),'! pes.subpes.NR.rxntype  7.1.2.Recombination-decomposition - termination'],
        [(('C2H4', 'O2'), ('C2H3', 'HO2'), (None,)),'! pes.subpes.NR.rxntype  9.1.2.H abstraction'],
        [(('C2H3OO',), ('CH2O', 'HCO'), (None,)),'! pes.subpes.NR.rxntype  8.1.1.Beta-scission'],
        [(('C2H3OO',), ('CH2O', 'H', 'CO'), (None,)),'! pes.subpes.NR.rxntype  8.1.1.Decomposition(lumped)'],
        [(('C2H3', 'O2'), ('C2H3OO',), (None,)),'! pes.subpes.NR.rxntype  8.1.2.Recombination O2'],
        [(('HOCH2CO',), ('CH2OH', 'CO'), (None,)),'! pes.subpes.NR.rxntype  8.2.1.Decomposition'],
        [(('C2H3', 'O2'), ('C2H2', 'HO2'), (None,)),'! pes.subpes.NR.rxntype  8.3.2.Recombination-decomposition - termination'],
        [(('C2H3', 'O2'), ('CH2CO', 'OH'), (None,)),'! pes.subpes.NR.rxntype  8.4.2.Recombination-decomposition - termination'],
        [(('C2H3', 'O2'), ('CHOCHO', 'H'), (None,)),'! pes.subpes.NR.rxntype  8.5.2.Recombination-decomposition - termination'],
        [(('C2H3', 'O2'), ('CH2O', 'H', 'CO'), (None,)),'! pes.subpes.NR.rxntype  8.6.2.Recombination-decomposition(lumped) - termination'],
        [(('C2H5O2',), ('C2H4', 'HO2'), (None,)),'! pes.subpes.NR.rxntype  10.1.1.Beta-scission +HO2'],
        [(('C2H4', 'CH3'), ('C2H3', 'CH4'), (None,)),'! pes.subpes.NR.rxntype  11.1.2.H abstraction'],
        [(('C3H4-A', 'O'), ('C2H4', 'CO'), (None,)),'! pes.subpes.NR.rxntype  12.1.2.Addition-decomposition - termination']
    ]
    results = dict(results)
    # Read the mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_C2H4_mech.dat')
    sort_path = os.path.join(CWD, 'data', 'sort.dat')

    spc_str, mech_str, sort_str = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism
    isolate_spc, sort_lst, _ = mparser.parse_sort(sort_str)
    param_dct_sort, _, cmts_dct, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    for rxn in param_dct_sort.keys():
        assert cmts_dct[rxn]['cmts_inline'] == results[rxn]

    # now sort with the singlespecies option active
    results = [
        [(('C2H4',), ('H2', 'H2CC'), ('(+M)',)),
         '! pes.subpes.NR.rxntype  1.1.1.Decomposition'],
        [(('C2H3', 'H'), ('C2H4',), ('(+M)',)),
         '! pes.subpes.NR.rxntype  1.1.2.Recombination H'],
        [(('C2H4', 'H'), ('C2H5',), ('(+M)',)),
         '! pes.subpes.NR.rxntype  2.1.2.Addition H'],
        [(('C2H4', 'H'), ('C2H3', 'H2'), (None,)),
         '! pes.subpes.NR.rxntype  2.2.2.H abstraction'],
        [(('CH2(S)', 'CH3'), ('C2H4', 'H'), (None,)),
         '! pes.subpes.NR.rxntype  2.3.2.Addition-decomposition - propagation'],
        [(('C2H5', 'H'), ('C2H4', 'H2'), (None,)),
         '! pes.subpes.NR.rxntype  3.2.2.Recombination-decomposition - termination'],
        [(('C2H4', 'O'), ('CH3', 'HCO'), (None,)),
         '! pes.subpes.NR.rxntype  5.1.2.Addition-decomposition - branching'],
        [(('C2H4', 'OH'), ('PC2H4OH',), (None,)),
         '! pes.subpes.NR.rxntype  6.1.2.Addition OH'],
        [(('C2H5', 'OH'), ('C2H4', 'H2O'), (None,)),
         '! pes.subpes.NR.rxntype  7.1.2.Recombination-decomposition - termination'],
        [(('C2H4', 'O2'), ('C2H3', 'HO2'), (None,)),
         '! pes.subpes.NR.rxntype  9.1.2.H abstraction'],
        [(('C2H5O2',), ('C2H4', 'HO2'), (None,)),
         '! pes.subpes.NR.rxntype  10.1.1.Beta-scission +HO2'],
        [(('C2H4', 'CH3'), ('C2H3', 'CH4'), (None,)),
         '! pes.subpes.NR.rxntype  11.1.2.H abstraction'],
        [(('C3H4-A', 'O'), ('C2H4', 'CO'), (None,)),
         '! pes.subpes.NR.rxntype  12.1.2.Addition-decomposition - termination']
    ]
    results = dict(results)
    sort_path = os.path.join(CWD, 'data', 'sort_singlespecies.dat')
    _, _, sort_str = _read_files(spc_path, mech_path, sort_path)
    # Sort mechanism
    isolate_spc, sort_lst, _ = mparser.parse_sort(sort_str)
    param_dct_sort, _, cmts_dct, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    for rxn in param_dct_sort.keys():
        assert cmts_dct[rxn]['cmts_inline'] == results[rxn]

def test__sortby_submech_deletelarge():
    """ test mechanalyzer.parser.sort

        sort by fuel submechanism: extract reactions of
        fuel, fuel radicals, R+O2, R+O4
        and exclude all reactions above a given stoichiometry

        rests show different options
    """

    results = [[(('OHV',), ('OH',), (None,)) ,'CORE'],
        [(('OHV','OH'), ('OH','OH'), (None,)) ,'CORE'],
        [(('OHV','O2'), ('OH','O2'), (None,)) ,'CORE'],
        [(('OHV','N2'), ('OH','N2'), (None,)) ,'CORE'],
        [(('OHV','H'), ('OH','H'), (None,)) ,'CORE'],
        [(('OHV','H2O'), ('OH','H2O'), (None,)) ,'CORE'],
        [(('OHV','H2'), ('OH','H2'), (None,)) ,'CORE'],
        [(('OHV','CO'), ('OH','CO'), (None,)) ,'CORE'],
        [(('OHV','CO2'), ('OH','CO2'), (None,)) ,'CORE'],
        [(('OHV','AR'), ('OH','AR'), (None,)) ,'CORE'],
        [(('OH','HO2'), ('H2O','O2'), (None,)) ,'CORE'],
        [(('OCHO','HO2'), ('HOCHO','O2'), (None,)) ,'CORE'],
        [(('OCHO','H2O2'), ('HOCHO','HO2'), (None,)) ,'CORE'],
        [(('OCH2O2H',), ('HOCH2O2',), (None,)) ,'CORE'],
        [(('O2','H'), ('O','OH'), (None,)) ,'CORE'],
        [(('O','O'), ('O2',), ('+M',)) ,'CORE'],
        [(('O','H'), ('OH',), ('+M',)) ,'CORE'],
        [(('O','H2O'), ('OH','OH'), (None,)) ,'CORE'],
        [(('HOCO',), ('CO2','H'), (None,)) ,'CORE'],
        [(('HOCO',), ('CO','OH'), (None,)) ,'CORE'],
        [(('HOCHO',), ('CO2','H2'), (None,)) ,'CORE'],
        [(('HOCHO',), ('CO','H2O'), (None,)) ,'CORE'],
        [(('HOCHO','OH'), ('H2O','CO2','H'), (None,)) ,'CORE'],
        [(('HOCHO','OH'), ('H2O','CO','OH'), (None,)) ,'CORE'],
        [(('HOCHO','O'), ('CO','OH','OH'), (None,)) ,'CORE'],
        [(('HOCHO','HO2'), ('H2O2','CO','OH'), (None,)) ,'CORE'],
        [(('HOCHO','H'), ('H2','CO2','H'), (None,)) ,'CORE'],
        [(('HOCHO','H'), ('H2','CO','OH'), (None,)) ,'CORE'],
        [(('HOCH2O2','HO2'), ('HOCH2O2H','O2'), (None,)) ,'CORE'],
        [(('HO2','O'), ('OH','O2'), (None,)) ,'CORE'],
        [(('HO2','HO2'), ('H2O2','O2'), (None,)) ,'CORE'],
        [(('HO2','H'), ('OH','OH'), (None,)) ,'CORE'],
        [(('HO2','H'), ('H2','O2'), (None,)) ,'CORE'],
        [(('HCOH','OH'), ('HCO','H2O'), (None,)) ,'CORE'],
        [(('HCOH','O'), ('CO2','H','H'), (None,)) ,'CORE'],
        [(('HCOH','O'), ('CO','OH','H'), (None,)) ,'CORE'],
        [(('HCOH','O2'), ('CO2','H2O'), (None,)) ,'CORE'],
        [(('HCOH','O2'), ('CO2','H','OH'), (None,)) ,'CORE'],
        [(('HCOH','H'), ('CH2O','H'), (None,)) ,'CORE'],
        [(('HCO',), ('H','CO'), ('+M',)) ,'CORE'],
        [(('HCO','OH'), ('CO','H2O'), (None,)) ,'CORE'],
        [(('HCO','O'), ('CO2','H'), (None,)) ,'CORE'],
        [(('HCO','O'), ('CO','OH'), (None,)) ,'CORE'],
        [(('HCO','O2'), ('O2CHO',), (None,)) ,'CORE'],
        [(('HCO','O2'), ('CO','HO2'), (None,)) ,'CORE'],
        [(('HCO','HO2'), ('CO2','H','OH'), (None,)) ,'CORE'],
        [(('HCO','HCO'), ('H2','CO','CO'), (None,)) ,'CORE'],
        [(('HCO','HCO'), ('CH2O','CO'), (None,)) ,'CORE'],
        [(('HCO','H'), ('CO','H2'), (None,)) ,'CORE'],
        [(('HCO','H'), ('CH2O',), ('(+M)',)) ,'CORE'],
        [(('H2O2',), ('OH','OH'), ('(+M)',)) ,'CORE'],
        [(('H2O2','OH'), ('H2O','HO2'), (None,)) ,'CORE'],
        [(('H2O2','O'), ('OH','HO2'), (None,)) ,'CORE'],
        [(('H2O2','H'), ('H2O','OH'), (None,)) ,'CORE'],
        [(('H2O2','H'), ('H2','HO2'), (None,)) ,'CORE'],
        [(('H2CC','O2'), ('HCO','HCO'), (None,)) ,'CORE'],
        [(('H2CC','H'), ('C2H2','H'), (None,)) ,'CORE'],
        [(('H2',), ('H','H'), ('+M',)) ,'CORE'],
        [(('H2','OH'), ('H','H2O'), (None,)) ,'CORE'],
        [(('H2','O'), ('H','OH'), (None,)) ,'CORE'],
        [(('H','OH'), ('H2O',), ('+M',)) ,'CORE'],
        [(('H','O'), ('OHV',), ('+M',)) ,'CORE'],
        [(('H','O2'), ('HO2',), ('(+M)',)) ,'CORE'],
        [(('H','CO2'), ('OCHO',), (None,)) ,'CORE'],
        [(('CO','OH'), ('CO2','H'), (None,)) ,'CORE'],
        [(('CO','O'), ('CO2',), ('(+M)',)) ,'CORE'],
        [(('CO','O2'), ('CO2','O'), (None,)) ,'CORE'],
        [(('CO','HO2'), ('CO2','OH'), (None,)) ,'CORE'],
        [(('CO','H2'), ('CH2O',), ('(+M)',)) ,'CORE'],
        [(('CHV',), ('CH',), (None,)) ,'CORE'],
        [(('CHV','O2'), ('CH','O2'), (None,)) ,'CORE'],
        [(('CHV','N2'), ('CH','N2'), (None,)) ,'CORE'],
        [(('CHV','H2O'), ('CH','H2O'), (None,)) ,'CORE'],
        [(('CHV','H2'), ('CH','H2'), (None,)) ,'CORE'],
        [(('CHV','CO'), ('CH','CO'), (None,)) ,'CORE'],
        [(('CHV','CO2'), ('CH','CO2'), (None,)) ,'CORE'],
        [(('CHV','AR'), ('CH','AR'), (None,)) ,'CORE'],
        [(('CH2O','OH'), ('HCO','H2O'), (None,)) ,'CORE'],
        [(('CH2O','OCHO'), ('HCO','HOCHO'), (None,)) ,'CORE'],
        [(('CH2O','O'), ('HCO','OH'), (None,)) ,'CORE'],
        [(('CH2O','O2'), ('HCO','HO2'), (None,)) ,'CORE'],
        [(('CH2O','HO2'), ('OCH2O2H',), (None,)) ,'CORE'],
        [(('CH2O','HO2'), ('HCO','H2O2'), (None,)) ,'CORE'],
        [(('CH2O','H'), ('HCO','H2'), (None,)) ,'CORE'],
        [(('CH2','OH'), ('CH','H2O'), (None,)) ,'CORE'],
        [(('CH2','O'), ('CO','H','H'), (None,)) ,'CORE'],
        [(('CH2','O2'), ('HCO','OH'), (None,)) ,'CORE'],
        [(('CH2','O2'), ('CO2','H','H'), (None,)) ,'CORE'],
        [(('CH2','H'), ('CH','H2'), (None,)) ,'CORE'],
        [(('CH2(S)','OH'), ('CH2O','H'), (None,)) ,'CORE'],
        [(('CH2(S)','O'), ('HCO','H'), (None,)) ,'CORE'],
        [(('CH2(S)','O'), ('CO','H2'), (None,)) ,'CORE'],
        [(('CH2(S)','O2'), ('H','OH','CO'), (None,)) ,'CORE'],
        [(('CH2(S)','O2'), ('CO','H2O'), (None,)) ,'CORE'],
        [(('CH2(S)','N2'), ('CH2','N2'), (None,)) ,'CORE'],
        [(('CH2(S)','H'), ('CH','H2'), (None,)) ,'CORE'],
        [(('CH2(S)','H2O'), ('CH2','H2O'), (None,)) ,'CORE'],
        [(('CH2(S)','CO'), ('CH2','CO'), (None,)) ,'CORE'],
        [(('CH2(S)','CO2'), ('CH2O','CO'), (None,)) ,'CORE'],
        [(('CH2(S)','CO2'), ('CH2','CO2'), (None,)) ,'CORE'],
        [(('CH2(S)','AR'), ('CH2','AR'), (None,)) ,'CORE'],
        [(('CH','OH'), ('HCO','H'), (None,)) ,'CORE'],
        [(('CH','O'), ('CO','H'), (None,)) ,'CORE'],
        [(('CH','O2'), ('HCO','O'), (None,)) ,'CORE'],
        [(('CH','O2'), ('CO','OHV'), (None,)) ,'CORE'],
        [(('CH','H'), ('C','H2'), (None,)) ,'CORE'],
        [(('CH','H2O'), ('H','CH2O'), (None,)) ,'CORE'],
        [(('CH','CO2'), ('HCO','CO'), (None,)) ,'CORE'],
        [(('C2H5','O2'), ('C2H4','HO2'), (None,)) ,'CORE'],
        [(('C2H5','H'), ('C2H4','H2'), (None,)) ,'CORE'],
        [(('C2H4',), ('H2','H2CC'), ('(+M)',)) ,'CORE'],
        [(('C2H4','OH'), ('C2H3','H2O'), (None,)) ,'CORE'],
        [(('C2H4','O2'), ('C2H3','HO2'), (None,)) ,'CORE'],
        [(('C2H4','H'), ('C2H5',), ('(+M)',)) ,'CORE'],
        [(('C2H4','H'), ('C2H3','H2'), (None,)) ,'CORE'],
        [(('C2H4','C2H4'), ('C2H5','C2H3'), (None,)) ,'CORE'],
        [(('C2H3','OH'), ('C2H2','H2O'), (None,)) ,'CORE'],
        [(('C2H3','O2'), ('CH2O','HCO'), (None,)) ,'CORE'],
        [(('C2H3','O2'), ('CH2O','H','CO'), (None,)) ,'CORE'],
        [(('C2H3','O2'), ('C2H2','HO2'), (None,)) ,'CORE'],
        [(('C2H3','H'), ('H2CC','H2'), (None,)) ,'CORE'],
        [(('C2H3','H'), ('C2H4',), ('(+M)',)) ,'CORE'],
        [(('C2H3','H'), ('C2H2','H2'), (None,)) ,'CORE'],
        [(('C2H3','CH2O'), ('C2H4','HCO'), (None,)) ,'CORE'],
        [(('C2H3','C2H3'), ('C2H2','C2H4'), (None,)) ,'CORE'],
        [(('C2H2',), ('H2CC',), ('(+M)',)) ,'CORE'],
        [(('C2H2','OH'), ('C2H','H2O'), (None,)) ,'CORE'],
        [(('C2H2','O'), ('CH2','CO'), (None,)) ,'CORE'],
        [(('C2H2','HCO'), ('C2H3','CO'), (None,)) ,'CORE'],
        [(('C2H2','H'), ('C2H3',), ('(+M)',)) ,'CORE'],
        [(('C2H','O'), ('CO','CHV'), (None,)) ,'CORE'],
        [(('C2H','O'), ('CO','CH'), (None,)) ,'CORE'],
        [(('C2H','O2'), ('HCO','CO'), (None,)) ,'CORE'],
        [(('C2H','O2'), ('CO2','CHV'), (None,)) ,'CORE'],
        [(('C2H','H'), ('C2H2',), ('(+M)',)) ,'CORE'],
        [(('C2H','H2'), ('H','C2H2'), (None,)) ,'CORE'],
        [(('C','OH'), ('CO','H'), (None,)) ,'CORE'],
        [(('C','O2'), ('CO','O'), (None,)) ,'CORE'],
        [(('C','H'), ('CHV',), ('+M',)) ,'CORE'],
        [(('OHV','CH4'), ('OH','CH4'), (None,)) ,'FUEL'],
        [(('IC4H8','CH3'), ('IC4H7-I1','CH4'), (None,)) ,'FUEL'],
        [(('IC4H8','CH3'), ('IC4H7','CH4'), (None,)) ,'FUEL'],
        [(('IC4H7-I1','H'), ('C3H4-P','CH4'), (None,)) ,'FUEL'],
        [(('IC4H7-I1','H'), ('C3H4-A','CH4'), (None,)) ,'FUEL'],
        [(('IC4H10','CH3'), ('TC4H9','CH4'), (None,)) ,'FUEL'],
        [(('IC4H10','CH3'), ('IC4H9','CH4'), (None,)) ,'FUEL'],
        [(('HOCHO','CH3'), ('CH4','CO','OH'), (None,)) ,'FUEL'],
        [(('HCO','CH3'), ('CH4','CO'), (None,)) ,'FUEL'],
        [(('CHV','CH4'), ('CH','CH4'), (None,)) ,'FUEL'],
        [(('CH4','OH'), ('CH3','H2O'), (None,)) ,'FUEL'],
        [(('CH4','O'), ('CH3','OH'), (None,)) ,'FUEL'],
        [(('CH4','HO2'), ('CH3','H2O2'), (None,)) ,'FUEL'],
        [(('CH4','H'), ('CH3','H2'), (None,)) ,'FUEL'],
        [(('CH4','CH3O2'), ('CH3','CH3O2H'), (None,)) ,'FUEL'],
        [(('CH4','CH2'), ('CH3','CH3'), (None,)) ,'FUEL'],
        [(('CH3OH','CH3'), ('CH3O','CH4'), (None,)) ,'FUEL'],
        [(('CH3OH','CH3'), ('CH2OH','CH4'), (None,)) ,'FUEL'],
        [(('CH3O','CH3'), ('CH2O','CH4'), (None,)) ,'FUEL'],
        [(('CH3','HO2'), ('CH4','O2'), (None,)) ,'FUEL'],
        [(('CH3','H'), ('CH4',), ('(+M)',)) ,'FUEL'],
        [(('CH2O','CH3'), ('HCO','CH4'), (None,)) ,'FUEL'],
        [(('CH','CH4'), ('C2H4','H'), (None,)) ,'FUEL'],
        [(('C8H142-6','CH3'), ('C8H132-6,SA','CH4'), (None,)) ,'FUEL'],
        [(('C8H142-6','CH3'), ('C8H132-6,PA','CH4'), (None,)) ,'FUEL'],
        [(('C8H141-5,3-4','CH3'), ('C8H131-5,3-4,TA','CH4'), (None,)) ,'FUEL'],
        [(('C8H141-5,3','CH3'), ('C8H131-5,3,TA','CH4'), (None,)) ,'FUEL'],
        [(('C8H141-5,3','CH3'), ('C8H131-5,3,SA','CH4'), (None,)) ,'FUEL'],
        [(('C8H141-5,3','CH3'), ('C8H131-5,3,PA','CH4'), (None,)) ,'FUEL'],
        [(('C4H8-2','CH3'), ('C4H72-2','CH4'), (None,)) ,'FUEL'],
        [(('C4H8-2','CH3'), ('C4H71-3','CH4'), (None,)) ,'FUEL'],
        [(('C4H8-1','CH3'), ('C4H71-4','CH4'), (None,)) ,'FUEL'],
        [(('C4H8-1','CH3'), ('C4H71-3','CH4'), (None,)) ,'FUEL'],
        [(('C4H8-1','CH3'), ('C4H71-2','CH4'), (None,)) ,'FUEL'],
        [(('C4H8-1','CH3'), ('C4H71-1','CH4'), (None,)) ,'FUEL'],
        [(('C4H612','CH3'), ('C4H5-I','CH4'), (None,)) ,'FUEL'],
        [(('C4H6-2','CH3'), ('C4H5-2','CH4'), (None,)) ,'FUEL'],
        [(('C4H6','CH3'), ('C4H5-N','CH4'), (None,)) ,'FUEL'],
        [(('C4H6','CH3'), ('C4H5-I','CH4'), (None,)) ,'FUEL'],
        [(('C4H4','CH3'), ('C4H3-N','CH4'), (None,)) ,'FUEL'],
        [(('C4H4','CH3'), ('C4H3-I','CH4'), (None,)) ,'FUEL'],
        [(('C4H10','CH3'), ('SC4H9','CH4'), (None,)) ,'FUEL'],
        [(('C4H10','CH3'), ('PC4H9','CH4'), (None,)) ,'FUEL'],
        [(('C3H8','CH3'), ('NC3H7','CH4'), (None,)) ,'FUEL'],
        [(('C3H8','CH3'), ('IC3H7','CH4'), (None,)) ,'FUEL'],
        [(('C3H6','CH3'), ('C3H5-T','CH4'), (None,)) ,'FUEL'],
        [(('C3H6','CH3'), ('C3H5-S','CH4'), (None,)) ,'FUEL'],
        [(('C3H6','CH3'), ('C3H5-A','CH4'), (None,)) ,'FUEL'],
        [(('C3H5-T','CH3'), ('C3H4-P','CH4'), (None,)) ,'FUEL'],
        [(('C3H5-S','CH3'), ('C3H4-P','CH4'), (None,)) ,'FUEL'],
        [(('C3H5-S','CH3'), ('C3H4-A','CH4'), (None,)) ,'FUEL'],
        [(('C3H5-A','CH3'), ('C3H4-A','CH4'), (None,)) ,'FUEL'],
        [(('C3H4-P','CH3'), ('C3H3','CH4'), (None,)) ,'FUEL'],
        [(('C3H4-A','CH3'), ('C3H3','CH4'), (None,)) ,'FUEL'],
        [(('C2H6','CH3'), ('C2H5','CH4'), (None,)) ,'FUEL'],
        [(('C2H5','CH3'), ('CH4','C2H4'), (None,)) ,'FUEL'],
        [(('C2H4','CH3'), ('C2H3','CH4'), (None,)) ,'FUEL'],
        [(('C2H3','CH3'), ('CH4','C2H2'), (None,)) ,'FUEL'],
        [(('IC4H8','CH3O'), ('IC4H7','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('IC4H10','CH3O'), ('TC4H9','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('IC4H10','CH3O'), ('IC4H9','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3OH',), ('CH2OH','H'), ('(+M)',)) ,'FUEL_ADD_O'],
        [(('CH3OH',), ('CH2(S)','H2O'), ('(+M)',)) ,'FUEL_ADD_O'],
        [(('CH3OH','OH'), ('CH3O','H2O'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3OH','OH'), ('CH2OH','H2O'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3OH','O'), ('CH3O','OH'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3OH','O'), ('CH2OH','OH'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3OH','O2'), ('CH3O','HO2'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3OH','O2'), ('CH2OH','HO2'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3OH','HO2'), ('CH3O','H2O2'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3OH','HO2'), ('CH2OH','H2O2'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3OH','HCO'), ('CH2OH','CH2O'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3OH','H'), ('CH3O','H2'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3OH','H'), ('CH2OH','H2'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3OH','CH3O'), ('CH2OH','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3OH','CH3O2'), ('CH2OH','CH3O2H'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3O2','OH'), ('CH3OH','O2'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3O2','CH3O2'), ('CH2O','CH3OH','O2'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3O','HCO'), ('CH3OH','CO'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3O','CH3O'), ('CH3OH','CH2O'), (None,)) ,'FUEL_ADD_O'],
        [(('CH2OH','HCO'), ('CH3OH','CO'), (None,)) ,'FUEL_ADD_O'],
        [(('CH2OH','CH3O'), ('CH2O','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('CH2OH','CH2OH'), ('CH2O','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('CH2O','CH3O'), ('HCO','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('C4H8-2','CH3O'), ('C4H71-3','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('C4H8-1','CH3O'), ('C4H71-4','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('C4H8-1','CH3O'), ('C4H71-3','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('C4H10','CH3O'), ('SC4H9','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('C4H10','CH3O'), ('PC4H9','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('C3H8','CH3O'), ('NC3H7','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('C3H8','CH3O'), ('IC3H7','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('C3H6','CH3O'), ('C3H5-A','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('C2H6','CH3O'), ('C2H5','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('C2H4','CH3O'), ('C2H3','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('IC4H8','CH3O2'), ('IC4H7','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('IC4H10','CH3O2'), ('TC4H9','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('IC4H10','CH3O2'), ('IC4H9','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('H2','CH3O2'), ('H','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('CH3O2H',), ('CH3O','OH'), (None,)) ,'FUEL_ADD_O2'],
        [(('CH3O2','HO2'), ('CH3O2H','O2'), (None,)) ,'FUEL_ADD_O2'],
        [(('CH3O2','H2O2'), ('CH3O2H','HO2'), (None,)) ,'FUEL_ADD_O2'],
        [(('CH2O','CH3O2'), ('HCO','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('C4H8-1','CH3O2'), ('C4H71-4','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('C4H8-1','CH3O2'), ('C4H71-3','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('C4H10','CH3O2'), ('SC4H9','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('C4H10','CH3O2'), ('PC4H9','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('C3H8','CH3O2'), ('NC3H7','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('C3H8','CH3O2'), ('IC3H7','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('C3H6','CH3O2'), ('C3H5-A','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('C3H4-P','CH3O2'), ('C3H3','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('C3H4-A','CH3O2'), ('C3H3','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('C2H6','CH3O2'), ('C2H5','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('C2H4','CH3O2'), ('C2H3','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('SC4H9',), ('C3H6','CH3'), (None,)) ,'FUEL_RAD'],
        [(('PC4H9',), ('C3H6','CH3'), (None,)) ,'FUEL_RAD'],
        [(('IC4H9',), ('C3H6','CH3'), (None,)) ,'FUEL_RAD'],
        [(('IC4H8',), ('C3H5-T','CH3'), (None,)) ,'FUEL_RAD'],
        [(('IC4H8','H'), ('C3H6','CH3'), (None,)) ,'FUEL_RAD'],
        [(('IC4H7-I1','CH3'), ('C3H4-P','C2H6'), (None,)) ,'FUEL_RAD'],
        [(('IC4H10',), ('CH3','IC3H7'), ('(+M)',)) ,'FUEL_RAD'],
        [(('IC3H7','H'), ('C2H5','CH3'), (None,)) ,'FUEL_RAD'],
        [(('CH3OH',), ('CH3','OH'), ('(+M)',)) ,'FUEL_RAD'],
        [(('CH3O2','CH3'), ('CH3O','CH3O'), (None,)) ,'FUEL_RAD'],
        [(('CH3','OH'), ('HCOH','H2'), (None,)) ,'FUEL_RAD'],
        [(('CH3','OH'), ('H','CH3O'), (None,)) ,'FUEL_RAD'],
        [(('CH3','OH'), ('CH2OH','H'), (None,)) ,'FUEL_RAD'],
        [(('CH3','OH'), ('CH2O','H2'), (None,)) ,'FUEL_RAD'],
        [(('CH3','OH'), ('CH2','H2O'), (None,)) ,'FUEL_RAD'],
        [(('CH3','OH'), ('CH2(S)','H2O'), (None,)) ,'FUEL_RAD'],
        [(('CH3','O'), ('CH2O','H'), (None,)) ,'FUEL_RAD'],
        [(('CH3','O2'), ('CH3O2',), ('(+M)',)) ,'FUEL_RAD'],
        [(('CH3','O2'), ('CH3O','O'), (None,)) ,'FUEL_RAD'],
        [(('CH3','O2'), ('CH2O','OH'), (None,)) ,'FUEL_RAD'],
        [(('CH3','HO2'), ('CH3O','OH'), (None,)) ,'FUEL_RAD'],
        [(('CH3','CH3'), ('H','C2H5'), (None,)) ,'FUEL_RAD'],
        [(('CH3','CH3'), ('C2H6',), ('(+M)',)) ,'FUEL_RAD'],
        [(('CH3','C2H2'), ('C3H5-S',), (None,)) ,'FUEL_RAD'],
        [(('CH2','H'), ('CH3',), ('(+M)',)) ,'FUEL_RAD'],
        [(('CH2(S)','H2'), ('CH3','H'), (None,)) ,'FUEL_RAD'],
        [(('CH2(S)','CH3'), ('C2H4','H'), (None,)) ,'FUEL_RAD'],
        [(('CH2(S)','C2H4'), ('C2H3','CH3'), (None,)) ,'FUEL_RAD'],
        [(('CC3H6',), ('C2H3','CH3'), (None,)) ,'FUEL_RAD'],
        [(('C4H8-2',), ('C3H5-A','CH3'), (None,)) ,'FUEL_RAD'],
        [(('C4H8-2','H'), ('C3H6','CH3'), (None,)) ,'FUEL_RAD'],
        [(('C4H8-1','H'), ('C3H6','CH3'), (None,)) ,'FUEL_RAD'],
        [(('C4H72-2',), ('C3H4-P','CH3'), (None,)) ,'FUEL_RAD'],
        [(('C4H71-2',), ('C3H4-A','CH3'), (None,)) ,'FUEL_RAD'],
        [(('C4H612','H'), ('C3H4-P','CH3'), (None,)) ,'FUEL_RAD'],
        [(('C4H612','H'), ('C3H4-A','CH3'), (None,)) ,'FUEL_RAD'],
        [(('C4H6-2','H'), ('CH3','C3H4-P'), (None,)) ,'FUEL_RAD'],
        [(('C4H6','H'), ('C3H4-P','CH3'), (None,)) ,'FUEL_RAD'],
        [(('C4H6','H'), ('C3H4-A','CH3'), (None,)) ,'FUEL_RAD'],
        [(('C4H5-I','H'), ('C3H3','CH3'), (None,)) ,'FUEL_RAD'],
        [(('C4H10',), ('NC3H7','CH3'), ('(+M)',)) ,'FUEL_RAD'],
        [(('C3H8',), ('CH3','C2H5'), ('(+M)',)) ,'FUEL_RAD'],
        [(('C3H6',), ('C2H3','CH3'), (None,)) ,'FUEL_RAD'],
        [(('C3H6','H'), ('C2H4','CH3'), (None,)) ,'FUEL_RAD'],
        [(('C3H5-T','H'), ('C2H3','CH3'), (None,)) ,'FUEL_RAD'],
        [(('C3H5-S','CH3'), ('C4H8-2',), ('(+M)',)) ,'FUEL_RAD'],
        [(('C3H5-A','CH3'), ('C4H8-1',), ('(+M)',)) ,'FUEL_RAD'],
        [(('C3H4-P','H'), ('CH3','C2H2'), (None,)) ,'FUEL_RAD'],
        [(('C3H4-A','H'), ('CH3','C2H2'), (None,)) ,'FUEL_RAD'],
        [(('C3H4-A','CH3'), ('IC4H7',), (None,)) ,'FUEL_RAD'],
        [(('C3H3','CH3'), ('C4H612',), ('(+M)',)) ,'FUEL_RAD'],
        [(('C2H6','CH2(S)'), ('C2H5','CH3'), (None,)) ,'FUEL_RAD'],
        [(('C2H5','C2H'), ('C3H3','CH3'), (None,)) ,'FUEL_RAD'],
        [(('C2H4','OH'), ('CH3','CH2O'), (None,)) ,'FUEL_RAD'],
        [(('C2H4','O'), ('CH3','HCO'), (None,)) ,'FUEL_RAD'],
        [(('C2H4','CH3'), ('NC3H7',), (None,)) ,'FUEL_RAD'],
        [(('C2H3','O2'), ('CO2','CH3'), (None,)) ,'FUEL_RAD'],
        [(('C2H3','CH3'), ('C3H6',), ('(+M)',)) ,'FUEL_RAD'],
        [(('C2H3','CH3'), ('C3H5-A','H'), (None,)) ,'FUEL_RAD'],
        [(('C2H2','OH'), ('CH3','CO'), (None,)) ,'FUEL_RAD'],
        [(('C2H2','CH3'), ('C3H5-T',), (None,)) ,'FUEL_RAD'],
        [(('C2H2','CH3'), ('C3H5-A',), (None,)) ,'FUEL_RAD'],
        [(('C2H','CH3'), ('C3H4-P',), (None,)) ,'FUEL_RAD'],
        [(('IC4H10','C2H5'), ('TC4H9','C2H6'), (None,)) ,'R_CH3'],
        [(('IC4H10','C2H5'), ('IC4H9','C2H6'), (None,)) ,'R_CH3'],
        [(('C4H10','C2H5'), ('SC4H9','C2H6'), (None,)) ,'R_CH3'],
        [(('C4H10','C2H5'), ('PC4H9','C2H6'), (None,)) ,'R_CH3'],
        [(('C3H8','C2H5'), ('NC3H7','C2H6'), (None,)) ,'R_CH3'],
        [(('C3H8','C2H5'), ('IC3H7','C2H6'), (None,)) ,'R_CH3'],
        [(('C3H6','C2H5'), ('C3H5-A','C2H6'), (None,)) ,'R_CH3'],
        [(('C3H5-A','C2H5'), ('C3H4-A','C2H6'), (None,)) ,'R_CH3'],
        [(('C2H6','OH'), ('C2H5','H2O'), (None,)) ,'R_CH3'],
        [(('C2H6','O'), ('C2H5','OH'), (None,)) ,'R_CH3'],
        [(('C2H6','O2'), ('C2H5','HO2'), (None,)) ,'R_CH3'],
        [(('C2H6','HO2'), ('C2H5','H2O2'), (None,)) ,'R_CH3'],
        [(('C2H6','H'), ('C2H5','H2'), (None,)) ,'R_CH3'],
        [(('C2H6','CH'), ('C2H5','CH2'), (None,)) ,'R_CH3'],
        [(('C2H5','H'), ('C2H6',), ('(+M)',)) ,'R_CH3'],
        [(('C2H5','C4H71-3'), ('C4H6','C2H6'), (None,)) ,'R_CH3'],
        [(('CH3O',), ('CH2O','H'), ('(+M)',)) ,'R_O'],
        [(('CH3O2','O'), ('CH3O','O2'), (None,)) ,'R_O'],
        [(('CH3O2','H'), ('CH3O','OH'), (None,)) ,'R_O'],
        [(('CH3O2','CH3O2'), ('O2','CH3O','CH3O'), (None,)) ,'R_O'],
        [(('CH3O','O2'), ('CH2O','HO2'), (None,)) ,'R_O'],
        [(('CH3O','HO2'), ('CH2O','H2O2'), (None,)) ,'R_O'],
        [(('CH3O','H'), ('CH2O','H2'), (None,)) ,'R_O'],
        [(('CH2OH','OH'), ('H2O','CH2O'), (None,)) ,'R_O'],
        [(('CH2OH','O'), ('OH','CH2O'), (None,)) ,'R_O'],
        [(('CH2OH','O2'), ('CH2O','HO2'), (None,)) ,'R_O'],
        [(('CH2OH','HO2'), ('HOCH2O','OH'), (None,)) ,'R_O'],
        [(('CH2OH','HO2'), ('CH2O','H2O2'), (None,)) ,'R_O'],
        [(('CH2OH','HCO'), ('CH2O','CH2O'), (None,)) ,'R_O'],
        [(('CH2OH','H'), ('CH2O','H2'), (None,)) ,'R_O'],
        [(('CH2O','H'), ('CH2OH',), ('(+M)',)) ,'R_O'],
        [(('C4H71-3','CH3O'), ('C4H8-1','CH2O'), (None,)) ,'R_O'],
        [(('C4H5-2','OH'), ('CH2OH','C3H3'), (None,)) ,'R_O'],
        [(('C2H3','O2'), ('CO','CH3O'), (None,)) ,'R_O'],
        [(('HOCH2O',), ('HOCHO','H'), (None,)) ,'R_O2'],
        [(('HOCH2O','OH'), ('HOCH2O2H',), (None,)) ,'R_O2'],
        [(('CH2O2H',), ('CH2O','OH'), (None,)) ,'R_O2'],
        [(('CH2O','OH'), ('HOCH2O',), (None,)) ,'R_O2'],
        [(('OCHO','OH'), ('HO2CHO',), (None,)) ,'R_O3-H'],
        [(('IC4H8','O2CHO'), ('IC4H7','HO2CHO'), (None,)) ,'R_O3-H'],
        [(('IC4H10','O2CHO'), ('TC4H9','HO2CHO'), (None,)) ,'R_O3-H'],
        [(('IC4H10','O2CHO'), ('IC4H9','HO2CHO'), (None,)) ,'R_O3-H'],
        [(('CH2O','O2CHO'), ('HCO','HO2CHO'), (None,)) ,'R_O3-H'],
        [(('C4H10','O2CHO'), ('SC4H9','HO2CHO'), (None,)) ,'R_O3-H'],
        [(('C4H10','O2CHO'), ('PC4H9','HO2CHO'), (None,)) ,'R_O3-H'],
        [(('C3H8','O2CHO'), ('NC3H7','HO2CHO'), (None,)) ,'R_O3-H'],
        [(('C3H8','O2CHO'), ('IC3H7','HO2CHO'), (None,)) ,'R_O3-H'],
        [(('TC4H9','O2'), ('IC4H8','HO2'), (None,)) ,'SUPFUEL'],
        [(('SC4H9',), ('PC4H9',), (None,)) ,'SUPFUEL'],
        [(('SC4H9',), ('C2H4','C2H5'), (None,)) ,'SUPFUEL'],
        [(('SC4H9','O2'), ('C4H8-2','HO2'), (None,)) ,'SUPFUEL'],
        [(('SC4H9','O2'), ('C4H8-1','HO2'), (None,)) ,'SUPFUEL'],
        [(('PC4H9',), ('C2H4','C2H5'), (None,)) ,'SUPFUEL'],
        [(('PC4H9','O2'), ('C4H8-1','HO2'), (None,)) ,'SUPFUEL'],
        [(('NC3H7','O2'), ('C3H6','HO2'), (None,)) ,'SUPFUEL'],
        [(('NC3H7','H'), ('C3H8',), (None,)) ,'SUPFUEL'],
        [(('IC4H9',), ('TC4H9',), (None,)) ,'SUPFUEL'],
        [(('IC4H9','O2'), ('IC4H8','HO2'), (None,)) ,'SUPFUEL'],
        [(('IC4H8',), ('IC4H7-I1','H'), (None,)) ,'SUPFUEL'],
        [(('IC4H8',), ('IC4H7','H'), (None,)) ,'SUPFUEL'],
        [(('IC4H8','OH'), ('IC4H7-I1','H2O'), (None,)) ,'SUPFUEL'],
        [(('IC4H8','OH'), ('IC4H7','H2O'), (None,)) ,'SUPFUEL'],
        [(('IC4H8','O'), ('IC4H7-I1','OH'), (None,)) ,'SUPFUEL'],
        [(('IC4H8','O'), ('IC4H7','OH'), (None,)) ,'SUPFUEL'],
        [(('IC4H8','O'), ('IC3H7','HCO'), (None,)) ,'SUPFUEL'],
        [(('IC4H8','O2'), ('IC4H7-I1','HO2'), (None,)) ,'SUPFUEL'],
        [(('IC4H8','O2'), ('IC4H7','HO2'), (None,)) ,'SUPFUEL'],
        [(('IC4H8','HO2'), ('IC4H7-I1','H2O2'), (None,)) ,'SUPFUEL'],
        [(('IC4H8','HO2'), ('IC4H7','H2O2'), (None,)) ,'SUPFUEL'],
        [(('IC4H8','H'), ('TC4H9',), (None,)) ,'SUPFUEL'],
        [(('IC4H8','H'), ('IC4H9',), (None,)) ,'SUPFUEL'],
        [(('IC4H8','H'), ('IC4H7-I1','H2'), (None,)) ,'SUPFUEL'],
        [(('IC4H8','H'), ('IC4H7','H2'), (None,)) ,'SUPFUEL'],
        [(('IC4H8','C3H5-T'), ('IC4H7','C3H6'), (None,)) ,'SUPFUEL'],
        [(('IC4H8','C3H5-S'), ('IC4H7','C3H6'), (None,)) ,'SUPFUEL'],
        [(('IC4H8','C3H5-A'), ('IC4H7','C3H6'), (None,)) ,'SUPFUEL'],
        [(('IC4H7',), ('IC4H7-I1',), (None,)) ,'SUPFUEL'],
        [(('IC4H7-I1','OH'), ('C3H6','HCO','H'), (None,)) ,'SUPFUEL'],
        [(('IC4H7-I1','O'), ('C3H6','HCO'), (None,)) ,'SUPFUEL'],
        [(('IC4H7-I1','HO2'), ('C3H6','HCO','OH'), (None,)) ,'SUPFUEL'],
        [(('IC4H7-I1','HCO'), ('IC4H8','CO'), (None,)) ,'SUPFUEL'],
        [(('IC4H10',), ('TC4H9','H'), (None,)) ,'SUPFUEL'],
        [(('IC4H10',), ('IC4H9','H'), (None,)) ,'SUPFUEL'],
        [(('IC4H10','OH'), ('TC4H9','H2O'), (None,)) ,'SUPFUEL'],
        [(('IC4H10','OH'), ('IC4H9','H2O'), (None,)) ,'SUPFUEL'],
        [(('IC4H10','O'), ('TC4H9','OH'), (None,)) ,'SUPFUEL'],
        [(('IC4H10','O'), ('IC4H9','OH'), (None,)) ,'SUPFUEL'],
        [(('IC4H10','O2'), ('TC4H9','HO2'), (None,)) ,'SUPFUEL'],
        [(('IC4H10','O2'), ('IC4H9','HO2'), (None,)) ,'SUPFUEL'],
        [(('IC4H10','IC4H9'), ('TC4H9','IC4H10'), (None,)) ,'SUPFUEL'],
        [(('IC4H10','HO2'), ('TC4H9','H2O2'), (None,)) ,'SUPFUEL'],
        [(('IC4H10','HO2'), ('IC4H9','H2O2'), (None,)) ,'SUPFUEL'],
        [(('IC4H10','H'), ('TC4H9','H2'), (None,)) ,'SUPFUEL'],
        [(('IC4H10','H'), ('IC4H9','H2'), (None,)) ,'SUPFUEL'],
        [(('IC3H7','OH'), ('C3H6','H2O'), (None,)) ,'SUPFUEL'],
        [(('IC3H7','H'), ('C3H8',), (None,)) ,'SUPFUEL'],
        [(('H2CCC(S)','O2'), ('CO2','C2H2'), (None,)) ,'SUPFUEL'],
        [(('H2CC','C2H4'), ('C4H6',), (None,)) ,'SUPFUEL'],
        [(('H2CC','C2H2'), ('C4H4',), ('(+M)',)) ,'SUPFUEL'],
        [(('H','C4H71-3'), ('C4H6','H2'), (None,)) ,'SUPFUEL'],
        [(('CH2(S)','C2H4'), ('CC3H6',), (None,)) ,'SUPFUEL'],
        [(('CH2(S)','C2H4'), ('C3H6',), (None,)) ,'SUPFUEL'],
        [(('CH2(S)','C2H4'), ('C3H5-A','H'), (None,)) ,'SUPFUEL'],
        [(('CC3H6',), ('C3H5-A','H'), (None,)) ,'SUPFUEL'],
        [(('CC3H4',), ('C3H4-P',), (None,)) ,'SUPFUEL'],
        [(('CC3H4',), ('C3H4-A',), (None,)) ,'SUPFUEL'],
        [(('C8H142-6','OH'), ('C8H132-6,SA','H2O'), (None,)) ,'SUPFUEL'],
        [(('C8H142-6','OH'), ('C8H132-6,PA','H2O'), (None,)) ,'SUPFUEL'],
        [(('C8H142-6','O'), ('C8H132-6,SA','OH'), (None,)) ,'SUPFUEL'],
        [(('C8H142-6','O'), ('C8H132-6,PA','OH'), (None,)) ,'SUPFUEL'],
        [(('C8H142-6','O2'), ('C8H132-6,SA','HO2'), (None,)) ,'SUPFUEL'],
        [(('C8H142-6','O2'), ('C8H132-6,PA','HO2'), (None,)) ,'SUPFUEL'],
        [(('C8H142-6','HO2'), ('C8H132-6,SA','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C8H142-6','HO2'), ('C8H132-6,PA','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C8H142-6','H'), ('C8H132-6,SA','H2'), (None,)) ,'SUPFUEL'],
        [(('C8H142-6','H'), ('C8H132-6,PA','H2'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3-4','OH'), ('C8H131-5,3-4,TA','H2O'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3-4','O'), ('C8H131-5,3-4,TA','OH'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3-4','O2'), ('C8H131-5,3-4,TA','HO2'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3-4','HO2'), ('C8H131-5,3-4,TA','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3-4','H'), ('C8H131-5,3-4,TA','H2'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3','OH'), ('C8H131-5,3,TA','H2O'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3','OH'), ('C8H131-5,3,SA','H2O'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3','OH'), ('C8H131-5,3,PA','H2O'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3','O'), ('C8H131-5,3,TA','OH'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3','O'), ('C8H131-5,3,SA','OH'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3','O'), ('C8H131-5,3,PA','OH'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3','O2'), ('C8H131-5,3,TA','HO2'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3','O2'), ('C8H131-5,3,SA','HO2'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3','O2'), ('C8H131-5,3,PA','HO2'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3','HO2'), ('C8H131-5,3,TA','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3','HO2'), ('C8H131-5,3,SA','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3','HO2'), ('C8H131-5,3,PA','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3','H'), ('C8H131-5,3,TA','H2'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3','H'), ('C8H131-5,3,SA','H2'), (None,)) ,'SUPFUEL'],
        [(('C8H141-5,3','H'), ('C8H131-5,3,PA','H2'), (None,)) ,'SUPFUEL'],
        [(('C6H101-3,3',), ('C2H3','C4H72-2'), (None,)) ,'SUPFUEL'],
        [(('C6H101-3,3','C2H3'), ('C8H131-5,3-4,TA',), (None,)) ,'SUPFUEL'],
        [(('C4H8-2',), ('H','C4H71-3'), (None,)) ,'SUPFUEL'],
        [(('C4H8-2','OH'), ('C4H72-2','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H8-2','OH'), ('C4H71-3','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H8-2','O'), ('C4H72-2','OH'), (None,)) ,'SUPFUEL'],
        [(('C4H8-2','O'), ('C4H71-3','OH'), (None,)) ,'SUPFUEL'],
        [(('C4H8-2','O2'), ('C4H72-2','HO2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-2','O2'), ('C4H71-3','HO2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-2','HO2'), ('C4H72-2','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-2','HO2'), ('C4H71-3','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-2','H'), ('SC4H9',), (None,)) ,'SUPFUEL'],
        [(('C4H8-2','H'), ('PC4H9',), (None,)) ,'SUPFUEL'],
        [(('C4H8-2','H'), ('C4H72-2','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-2','H'), ('C4H71-3','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-2','H'), ('C2H4','C2H5'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','OH'), ('C4H71-4','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','OH'), ('C4H71-3','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','OH'), ('C4H71-2','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','OH'), ('C4H71-1','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','O'), ('NC3H7','HCO'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','O'), ('C4H71-4','OH'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','O'), ('C4H71-3','OH'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','O'), ('C4H71-2','OH'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','O'), ('C4H71-1','OH'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','O2'), ('C4H71-4','HO2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','O2'), ('C4H71-3','HO2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','O2'), ('C4H71-2','HO2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','O2'), ('C4H71-1','HO2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','HO2'), ('C4H71-4','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','HO2'), ('C4H71-3','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','HO2'), ('C4H71-2','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','HO2'), ('C4H71-1','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','H'), ('SC4H9',), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','H'), ('PC4H9',), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','H'), ('C4H8-2','H'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','H'), ('C4H71-4','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','H'), ('C4H71-3','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','H'), ('C4H71-2','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','H'), ('C4H71-1','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','H'), ('C2H4','C2H5'), (None,)) ,'SUPFUEL'],
        [(('C4H8-1','C3H5-A'), ('C4H71-3','C3H6'), (None,)) ,'SUPFUEL'],
        [(('C4H72-2',), ('C4H612','H'), (None,)) ,'SUPFUEL'],
        [(('C4H72-2',), ('C4H6-2','H'), (None,)) ,'SUPFUEL'],
        [(('C4H71-4',), ('C4H6','H'), (None,)) ,'SUPFUEL'],
        [(('C4H71-4',), ('C2H4','C2H3'), (None,)) ,'SUPFUEL'],
        [(('C4H71-4','O2'), ('C4H6','HO2'), (None,)) ,'SUPFUEL'],
        [(('C4H71-4','H'), ('C4H8-1',), ('(+M)',)) ,'SUPFUEL'],
        [(('C4H71-3',), ('C4H72-2',), (None,)) ,'SUPFUEL'],
        [(('C4H71-3',), ('C4H71-4',), (None,)) ,'SUPFUEL'],
        [(('C4H71-3',), ('C4H612','H'), (None,)) ,'SUPFUEL'],
        [(('C4H71-3',), ('C4H6','H'), (None,)) ,'SUPFUEL'],
        [(('C4H71-3','O2'), ('C4H6','HO2'), (None,)) ,'SUPFUEL'],
        [(('C4H71-3','H'), ('C4H8-1',), (None,)) ,'SUPFUEL'],
        [(('C4H71-3','C4H71-3'), ('C8H142-6',), (None,)) ,'SUPFUEL'],
        [(('C4H71-3','C4H71-3'), ('C8H141-5,3-4',), (None,)) ,'SUPFUEL'],
        [(('C4H71-3','C4H71-3'), ('C8H141-5,3',), (None,)) ,'SUPFUEL'],
        [(('C4H71-3','C2H5'), ('C4H8-1','C2H4'), (None,)) ,'SUPFUEL'],
        [(('C4H71-2',), ('C4H612','H'), (None,)) ,'SUPFUEL'],
        [(('C4H71-2',), ('C4H6-1','H'), (None,)) ,'SUPFUEL'],
        [(('C4H71-1',), ('C4H6-1','H'), (None,)) ,'SUPFUEL'],
        [(('C4H71-1',), ('C2H5','C2H2'), (None,)) ,'SUPFUEL'],
        [(('C4H6',), ('C4H5-N','H'), (None,)) ,'SUPFUEL'],
        [(('C4H6',), ('C4H5-I','H'), (None,)) ,'SUPFUEL'],
        [(('C4H6',), ('C4H4','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H612',), ('C4H6',), (None,)) ,'SUPFUEL'],
        [(('C4H612',), ('C4H5-I','H'), (None,)) ,'SUPFUEL'],
        [(('C4H612','OH'), ('C4H5-I','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H612','O'), ('C4H5-I','OH'), (None,)) ,'SUPFUEL'],
        [(('C4H612','H'), ('C4H6','H'), (None,)) ,'SUPFUEL'],
        [(('C4H612','H'), ('C4H5-I','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H6-2',), ('H','C4H5-2'), (None,)) ,'SUPFUEL'],
        [(('C4H6-2',), ('C4H612',), (None,)) ,'SUPFUEL'],
        [(('C4H6-2',), ('C4H6',), (None,)) ,'SUPFUEL'],
        [(('C4H6-2','H'), ('C4H612','H'), (None,)) ,'SUPFUEL'],
        [(('C4H6-2','H'), ('C4H5-2','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H6','OH'), ('C4H5-N','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H6','OH'), ('C4H5-I','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H6','OH'), ('C3H5-A','CH2O'), (None,)) ,'SUPFUEL'],
        [(('C4H6','O'), ('C4H5-N','OH'), (None,)) ,'SUPFUEL'],
        [(('C4H6','O'), ('C4H5-I','OH'), (None,)) ,'SUPFUEL'],
        [(('C4H6','H'), ('C4H5-N','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H6','H'), ('C4H5-I','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H6','H'), ('C2H4','C2H3'), (None,)) ,'SUPFUEL'],
        [(('C4H6','C4H71-3'), ('C8H132-6,PA',), (None,)) ,'SUPFUEL'],
        [(('C4H6','C4H71-3'), ('C8H131-5,3,PA',), (None,)) ,'SUPFUEL'],
        [(('C4H6','C3H5-A'), ('C4H5-N','C3H6'), (None,)) ,'SUPFUEL'],
        [(('C4H6','C3H5-A'), ('C4H5-I','C3H6'), (None,)) ,'SUPFUEL'],
        [(('C4H6','C3H3'), ('C4H5-N','C3H4-A'), (None,)) ,'SUPFUEL'],
        [(('C4H6','C3H3'), ('C4H5-I','C3H4-A'), (None,)) ,'SUPFUEL'],
        [(('C4H6','C2H3'), ('C4H5-N','C2H4'), (None,)) ,'SUPFUEL'],
        [(('C4H6','C2H3'), ('C4H5-I','C2H4'), (None,)) ,'SUPFUEL'],
        [(('C4H5-N',), ('C4H5-I',), (None,)) ,'SUPFUEL'],
        [(('C4H5-N','OH'), ('C4H4','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H5-N','HO2'), ('C4H6','O2'), (None,)) ,'SUPFUEL'],
        [(('C4H5-N','HCO'), ('C4H6','CO'), (None,)) ,'SUPFUEL'],
        [(('C4H5-N','H'), ('C4H5-I','H'), (None,)) ,'SUPFUEL'],
        [(('C4H5-N','H'), ('C4H4','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H5-N','H2O2'), ('C4H6','HO2'), (None,)) ,'SUPFUEL'],
        [(('C4H5-I','OH'), ('C4H4','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H5-I','HO2'), ('C4H6','O2'), (None,)) ,'SUPFUEL'],
        [(('C4H5-I','HCO'), ('C4H6','CO'), (None,)) ,'SUPFUEL'],
        [(('C4H5-I','H'), ('C4H4','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H5-I','H2O2'), ('C4H6','HO2'), (None,)) ,'SUPFUEL'],
        [(('C4H5-2',), ('C4H5-I',), (None,)) ,'SUPFUEL'],
        [(('C4H5-2','O'), ('CH2O','C3H3'), (None,)) ,'SUPFUEL'],
        [(('C4H5-2','H'), ('C4H5-I','H'), (None,)) ,'SUPFUEL'],
        [(('C4H4','OH'), ('CH2O','C3H3'), (None,)) ,'SUPFUEL'],
        [(('C4H4','OH'), ('C4H3-N','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H4','OH'), ('C4H3-I','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H4','O'), ('C3H3','HCO'), (None,)) ,'SUPFUEL'],
        [(('C4H4','H'), ('C4H5-N',), (None,)) ,'SUPFUEL'],
        [(('C4H4','H'), ('C4H5-I',), (None,)) ,'SUPFUEL'],
        [(('C4H4','H'), ('C4H3-N','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H4','H'), ('C4H3-I','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H3-N',), ('C4H3-I',), (None,)) ,'SUPFUEL'],
        [(('C4H3-N','OH'), ('C4H2','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H3-N','H'), ('C4H4',), (None,)) ,'SUPFUEL'],
        [(('C4H3-N','H'), ('C4H3-I','H'), (None,)) ,'SUPFUEL'],
        [(('C4H3-N','H'), ('C4H2','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H3-N','H'), ('C2H2','H2CC'), (None,)) ,'SUPFUEL'],
        [(('C4H3-N','C2H3'), ('C3H3','C3H3'), (None,)) ,'SUPFUEL'],
        [(('C4H3-I','OH'), ('C4H2','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H3-I','H'), ('C4H4',), (None,)) ,'SUPFUEL'],
        [(('C4H3-I','H'), ('C4H2','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H3-I','H'), ('C2H2','H2CC'), (None,)) ,'SUPFUEL'],
        [(('C4H3-I','CH2'), ('C3H4-A','C2H'), (None,)) ,'SUPFUEL'],
        [(('C4H2','OH'), ('CO','C3H3'), (None,)) ,'SUPFUEL'],
        [(('C4H2','H'), ('C4H3-N',), (None,)) ,'SUPFUEL'],
        [(('C4H2','H'), ('C4H3-I',), (None,)) ,'SUPFUEL'],
        [(('C4H10',), ('SC4H9','H'), (None,)) ,'SUPFUEL'],
        [(('C4H10',), ('PC4H9','H'), (None,)) ,'SUPFUEL'],
        [(('C4H10',), ('C2H5','C2H5'), ('(+M)',)) ,'SUPFUEL'],
        [(('C4H10','PC4H9'), ('SC4H9','C4H10'), (None,)) ,'SUPFUEL'],
        [(('C4H10','OH'), ('SC4H9','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H10','OH'), ('PC4H9','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H10','O'), ('SC4H9','OH'), (None,)) ,'SUPFUEL'],
        [(('C4H10','O'), ('PC4H9','OH'), (None,)) ,'SUPFUEL'],
        [(('C4H10','O2'), ('SC4H9','HO2'), (None,)) ,'SUPFUEL'],
        [(('C4H10','O2'), ('PC4H9','HO2'), (None,)) ,'SUPFUEL'],
        [(('C4H10','HO2'), ('SC4H9','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C4H10','HO2'), ('PC4H9','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C4H10','H'), ('SC4H9','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H10','H'), ('PC4H9','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H10','C3H5-A'), ('SC4H9','C3H6'), (None,)) ,'SUPFUEL'],
        [(('C4H10','C3H5-A'), ('PC4H9','C3H6'), (None,)) ,'SUPFUEL'],
        [(('C4H10','C2H3'), ('SC4H9','C2H4'), (None,)) ,'SUPFUEL'],
        [(('C4H10','C2H3'), ('PC4H9','C2H4'), (None,)) ,'SUPFUEL'],
        [(('C3H8','OH'), ('NC3H7','H2O'), (None,)) ,'SUPFUEL'],
        [(('C3H8','OH'), ('IC3H7','H2O'), (None,)) ,'SUPFUEL'],
        [(('C3H8','O'), ('NC3H7','OH'), (None,)) ,'SUPFUEL'],
        [(('C3H8','O'), ('IC3H7','OH'), (None,)) ,'SUPFUEL'],
        [(('C3H8','O2'), ('NC3H7','HO2'), (None,)) ,'SUPFUEL'],
        [(('C3H8','O2'), ('IC3H7','HO2'), (None,)) ,'SUPFUEL'],
        [(('C3H8','IC3H7'), ('NC3H7','C3H8'), (None,)) ,'SUPFUEL'],
        [(('C3H8','HO2'), ('NC3H7','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C3H8','HO2'), ('IC3H7','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C3H8','H'), ('NC3H7','H2'), (None,)) ,'SUPFUEL'],
        [(('C3H8','H'), ('IC3H7','H2'), (None,)) ,'SUPFUEL'],
        [(('C3H8','C3H5-A'), ('NC3H7','C3H6'), (None,)) ,'SUPFUEL'],
        [(('C3H8','C3H5-A'), ('IC3H7','C3H6'), (None,)) ,'SUPFUEL'],
        [(('C3H8','C2H3'), ('NC3H7','C2H4'), (None,)) ,'SUPFUEL'],
        [(('C3H8','C2H3'), ('IC3H7','C2H4'), (None,)) ,'SUPFUEL'],
        [(('C3H6',), ('CC3H6',), (None,)) ,'SUPFUEL'],
        [(('C3H6',), ('C3H5-S','H'), (None,)) ,'SUPFUEL'],
        [(('C3H6',), ('C3H5-A','H'), (None,)) ,'SUPFUEL'],
        [(('C3H6','OH'), ('C3H5-T','H2O'), (None,)) ,'SUPFUEL'],
        [(('C3H6','OH'), ('C3H5-S','H2O'), (None,)) ,'SUPFUEL'],
        [(('C3H6','OH'), ('C3H5-A','H2O'), (None,)) ,'SUPFUEL'],
        [(('C3H6','O'), ('C3H5-T','OH'), (None,)) ,'SUPFUEL'],
        [(('C3H6','O'), ('C3H5-S','OH'), (None,)) ,'SUPFUEL'],
        [(('C3H6','O'), ('C3H5-A','OH'), (None,)) ,'SUPFUEL'],
        [(('C3H6','O'), ('C2H5','HCO'), (None,)) ,'SUPFUEL'],
        [(('C3H6','O2'), ('C3H5-T','HO2'), (None,)) ,'SUPFUEL'],
        [(('C3H6','O2'), ('C3H5-S','HO2'), (None,)) ,'SUPFUEL'],
        [(('C3H6','O2'), ('C3H5-A','HO2'), (None,)) ,'SUPFUEL'],
        [(('C3H6','HO2'), ('IC3H7','O2'), (None,)) ,'SUPFUEL'],
        [(('C3H6','HO2'), ('C3H5-T','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C3H6','HO2'), ('C3H5-S','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C3H6','HO2'), ('C3H5-A','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C3H6','H'), ('NC3H7',), (None,)) ,'SUPFUEL'],
        [(('C3H6','H'), ('IC3H7',), (None,)) ,'SUPFUEL'],
        [(('C3H6','H'), ('C3H5-T','H2'), (None,)) ,'SUPFUEL'],
        [(('C3H6','H'), ('C3H5-S','H2'), (None,)) ,'SUPFUEL'],
        [(('C3H6','H'), ('C3H5-A','H2'), (None,)) ,'SUPFUEL'],
        [(('C3H5-T',), ('C3H5-S',), (None,)) ,'SUPFUEL'],
        [(('C3H5-T','O2'), ('C3H4-A','HO2'), (None,)) ,'SUPFUEL'],
        [(('C3H5-T','HCO'), ('C3H6','CO'), (None,)) ,'SUPFUEL'],
        [(('C3H5-T','H'), ('C3H6',), (None,)) ,'SUPFUEL'],
        [(('C3H5-T','H'), ('C3H5-A','H'), (None,)) ,'SUPFUEL'],
        [(('C3H5-T','H'), ('C3H4-P','H2'), (None,)) ,'SUPFUEL'],
        [(('C3H5-T','CH2O'), ('C3H6','HCO'), (None,)) ,'SUPFUEL'],
        [(('C3H5-S','OH'), ('C2H4','HCO','H'), (None,)) ,'SUPFUEL'],
        [(('C3H5-S','O'), ('C2H4','HCO'), (None,)) ,'SUPFUEL'],
        [(('C3H5-S','HO2'), ('C2H4','HCO','OH'), (None,)) ,'SUPFUEL'],
        [(('C3H5-S','HCO'), ('C3H6','CO'), (None,)) ,'SUPFUEL'],
        [(('C3H5-S','H'), ('C3H4-P','H2'), (None,)) ,'SUPFUEL'],
        [(('C3H5-S','H'), ('C3H4-A','H2'), (None,)) ,'SUPFUEL'],
        [(('C3H5-S','CH2O'), ('C3H6','HCO'), (None,)) ,'SUPFUEL'],
        [(('C3H5-A',), ('C3H5-T',), (None,)) ,'SUPFUEL'],
        [(('C3H5-A',), ('C3H5-S',), (None,)) ,'SUPFUEL'],
        [(('C3H5-A','OH'), ('C3H4-A','H2O'), (None,)) ,'SUPFUEL'],
        [(('C3H5-A','O2'), ('C3H4-A','HO2'), (None,)) ,'SUPFUEL'],
        [(('C3H5-A','HCO'), ('C3H6','CO'), (None,)) ,'SUPFUEL'],
        [(('C3H5-A','H'), ('C3H4-A','H2'), (None,)) ,'SUPFUEL'],
        [(('C3H5-A','C4H71-3'), ('C3H6','C4H6'), (None,)) ,'SUPFUEL'],
        [(('C3H5-A','C3H5-A'), ('C3H4-A','C3H6'), (None,)) ,'SUPFUEL'],
        [(('C3H5-A','C2H5'), ('C2H4','C3H6'), (None,)) ,'SUPFUEL'],
        [(('C3H5-A','C2H3'), ('C3H4-A','C2H4'), (None,)) ,'SUPFUEL'],
        [(('C3H4-P',), ('C3H3','H'), (None,)) ,'SUPFUEL'],
        [(('C3H4-P','OH'), ('C3H3','H2O'), (None,)) ,'SUPFUEL'],
        [(('C3H4-P','O'), ('C3H3','OH'), (None,)) ,'SUPFUEL'],
        [(('C3H4-P','O'), ('C2H4','CO'), (None,)) ,'SUPFUEL'],
        [(('C3H4-P','O'), ('C2H3','HCO'), (None,)) ,'SUPFUEL'],
        [(('C3H4-P','O2'), ('C3H3','HO2'), (None,)) ,'SUPFUEL'],
        [(('C3H4-P','HO2'), ('C3H3','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C3H4-P','HO2'), ('C2H4','CO','OH'), (None,)) ,'SUPFUEL'],
        [(('C3H4-P','H'), ('C3H5-T',), (None,)) ,'SUPFUEL'],
        [(('C3H4-P','H'), ('C3H5-S',), (None,)) ,'SUPFUEL'],
        [(('C3H4-P','H'), ('C3H5-A',), (None,)) ,'SUPFUEL'],
        [(('C3H4-P','H'), ('C3H3','H2'), (None,)) ,'SUPFUEL'],
        [(('C3H4-P','C3H5-A'), ('C3H3','C3H6'), (None,)) ,'SUPFUEL'],
        [(('C3H4-P','C3H3'), ('C3H4-A','C3H3'), (None,)) ,'SUPFUEL'],
        [(('C3H4-P','C2H'), ('C2H2','C3H3'), (None,)) ,'SUPFUEL'],
        [(('C3H4-P','C2H3'), ('C3H3','C2H4'), (None,)) ,'SUPFUEL'],
        [(('C3H4-A',), ('C3H4-P',), (None,)) ,'SUPFUEL'],
        [(('C3H4-A',), ('C3H3','H'), (None,)) ,'SUPFUEL'],
        [(('C3H4-A','OH'), ('C3H3','H2O'), (None,)) ,'SUPFUEL'],
        [(('C3H4-A','O'), ('C2H4','CO'), (None,)) ,'SUPFUEL'],
        [(('C3H4-A','O'), ('C2H2','CH2O'), (None,)) ,'SUPFUEL'],
        [(('C3H4-A','O2'), ('C3H3','HO2'), (None,)) ,'SUPFUEL'],
        [(('C3H4-A','HO2'), ('C3H3','H2O2'), (None,)) ,'SUPFUEL'],
        [(('C3H4-A','HO2'), ('C2H4','CO','OH'), (None,)) ,'SUPFUEL'],
        [(('C3H4-A','H'), ('C3H5-T',), (None,)) ,'SUPFUEL'],
        [(('C3H4-A','H'), ('C3H5-S',), (None,)) ,'SUPFUEL'],
        [(('C3H4-A','H'), ('C3H5-A',), (None,)) ,'SUPFUEL'],
        [(('C3H4-A','H'), ('C3H4-P','H'), (None,)) ,'SUPFUEL'],
        [(('C3H4-A','H'), ('C3H3','H2'), (None,)) ,'SUPFUEL'],
        [(('C3H4-A','C3H5-A'), ('C3H3','C3H6'), (None,)) ,'SUPFUEL'],
        [(('C3H4-A','C3H4-A'), ('C3H5-A','C3H3'), (None,)) ,'SUPFUEL'],
        [(('C3H4-A','C2H'), ('C2H2','C3H3'), (None,)) ,'SUPFUEL'],
        [(('C3H3','OH'), ('H2CCC(S)','H2O'), (None,)) ,'SUPFUEL'],
        [(('C3H3','OH'), ('CH2O','C2H2'), (None,)) ,'SUPFUEL'],
        [(('C3H3','OH'), ('C3H2','H2O'), (None,)) ,'SUPFUEL'],
        [(('C3H3','OH'), ('C3H2(S)','H2O'), (None,)) ,'SUPFUEL'],
        [(('C3H3','OH'), ('C2H4','CO'), (None,)) ,'SUPFUEL'],
        [(('C3H3','OH'), ('C2H3','HCO'), (None,)) ,'SUPFUEL'],
        [(('C3H3','O'), ('CH2O','C2H'), (None,)) ,'SUPFUEL'],
        [(('C3H3','HO2'), ('OH','CO','C2H3'), (None,)) ,'SUPFUEL'],
        [(('C3H3','HCO'), ('C3H4-P','CO'), (None,)) ,'SUPFUEL'],
        [(('C3H3','HCO'), ('C3H4-A','CO'), (None,)) ,'SUPFUEL'],
        [(('C3H3','H'), ('H2CCC(S)','H2'), (None,)) ,'SUPFUEL'],
        [(('C3H3','H'), ('CC3H4',), (None,)) ,'SUPFUEL'],
        [(('C3H3','H'), ('C3H2C','H2'), (None,)) ,'SUPFUEL'],
        [(('C3H3','H'), ('C3H2','H2'), (None,)) ,'SUPFUEL'],
        [(('C3H3','H'), ('C3H2(S)','H2'), (None,)) ,'SUPFUEL'],
        [(('C3H3','CH'), ('C4H3-N','H'), (None,)) ,'SUPFUEL'],
        [(('C3H3','CH'), ('C4H3-I','H'), (None,)) ,'SUPFUEL'],
        [(('C3H3','CH2'), ('C4H4','H'), (None,)) ,'SUPFUEL'],
        [(('C3H2C','O2'), ('C2H2','CO2'), (None,)) ,'SUPFUEL'],
        [(('C3H2(S)',), ('C3H2',), ('+M',)) ,'SUPFUEL'],
        [(('C3H2(S)','H'), ('H2CCC(S)','H'), (None,)) ,'SUPFUEL'],
        [(('C2H5','C2H3'), ('C4H8-1',), ('(+M)',)) ,'SUPFUEL'],
        [(('C2H3','C4H71-3'), ('C2H4','C4H6'), (None,)) ,'SUPFUEL'],
        [(('C2H3','C2H3'), ('C4H6',), (None,)) ,'SUPFUEL'],
        [(('C2H3','C2H3'), ('C4H5-N','H'), (None,)) ,'SUPFUEL'],
        [(('C2H3','C2H3'), ('C4H5-I','H'), (None,)) ,'SUPFUEL'],
        [(('C2H3','C2H2'), ('C4H5-N',), (None,)) ,'SUPFUEL'],
        [(('C2H3','C2H2'), ('C4H5-I',), (None,)) ,'SUPFUEL'],
        [(('C2H3','C2H2'), ('C4H4','H'), (None,)) ,'SUPFUEL'],
        [(('C2H2','CH2'), ('C3H3','H'), (None,)) ,'SUPFUEL'],
        [(('C2H2','CH2(S)'), ('C3H3','H'), (None,)) ,'SUPFUEL'],
        [(('C2H2','C2H'), ('C4H3-N',), ('(+M)',)) ,'SUPFUEL'],
        [(('C2H2','C2H'), ('C4H3-I',), ('(+M)',)) ,'SUPFUEL'],
        [(('C2H2','C2H'), ('C4H2','H'), (None,)) ,'SUPFUEL'],]
    results = dict(results)

    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'heptane_cut_species.csv')
    mech_path = os.path.join(CWD, 'data', 'heptane_cut_mech.txt')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort with headers for species subset
    isolate_spc = ['CH4', 'deleteabove C2O'] # delete all with more than 2 carbons and one oxygen. more sophisticated: delete low temperature of heavier species
    sort_lst = ['submech_deletelarge', 0]

    param_dct_sort, _, cmts_dct, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    check_results = {}
    for rxn in param_dct_sort.keys():
        check_results[rxn] =  cmts_dct[rxn]['cmts_inline'].split('submech_deletelarge')[1].strip()

    check_results = dict(check_results)
    for key, val in results.items():
        assert val == check_results[key]

############ COMMENTED TO SPEED UP TESTS
    # option 2: do not sort according to submech_deletelarge, but e.g., according to pes. #
    # the list of reactions should be the same though
    # isolate_spc = ['CH4', 'deleteabove C2O']
    # sort_lst = ['pes', 0]
    # param_dct_sort2, _, _, _, _ = sorter.sorted_mech(
    #     spc_str, mech_str, isolate_spc, sort_lst)

    # assert all(key in param_dct_sort.keys() for key in param_dct_sort2.keys())
  ############ COMMENTED TO SPEED UP TESTS
    # option 2b: does not specify the stoichiometry, but classifies according to submech_deletelarge.
    # it will give the same result as +1C, +2O, which would be 'deleteabove C2H4O2' i.e., the default subfuel options.
    examples_keys_tocheck = [[(('SC2H4OH',), ('PC2H4OH',), (None,)) ,'CORE'],
        [(('SC2H2OH',), ('CH2CO','H'), (None,)) ,'CORE'],
        [(('SC2H2OH','O2'), ('CH2CO','HO2'), (None,)) ,'CORE'],
        [(('HCCO','H'), ('CH2(S)','CO'), (None,)) ,'CORE'],
        [(('H2CC','OH'), ('CH2CO','H'), (None,)) ,'CORE'],
        [(('CH3OCH3','H'), ('CH3OCH2','H2'), (None,)) ,'CORE'],
        [(('CH3OCH2','O2'), ('CH2O','CH2O','OH'), (None,)) ,'CORE'],
        [(('CH3CO','H'), ('CH2CO','H2'), (None,)) ,'CORE'],
        [(('CH3CHO','OH'), ('CH3CO','H2O'), (None,)) ,'CORE'],
        [(('CH2CO','OH'), ('HCCO','H2O'), (None,)) ,'CORE'],
        [(('CH2CHO','O2'), ('CH2CO','HO2'), (None,)) ,'CORE'],
        [(('CH2','CO'), ('CH2CO',), ('(+M)',)) ,'CORE'],
        [(('CH','CH2O'), ('H','CH2CO'), (None,)) ,'CORE'],
        [(('C2H5OH',), ('C2H5','OH'), (None,)) ,'CORE'],
        [(('C2H5O','O2'), ('CH3CHO','HO2'), (None,)) ,'CORE'],
        [(('C2H5','O'), ('CH3CHO','H'), (None,)) ,'CORE'],
        [(('C2H4O1-2',), ('CH3CHO',), (None,)) ,'CORE'],
        [(('C2H4','OH'), ('PC2H4OH',), (None,)) ,'CORE'],
        [(('C2H3OH',), ('CH3CHO',), (None,)) ,'CORE'],
        [(('C2H3O1-2',), ('CH3CO',), (None,)) ,'CORE'],
        [(('C2H3','O2'), ('CHCHO','OH'), (None,)) ,'CORE'],
        [(('C2H2','HO2'), ('CH2CO','OH'), (None,)) ,'CORE'],
        [(('C2H','OH'), ('H','HCCO'), (None,)) ,'CORE'],
        [(('TC3H6CHO','CH3'), ('IC3H5CHO','CH4'), (None,)) ,'FUEL'],
        [(('SC4H7OH-I','CH3'), ('IC4H6OH','CH4'), (None,)) ,'FUEL'],
        [(('C2H4O1-2','CH3'), ('C2H3O1-2','CH4'), (None,)) ,'FUEL'],
        [(('C2H3OH','CH3'), ('CH2CHO','CH4'), (None,)) ,'FUEL'],
        [(('C2H3CHOCH2','CH3'), ('C2H3','CH2CO','CH4'), (None,)) ,'FUEL'],
        [(('C2H3CHO','CH3'), ('C2H3CO','CH4'), (None,)) ,'FUEL'],
        [(('SC4H7OH-I','CH3O'), ('IC4H6OH','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('SC3H5OH','CH3O'), ('C2H3CHO','H','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3OCH3','CH3O'), ('CH3OCH2','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('CH3COCH3','CH3O'), ('CH3COCH2','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('C2H5COCH3','CH3O'), ('CH3CHCOCH3','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('C2H5CHO','CH3O'), ('C2H5CO','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('C2H4O1-2','CH3O'), ('C2H3O1-2','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('C2H3CHO','CH3O'), ('C2H3CO','CH3OH'), (None,)) ,'FUEL_ADD_O'],
        [(('SC4H7OH-I','CH3O2'), ('IC4H6OH','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('SC3H5OH','CH3O2'), ('C2H3CHO','H','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('C4H8O2-3','CH3O2'), ('CH2O','C3H5-A','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('C2H5OH','CH3O2'), ('C2H5O','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('C2H5COCH3','CH3O2'), ('CH3CHCOCH3','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('C2H3CHOCH2','CH3O2'), ('C2H3','CH2CO','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('C2H3CHO','CH3O2'), ('C2H3CO','CH3O2H'), (None,)) ,'FUEL_ADD_O2'],
        [(('SC3H5OH','CH3'), ('SC4H8OH-3',), (None,)) ,'FUEL_RAD'],
        [(('SC3H5CHO','H'), ('CH3','C2H3CHO'), (None,)) ,'FUEL_RAD'],
        [(('SC3H4OH',), ('CH2CO','CH3'), (None,)) ,'FUEL_RAD'],
        [(('IC4H8','OH'), ('SC3H5OH','CH3'), (None,)) ,'FUEL_RAD'],
        [(('IC4H8','OH'), ('IC3H5OH','CH3'), (None,)) ,'FUEL_RAD'],
        [(('IC4H8','OH'), ('CH3COCH3','CH3'), (None,)) ,'FUEL_RAD'],
        [(('C2H3OH','CH3'), ('C3H6OH1-1',), (None,)) ,'FUEL_RAD'],
        [(('C2H3CHO','CH3'), ('IC3H6CHO',), (None,)) ,'FUEL_RAD'],
        [(('C2H3CHO','CH3'), ('C3H6CHO-3',), (None,)) ,'FUEL_RAD'],
        [(('C2H5OH','C2H5'), ('SC2H4OH','C2H6'), (None,)) ,'R_CH3'],
        [(('C2H5OH','C2H5'), ('PC2H4OH','C2H6'), (None,)) ,'R_CH3'],
        [(('C2H5COCH3','C2H5'), ('CH3CHCOCH3','C2H6'), (None,)) ,'R_CH3'],
        [(('C2H5COCH3','C2H5'), ('CH2CH2COCH3','C2H6'), (None,)) ,'R_CH3'],
        [(('C2H5COCH3','C2H5'), ('C2H5COCH2','C2H6'), (None,)) ,'R_CH3'],
        [(('C2H5CHO','C2H5'), ('C2H5CO','C2H6'), (None,)) ,'R_CH3'],
        [(('TC4H9','CH3O2'), ('TC4H9O','CH3O'), (None,)) ,'R_O'],
        [(('SC4H9','CH3O2'), ('CH3O','SC4H9O'), (None,)) ,'R_O'],
        [(('C4H71-3','CH3O2'), ('C4H7O2-1','CH3O'), (None,)) ,'R_O'],
        [(('C4H71-3','CH3O2'), ('C4H71-O','CH3O'), (None,)) ,'R_O'],
        [(('C3H6','CH2OH'), ('PC4H8OH-4',), (None,)) ,'R_O'],
        [(('C3H6','CH2OH'), ('PC4H8OH-3',), (None,)) ,'R_O'],
        [(('C3H5-A','CH3O2'), ('C3H5O','CH3O'), (None,)) ,'R_O'],
        [(('C3H4-A','CH2OH'), ('IC4H6OH',), (None,)) ,'R_O'],
        [(('C2H5','CH3O2'), ('C2H5O','CH3O'), (None,)) ,'R_O'],
        [(('C2H4','CH3O2'), ('C2H4O1-2','CH3O'), (None,)) ,'R_O'],
        [(('C2H3COCH3','CH3O2'), ('CH2CHO','CH3CO','CH3O'), (None,)) ,'R_O'],
        [(('CH3OCH3','O2CHO'), ('CH3OCH2','HO2CHO'), (None,)) ,'R_O3-H'],
        [(('TC4H9O','O2'), ('IC4H8O','HO2'), (None,)) ,'SUPFUEL'],
        [(('TC3H6CHO','H2'), ('IC3H7CHO','H'), (None,)) ,'SUPFUEL'],
        [(('TC3H6CHO','CH2O'), ('IC3H7CHO','HCO'), (None,)) ,'SUPFUEL'],
        [(('SC4H9','HO2'), ('SC4H9O','OH'), (None,)) ,'SUPFUEL'],
        [(('SC4H8OH-1',), ('PC4H8OH-4',), (None,)) ,'SUPFUEL'],
        [(('SC4H8OH-1',), ('PC4H8OH-3',), (None,)) ,'SUPFUEL'],
        [(('SC3H5OCH2-1',), ('C3H6CHO-2',), (None,)) ,'SUPFUEL'],
        [(('SC3H5OCH2-1',), ('C3H6','HCO'), (None,)) ,'SUPFUEL'],
        [(('PC3H4OH-1','O2'), ('CH3CHCO','HO2'), (None,)) ,'SUPFUEL'],
        [(('NC3H7CO',), ('NC3H7','CO'), (None,)) ,'SUPFUEL'],
        [(('NC3H7CHO','OH'), ('NC3H7CO','H2O'), (None,)) ,'SUPFUEL'],
        [(('NC3H7','HO2'), ('NC3H7O','OH'), (None,)) ,'SUPFUEL'],
        [(('NC3H7','CH2O'), ('PC4H9O',), (None,)) ,'SUPFUEL'],
        [(('IC4H9O','OH'), ('IC3H7CHO','H2O'), (None,)) ,'SUPFUEL'],
        [(('IC4H8O',), ('IC3H7CHO',), (None,)) ,'SUPFUEL'],
        [(('IC4H8','OH'), ('IC4H7OH','H'), (None,)) ,'SUPFUEL'],
        [(('IC4H7O',), ('IC3H5CHO','H'), (None,)) ,'SUPFUEL'],
        [(('IC4H7O',), ('C3H6','HCO'), (None,)) ,'SUPFUEL'],
        [(('IC4H7O',), ('C3H5-T','CH2O'), (None,)) ,'SUPFUEL'],
        [(('IC4H7O','OH'), ('IC3H5CHO','H2O'), (None,)) ,'SUPFUEL'],
        [(('IC4H7O','O'), ('IC3H5CHO','OH'), (None,)) ,'SUPFUEL'],
        [(('IC4H7-I1','O2'), ('CH3COCH3','HCO'), (None,)) ,'SUPFUEL'],
        [(('IC4H7','OH'), ('IC4H7OH',), (None,)) ,'SUPFUEL'],
        [(('IC4H7','O'), ('IC3H5CHO','H'), (None,)) ,'SUPFUEL'],
        [(('IC4H7','HO2'), ('IC4H7O','OH'), (None,)) ,'SUPFUEL'],
        [(('IC4H7','HO2'), ('IC3H5CHO','H2O'), (None,)) ,'SUPFUEL'],
        [(('IC4H6OH','IC4H8'), ('IC4H7OH','IC4H7'), (None,)) ,'SUPFUEL'],
        [(('IC3H7CHO','H'), ('IC3H7CO','H2'), (None,)) ,'SUPFUEL'],
        [(('IC3H7','O'), ('CH3COCH3','H'), (None,)) ,'SUPFUEL'],
        [(('IC3H7','HO2'), ('IC3H7O','OH'), (None,)) ,'SUPFUEL'],
        [(('IC3H5CO',), ('C3H5-T','CO'), (None,)) ,'SUPFUEL'],
        [(('IC3H5CHO','OH'), ('IC3H5CO','H2O'), (None,)) ,'SUPFUEL'],
        [(('H2C4O','OH'), ('CH2CO','HCCO'), (None,)) ,'SUPFUEL'],
        [(('H2C4O','H'), ('C2H2','HCCO'), (None,)) ,'SUPFUEL'],
        [(('CH3COCH3','OH'), ('CH3COCH2','H2O'), (None,)) ,'SUPFUEL'],
        [(('CH3CHCHO',), ('CH3CHCO','H'), (None,)) ,'SUPFUEL'],
        [(('CH3CHCHO',), ('C2H3CHO','H'), (None,)) ,'SUPFUEL'],
        [(('CH3CHCHO','H2'), ('C2H5CHO','H'), (None,)) ,'SUPFUEL'],
        [(('CH2O','IC3H7'), ('IC4H9O',), (None,)) ,'SUPFUEL'],
        [(('CH2CHOCH2',), ('CH2CH2CHO',), (None,)) ,'SUPFUEL'],
        [(('CH2CH2COCH3',), ('C2H3','CH3CHO'), (None,)) ,'SUPFUEL'],
        [(('CH2CCH2OH','H'), ('C3H5OH',), (None,)) ,'SUPFUEL'],
        [(('CH2CCH2OH','H2O2'), ('C3H5OH','HO2'), (None,)) ,'SUPFUEL'],
        [(('CC4H8O','OH'), ('CH2O','C3H5-A','H2O'), (None,)) ,'SUPFUEL'],
        [(('CC4H8O','O'), ('CH2O','C3H5-A','OH'), (None,)) ,'SUPFUEL'],
        [(('CC4H8O','HO2'), ('CH2O','C3H5-A','H2O2'), (None,)) ,'SUPFUEL'],
        [(('CC4H8O','H'), ('CH2O','C3H5-A','H2'), (None,)) ,'SUPFUEL'],
        [(('C8H132-6,SA','HO2'), ('C8H132-6,SAO','OH'), (None,)) ,'SUPFUEL'],
        [(('C8H132-6,PA','HO2'), ('C8H132-6,PAO','OH'), (None,)) ,'SUPFUEL'],
        [(('C8H131-5,3-4,TA','HO2'), ('C8H131-5,3-4,TAO','OH'), (None,)) ,'SUPFUEL'],
        [(('C8H131-5,3,TA','HO2'), ('C8H131-5,3,TAO','OH'), (None,)) ,'SUPFUEL'],
        [(('C8H131-5,3,SA','HO2'), ('C8H131-5,3,SAO','OH'), (None,)) ,'SUPFUEL'],
        [(('C8H131-5,3,PA','HO2'), ('C8H131-5,3,PAO','OH'), (None,)) ,'SUPFUEL'],
        [(('C7H111-5,3,6P','CH2O'), ('C8H131-5,3,PAO',), (None,)) ,'SUPFUEL'],
        [(('C7H111-5,1P','CH2O'), ('C8H132-6,PAO',), (None,)) ,'SUPFUEL'],
        [(('C4H8O2-3','OH'), ('CH2O','C3H5-A','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H8O1-2','H'), ('CH2O','C3H5-A','H2'), (None,)) ,'SUPFUEL'],
        [(('C4H8-2','OH'), ('SC4H8OH-3',), (None,)) ,'SUPFUEL'],
        [(('C4H8-2','O'), ('CH2CO','C2H5','H'), (None,)) ,'SUPFUEL'],
        [(('C4H8-2','O'), ('C2H5CHCO','H','H'), (None,)) ,'SUPFUEL'],
        [(('C4H7O2-1',), ('C3H6','HCO'), (None,)) ,'SUPFUEL'],
        [(('C4H7O2-1',), ('C3H5-S','CH2O'), (None,)) ,'SUPFUEL'],
        [(('C4H72-2OH','OH'), ('C4H63,1-3OH','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H71-O',), ('C2H4','CH3CO'), (None,)) ,'SUPFUEL'],
        [(('C4H71-4OH','OH'), ('C4H64,2-1OH','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H71-4','HO2'), ('C4H7O1-4','OH'), (None,)) ,'SUPFUEL'],
        [(('C4H71-3','C2H3COCH3'), ('C8H131-5,3,TAO',), (None,)) ,'SUPFUEL'],
        [(('C4H71-2OH','OH'), ('C4H63,1-2OH','H2O'), (None,)) ,'SUPFUEL'],
        [(('C4H2','OH'), ('H2C4O','H'), (None,)) ,'SUPFUEL'],
        [(('C4H10','C2H5O'), ('SC4H9','C2H5OH'), (None,)) ,'SUPFUEL'],
        [(('C3H6O1-3',), ('C2H4','CH2O'), (None,)) ,'SUPFUEL'],
        [(('C3H3','O2'), ('CH2CO','HCO'), (None,)) ,'SUPFUEL'],
        [(('C2HCHO',), ('C2H2','CO'), (None,)) ,'SUPFUEL'],
        [(('C2H5COCH3','OH'), ('CH3CHCOCH3','H2O'), (None,)) ,'SUPFUEL'],
        [(('C2H5COCH2',), ('CH2CO','C2H5'), (None,)) ,'SUPFUEL'],
        [(('C2H5CHO','OH'), ('C2H5CO','H2O'), (None,)) ,'SUPFUEL'],
        [(('C2H3CHOCH2',), ('C4H6O23',), (None,)) ,'SUPFUEL'],
        [(('C2H3CHO','O'), ('C2H3CO','OH'), (None,)) ,'SUPFUEL'],
        [(('C2H','CH2O'), ('C3H3O',), (None,)) ,'SUPFUEL'],
        [(('AC3H5OCH2',), ('C2H3COCH3','H'), (None,)) ,'SUPFUEL'],
]
    examples_keys_tocheck = dict(examples_keys_tocheck)

    isolate_spc = ['CH4'] # equivalent to ['CH4', 'deleteabove C2H4O2']
    sort_lst = ['submech_deletelarge', 0]
    param_dct_sort2b, _, cmts_dct, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    assert all(key in param_dct_sort2b.keys() for key in examples_keys_tocheck.keys())
    assert all(key not in param_dct_sort.keys() for key in examples_keys_tocheck.keys())

    # option 4: same deleteabove C2H4O2 as option2b, but without fuel specification and submech option:
    # it will exclude some of the reactions of the fuel similar to deleteabove C2O, but fewer
    checks_in_4_and_2b_but_not_in_1 = [(('CH', 'CO'), ('HCCO',), ('+M',)),
        (('HCCO', 'H'), ('CH2(S)', 'CO'), (None,)),
        (('CH2', 'CO'), ('CH2CO',), ('(+M)',)),
        (('CH2CO', 'H'), ('HCCO', 'H2'), (None,)),
        (('CH2CO', 'H'), ('CH3', 'CO'), (None,)),
        (('CH2CHO',), ('CH3', 'CO'), ('(+M)',)),
        (('CH2CHO',), ('CH2CO', 'H'), ('(+M)',)),
        (('CH', 'CH2O'), ('H', 'CH2CO'), (None,)),
        (('C2H3O1-2',), ('CH3CO',), (None,)),
        (('C2H3O1-2',), ('CH2CHO',), (None,)),
        (('C2H2', 'OH'), ('HCCOH', 'H'), (None,)),
        (('C2H2', 'OH'), ('CH2CO', 'H'), (None,)),
        (('C2H2', 'OH'), ('C2H2OH',), (None,)),
        (('CH3CO', 'H'), ('CH2CO', 'H2'), (None,)),
        (('CH3CHO',), ('CH4', 'CO'), ('(+M)',)),
        (('CH3CHO',), ('CH3', 'HCO'), ('(+M)',)),
        (('C2H4O1-2',), ('CH3CHO',), (None,)),
        (('C2H4O1-2',), ('CH3', 'HCO'), (None,)),
        (('CH3CHO', 'OH'), ('CH3CO', 'H2O'), (None,)),
        (('CH3CHO', 'OH'), ('CH3', 'HOCHO'), (None,)),
        (('CH3OCH2', 'O2'), ('CH2O', 'CH2O', 'OH'), (None,)),
        (('CH3CHO', 'HO2'), ('CH3CO', 'H2O2'), (None,)),
        (('C2H5O', 'O2'), ('CH3CHO', 'HO2'), (None,)),
        (('C2H4O1-2', 'HO2'), ('C2H3O1-2', 'H2O2'), (None,)),
        (('CH3CHCHO',), ('CH3CHCO', 'H'), (None,)),
        (('CH3CHCHO',), ('C2H3CHO', 'H'), (None,)),
        (('C3H6', 'OH'), ('C3H5OH', 'H'), (None,)),
        (('C3H6', 'OH'), ('C2H3OH', 'CH3'), (None,)),
        (('C3H5OH', 'H'), ('CH2CCH2OH', 'H2'), (None,)),
        (('C2H5CHO', 'H'), ('NC3H7O',), (None,)),
        (('C3H5-A', 'O2'), ('CH3CO', 'CH2O'), (None,)),
        (('C3H5-A', 'O2'), ('C2H3CHO', 'OH'), (None,)),
        (('C3H4-A', 'HO2'), ('CH2CO', 'CH2', 'OH'), (None,)),
        (('CH3OCH3', 'O2CHO'), ('CH3OCH2', 'HO2CHO'), (None,)),
        (('H2C4O', 'H'), ('C2H2', 'HCCO'), (None,)),
        (('C2H3CHOCH2',), ('C4H6O23',), (None,)),
        (('SC3H5OCH2-1',), ('SC3H5CHO', 'H'), (None,)),
        (('IC3H7', 'CO'), ('IC3H7CO',), (None,)),
        (('C4H7O2-1',), ('SC3H5OCH2-1',), (None,)),
        (('C3H6CHO-2',), ('SC3H5CHO', 'H'), (None,)),
        (('C2H5CHCO', 'H'), ('NC3H7', 'CO'), (None,)),
        (('C2H5CHCO', 'H'), ('C3H6CHO-3',), (None,)),
        (('CH2CCH2OH', 'CH3'), ('IC4H7OH',), (None,)),
        (('IC4H8', 'OH'), ('SC4H7OH-I', 'H'), (None,)),
        (('C4H8O1-3', 'H'), ('CH2O', 'C3H5-A', 'H2'), (None,)),
        (('C3H6', 'CH2OH'), ('PC4H8OH-3',), (None,)),
        (('C3H5OH', 'CH3'), ('CH2CCH2OH', 'CH4'), (None,)),
        (('C2H5COCH3', 'H'), ('SC4H9O',), (None,)),
        (('C4H3-I', 'O2'), ('HCCO', 'CH2CO'), (None,)),
        (('IC4H7-I1', 'O2'), ('IC3H5CHO', 'OH'), (None,)),
        (('C4H71-3', 'HO2'), ('C4H7O2-1', 'OH'), (None,)),
        (('SC3H5OH', 'CH3O'), ('C2H3CHO', 'H', 'CH3OH'), (None,)),
        (('IC3H7CHO', 'OH'), ('TC3H6CHO', 'H2O'), (None,)),
        (('C4H72-2OH', 'OH'), ('C4H63,1-3OH', 'H2O'), (None,)),
        (('C2H3COCH3', 'HO2'), ('CH2CO', 'C2H3', 'H2O2'), (None,)),
        (('NC3H7CHO', 'HO2'), ('NC3H7CO', 'H2O2'), (None,)),
        (('IC4H9O', 'O2'), ('IC3H7CHO', 'HO2'), (None,)),
        (('CC4H8O', 'HO2'), ('CH2O', 'C3H5-A', 'H2O2'), (None,)),
        (('C4H71-1OH', 'HO2'), ('C4H63,1-1OH', 'H2O2'), (None,)),
        (('C3H6O1-3', 'CH3O2'), ('CH2O', 'C2H3', 'CH3O2H'), (None,)),
        (('C2H3CHOCH2', 'CH3'), ('C2H3', 'CH2CO', 'CH4'), (None,)),
        (('TC3H6CHO', 'CH3'), ('IC3H5CHO', 'CH4'), (None,)),
        (('IC4H7O', 'CH3'), ('IC3H5CHO', 'CH4'), (None,)),
        (('SC4H7OH-I', 'CH3'), ('IC4H6OH', 'CH4'), (None,)),
        (('NC3H7CHO', 'CH3'), ('NC3H7CO', 'CH4'), (None,)),
        (('IC4H8O', 'CH3'), ('IC3H6CHO', 'CH4'), (None,)),
        (('C2H5COCH3', 'CH3'), ('CH3CHCOCH3', 'CH4'), (None,)),
        (('C2H5COCH3', 'CH3'), ('CH2CH2COCH3', 'CH4'), (None,)),
        (('C2H5COCH3', 'CH3'), ('C2H5COCH2', 'CH4'), (None,)),
        (('C2H5CHO', 'C2H5'), ('C2H5CO', 'C2H6'), (None,)),
        (('IC4H9O', 'CH3'), ('IC3H7CHO', 'CH4'), (None,)),
        (('TC3H6CHO', 'CH2O'), ('IC3H7CHO', 'HCO'), (None,)),
        (('SC3H5CHO', 'C2H3'), ('CH2CHCHCHO', 'C2H4'), (None,)),
        (('C2H5COCH3', 'C2H3'), ('CH3CHCOCH3', 'C2H4'), (None,)),
        (('IC4H6OH', 'IC4H8'), ('IC4H7OH', 'IC4H7'), (None,)),
        (('C8H131-5,3,PA', 'HO2'), ('C8H131-5,3,PAO', 'OH'), (None,)),
        ]
    isolate_spc = ['deleteabove C2H4O2']
    sort_lst = ['pes', 0]
    param_dct_sort4, _, _, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    for key in param_dct_sort4.keys(): # compare with less restrictive
        assert key in param_dct_sort2b.keys()
    for key in param_dct_sort.keys(): #compare with more restrictive
        assert key in param_dct_sort4.keys()

    assert all(key in param_dct_sort2b.keys() and key not in param_dct_sort.keys() for key in checks_in_4_and_2b_but_not_in_1)
    assert all(key in param_dct_sort4.keys() and key not in param_dct_sort.keys() for key in checks_in_4_and_2b_but_not_in_1)

def test__sortby_submech_keepsubfuel():
    """ test mechanalyzer.parser.sort

        sort by fuel submechanism: extract reactions of
        fuel, fuel radicals, R+O2, R+O4
        and also the relative submech
    """
    must_be_in_results = [
        [(('C2H4',), ('H2', 'H2CC'), ('(+M)',)), 'FUEL'],
        [(('CH2(S)', 'C2H4'), ('CC3H6',), (None,)), 'FUEL'],
        [(('C2H4', 'CH3O'), ('C2H3', 'CH3OH'), (None,)), 'FUEL'],
        [(('C3H5-A', 'C2H5'), ('C2H4', 'C3H6'), (None,)), 'FUEL'],
        [(('C3H8', 'O2'), ('IC3H7', 'HO2'), (None,)), 'FUEL_ADD_CH3'],
        [(('C2H5CHCO', 'OH'), ('NC3H7', 'CO2'), (None,)), 'FUEL_ADD_CH3'],
        [(('C2H5', 'O2'), ('C2H4O1-2', 'OH'), (None,)), 'FUEL_ADD_H'],
        [(('C4H71-1',), ('C2H5', 'C2H2'), (None,)), 'FUEL_ADD_H'],
        [(('C2H3OH', 'H'), ('PC2H4OH',), (None,)), 'FUEL_ADD_O'],
        [(('C2H3OH', 'HO2'), ('CH3CHO', 'HO2'), (None,)), 'FUEL_ADD_O'],
        [(('C4H6', 'O'), ('C2H2', 'C2H4O1-2'), (None,)), 'FUEL_ADD_O'],
        [(('CH3OCHO', 'O2'), ('CH2OCHO', 'HO2'), (None,)), 'FUEL_ADD_O2'],
        [(('CH3', 'CH2O'), ('C2H5O',), (None,)), 'FUEL_ADD_OH'],
        [(('C2H5OH', 'O2'), ('SC2H4OH', 'HO2'), (None,)), 'FUEL_ADD_OH'],
        [(('CH3OCH3', 'CH3O2'), ('CH3OCH2', 'CH3O2H'), (None,)), 'FUEL_ADD_OH'],
        [(('C2H3', 'CH3'), ('CH4', 'C2H2'), (None,)), 'FUEL_RAD'],
        [(('C4H71-O',), ('C2H3', 'CH3CHO'), (None,)), 'FUEL_RAD'],
        [(('C3H6', 'OH'), ('IC3H5OH', 'H'), (None,)), 'R_CH3'],
        [(('C4H8-2', 'H'), ('C3H6', 'CH3'), (None,)), 'R_CH3'],
        [(('C2H5CHCO', 'O'), ('C3H6', 'CO2'), (None,)), 'R_CH3'],
        [(('SC2H2OH', 'O2'), ('CH2CO', 'HO2'), (None,)), 'R_O'],
        [(('C2H3OO',), ('CH2CO', 'OH'), (None,)), 'R_O2'],
        [(('O', 'O'), ('O2',), ('+M',)), 'SUBFUEL'],
        [(('CH', 'H'), ('C', 'H2'), (None,)), 'SUBFUEL'],
        [(('CH3O',), ('CH2O', 'H'), ('(+M)',)), 'SUBFUEL'],
        [(('CH4', 'O'), ('CH3', 'OH'), (None,)), 'SUBFUEL'],
        [(('CH2(S)', 'O2'), ('CO', 'H2O'), (None,)), 'SUBFUEL'],
        [(('CH3', 'HO2'), ('CH3O', 'OH'), (None,)), 'SUBFUEL'],
        [(('CH2O', 'HO2'), ('OCH2O2H',), (None,)), 'SUBFUEL'],
        [(('CH3O2', 'H2O2'), ('CH3O2H', 'HO2'), (None,)), 'SUBFUEL'],
        [(('C2H2', 'OH'), ('HCCOH', 'H'), (None,)), 'SUBFUEL'],
        [(('CH2(S)', 'CO2'), ('CH2O', 'CO'), (None,)), 'SUBFUEL'],
        [(('CH3O', 'CH3O'), ('CH3OH', 'CH2O'), (None,)), 'SUBFUEL'],
        [(('C2H', 'CH3'), ('C3H4-P',), (None,)), 'SUBFUEL'],
        [(('C3H2C', 'O2'), ('C2H2', 'CO2'), (None,)), 'SUBFUEL'],
    ]
    must_be_in_results = dict(must_be_in_results)
    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'heptane_cut_species.csv')
    mech_path = os.path.join(CWD, 'data', 'heptane_cut_mech.txt')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort with headers for species subset
    isolate_spc = ['C2H4', 'keepbelow C2H6O2']
    sort_lst = ['submech_keepsubfuel', 0]

    param_dct_sort, _, cmts_dct, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    check_results = {}
    for rxn in param_dct_sort.keys():
        check_results[rxn] =  cmts_dct[rxn]['cmts_inline'].split('submech_keepsubfuel')[1].strip()

    check_results = dict(check_results)
    for key, val in must_be_in_results.items():
        assert val == check_results[key]

############ COMMENTED TO SPEED UP TESTS
    # option 2: do not sort according to submech_keepsubfuel, but e.g., according to pes. #
    # the list of reactions should be the same though
    # isolate_spc = ['C2H4', 'keepbelow C2H6O2']
    # sort_lst = ['pes', 0]
    # param_dct_sort2, _, _, _, _ = sorter.sorted_mech(
    #     spc_str, mech_str, isolate_spc, sort_lst)

    # for key in param_dct_sort2.keys():
    #     assert key in param_dct_sort

############ COMMENTED TO SPEED UP TESTS
    # option 2b: does not specify the keepbelow, but classifies according to submech_keepsubfuel.
    # it will give the same result as +2H, +2O are the default subfuel options.
    # isolate_spc = ['C2H4']
    # sort_lst = ['submech_keepsubfuel', 0]
    # param_dct_sort2b, _, _, _, _ = sorter.sorted_mech(
    #     spc_str, mech_str, isolate_spc, sort_lst)

    # for key in param_dct_sort2b.keys():
    #     assert key in param_dct_sort

    # option 3: only specify the stoichiometry, so does NOT filter according to the fuel.
    # all rxns will be classified as "SUBFUEL". there won't be as many reactions as those kept above,
    #because many are those of the fuel submech.
    # keys to check here will not be found in the example below but only in that above
    examples_keys_tocheck = [
        (('C3H8', 'C2H3'), ('IC3H7', 'C2H4'), (None,)),
        (('C3H5-A', 'C2H5'), ('C2H4', 'C3H6'), (None,)),
        (('NC3H7CO',), ('NC3H7', 'CO'), (None,)),
        (('NC3H7', 'O2'), ('C3H6', 'HO2'), (None,)),
        (('IC4H8', 'O'), ('IC3H7', 'HCO'), (None,)),
        (('IC4H10',), ('CH3', 'IC3H7'), ('(+M)',)),
        (('IC3H7', 'OH'), ('C3H6', 'H2O'), (None,)),
        (('IC3H7', 'O'), ('CH3COCH3', 'H'), (None,)),
        (('IC3H7', 'O2'), ('IC3H7O2',), (None,)),
        (('IC3H7', 'HO2'), ('IC3H7O', 'OH'), (None,)),
        (('IC3H7', 'HCO'), ('IC3H7CHO',), (None,)),
        (('IC3H6CO', 'OH'), ('IC3H7', 'CO2'), (None,)),
        (('CH2O', 'IC3H7'), ('IC4H9O',), (None,)),
        (('C4H8-1', 'O'), ('NC3H7', 'HCO'), (None,)),
        (('C4H10',), ('NC3H7', 'CH3'), ('(+M)',)),
        (('C3H8', 'OH'), ('NC3H7', 'H2O'), (None,)),
        (('C3H8', 'OH'), ('IC3H7', 'H2O'), (None,)),
        (('C3H8', 'O'), ('NC3H7', 'OH'), (None,)),
        (('C3H8', 'CH3O'), ('NC3H7', 'CH3OH'), (None,)),
        (('C3H8', 'CH3O2'), ('NC3H7', 'CH3O2H'), (None,)),
        (('C3H8', 'CH3'), ('NC3H7', 'CH4'), (None,)),
        (('C3H8', 'C2H5O2'), ('IC3H7', 'C2H5O2H'), (None,)),
        (('C3H6', 'HO2'), ('IC3H7', 'O2'), (None,)),
        (('C3H6', 'H'), ('NC3H7',), (None,)),
        (('C2H5CHCO', 'OH'), ('NC3H7', 'CO2'), (None,)),
        (('SC3H5CHO',), ('C3H6', 'CO'), (None,)),
        (('SC3H5CHO', 'H'), ('C3H6', 'HCO'), (None,)),
        (('PC4H9',), ('C3H6', 'CH3'), (None,)),
        (('NC3H7O2',), ('C3H6', 'HO2'), (None,)),
        (('IC4H7-I1', 'OH'), ('C3H6', 'HCO', 'H'), (None,)),
        (('IC4H7-I1', 'O'), ('C3H6', 'HCO'), (None,)),
        (('IC4H7-I1', 'HO2'), ('C3H6', 'HCO', 'OH'), (None,)),
        (('IC3H7O2',), ('C3H6', 'HO2'), (None,)),
        (('C3H6OOH1-2',), ('C3H6', 'HO2'), (None,)),
        (('C3H6CHO-2',), ('C3H6', 'HCO'), (None,)),
        (('C3H6',), ('CC3H6',), (None,)),
        (('C3H6',), ('C3H5-S', 'H'), (None,)),
        (('C3H6',), ('C3H5-A', 'H'), (None,)),
        (('C3H6', 'OH'), ('SC3H5OH', 'H'), (None,)),
        (('C3H6', 'O'), ('CH3CHCO', 'H', 'H'), (None,)),
        (('C3H6', 'O'), ('C3H5-T', 'OH'), (None,)),
        (('C3H6', 'O2'), ('C3H5-T', 'HO2'), (None,)),
        (('C3H6', 'O2'), ('C3H5-S', 'HO2'), (None,)),
        (('C3H6', 'HO2'), ('C3H6OOH2-1',), (None,)),
        (('C3H6', 'HO2'), ('C3H5-A', 'H2O2'), (None,)),
        (('C3H6', 'H'), ('C3H5-A', 'H2'), (None,)),
        (('C3H6', 'CH3O'), ('C3H5-A', 'CH3OH'), (None,)),
        (('C3H6', 'CH3O2'), ('C3H5-A', 'CH3O2H'), (None,)),
        (('C3H6', 'CH3'), ('C3H5-A', 'CH4'), (None,)),
        (('C3H6', 'CH2OH'), ('PC4H8OH-4',), (None,)),
        (('C3H6', 'CH2OH'), ('PC4H8OH-3',), (None,)),
        (('C3H6', 'C2H5O2'), ('C3H5-A', 'C2H5O2H'), (None,)),
        (('C3H5-T', 'CH2O'), ('C3H6', 'HCO'), (None,)),
        (('C3H5-S', 'HCO'), ('C3H6', 'CO'), (None,)),
        (('C3H5-A', 'HCO'), ('C3H6', 'CO'), (None,)),
        (('C2H5CHCO', 'O'), ('C3H6', 'CO2'), (None,)),
        (('AC3H5OCH2',), ('C3H6', 'HCO'), (None,)),
    ]

    isolate_spc = ['keepbelow C2H6O2']
    sort_lst = ['submech_keepsubfuel', 0]
    param_dct_sort3, _, _, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    assert param_dct_sort3 != param_dct_sort
    for key in examples_keys_tocheck:
        assert key in param_dct_sort.keys()
        assert key not in param_dct_sort3.keys()

############ COMMENTED TO SPEED UP TESTS
    # option 4: as above, but does not specify classification
    # isolate_spc = ['keepbelow C2H6O2']
    # sort_lst = ['pes', 0]
    # param_dct_sort4, _, _, _, _ = sorter.sorted_mech(
    #     spc_str, mech_str, isolate_spc, sort_lst)

    # for key in param_dct_sort4.keys():
    #     assert key in param_dct_sort3

def test__sortby_submech_class():
    """ test mechanalyzer.parser.sort

        sort by fuel submechanism: extract reactions of
        fuel, fuel radicals, R+O2, R+O4
        then order by subpes and broad class
    """
    results = [
        [(('IC8',), ('NEOC5H11', 'IC3H7'), (None,)),
         '  FUEL.2.1.Bond fission'],
        [(('IC8', 'O2'), ('IC8-1R', 'HO2'), (None,)),
         '  FUEL.4.2.H abstraction'],
        [(('IC8OOH1',), ('IC8-1OR', 'OH'), (None,)),
         '  FUEL_ADD_O2.4.1.Bond fission +OH'],
        [(('IC8-1O2R', 'HO2'), ('IC8OOH1', 'O2'), (None,)),
         '  FUEL_ADD_O2.6.1.Recombination-decomposition - termination'],
        [(('IC8-1O2R', 'H2O2'), ('IC8OOH1', 'HO2'), (None,)),
         '  FUEL_ADD_O2.7.1.H abstraction'],
        [(('IC8-1R',), ('IC8-5R',), (None,)),
         '  FUEL_RAD.1.1.Isomerization'],
        [(('IC8-1R', 'O2'), ('IC8-1O2R',), (None,)),
         '  FUEL_RAD.3.1.Recombination O2'],
        [(('IC8-3R', 'O2'), ('IC8D3', 'HO2'), (None,)),
         '  FUEL_RAD.3.3.Recombination-decomposition - termination'],
        [(('IC8-1R', 'CH3O2'), ('IC8-1OR', 'CH3O'), (None,)),
         '  FUEL_RAD.8.1.Recombination-decomposition - propagation'],
        [(('IC8OOH1-1AR',), ('IC4H7OOH', 'IC4H9'), (None,)),
         '  R_O2.3.1.Beta-scission'],
        [(('IC8OOH1-1AR',), ('IC8O1-1A', 'OH'), (None,)),
         '  R_O2.3.1.Beta-scission +OH'],
        [(('IC8OOH1-1AR',), ('CH2O', 'I24C7D1', 'OH'), (None,)),
         '  R_O2.3.1.Decomposition(lumped)'],
        [(('IC8-1O2R',), ('IC8OOH1-1AR',), (None,)),
         '  R_O2.3.1.Isomerization'],
        [(('IC8-3O2R',), ('IC8D3', 'HO2'), (None,)),
         '  R_O2.3.2.Beta-scission +HO2'],
        [(('IC8OOH1-1AR', 'O2'), ('IC8OOH1-1AO2R',), (None,)),
         '  R_O2.5.1.Recombination O2'],
    ]

    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_IC8_red_mech.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort with headers for species subset
    isolate_spc = ['IC8']
    sort_lst = ['submech', 'subpes', 'rxn_class_broad', 0]

    param_dct_sort, _, cmts_dct, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    sorted_results = []
    for rxn in param_dct_sort.keys():
        sorted_results.append(
            [rxn, cmts_dct[rxn]['cmts_inline'].split('rxntype')[1]])

    assert results == sorted_results


def test__sortby_submech_prompt():
    """ test mechanalyzer.parser.sort

        sort by prompt reactions identified
        based on radical type
    """
    results = [
        [(('C4H72-2',), ('C4H612', 'H'), (None,)), '111.1.1.RAD_DECO_C4H72-2'],
        [(('C4H71-4', 'H'), ('C4H8-1',), ('(+M)',)), '112.1.1.RAD_GEN_C4H71-4'],
        [(('C4H8-1', 'H'), ('C4H71-3', 'H2'), (None,)), '113.8.25.RAD_GEN_C4H71-3'],
        [(('C4H71-3', 'O'), ('C2H3CHO', 'CH3'), (None,)), '120.23.68.RAD_GEN_C4H71-3'],
        [(('C4H8-1', 'O'), ('C4H71-3', 'OH'), (None,)), '121.5.11.RAD_GEN_C4H71-3'],
        [(('C4H8-1', 'OH'), ('C4H71-3', 'H2O'), (None,)), '122.21.38.RAD_GEN_C4H71-3'],
        [(('C4H71-4O2',), ('C4H61-3OOH4',), (None,)), '128.2.19.'],
        [(('C4H71-3OOH',), ('CH3CHO', 'C2H3', 'OH'), (None,)), '129.3.6.'],
        [(('C4H8-1', 'HO2'), ('C4H71-3', 'H2O2'), (None,)), '130.22.61.RAD_GEN_C4H71-3'],
        [(('C4H8-1', 'CH3'), ('C4H71-3', 'CH4'), (None,)), '153.6.6.RAD_GEN_C4H71-3'],
        [(('C4H71-3', 'CH3O'), ('C4H8-1', 'CH2O'), (None,)), '160.1.1.RAD_GEN_C4H71-3'],
        [(('C4H8-1', 'CH3O'), ('C4H71-3', 'CH3OH'), (None,)), '161.5.5.RAD_GEN_C4H71-3'],
        [(('C4H71-3OOCH3',), ('C4H71-O', 'CH3O'), (None,)), '166.1.1.'],
        [(('C4H8-1', 'CH3O2'), ('C4H71-3', 'CH3O2H'), (None,)), '167.7.7.RAD_GEN_C4H71-3'],
        [(('C6H101-3,3',), ('C2H3', 'C4H72-2'), (None,)), '183.1.1.RAD_GEN_C4H72-2'],
        [(('C2H5', 'C4H71-3'), ('C4H6', 'C2H6'), (None,)), '184.1.1.RAD_GEN_C4H71-3'],
        [(('C4H71-3', 'C2H5O2'), ('C4H71-O', 'C2H5O'), (None,)), '192.1.1.RAD_GEN_C4H71-3'],
        [(('C4H8-1', 'C2H5O2'), ('C4H71-3', 'C2H5O2H'), (None,)), '193.2.2.RAD_GEN_C4H71-3'],
        [(('C4H8-1', 'CH3CO3'), ('C4H71-3', 'CH3CO3H'), (None,)), '196.2.2.RAD_GEN_C4H71-3'],
        [(('C3H5-A', 'C4H71-3'), ('C3H6', 'C4H6'), (None,)), '206.1.1.RAD_GEN_C4H71-3'],
        [(('C4H8-1', 'C3H5-A'), ('C4H71-3', 'C3H6'), (None,)), '207.1.1.RAD_GEN_C4H71-3'],
        [(('IC3H7O2', 'C4H71-3'), ('IC3H7O', 'C4H71-O'), (None,)), '212.1.1.RAD_GEN_C4H71-3'],
        [(('C4H8-1', 'IC3H7O2'), ('C4H71-3', 'IC3H7O2H'), (None,)), '213.3.3.RAD_GEN_C4H71-3'],
        [(('C4H6', 'C4H71-3'), ('C8H131-5,3,PA',), (None,)), '221.1.1.RAD_GEN_C4H71-3'],
        [(('C4H71-3', 'C4H71-3'), ('C8H141-5,3',), (None,)), '222.1.1.RAD_GEN_C4H71-3'],
        [(('C4H71-3', 'C2H3COCH3'), ('C8H131-5,3,TAO',), (None,)), '225.1.1.RAD_GEN_C4H71-3'],
        [(('C4H71-4O2', 'C4H71-3'), ('C4H7O1-4', 'C4H71-O'), (None,)), '228.2.4.RAD_GEN_C4H71-3'],
        [(('IC4H9O2', 'C4H71-3'), ('IC4H9O', 'C4H71-O'), (None,)), '230.1.1.RAD_GEN_C4H71-3'],
        [(('IC4H9O2', 'C4H8-1'), ('IC4H9O2H', 'C4H71-3'), (None,)), '231.1.1.RAD_GEN_C4H71-3'],
    ]

    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'heptane_cut_species.csv')
    mech_path = os.path.join(CWD, 'data', 'heptane_cut_mech.txt')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort with headers for species subset
    isolate_spc = ['C4H71-3', 'C4H71-4', 'C4H72-2']
    sort_lst = ['submech_prompt', 0]

    param_dct_sort, _, cmts_dct, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    print('Sort by submech_prompt: check only 1st channel of each PES for simplicity')
    sorted_results = []
    pess = []
    for rxn in param_dct_sort.keys():
        cmt = cmts_dct[rxn]['cmts_inline'].split('submech_prompt')[1].strip()
        pes = cmt.split('.')[0]
        if pes not in pess:
            sorted_results.append(
                [rxn, cmt])
            pess.append(pes)

    assert results == sorted_results


def test__filter_pesgroups():
    """ test mechanalyzer.parser.sort

        sort by fuel submechanism: extract reactions of
        fuel, fuel radicals, R+O2, R+O4
        and also the relative submech
        then order by subpes and broad class
    """
    results = [
        {'grp': 1, 'idxs': ['113:8', '113:9', '113:11', '111:1'], 'peds': [['C4H8-1+H=C4H71-3+H2'], ['C4H8-1+H=C4H71-4+H2'], [
            'C4H8-2+H=C4H72-2+H2'], []], 'hot': [[], [], [], ['C4H71-3', 'C4H71-4', 'C4H72-2']], 'modeltype': 'rovib_dos'},
        {'grp': 2, 'idxs': ['121:6', '111:1'], 'peds': [['C4H8-1+O=C4H71-4+OH'], []], 'hot': [[], ['C4H71-4']], 'modeltype': 'rovib_dos'},
        {'grp': 3, 'idxs': ['122:21', '122:22', '122:23', '122:24', '111:1'], 'peds': [['C4H8-1+OH=C4H71-3+H2O'], ['C4H8-1+OH=C4H71-4+H2O'], [
            'C4H8-2+OH=C4H71-3+H2O'], ['C4H8-2+OH=C4H72-2+H2O'], []], 'hot': [[], [], [], [], ['C4H71-3', 'C4H71-4', 'C4H72-2']],
         'modeltype': 'rovib_dos'},
        {'grp': 4, 'idxs': ['153:7', '153:9', '111:1'], 'peds': [['C4H8-1+CH3=C4H71-4+CH4'], ['C4H8-2+CH3=C4H72-2+CH4'], []],
         'hot': [[], [], ['C4H71-4', 'C4H72-2']], 'modeltype': 'rovib_dos'},
        {'grp': 5, 'idxs': ['161:6', '111:1'], 'peds': [['C4H8-1+CH3O=C4H71-4+CH3OH'], []], 'hot': [[], ['C4H71-4']], 'modeltype': 'rovib_dos'}
        ]

    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'heptane_cut_species.csv')
    mech_path = os.path.join(CWD, 'data', 'heptane_cut_mech.txt')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)
    therm_str = pathtools.read_file(os.path.join(CWD, 'data'), 'therm.dat')
    spc_therm_dct = ckin_parser.parse_spc_therm_dct(therm_str, [300, 1000, 1500, 2000])

    # Sort with headers for species subset
    isolate_spc = ['C4H71-3', 'C4H71-4', 'C4H72-2']
    sort_lst = ['submech_prompt', 0]

    _, _, _, pes_groups, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst, spc_therm_dct=spc_therm_dct, dct_flt_grps={'DH':30., 'lookforpromptchains': 0}) #

    assert results == pes_groups

# Helper function


def _read_files(spc_path, mech_path, sort_path):
    """ read file names
    """

    spc_str, mech_str, sort_str = '', '', ''

    if spc_path is not None:
        with open(spc_path, encoding='utf-8') as fobj:
            spc_str = fobj.read()
    if mech_path is not None:
        with open(mech_path, encoding='utf-8') as fobj:
            mech_str = fobj.read()
    if sort_path is not None:
        with open(sort_path, encoding='utf-8') as fobj:
            sort_str = fobj.read()

    return spc_str, mech_str, sort_str


if __name__ == '__main__':
    test__sortby_submech_deletelarge() # 580 s
    test__sort_submech() # 70 s
    test__sortby_submech_keepsubfuel() # 713 s
    test__sortby_submech_prompt() # 136 s
    test__filter_pesgroups() #165 s
    test__sortby_submech_class() #34 s




