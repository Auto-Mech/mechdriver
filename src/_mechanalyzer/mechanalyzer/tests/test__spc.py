""" species dict function tests
"""

import os
import ioformat
import mechanalyzer
from ioformat import pathtools

# Set paths
PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')

# Set information for running tests
SMILES = ('CC', 'CC(O)Cl')

HEADERS = ('smiles', 'inchi', 'mult', 'charge')


def test__csv_io():
    """ test mechanalyzer.parser.spc.dct
        test mechanalyzer.parser.spc.dct_to_csv_str
    """

    ref_csv_str = ioformat.pathtools.read_file(DAT_PATH, 'spc3.csv')

    spc_dct = mechanalyzer.parser.spc.build_spc_dct(ref_csv_str, 'csv')
    csv_str = mechanalyzer.parser.spc.csv_string(spc_dct, HEADERS)

    assert ref_csv_str == csv_str


# modify/add functionality
def test__mod_spc_dct_atomcount():
    """ test mechanalyzer.parser.new_spc.reorder_by_atomcount
    """

    ref_spc_dct = {
        'O2': {'inchi': 'InChI=1S/O2/c1-2'},
        'H2': {'inchi': 'InChI=1S/H2/h1H'},
        'OH': {'inchi': 'InChI=1S/HO/h1H'},
        'CH4': {'inchi': 'InChI=1S/CH4/h1H4'},
        'CH3OH': {'inchi': 'InChI=1S/CH4O/c1-2/h2H,1H3'},
        'CH3SH': {'inchi': 'InChI=1S/CH4S/c1-2/h2H,1H3'},
        'C2H6': {'inchi': 'InChI=1S/C2H6/c1-2/h1-2H3'},
        'C2H5NH2': {'inchi': 'InChI=1S/C2H6/c1-2/h1-2H3'},
        'C2H5OH': {'inchi': 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'},
        'C3H8': {'inchi': 'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3'}
    }

    spc_dct = {
        'C2H6': {'inchi': 'InChI=1S/C2H6/c1-2/h1-2H3'},
        'O2': {'inchi': 'InChI=1S/O2/c1-2'},
        'C3H8': {'inchi': 'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3'},
        'CH3OH': {'inchi': 'InChI=1S/CH4O/c1-2/h2H,1H3'},
        'C2H5OH': {'inchi': 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'},
        'CH4': {'inchi': 'InChI=1S/CH4/h1H4'},
        'H2': {'inchi': 'InChI=1S/H2/h1H'},
        'CH3SH': {'inchi': 'InChI=1S/CH4S/c1-2/h2H,1H3'},
        'OH': {'inchi': 'InChI=1S/HO/h1H'},
        'C2H5NH2': {'inchi': 'InChI=1S/C2H6/c1-2/h1-2H3'},
    }

    spc_dct = mechanalyzer.parser.spc.reorder_by_atomcount(ref_spc_dct)
    assert ref_spc_dct == spc_dct


def test__mod_spc_dct_hof_basis():
    """ test mechanalyzer.parser.new_spc.add_heat_of_formation_basis
    """

    spc_dct = {
        'C4H9OH': {'inchi': 'InChI=1S/C4H10O/c1-2-3-4-5/h5H,2-4H2,1H3', 
                   'canon_enant_ich': 'InChI=1S/C4H10O/c1-2-3-4-5/h5H,2-4H2,1H3'},
        'CH4': {'inchi': 'InChI=1S/CH4/h1H4', 
                'canon_enant_ich': 'InChI=1S/CH4/h1H4'},
        'H2O': {'inchi': 'InChI=1S/H2O/h1H2', 
                'canon_enant_ich': 'InChI=1S/H2O/h1H2'},
    }

    ref_spc_dct = {
        'C4H9OH': {
            'inchi': 'InChI=1S/C4H10O/c1-2-3-4-5/h5H,2-4H2,1H3',
            'canon_enant_ich': 'InChI=1S/C4H10O/c1-2-3-4-5/h5H,2-4H2,1H3'},
        'CH4': {
            'inchi': 'InChI=1S/CH4/h1H4',
            'canon_enant_ich': 'InChI=1S/CH4/h1H4'},
        'H2O': {
            'inchi': 'InChI=1S/H2O/h1H2',
            'canon_enant_ich': 'InChI=1S/H2O/h1H2'},
        'cbh0_[H][H]': {
            'smiles': '[H][H]',
            'inchi': 'InChI=1S/H2/h1H',
            'canon_enant_ich': 'InChI=1S/H2/h1H',
            'inchikey': 'UFHFLCQGNIYNRP-UHFFFAOYSA-N',
            'charge': 0, 'mult': 1, 
	        'fml': {'H', 2},
            'mc_nsamp': (True, 3, 1, 3, 100, 12),
            'hind_inc': 0.5235987755982988,
            'hbond_cutoffs': (4.55, 1.92)},
        'cbh1_CO': {
            'smiles': 'CO',
            'inchi': 'InChI=1S/CH4O/c1-2/h2H,1H3',
            'canon_enant_ich': 'InChI=1S/CH4O/c1-2/h2H,1H3', 
            'inchikey': 'OKKJLVBELUTLKV-UHFFFAOYSA-N',
            'charge': 0, 'mult': 1, 
	        'fml': {'C': 1, 'H': 4, 'O': 1}, 
            'mc_nsamp': (True, 3, 1, 3, 100, 12),
            'hind_inc': 0.5235987755982988,
            'hbond_cutoffs': (4.55, 1.92)},
        'cbh1_CC': {
            'smiles': 'CC',
            'inchi': 'InChI=1S/C2H6/c1-2/h1-2H3',
	        'canon_enant_ich': 'InChI=1S/C2H6/c1-2/h1-2H3', 
            'inchikey': 'OTMSDBZUPAUEDD-UHFFFAOYSA-N',
            'charge': 0, 'mult': 1,
	        'fml': {'C': 2, 'H': 6}, 
            'mc_nsamp': (True, 3, 1, 3, 100, 12),
            'hind_inc': 0.5235987755982988,
            'hbond_cutoffs': (4.55, 1.92)},
        'cbh2_CCC': {
            'smiles': 'CCC',
            'inchi': 'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3',
	        'canon_enant_ich': 'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3', 
            'inchikey': 'ATUOYWHBWRKTHZ-UHFFFAOYSA-N',
            'charge': 0, 'mult': 1,
            'fml': {'C': 3, 'H': 8}, 
            'mc_nsamp': (True, 3, 1, 3, 100, 12),
            'hind_inc': 0.5235987755982988,
            'hbond_cutoffs': (4.55, 1.92)},
        'cbh2_CCO': {
            'smiles': 'CCO',
            'inchi': 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
	        'canon_enant_ich': 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', 
            'inchikey': 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N',
            'charge': 0, 'mult': 1,
	        'fml': {'C': 2, 'H': 6, 'O': 1}, 
            'mc_nsamp': (True, 3, 1, 3, 100, 12),
            'hind_inc': 0.523598775598298,
            'hbond_cutoffs': (4.55, 1.92)}
    }

    spc_dct = mechanalyzer.parser.spc.add_heat_of_formation_basis(
        spc_dct, ref_schemes=('cbh0', 'cbh1', 'cbh2'))
    spc_dct2 = mechanalyzer.parser.spc.add_heat_of_formation_basis(
        spc_dct, ref_schemes=('cbh0', 'cbh1', 'cbh2'))

    assert set(ref_spc_dct.keys()) == set(spc_dct.keys())
    assert set(ref_spc_dct.keys()) == set(spc_dct2.keys())
    for name, ref_dct in ref_spc_dct.items():
        assert set(ref_dct.keys()) == set(spc_dct[name].keys())
        assert set(ref_dct.keys()) == set(spc_dct2[name].keys())


def test__mod_spc_dct_stereo():
    """ test mechanalyzer.parser.new_spc.stereochemical_spc_dct
    """

    spc_dct = {
        'CC(O)Cl': {'inchi': 'InChI=1S/C2H5ClO/c1-2(3)4/h2,4H,1H3'},
        'CC=CC': {'inchi': 'InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+'},
        'OC=CN': {'inchi': 'InChI=1S/C2H5NO/c3-1-2-4/h1-2,4H,3H2'},
        'CC': {'inchi': 'InChI=1S/C2H6/c1-2/h1-2H3'},
    }

    ref_spc_dct = {
        'CC(O)Cl': {
            'inchi': 'InChI=1S/C2H5ClO/c1-2(3)4/h2,4H,1H3/t2-/m0/s1',
            'inchikey': 'KJESGYZFVCIMDE-REOHCLBHSA-N'
        },
        'CC=CC': {
            'inchi': 'InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+',
            'inchikey': 'IAQRGUVFOMOMEM-ONEGZZNKSA-N'
        },
        'OC=CN': {
            'inchi': 'InChI=1S/C2H5NO/c3-1-2-4/h1-2,4H,3H2/b2-1-',
            'inchikey': 'UEVZOFFYIPNZNW-UPHRSURJSA-N'
        },
        'CC': {
            'inchi': 'InChI=1S/C2H6/c1-2/h1-2H3',
            'inchikey': 'OTMSDBZUPAUEDD-UHFFFAOYSA-N'
        }
    }

    spc_dct = mechanalyzer.parser.spc.stereochemical_spc_dct(
        ref_spc_dct, all_stereo=False)

    assert ref_spc_dct == spc_dct

def test_lumped_species():
    """ generate species name of type SPECIES-LUMPED-X if multiple entries are found
    """
    lumped_list = [
        'C10H10-LUMPED-1',
        'C10H6(C2H)2-LUMPED-1',
        'C10H7CH3-LUMPED-1',
        'C14H10-LUMPED-1',
        'C18H12-LUMPED-1',
        'C18H12-LUMPED-2',
        'C18H12-LUMPED-3',
        'C2H4O2-LUMPED-1',
        'C4H3-LUMPED-1',
        'C4H4-LUMPED-1',
        'C4H5-LUMPED-1',
        'C4H5-LUMPED-2',
        'C4H5-LUMPED-3',
        'C4H5-LUMPED-4',
        'C4H5-LUMPED-5',
        'C5H5CH3-LUMPED-1',
        'C5H5CH3-LUMPED-2',
        'C5H5O-LUMPED-1',
        'C5H5O-LUMPED-2',
        'C5H5OH-LUMPED-1',
        'C5H5OH-LUMPED-2',
        'C6H11-LUMPED-1',
        'C6H4O2-LUMPED-1',
        'C6H5C3H2-LUMPED-1',
        'C6H5C3H3-LUMPED-1',
        'C9H6CH2-LUMPED-1',
        'C9H6CH3-LUMPED-1',
        'C9H6CH3-LUMPED-2',
        'C9H6CH3-LUMPED-3',
        'C9H7CH3-LUMPED-1',
        'C9H7CH3-LUMPED-2',
        'C9H7CH3-LUMPED-3',
        'CH3C6H4-LUMPED-1',
        'CH3C6H4-LUMPED-2',
        'CRESOL-LUMPED-1',
        'CRESOL-LUMPED-2',
        'CYC5H7-LUMPED-1',
        'FC10H10-LUMPED-1',
        'FC10H10-LUMPED-2',
        'GUAIACOL-LUMPED-1',
        'GUAIACOL-LUMPED-2',
        'LC5H8-LUMPED-1',
        'NC5H10-LUMPED-1',
        'RXYLENE-LUMPED-1',
        'SALICALD-LUMPED-1',
        'XYLENE-LUMPED-1',
        'XYLENE-LUMPED-2',
    ]
    
    spc_dct = mechanalyzer.parser.new_spc.load_mech_spc_dct(
        'SpeciesCRECK.csv', DAT_PATH, quotechar='"', canon_ent=False, allow_multiple=True)
    assert all(spc in spc_dct.keys() for spc in lumped_list)


def test_parse_fct_grp():
    """ read functional group dictionary, parse it, convert to dictionary, and back
    """
    groups_check = {"ACENAPH":{'A1-M': 2},
        "BENZINDENE":{'C5-M': 1, 'A1-M': 2},
        "BENZOFLUORENE":{'C5-M': 1, 'A1-M': 3},
        "BENZOPYRENE":{'A1-M': 5},
        "BIN1A":{'A1-M': 3},
        "BIN1B":{'A1-M': 5},
        "BINAPH":{'A1-M': 4},
        "BIPHENYL":{'A1-M': 2},
        "BZCOOH":{'A1-M': 1},
        "BZFUR":{'FUR-M': 1, 'A1-M': 1},
        "C10H10":{'A1-M': 1},
        "C10H10-LUMPED-1":{'A1-M': 1},
        "C10H6(C2H)2":{'A1-M': 2, 'A1,C2H-M': 2},
        "C10H6(C2H)2-LUMPED-1":{'A1-M': 2, 'A1,C2H-M': 2},
        "C10H6CH3":{'A1-R': 1, 'A1,CH3-R': 1},
        "C10H7":{'A1-R': 1},
        "C10H7C2H":{'A1-M': 2, 'A1,C2H-M': 1},
        "C10H7C3H5":{'A1-M': 2, 'A1,C2H3-M': 1},
        "C10H7C6H5":{'A1-M': 3},
        "C10H7CH2":{'A1CH2-RSR': 1},
        "C10H7CH2C6H57H7":{'A1-M': 3},
        "C10H7CH3":{'A1-M': 2, 'A1,CH3-M': 1},
        "C10H7CH3-LUMPED-1":{'A1-M': 2, 'A1,CH3-M': 1},
        "C10H7CHO":{'A1-M': 2, 'A1,CHO-M': 1},
        "C10H7O":{'A1O-RSR': 1},
        "C10H7OH":{'A1-M': 2, 'A1,OH-M': 1},
        "C10H8":{'A1-M': 2},
        "C12H8":{'A1-M': 2, 'A1,C2H3-M': 1},
        "C12H9":{'A1-R': 1},
        "C13H8CH2":{'C5CH2-M': 1, 'A1-M': 2},
        "C14H10":{'A1-M': 3},
        "C14H10-LUMPED-1":{'A1-M': 3},
        "C14H12":{'A1-M': 2},
        "C14H9":{'A1-R': 1},
        "C14H9CH3":{'A1-M': 3, 'A1,CH3-M': 1},
        "C15H10":{'C5-M': 1, 'A1-M': 3},
        "C16H10":{'A1-M': 4},
        "C16H9":{'A1-R': 1},
        "C18H10":{'A1-M': 4},
        "C18H12":{'A1-M': 4},
        "C18H12-LUMPED-1":{'A1-M': 4},
        "C18H12-LUMPED-2":{'A1-M': 4},
        "C18H12-LUMPED-3":{'A1-M': 4},
        "C18H14":{'C5-M': 2, 'A1-M': 2},
        "C18H9":{'A1-R': 1},
        "C2HC6H4C6H4C2H":{'A1-M': 2, 'A1,C2H-M': 2},
        "C5H4CH2":{'C5CH2-M': 1},
        "C5H4O":{'C5O-M': 1},
        "C5H4O2":{'FUR-M': 1},
        "C5H4OH":{'C5-RSR': 1},
        "C5H5":{'C5-RSR': 1},
        "C5H5CH3":{'C5-M': 1, 'C5,CH3-M': 1},
        "C5H5CH3-LUMPED-1":{'C5-M': 1, 'C5,CH3-M': 1},
        "C5H5CH3-LUMPED-2":{'C5-M': 1, 'C5,CH3-M': 1},
        "C5H5O":{'C5-M': 1},
        "C5H5O-LUMPED-1":{'C5O-RSR': 1},
        "C5H5O-LUMPED-2":{'C5O-RSR': 1},
        "C5H5OH":{'C5-M': 1, 'C5,OH-M': 1},
        "C5H5OH-LUMPED-1":{'C5-M': 1, 'C5,OH-M': 1},
        "C5H5OH-LUMPED-2":{'C5-M': 1, 'C5,OH-M': 1},
        "C5H6":{'C5-M': 1},
        "C6H4O2":{'A1-M': 1},
        "C6H4O2-LUMPED-1":{'A1-M': 1},
        "C6H4OH":{'A1-R': 1, 'A1,OH-R': 1},
        "C6H5":{'A1-R': 1},
        "C6H5C2H":{'A1-M': 1, 'A1,C2H-M': 1},
        "C6H5C2H3":{'A1-M': 1, 'A1,C2H3-M': 1},
        "C6H5C2H4C6H5":{'A1-M': 2},
        "C6H5C2H5":{'A1-M': 1},
        "C6H5C3H2":{'A1CH2-RSR': 1},
        "C6H5C3H2-LUMPED-1":{'A1CH2-RSR': 1},
        "C6H5C3H3":{'A1-M': 1, 'A1,C3.DD-M': 1},
        "C6H5C3H3-LUMPED-1":{'A1-M': 1, 'A1,C3.ST-M': 1},
        "C6H5C3H7":{'A1-M': 1},
        "C6H5C4H5":{'A1-M': 1, 'A1,C2H3-M': 1},
        "C6H5C4H9":{'A1-M': 1},
        "C6H5CCC6H5":{'A1-M': 2, 'A1,C2H-M': 1},
        "C6H5CCCH3":{'A1-M': 1, 'A1,C2H-M': 1, 'A1,C3.ST-M': 1},
        "C6H5CH2C6H5":{'A1-M': 2},
        "C6H5CH2OH":{'A1-M': 1},
        "C6H5CHCH3":{'A1CH2-RSR': 1},
        "C6H5CHO":{'A1-M': 1, 'A1,CHO-M': 1},
        "C6H5CO":{'A1CH2-RSR': 1},
        "C6H5O":{'A1O-RSR': 1},
        "C6H5OCH3":{'A1-M': 1, 'A1,OCH3-M': 1},
        "C6H5OH":{'A1-M': 1, 'A1,OH-M': 1},
        "C6H6":{'A1-M': 1},
        "C7H5":{'C5-RSR': 1},
        "C7H7":{'A1CH2-RSR': 1},
        "C7H8":{'A1-M': 1, 'A1,CH3-M': 1},
        "C9H6CH2":{'C5CH2-M': 1, 'A1-M': 1},
        "C9H6CH2-LUMPED-1":{'C5CH2-M': 1, 'A1-M': 1},
        "C9H6CH3":{'C5-RSR': 1, 'C5,CH3-RSR': 1},
        "C9H6CH3-LUMPED-3":{'A1CH2-RSR': 2},
        "C9H6O":{'C5O-M': 1, 'A1-M': 1},
        "C9H6O-M":{'C5O-M': 1, 'A1-M': 1},
        "C9H6OH":{'C5-RSR': 1},
        "C9H6OH-M":{'C5-RSR': 1},
        "C9H7CH3":{'C5-M': 1, 'A1-M': 1, 'C5,CH3-M': 1},
        "C9H7CH3-LUMPED-1":{'A1-M': 1, 'A1,C2H3-M': 1},
        "C9H7CH3-LUMPED-2":{'C5-M': 1, 'A1-M': 1, 'C5,CH3-M': 1},
        "C9H7CH3-LUMPED-3":{'C5-M': 1, 'A1-M': 1, 'C5,CH3-M': 1},
        "C9H7O":{'C5O-RSR': 1},
        "C9H7O-M":{'C5O-RSR': 1},
        "C9H7OH":{'C5-M': 1, 'A1-M': 1, 'C5,OH-M': 1},
        "C9H7OH-M":{'C5-M': 1, 'A1-M': 1, 'C5,OH-M': 1},
        "CATECHOL":{'A1-M': 1, 'A1,OH-M': 2, 'A1,OH,OH-M': 1},
        "CH3C10H6O":{'A1O-RSR': 1},
        "CH3C10H6OH":{'A1-M': 2, 'A1,CH3-M': 1, 'A1,OH-M': 1, 'A1,OH,CH3-M': 1},
        "CH3C6H4":{'A1-R': 1, 'A1,CH3-R': 1},
        "CH3C6H4-LUMPED-1":{'A1-R': 1, 'A1,CH3-R': 1},
        "CH3C6H4-LUMPED-2":{'A1-R': 1, 'A1,CH3-R': 1},
        "CRESOL":{'A1-M': 1, 'A1,CH3-M': 1, 'A1,OH-M': 1, 'A1,OH,CH3-M': 1},
        "CRESOL-LUMPED-1":{'A1-M': 1, 'A1,CH3-M': 1, 'A1,OH-M': 1, 'A1,OH,CH3-M': 1},
        "CRESOL-LUMPED-2":{'A1-M': 1, 'A1,CH3-M': 1, 'A1,OH-M': 1, 'A1,OH,CH3-M': 1},
        "CYC5H7-LUMPED-1":{'C5H2-RSR': 1},
        "DIBZFUR":{'FUR-M': 1, 'A1-M': 2},
        "FC10H10":{'C5-M': 2},
        "FC10H10-LUMPED-1":{'C5-M': 2},
        "FC10H10-LUMPED-2":{'C5-M': 2},
        "FLUORANTHENE":{'A1-M': 3},
        "FLUORENE":{'C5-M': 1, 'A1-M': 2},
        "FURAN":{'FUR-M': 1},
        "GUAIACOL":{'A1-M': 1, 'A1,OH-M': 1, 'A1,OH,OCH3-M': 1, 'A1,OCH3-M': 1},
        "GUAIACOL-LUMPED-1":{'A1-M': 1, 'A1,OH-M': 1, 'A1,OH,OCH3-M': 1, 'A1,OCH3-M': 1},
        "GUAIACOL-LUMPED-2":{'A1-M': 1, 'A1,OH-M': 1, 'A1,OH,OCH3-M': 1, 'A1,OCH3-M': 1},
        "HOC6H4CH2":{'A1CH2-RSR': 1},
        "INDANE":{'A1-M': 1},
        "INDENE":{'C5-M': 1, 'A1-M': 1},
        "INDENYL":{'C5-RSR': 1},
        "METHYLCHRYSENE":{'A1-M': 4, 'A1,CH3-M': 1},
        "OC6H4CH2":{'A1-M': 1},
        "OC6H4CH3":{'A1O-RSR': 1},
        "OC6H4OH":{'A1O-RSR': 1, 'A1O,OH-RSR': 1},
        "PHENALENE":{'A1-M': 2},
        "RBBENZ":{'A1CH2-RSR': 1},
        "RTETRALIN":{'A1CH2-RSR': 1},
        "RXYLENE":{'A1CH2-RSR': 1},
        "RXYLENE-LUMPED-1":{'A1CH2-RSR': 1},
        "SALICALD":{'A1-M': 1, 'A1,OH-M': 1, 'A1,OH,CHO-M': 1, 'A1,CHO-M': 1},
        "SALICALD-LUMPED-1":{'A1-M': 1, 'A1,OH-M': 1, 'A1,OH,CHO-M': 1, 'A1,CHO-M': 1},
        "STILB":{'A1-M': 2, 'A1,C2H3-M': 1},
        "TETRALIN":{'A1-M': 1},
        "TMBENZ":{'A1-M': 1, 'A1,CH3-M': 3},
        "XYLENE":{'A1-M': 1, 'A1,CH3-M': 2},
        "XYLENE-LUMPED-1":{'A1-M': 1, 'A1,CH3-M': 2},
        "XYLENE-LUMPED-2":{'A1-M': 1, 'A1,CH3-M': 2},
        "m-TERPH":{'A1-M': 3}}
    
    spc_dct = mechanalyzer.parser.new_spc.load_mech_spc_dct(
        'SpeciesCRECK_wfctgrp.csv', DAT_PATH, quotechar="'", canon_ent=False, allow_multiple=True)

    for spc, vals in spc_dct.items():
        if vals['fct_grp'] != {}:
                assert vals['fct_grp'] == groups_check[spc]
    headers = ['smiles', 'inchi', 'fct_grp', 'fml']
    sortd_csv_str = mechanalyzer.parser.spc.csv_string(spc_dct, headers)
    new_spc_str = pathtools.read_file(
        DAT_PATH, 'SpeciesCRECK_wfctgrp.csv', print_debug=True)

    assert sortd_csv_str == new_spc_str
    
if __name__ == '__main__':
    test_parse_fct_grp()
    test_lumped_species()
    test__csv_io()
    test__mod_spc_dct_atomcount()
    test__mod_spc_dct_hof_basis()
    test__mod_spc_dct_stereo()
