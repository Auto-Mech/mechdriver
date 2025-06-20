""" Tests the writer functions in mechanism.py pertaining to species
"""

import difflib
from chemkin_io.writer import mechanism

# Define mechanism objects
MECH_SPC_DCT = {
    'O': {'smiles': 'smiles_1',
          'inchi': 'inchi_1',
          'charge': '',
          'mult': '',
          'fml': {'O': 1},
          'sens': ''},
    'H': {'smiles': 'smiles_2',
          'inchi': 'inchi_2',
          'charge': '',
          'mult': '',
          'fml': {'H': 1},
          'sens': ''}
}

SPC_STR = (
    'SPECIES\n\n'
    'O     ! SMILES: smiles_1         InChI: inchi_1  \n'  # note: dumb spaces
    'H     ! SMILES: smiles_2         InChI: inchi_2  \n\n'
    'END\n\n\n'
)

ELEM_STR = (
    'ELEMENTS\n\n'
    'H\n'
    'O\n\n'
    'END\n\n\n'
)


def test_spc():
    """ Tests the species block writer
    """
    out_str = mechanism.species_block(MECH_SPC_DCT)
    assert out_str == SPC_STR


def test_elem():
    """ Tests the elements block writer
    """
    out_str = mechanism.elements_block(MECH_SPC_DCT)
    assert out_str == ELEM_STR


def test_both():
    """ Tests both together
    """

    out_str = mechanism.write_chemkin_file(mech_spc_dct=MECH_SPC_DCT)
    assert out_str == (ELEM_STR + SPC_STR)


if __name__ == '__main__':
    test_spc()
    test_elem()
    test_both()
