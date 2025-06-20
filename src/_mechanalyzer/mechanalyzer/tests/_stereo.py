""" test the stereochem
"""

import os
import chemkin_io
import mechanalyzer


# Get mechanism information
def _read_file(file_name):
    with open(file_name, encoding='utf8', errors='ignore') as file_obj:
        file_str = file_obj.read()
    return file_str


PATH = os.path.dirname(os.path.realpath(__file__))
DATA_PATH = os.path.join(PATH, 'data')
MECH_NAME = 'ste.txt'
CSV_NAME = 'ste2.csv'

MECH_STR = _read_file(
    os.path.join(DATA_PATH, MECH_NAME))
CSV_STR = _read_file(
    os.path.join(DATA_PATH, CSV_NAME))

# Get the species
SPC_DCT = mechanalyzer.parser.spc.build_spc_dct(CSV_STR, 'csv')

# Get the reactions
RXN_BLOCK = chemkin_io.parser.mechanism.reaction_block(
    MECH_STR)
RXN_DCT = chemkin_io.parser.reaction.param_dct(
    RXN_BLOCK, 'cal/mole', 'moles')


def test__():
    """ test
    """

    ste_rxn_dct, ste_spc_dct = mechanalyzer.builder.expand_mech_stereo(
        RXN_DCT, SPC_DCT)
    print('\nrxn dct')
    for rxn, pars in ste_rxn_dct.items():
        print(rxn)
        print(pars)
    print('\nspc dct')
    for spc in ste_spc_dct:
        print(spc, ste_spc_dct[spc]['inchi'])


if __name__ == '__main__':
    test__()
