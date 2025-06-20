""" test mechanalyzer.parser.pes
"""

import os

from ioformat import pathtools
from chemkin_io import parser as cki_parser
from mechanalyzer.parser.mech import parse_classtype
from mechanalyzer.builder.creckclass import build_creckclass_fromdct

CWD = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(CWD, 'data', 'osclass')


def test__build_creckclass():
    """generate creck class dataframe
    """
    ###### parsing operations
    # read external inputs
    class_str = pathtools.read_file(DAT_PATH, 'rxn_class_type.txt', remove_comments = '#')   
    ckin_str = pathtools.read_file(DAT_PATH, 'pah_block.dat')
    # parse classtype
    classtype_dct = parse_classtype(class_str)
    # parse reactions and assign class
    rxn_block_comments = cki_parser.mechanism.reaction_block(ckin_str, remove_comments=False)
    rxn_creckclass_dct = cki_parser.reaction.get_rxn_osclass_dct(rxn_block_comments)
    ####### generate creckclass dictionary
    creckclass_dct = build_creckclass_fromdct(rxn_creckclass_dct)
    print(creckclass_dct)
    ######## again but with classtype assigned
    creckclass_dct = build_creckclass_fromdct(rxn_creckclass_dct, classtype_dct=classtype_dct)
    
if __name__ == '__main__':
    test__build_creckclass()

    
    
