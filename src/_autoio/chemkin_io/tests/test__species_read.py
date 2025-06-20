""" tests chemkin_io.writer.mechanism.species_block
"""

from chemkin_io.parser.species import names as parser


SPC_NAMES_STR = 'OH \nHO2 \nC3H8 \nN2O'
SPC_NAMES_TUPLE = ('OH', 'HO2', 'C3H8', 'N2O')


def test__read_spc_names():
    """ Tests the parsing of species names
    """
    spc_tuple = parser(SPC_NAMES_STR)
    assert spc_tuple == SPC_NAMES_TUPLE

if __name__ == '__main__':
    test__read_spc_names()
