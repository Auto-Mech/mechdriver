""" Tests parsing of mechanism blocks
"""

from chemkin_io.parser import mechanism

SPECIES_BLOCK_GOOD = 'SPECIES\nCH4 H2 O2 N2\nEND'
SPECIES_BLOCK_BAD = 'SPECIES\nCH4 H2 O2 N2\n'
REACTION_BLOCK_GOOD = 'REACTIONS\nCH4+H=CH3+H2\nEND'
REACTION_BLOCK_BAD = 'REACTIONS\nCH4+H=CH3+H2\n'
THERMO_BLOCK_GOOD = 'THERMO\nC7H15OOH-1\nEND'
THERMO_BLOCK_BAD = 'THERMO\nC7H15OOH-1\n'
ELEMENT_BLOCK_GOOD = 'ELEMENTS\nC H O N\nEND'
ELEMENT_BLOCK_BAD = 'ELEMENTS\nC H O N\n'

REACTION_BLOCK_COMMENTS = "REACTIONS\n\
C4H7ORvE4fmAB0 = C4H7O4H74fm1                                             1.000E+00     0.000        0  ! pes.subpes.channel  1.1.1\n\n \
C4H7ORvE4fmAB0 = C4H7O-kSV4fm                                             1.000E+00     0.000        0  ! pes.subpes.channel  1.1.2\n\n \
C4H7ORvE4fmAB0 = C4H6O-RvErx51 + H-TcYTcY                                 1.000E+00     0.000        0  ! pes.subpes.channel  1.1.3\n\n \
C4H7ORvE4fmAA0 = C4H7O4H74fm0                                             1.000E+00     0.000        0  ! pes.subpes.channel  1.1.4\n\n \
C4H7ORvE4fmAA0 = C4H7O-kSV4fm                                             1.000E+00     0.000        0  ! pes.subpes.channel  1.1.5\n\n \
C4H7ORvE4fmAA0 = C4H6O-RvErx50 + H-TcYTcY                                 1.000E+00     0.000        0  ! pes.subpes.channel  1.1.6\n\n \
C4H7O4H74fm0 = C3H4OALAD-Wv9FbZ + CH3                                     1.000E+00     0.000        0  ! pes.subpes.channel  1.1.7\n\n \
C4H7O4H74fm0 = C2H4OALD-UPQWKw + C2H3ALK-S58hH1                           1.000E+00     0.000        0  ! pes.subpes.channel  1.1.8\n\n \
C4H7O-kSV4fm = C2H4OALD-UPQWKw + C2H3ALK-S58hH1                           1.000E+00     0.000        0  ! pes.subpes.channel  1.1.9\n\n \
C4H8ORvEsWvAA0 + OH = C4H7ORvE4fmAA0 + H2O                                1.000E+00     0.000        0  ! pes.subpes.channel  2.1.1\n\n \
C4H8ORvEsWvAB + OH = C4H7ORvE4fmAB0 + H2O                                 1.000E+00     0.000        0  ! pes.subpes.channel  2.2.2\n\n \
C4H8ORvEsWvAA0 + HO2-S580KW = C4H7ORvE4fmAA0 + H2O2-S58pAY                 1.000E+00     0.000        0  ! pes.subpes.channel  3.1.1\n\n \
C4H8ORvEsWvAA0 + CH3 = C4H7ORvE4fmAA0 + CH4                                1.000E+00     0.000        0  ! pes.subpes.channel  4.1.1\n\n \
C4H8ORvEsWvAA0 + CH3O2RO2-2LTcwB = C4H7ORvE4fmAA0 + CH4O2QOOH-2LTWKw       1.000E+00     0.000        0  ! pes.subpes.channel  5.1.1\n\n \
C4H8ORvEsWvAA0 + CH3O-S58cwB = C4H7ORvE4fmAA0 + CH4O-S58WKw                1.000E+00     0.000        0  ! pes.subpes.channel  6.1.1\n\n \
C4H8ORvEsWvAA0 + Cl = C4H7ORvE4fmAA0 + HCl                                 1.000E+00     0.000        0  ! pes.subpes.channel  7.1.1\n\n \
END"

REACTION_BLOCK_OKCOMM = "REACTIONS\n\
C4H7ORvE4fmAB0 = C4H7O4H74fm1                                             1.000E+00     0.000        0  # pes.subpes.channel  1.1.1\n \
C4H7ORvE4fmAB0 = C4H7O-kSV4fm                                             1.000E+00     0.000        0  # pes.subpes.channel  1.1.2\n \
C4H7ORvE4fmAB0 = C4H6O-RvErx51 + H-TcYTcY                                 1.000E+00     0.000        0  # pes.subpes.channel  1.1.3\n \
C4H7ORvE4fmAA0 = C4H7O4H74fm0                                             1.000E+00     0.000        0  # pes.subpes.channel  1.1.4\n \
C4H7ORvE4fmAA0 = C4H7O-kSV4fm                                             1.000E+00     0.000        0  # pes.subpes.channel  1.1.5\n \
C4H7ORvE4fmAA0 = C4H6O-RvErx50 + H-TcYTcY                                 1.000E+00     0.000        0  # pes.subpes.channel  1.1.6\n \
C4H7O4H74fm0 = C3H4OALAD-Wv9FbZ + CH3                                     1.000E+00     0.000        0  # pes.subpes.channel  1.1.7\n \
C4H7O4H74fm0 = C2H4OALD-UPQWKw + C2H3ALK-S58hH1                           1.000E+00     0.000        0  # pes.subpes.channel  1.1.8\n \
C4H7O-kSV4fm = C2H4OALD-UPQWKw + C2H3ALK-S58hH1                           1.000E+00     0.000        0  # pes.subpes.channel  1.1.9\n \
C4H8ORvEsWvAA0 + OH = C4H7ORvE4fmAA0 + H2O                                1.000E+00     0.000        0  # pes.subpes.channel  2.1.1\n \
C4H8ORvEsWvAB + OH = C4H7ORvE4fmAB0 + H2O                                 1.000E+00     0.000        0  # pes.subpes.channel  2.2.2\n \
C4H8ORvEsWvAA0 + HO2-S580KW = C4H7ORvE4fmAA0 + H2O2-S58pAY                 1.000E+00     0.000        0  # pes.subpes.channel  3.1.1\n \
C4H8ORvEsWvAA0 + CH3 = C4H7ORvE4fmAA0 + CH4                                1.000E+00     0.000        0  # pes.subpes.channel  4.1.1\n \
C4H8ORvEsWvAA0 + CH3O2RO2-2LTcwB = C4H7ORvE4fmAA0 + CH4O2QOOH-2LTWKw       1.000E+00     0.000        0  # pes.subpes.channel  5.1.1\n \
C4H8ORvEsWvAA0 + CH3O-S58cwB = C4H7ORvE4fmAA0 + CH4O-S58WKw                1.000E+00     0.000        0  # pes.subpes.channel  6.1.1\n \
C4H8ORvEsWvAA0 + Cl = C4H7ORvE4fmAA0 + HCl                                 1.000E+00     0.000        0  # pes.subpes.channel  7.1.1\n \
END"

def test_species_block():
    """ Tests the species_block function
    """
    species_block_good = mechanism.species_block(SPECIES_BLOCK_GOOD)
    species_block_bad = mechanism.species_block(SPECIES_BLOCK_BAD)
    assert species_block_good == '\nCH4 H2 O2 N2'
    assert species_block_bad is None


def test_reaction_block():
    """ Tests the reaction_block function
    """
    reaction_block_good = mechanism.reaction_block(REACTION_BLOCK_GOOD)
    reaction_block_bad = mechanism.reaction_block(REACTION_BLOCK_BAD)
    assert reaction_block_good == '\nCH4+H=CH3+H2'
    assert reaction_block_bad is None

def test_reaction_block_comments():
    """ Tests the reaction_block function with comments
    """
    reaction_block_comments = mechanism.reaction_block(REACTION_BLOCK_COMMENTS, remove_comments=False)
    # check that first line is correct
    assert reaction_block_comments.split('\n')[1] == 'C4H7ORvE4fmAB0 = C4H7O4H74fm1                                             1.000E+00     0.000        0  ! pes.subpes.channel  1.1.1'
    reaction_block_comments = mechanism.reaction_block(REACTION_BLOCK_OKCOMM, remove_comments=True) #will not remove hashtags
    assert reaction_block_comments.split('\n')[1] == 'C4H7ORvE4fmAB0 = C4H7O4H74fm1                                             1.000E+00     0.000        0  # pes.subpes.channel  1.1.1'
    
def test_thermo_block():
    """ Tests the thermo_block function
    """
    thermo_block_good = mechanism.thermo_block(THERMO_BLOCK_GOOD)
    thermo_block_bad = mechanism.thermo_block(THERMO_BLOCK_BAD)
    assert thermo_block_good == '\nC7H15OOH-1'
    assert thermo_block_bad is None


def test_element_block():
    """ Tests the element_block function
    """
    element_block_good = mechanism.element_block(ELEMENT_BLOCK_GOOD)
    element_block_bad = mechanism.element_block(ELEMENT_BLOCK_BAD)
    assert element_block_good == '\nC H O N'
    assert element_block_bad is None


if __name__ == '__main__':
    test_species_block()
    test_reaction_block()
    test_reaction_block_comments()
    test_thermo_block()
    test_element_block()
