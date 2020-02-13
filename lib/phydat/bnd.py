""" Dictionary of standard bond lengths
"""

from lib.phydat.phycon import ANG2BOHR

LEN_DCT = {
    ('C', 'C'): 1.54 * ANG2BOHR,
    ('C', 'H'): 1.09 * ANG2BOHR,
    ('H', 'H'): 0.74 * ANG2BOHR,
    ('N', 'N'): 1.45 * ANG2BOHR,
    ('O', 'O'): 1.48 * ANG2BOHR,
    ('C', 'N'): 1.47 * ANG2BOHR,
    ('C', 'O'): 1.43 * ANG2BOHR,
    ('H', 'O'): 0.95 * ANG2BOHR,
}
