"""
Library of physical constants
"""

from qcelemental import constants as qcc


# Unit Conversion Factors
ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')
DEG2RAD = qcc.conversion_factor('degree', 'radian')
WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

# Physical Constants
NAVO = 6.02214076e23
