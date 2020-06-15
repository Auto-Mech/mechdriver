"""
Library of physical constants
"""

from qcelemental import constants as qcc


# Unit Conversion Factors
ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')
BOHR2ANG = qcc.conversion_factor('bohr', 'angstrom')
DEG2RAD = qcc.conversion_factor('degree', 'radian')
RAD2DEG = qcc.conversion_factor('radian', 'degree')
WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
KCAL2WAVEN = qcc.conversion_factor('kcal/mol', 'wavenumber')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')
KCAL2EH = qcc.conversion_factor('kcal/mol', 'hartree')
WAVEN2EH = qcc.conversion_factor('wavenumber', 'hartree')
EH2WAVEN = qcc.conversion_factor('hartree', 'wavenumber')

# Physical Constants
NAVO = 6.02214076e23
