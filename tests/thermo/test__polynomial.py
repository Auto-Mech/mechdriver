"""
Tests for interacting with the polynomial
"""

import os
from thermo import nasapoly


# Input data
FORMULA = 'CH4'
PAC99_OUTFILE_NAME = os.path.join(os.getcwd(), 'run', FORMULA+'.c97')


def test__polynomial():
    """ reads the polynomial from pac99 and converts it
    """

    # Open the output from pac99
    with open(PAC99_OUTFILE_NAME, 'r') as pac99_file:
        pac99_str = pac99_file.read()

    # Get the pac99 polynomial
    pac99_poly_str = nasapoly.get_pac99_polynomial(pac99_str)
    print('\nPAC99 Polynomial:')
    print(pac99_poly_str)

    # Convert the pac99 polynomial to chemkin polynomial
    chemkin_poly_str = nasapoly.convert_pac_to_chemkin(pac99_poly_str)
    print('\nCHEMKIN Polynomial:')
    print(chemkin_poly_str)


if __name__ == '__main__':
    test__polynomial()
