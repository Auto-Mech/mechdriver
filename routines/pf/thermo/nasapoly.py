"""
 Generates NASA Polynomial from MESS+THERMP+PAC99 outputs
"""

from . import util


def get_pac99_polynomial(output_string):
    """ read in the polyn
    """

    lines = output_string.splitlines()
    pac_polynomial = '\n'.join([lines[i] for i in range(11)])

    return pac_polynomial


def convert_pac_to_chemkin(name, atom_dict, comment_str, pac_poly_str):
    """ convert the polynimal from pac format to chemkin polynomial
    """

    # Parse the lines of the pac string containing the desired coefficients
    las = [0.0 for i in range(7)]
    has = [0.0 for i in range(7)]
    las[0:5] = util.parse_line16(pac_poly_str.splitlines()[6][0:80])
    las[5:7] = util.parse_line16(pac_poly_str.splitlines()[7][48:80])
    has[0:5] = util.parse_line16(pac_poly_str.splitlines()[9][0:80])
    has[5:7] = util.parse_line16(pac_poly_str.splitlines()[10][48:80])

    if 'H' in atom_dict:
        num_h = atom_dict['H']
    else:
        num_h = 0
    if 'C' in atom_dict:
        num_c = atom_dict['C']
    else:
        num_c = 0
    if 'N' in atom_dict:
        num_n = atom_dict['N']
    else:
        num_n = 0
    if 'O' in atom_dict:
        num_o = atom_dict['O']
    else:
        num_o = 0
    # Build a string for the NASA polynomial in ChemKin format
    line1 = "%s        H%3d C%3d O%3d N%3d G%9.1F%10.1F%9.1F      1\n" % (
        name.ljust(16)[0:16],
        num_h, num_c, num_o, num_n,
        200.0, 3000.0, 1000.0)
    line2 = "% 15.8E% 15.8E% 15.8E% 15.8E% 15.8E    2\n" % (
        has[0], has[1], has[2], has[3], has[4])
    line3 = "% 15.8E% 15.8E% 15.8E% 15.8E% 15.8E    3\n" % (
        has[5], has[6], las[0], las[1], las[2])
    line4 = "% 15.8E% 15.8E% 15.8E% 15.8E                   4\n" % (
        las[3], las[4], las[5], las[6])
    poly_str = comment_str + line1 + line2 + line3 + line4

    return poly_str
