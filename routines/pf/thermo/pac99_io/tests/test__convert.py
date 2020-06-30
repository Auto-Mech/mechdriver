""" Test pac99_io._convert.pac2ckin_poly
"""

from pac99_io import pac2ckin_poly


PAC99_POLY = """

"""


def test__convert():
    """ test pac99_io._convert.pac2ckin_poly
    """

    ckin_poly = pac2ckin_poly(name, atom_dict, comment_str, pac_poly_str):


    ref_ckin_poly = """

    """

    assert ckin_poly == ref_ckin_poly


if __name__ == '__main__':
    test__convert()
