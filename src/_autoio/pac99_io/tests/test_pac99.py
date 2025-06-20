""" Test pac99_io.reader
"""

import os
from ioformat import pathtools
import pac99_io


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')


def test__polynomial():
    """ test pac99_io.reader.nasa_polynomial
        test pac99_io._convert.pac3ckin_poly
    """

    # Read and Convert C2H4N2O6
    out_str = pathtools.read_file(DAT_PATH, 'C2H4N2O6.c97')
    pac99_poly = pac99_io.reader.nasa_polynomial(out_str)
    assert pac99_poly == (
        pathtools.read_file(DAT_PATH, 'C2H4N2O6.paccpoly').rstrip())

    name = 'C2H4N2O6'
    atom_dct = {'C': 2, 'H': 4, 'N': 2, 'O': 6}
    ckin_poly = pac99_io.pac2ckin_poly(name, atom_dct, pac99_poly)
    assert ckin_poly == pathtools.read_file(DAT_PATH, 'C2H4N2O6.ckinpoly')

    # Read and Convert CO
    out_str = pathtools.read_file(DAT_PATH, 'CO.c97')
    pac99_poly = pac99_io.reader.nasa_polynomial(out_str)
    assert pac99_poly == (
        pathtools.read_file(DAT_PATH, 'CO.paccpoly')).rstrip()

    name = 'CO'
    atom_dct = {'C': 1, 'O': 1}
    ckin_poly = pac99_io.pac2ckin_poly(name, atom_dct, pac99_poly)
    assert ckin_poly == pathtools.read_file(DAT_PATH, 'CO.ckinpoly')

    # Read and Convert H2
    out_str = pathtools.read_file(DAT_PATH, 'H2.c97')
    pac99_poly = pac99_io.reader.nasa_polynomial(out_str)
    assert pac99_poly == (
        pathtools.read_file(DAT_PATH, 'H2.paccpoly')).rstrip()

    name = 'H2'
    atom_dct = {'H': 2}
    ckin_poly = pac99_io.pac2ckin_poly(name, atom_dct, pac99_poly)
    print(ckin_poly)
    assert ckin_poly == pathtools.read_file(DAT_PATH, 'H2.ckinpoly')
