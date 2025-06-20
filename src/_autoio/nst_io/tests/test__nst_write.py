""" test intder_io.writer
"""

import os
from ioformat import pathtools
import nst_io.writer


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')

JOBTYPE = 'OPTGM'
ZERO_ENE = -184.17192618
GEO = (('N',  (0.4861852503, -0.5279763954, 0.0000000000)),
       ('N',  (2.4573988480, 0.3033407113, 0.0000000000)),
       ('O',  (-2.6043448708, 0.1832871600, 0.0000000000)))


def test__input():
    """ test nst_io.writer.input_file
    """

    inp_str = nst_io.writer.input_file(JOBTYPE, GEO, ZERO_ENE)
    assert inp_str == pathtools.read_file(DAT_PATH, 'input.dat')


if __name__ == '__main__':
    test__input()
