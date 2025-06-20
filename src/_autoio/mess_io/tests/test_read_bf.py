""" test mess_io.reader.ped
"""

import os
import numpy as np
from ioformat import pathtools
import mess_io


PATH = os.path.dirname(os.path.realpath(__file__))
OPATH = os.path.join(PATH, 'data', 'out')
OFILE = pathtools.read_file(OPATH, 'bf.out')

def test_get_bf():
    bf_dct = mess_io.reader.bf.get_bf(OFILE)
    assert np.isclose(bf_dct[1.0][0][4], 900)
    assert np.isclose(bf_dct[1.0][1][4], 0.463)


if __name__ == '__main__':
    test_get_bf()
