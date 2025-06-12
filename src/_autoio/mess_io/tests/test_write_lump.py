""" test mess_io.writer._lump
"""

import os
from ioformat import pathtools
import mess_io.writer


PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH = os.path.join(PATH, 'data', 'inp')

WELL_MERGE_LST = (
    ('W2', 'W1'),
    ('W5', 'W4', 'W3'),
    ('W9', 'W8', 'W7', 'W6')
)


def test__():
    """ test mess_io.writer._lump.well_lump_scheme
    """

    ref_lump1_str = pathtools.read_file(INP_PATH, 'well_lump1.inp')
    ref_lump2_str = pathtools.read_file(INP_PATH, 'well_lump2.inp')

    lump1_str = mess_io.writer.well_lump_scheme(WELL_MERGE_LST)
    lump2_str = mess_io.writer.well_lump_scheme(
        WELL_MERGE_LST, separator='-')

    assert lump1_str == ref_lump1_str.rstrip()
    assert lump2_str == ref_lump2_str.rstrip()

if __name__ == '__main__':
    test__()