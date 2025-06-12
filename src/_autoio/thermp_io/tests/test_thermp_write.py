""" Test writing files for ThermP
"""

import os
from ioformat import pathtools
import thermp_io


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')


TEMPS = (100.0, 200.0, 300.0, 400.0, 500.0,
         600.0, 700.0, 800.0, 900.0, 1000.0)
NTEMPS = len(TEMPS)
FORMULA = 'CH4'
H0FORM = -0.02535174495
ENTHALPYT = 0.
BREAKT = 1000.


def test__input_file():
    """ test thermp_io.writer.input_file
    """

    inp_str = thermp_io.writer.input_file(
        NTEMPS,
        formula=FORMULA,
        delta_h=H0FORM,
        enthalpy_temp=ENTHALPYT,
        break_temp=BREAKT)
    assert inp_str == pathtools.read_file(DAT_PATH, 'thermp.inp')
