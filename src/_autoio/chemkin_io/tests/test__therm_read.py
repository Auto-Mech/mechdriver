""" test chemkin_io.writer.mechanism.thermo_block
"""

import numpy
from chemkin_io.parser.thermo import create_spc_nasa7_dct as parser

THERM_STR = ( 
    'THERMO\n'
    '200.00    1000.00   5000.000\n\n'
    'O2                RUS 89O   2               G    200.00   6000.00 1000.00      1\n'
    ' 2.54363697E+00-2.73162486E-05-4.19029520E-09 4.95481845E-12-4.79553694E-16    2\n'
    ' 2.92260120E+04 4.92229457E+00 3.16826710E+00-3.27931884E-03 6.64306396E-06    3\n'
    '-6.12806624E-09 2.11265971E-12 2.91222592E+04 2.05193346E+00                   4\n'
    'END\n\n\n'
)


def test_read():
    """ Tests the chemkin_io parsing for thermo
    """

    spc_nasa7_dct = parser(THERM_STR)
    nasa7_params = spc_nasa7_dct['O2']
    assert nasa7_params[0] == 'RUS 89'
    assert nasa7_params[1] == 'O   2               '
    assert nasa7_params[2] == 'G'
    assert numpy.allclose(nasa7_params[3], [200.0, 6000.0, 1000.0])
    hight = nasa7_params[4][0]
    lowt = nasa7_params[4][1]
    ref_hight = [2.54363697, -2.73162486e-05, -4.1902952e-09, 4.95481845e-12,
                 -4.79553694e-16, 29226.012, 4.92229457] 
    ref_lowt = [3.1682671, -0.00327931884, 6.64306396e-06, -6.12806624e-09,
                2.11265971e-12, 29122.2592, 2.05193346]
    assert numpy.allclose(hight, ref_hight)
    assert numpy.allclose(lowt, ref_lowt)


if __name__ == '__main__':
    test_read()
