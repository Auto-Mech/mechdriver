"""
Tests running thermochemistry codes
"""

import os
from thermo import runner


# Input data
FORMULA = 'CH4'
DELTAH = -15.91
ENTHALPYT = 0.0
BREAKT = 1000.0
THERMP_FILE_NAME = 'thermp.dat'
PF_FILE_NAME = 'pf.out'
DATA_PATH = os.path.join(os.getcwd(), 'data')
RUN_PATH = os.path.join(os.getcwd(), 'run')


def test__run():
    """ run thermp and pac99
    """

    # Go into run directory (needed?)
    os.chdir(RUN_PATH)

    # Write thermp input file
    runner.write_thermp_input(FORMULA, DELTAH,
                              enthalpyT=ENTHALPYT, breakT=BREAKT,
                              thermp_file_name=THERMP_FILE_NAME)

    # Run thermp
    runner.run_thermp(DATA_PATH, RUN_PATH,
                      thermp_file_name=THERMP_FILE_NAME,
                      pf_file_name=PF_FILE_NAME)

    # Run pac99
    runner.run_pac99(RUN_PATH, FORMULA)


if __name__ == '__main__':
    test__run()
