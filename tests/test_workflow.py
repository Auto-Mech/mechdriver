""" Runs automech instancs for tests
"""

import os
# import tempfile
import numpy
from _util import run_mechdriver
from _util import chk_therm
from _util import chk_rates


# Set path where test input files and output data comparison exist
PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')

# Paths to the current input directory
CWD_INP_DIR = os.path.join(PATH, 'inp')

# Paths to the temp directory where tests will be run
# TMP_DIR = tempfile.mkdtemp()
TMP_DIR = os.path.join(os.getcwd(), 'tmp')
TMP_INP_DIR = os.path.join(TMP_DIR, 'inp')
TMP_RUN_DIR = os.path.join(TMP_DIR, 'run')
TMP_SAVE_DIR = os.path.join(TMP_DIR, 'save')
print(TMP_DIR)

# Set conditions
TEMPS = numpy.array([500.0, 1000.0, 1500.0, 2000.0])
PRESSURES = (1.0, 'high')


def test__rrho():
    """ Run es, thermo, and rates for PES; standard run
    """
    print('RUN test__rrho')
    run_mechdriver('run_c2h6_h_rrho.temp',
                   TMP_DIR,
                   TMP_INP_DIR, CWD_INP_DIR,
                   TMP_RUN_DIR, TMP_SAVE_DIR)
    # chk_therm('')
    # chk_rates('')


def test__1dhrfa():
    """ Run es, thermo, and rates for PES; standard run
    """

    print('RUN test__1dhrfa')
    run_mechdriver('run_c2h6_h_1dhrfa.temp',
                   TMP_DIR,
                   TMP_INP_DIR, CWD_INP_DIR,
                   TMP_RUN_DIR, TMP_SAVE_DIR)
    # chk_therm('c2h6_h_1dhrfa_therm.ckin', 'all_therm.ckin',
    #           DAT_PATH, TMP_DIR,
    #           TMP_INP_DIR,
    #           TEMPS)
    # chk_rates('c2h6_h_1dhrfa_rate.ckin', 'C2H7.ckin',
    #           DAT_PATH,
    #           TMP_DIR,  # TMP_INP_DIR,
    #           PRESSURES, TEMPS)


def __etoh():
    """ Run es, thermo, for EtOH with different rotor types
        need a species that uses theory methods scaling
    """
    run_mechdriver('run_c2h5oh_full.temp',
                   TMP_DIR,
                   TMP_INP_DIR, CWD_INP_DIR,
                   TMP_RUN_DIR, TMP_SAVE_DIR)
    # chk_therm('')
    # chk_rates('')


def __instab():
    """ Run es, thermo, and rates for PES with instabilities
    """
    run_mechdriver('run_ch2ooh_rrho.temp',
                   TMP_DIR,
                   TMP_INP_DIR, CWD_INP_DIR,
                   TMP_RUN_DIR, TMP_SAVE_DIR)


def test__radrad():
    """ Run es, thermo, and rates for PES with instabilities
    """

    print('RUN test__radrad')
    run_mechdriver('run_c2h5_h_1dhrfa.temp',
                   TMP_DIR,
                   TMP_INP_DIR, CWD_INP_DIR,
                   TMP_RUN_DIR, TMP_SAVE_DIR)
    # later: switch c2h5 to vtst/vrctst and ch3_h to pst


if __name__ == '__main__':
    test__rrho()
    test__1dhrfa()
    test__radrad()
