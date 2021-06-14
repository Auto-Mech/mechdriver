""" Runs automech instancs for tests
"""

import os
import shutil
import tempfile
import subprocess


# Set path where test input file exist
PATH = os.path.dirname(os.path.realpath(__file__))
CWD_INP_DIR = os.path.join(PATH, 'inp')

# Set paths where tests will run
# TMP_DIR = tempfile.mkdtemp()
TMP_DIR = os.path.join(os.getcwd(), 'tmp')
TMP_INP_DIR = os.path.join(TMP_DIR, 'inp')
TMP_RUN_DIR = os.path.join(TMP_DIR, 'run')
TMP_SAVE_DIR = os.path.join(TMP_DIR, 'save')
print(TMP_DIR)

# Set command line
EXE_PATH = os.path.join(PATH, '../bin/automech.py')
CMD_LINE = 'python -u {0} {1} & disown %1'.format(EXE_PATH, TMP_DIR)


# Test functions
def test__rrho():
    """ Run es, thermo, and rates for PES; standard run
    """
    _run('run_c2h6_h_rrho.temp')


def test__1dhr():
    """ Run es, thermo, and rates for PES; standard run
    """
    _run('run_c2h6_h_1dhr.temp')


def __etoh():
    """ Run es, thermo, for EtOH with different rotor types

        need a species that uses theory methods scaling
    """
    _run('run_c2h5oh_full.temp')


def __instab():
    """ Run es, thermo, and rates for PES with instabilities
    """
    _run('run_ch2ooh_rrho.temp')


def test__radrad():
    """ Run es, thermo, and rates for PES with instabilities
    """
    _run('run_c2h5_h_1dhrfa.temp')
    # later: switch c2h5 to vtst/vrctst and ch3_h to pst


def test__therm_basic():
    """ Run minimal tasks for simple thermo calc
    """
    _run('run_therm_basic.temp')


def __trans():
    """ Run minimal tasks to generate ckin transport
    """
    _run('run_trans.temp')


def __proc():
    """ Run minimal tasks to generate and produce output
    """
    _run('run_proc.temp')


# Helper functions to run a single instance of MechDriver
def _run(run_template):
    """ test automech.py
    """
    # Copy input to tmp directory and replace it
    if not os.path.exists(TMP_INP_DIR):
        shutil.copytree(CWD_INP_DIR, TMP_INP_DIR)
    _fill_template_and_write_file(run_template, 'run.dat')

    logfile = open('{0}/run.log'.format(TMP_DIR), 'w')
    subprocess.call(CMD_LINE.split(), stdout=logfile, stderr=logfile)


def _fill_template_and_write_file(templatefile, inpfile):
    """ Read the run.dat and replace run_prefix and save_prefix
    """

    # Set names of template and input for the calculation
    inp_file = os.path.join(TMP_INP_DIR, templatefile)
    new_inp_file = os.path.join(TMP_INP_DIR, inpfile)

    # Read template and fill with the run and save prefix
    with open(inp_file, 'r') as fobj:
        inp_str = fobj.read()
    new_inp_str = inp_str.format(TMP_RUN_DIR, TMP_SAVE_DIR)

    # Write the run.dat and models.dat for the calculation
    with open(new_inp_file, 'w') as fobj:
        fobj.write(new_inp_str)


if __name__ == '__main__':
    test__rrho()
    # test__1dhr()
    # test__etoh()
    # test__instab()
    # test__radrad()
    # test__therm_basic()
    # test__proc()
    # test__proc()
