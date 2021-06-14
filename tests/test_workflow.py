""" Runs automech instancs for tests
"""

import os
import shutil
import tempfile
import subprocess
import numpy
import mechanalyzer
import chemkin_io
from ioformat import read_text_file


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


def test__etoh():
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
    # _run('run_therm_basic.temp')
    _chk_therm('all.ckin')


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


# Checker functions for assessing the proper output
def _chk_therm(therm_dat_file):
    """ Read Check the
    """

    # Read the data in the thermo and rate CKIN files
    ckin_path = os.path.join(TMP_DIR, 'CKIN')

    therm_calc_str = read_text_file([ckin_path], 'all.ckin')
    therm_dat_str = read_text_file(['data'], therm_dat_file, path=PATH)

    nasa7_calc = chemkin_io.parser.thermo.create_spc_nasa7_dct(therm_calc_str)
    nasa7_dat = chemkin_io.parser.thermo.create_spc_nasa7_dct(therm_dat_str)

    thm_calc = mechanalyzer.calculator.thermo.create_spc_thermo_dct(
        nasa7_calc, (500., 1000., 1500., 2000.))
    thm_dat = mechanalyzer.calculator.thermo.create_spc_thermo_dct(
        nasa7_dat, (500., 1000., 1500., 2000.))

    spc_str = read_text_file([], 'species.csv', path=TMP_INP_DIR)
    spc_ident_dct = mechanalyzer.parser.spc.build_spc_dct(spc_str, 'csv')

    thm_dct = mechanalyzer.calculator.compare.get_aligned_spc_thermo_dct(
        [thm_dat, thm_calc], [spc_ident_dct, spc_ident_dct])

    for thm_data in thm_dct.values():
        # 0 is data/therm and 1 is calc'd therm
        assert numpy.allclose(thm_data[0][0], thm_data[1][0])
        assert numpy.allclose(thm_data[0][1], thm_data[1][1])
        assert numpy.allclose(thm_data[0][2], thm_data[1][2])
        assert numpy.allclose(thm_data[0][3], thm_data[1][3])


if __name__ == '__main__':
    test__rrho()
    # test__1dhr()
    # test__etoh()
    # test__instab()
    # test__radrad()
    # test__therm_basic()
    # test__proc()
    # test__proc()
