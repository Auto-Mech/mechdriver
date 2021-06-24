""" Runs automech instancs for tests
"""

import os
import shutil
import tempfile
import subprocess
import numpy
import chemkin_io
from ioformat import pathtools
import mechanalyzer
import ratefit


# Set path where test input files and output data comparison exist
PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')

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

# Set conditions
TEMPS = numpy.array([500.0, 1000.0, 1500.0, 2000.0])
PRESSURES = (1.0, 'high')


# Test functions
def __rrho():
    """ Run es, thermo, and rates for PES; standard run
    """
    _run('run_c2h6_h_rrho.temp')
    # _chk_therm('')
    # _chk_rates('')


def test__1dhrfa():
    """ Run es, thermo, and rates for PES; standard run
    """
    _run('run_c2h6_h_1dhrfa.temp')
    _chk_therm('c2h6_h_1dhrfa_therm.ckin', 'all_therm.ckin')
    _chk_rates('c2h6_h_1dhrfa_rate.ckin', 'C2H7.ckin')


def __etoh():
    """ Run es, thermo, for EtOH with different rotor types
        need a species that uses theory methods scaling
    """
    _run('run_c2h5oh_full.temp')
    # _chk_therm('')
    # _chk_rates('')


def __instab():
    """ Run es, thermo, and rates for PES with instabilities
    """
    _run('run_ch2ooh_rrho.temp')


def test__radrad():
    """ Run es, thermo, and rates for PES with instabilities
    """
    _run('run_c2h5_h_1dhrfa.temp')
    # later: switch c2h5 to vtst/vrctst and ch3_h to pst


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


# Helper functions to ensure consistent numerical results
def _chk_therm(therm_dat_file, therm_calc_file):
    """ Read Check the
    """

    # Read the data in the thermo and rate CKIN files
    ckin_path = os.path.join(TMP_DIR, 'CKIN')

    therm_calc_str = pathtools.read_file([ckin_path], therm_calc_file)
    therm_dat_str = pathtools.read_file(DAT_PATH, therm_dat_file)

    nasa7_calc = chemkin_io.parser.thermo.create_spc_nasa7_dct(therm_calc_str)
    nasa7_dat = chemkin_io.parser.thermo.create_spc_nasa7_dct(therm_dat_str)

    thm_calc = mechanalyzer.calculator.thermo.create_spc_thermo_dct(
        nasa7_calc, TEMPS)
    thm_dat = mechanalyzer.calculator.thermo.create_spc_thermo_dct(
        nasa7_dat, TEMPS)

    spc_str = pathtools.read_file([], 'species.csv', path=TMP_INP_DIR)
    spc_ident_dct = mechanalyzer.parser.spc.build_spc_dct(spc_str, 'csv')

    thm_dct = mechanalyzer.calculator.compare.get_aligned_spc_thermo_dct(
        [thm_dat, thm_calc], [spc_ident_dct, spc_ident_dct])

    for thm_data in thm_dct.values():
        # 0 is data/therm and 1 is calc'd therm
        assert _assess(thm_data[0][0], thm_data[1][0], thresh=3.0)
        assert _assess(thm_data[0][1], thm_data[1][1], thresh=3.0)
        assert _assess(thm_data[0][2], thm_data[1][2], thresh=3.0)
        assert _assess(thm_data[0][3], thm_data[1][3], thresh=3.0)


def _chk_rates(rates_dat_file, rates_calc_file):
    """ Read Check the
    """

    # Read the data in the thermo and rate CKIN files
    ckin_path = os.path.join(TMP_DIR, 'CKIN')

    rates_calc_str = pathtools.read_file([ckin_path], rates_calc_file)
    rates_dat_str = pathtools.read_file(DAT_PATH, rates_dat_file)

    rxn_str_calc = chemkin_io.parser.mechanism.reaction_block(rates_calc_str)
    rxn_str_dat = chemkin_io.parser.mechanism.reaction_block(rates_dat_str)
    units = chemkin_io.parser.mechanism.reaction_units(rates_dat_str)

    par_dct_calc = chemkin_io.parser.reaction.param_dct(rxn_str_calc, *units)
    par_dct_dat = chemkin_io.parser.reaction.param_dct(rxn_str_dat, *units)

    rxn_ktp_dct_calc = mechanalyzer.calculator.rates.eval_rxn_param_dct(
        par_dct_calc, PRESSURES, TEMPS)
    rxn_ktp_dct_dat = mechanalyzer.calculator.rates.eval_rxn_param_dct(
        par_dct_dat, PRESSURES, TEMPS)

    # spc_str = pathtools.read_file([], 'species.csv', path=TMP_INP_DIR)
    # spc_ident_dct = mechanalyzer.parser.spc.build_spc_dct(spc_str, 'csv')

    # ktp_dct = mechanalyzer.calculator.compare.get_aligned_rxn_ktp_dct(
    #     [ktp_dct_dat, ktp_dct_calc], [],
    #     [spc_ident_dct, spc_ident_dct],
    #     (500., 1000., 1500., 2000.))
    # print(ktp_dct)

    for rxn, ktp_dct_calc in rxn_ktp_dct_calc.items():
        for pressure in ktp_dct_calc:
            calc = ratefit.ktpdct.read(
                rxn_ktp_dct_calc, rxn, pressure, 'rates')
            dat = ratefit.ktpdct.read(
                rxn_ktp_dct_dat, rxn, pressure, 'rates')
            assert _assess(dat, calc, thresh=3.0)


def _assess(dat1, dat2, thresh):
    """ Assess if two sets of values are within some numerical threshold
    """

    cond = True
    for val1, val2 in zip(dat1, dat2):
        cond = bool((abs(val1 - val2) / val1) * 100.0 < thresh)

    return cond


if __name__ == '__main__':
    # test__rrho()
    test__1dhrfa()
    # test__etoh()
    # test__instab()
    # test__radrad()
    # test__therm_basic()
    # test__proc()
    # test__proc()
