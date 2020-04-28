"""
  Run Programs for thermo param calculations
"""

import os
import sys
import subprocess
import shutil
import numpy
import automol
import autofile
import mess_io
import thermp_io
from lib import filesys
from lib.submission import run_script
from lib.submission import DEFAULT_SCRIPT_DCT


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
CUR_PATH = os.path.dirname(os.path.realpath(__file__))


# MESSPF
def run_pf(pf_path, pf_script_str=DEFAULT_SCRIPT_DCT['messpf']):
    """ run messpf
    """
    run_script(pf_script_str, pf_path)


def read_messpf_temps(pf_path):
    """ Obtain the temperatures from the MESSPF file
    """

    # Read MESSPF file
    messpf_file = os.path.join(pf_path, 'pf.dat')
    with open(messpf_file, 'r') as pffile:
        output_string = pffile.read()

    # Obtain the temperatures, remove the 298.2 value
    temps, _, _, _ = mess_io.reader.pfs.partition_fxn(output_string)
    temps = [temp for temp in temps if not numpy.isclose(temp, 298.2)]

    return temps


# THERMP
def write_thermp_inp(spc_dct_i, temps, thermp_file_name='thermp.dat'):
    """ write the thermp input file
    """
    ich = spc_dct_i['ich']
    h0form = spc_dct_i['Hfs'][0]
    formula = automol.inchi.formula(ich)

    # Write thermp input file
    enthalpyt = 0.
    breakt = 1000.
    thermp_str = thermp_io.writer.thermp_input(
        ntemps=len(temps),
        formula=formula,
        delta_h=h0form,
        enthalpy_temp=enthalpyt,
        break_temp=breakt,
        thermp_file_name=thermp_file_name)

    # Write the file
    with open(thermp_file_name, 'w') as thermp_file:
        thermp_file.write(thermp_str)


def run_thermp(pf_path, thermp_path,
               thermp_file_name='thermp.dat', pf_file_name='pf.dat'):
    """
    Runs thermp.exe
    Requires thermp input file to be present
    partition function (pf) output file and
    """

    # Set full paths to files
    thermp_file = os.path.join(thermp_path, thermp_file_name)
    pf_outfile = os.path.join(pf_path, pf_file_name)

    # Copy MESSPF output file to THERMP run dir and rename to pf.dat
    pf_datfile = os.path.join(thermp_path, 'pf.dat')
    shutil.copyfile(pf_outfile, pf_datfile)

    # Check for the existance of ThermP input and PF output
    assert os.path.exists(thermp_file)
    assert os.path.exists(pf_outfile)

    # Run thermp
    subprocess.check_call(['thermp', thermp_file])

    # Read ene from file
    with open('thermp.out', 'r') as thermfile:
        lines = thermfile.readlines()
    line = lines[-1]
    hf298k = line.split()[-1]
    return hf298k


# PAC99
def run_pac(spc_dct_i, nasa_path):
    """
    Run pac99 for a given species name (formula)
    https://www.grc.nasa.gov/WWW/CEAWeb/readme_pac99.htm
    requires formula+'i97' and new.groups files
    """
    ich = spc_dct_i['ich']
    formula = automol.inchi.formula(ich)

    # Run pac99
    # Set file names for pac99
    i97_file = os.path.join(nasa_path, formula + '.i97')
    newgroups_file = os.path.join(nasa_path, 'new.groups')
    newgroups_ref = os.path.join(CUR_PATH, 'new.groups')

    # Copy new.groups file from thermo src dir to run dir
    shutil.copyfile(newgroups_ref, newgroups_file)

    # Check for the existance of pac99 files
    assert os.path.exists(i97_file)
    assert os.path.exists(newgroups_file)

    # Run pac99
    proc = subprocess.Popen('pac99', stdin=subprocess.PIPE)
    proc.communicate(bytes(formula, 'utf-8'))

    # Check to see if pac99 does not have error message
    with open(os.path.join(nasa_path, formula+'.o97'), 'r') as pac99_file:
        pac99_out_str = pac99_file.read()
    if 'INSUFFICIENT DATA' in pac99_out_str:
        print('*ERROR: PAC99 fit failed, maybe increase temperature ranges?')
        sys.exit()
    else:
        # Read the pac99 polynomial
        with open(os.path.join(nasa_path, formula+'.c97'), 'r') as pac99_file:
            pac99_str = pac99_file.read()
        if not pac99_str:
            print('No polynomial produced from PAC99 fits, check for errors')
            sys.exit()

    return pac99_str


# PATH CONTROL
def get_starting_path():
    """ get original working directory
    """
    starting_path = os.getcwd()
    return starting_path


def go_to_path(path):
    """ change directory to path and return the original working directory
    """
    os.chdir(path)


def return_to_path(path):
    """ change directory to starting path
    """
    os.chdir(path)


def prepare_path(path, loc):
    """ change directory to starting path, return chemkin path
    """
    new_path = os.path.join(path, loc)
    return new_path


def get_thermo_paths(spc_save_path, spc_info, har_level):
    """ set up the path for saving the pf input and output
    currently using the harmonic theory directory for this because
    there is no obvious place to save this information for a random
    assortment of har_level, tors_level, vpt2_level
    """
    har_levelp = filesys.inf.mod_orb_restrict(
        spc_info, har_level)

    thy_save_fs = autofile.fs.theory(spc_save_path)
    thy_save_fs[-1].create(har_levelp[1:4])
    thy_save_path = thy_save_fs[-1].path(har_levelp[1:4])
    bld_locs = ['PF', 0]
    bld_save_fs = autofile.fs.build(thy_save_path)
    bld_save_fs[-1].create(bld_locs)
    pf_path = bld_save_fs[-1].path(bld_locs)

    # prepare NASA polynomials
    bld_locs = ['NASA_POLY', 0]
    bld_save_fs[-1].create(bld_locs)
    nasa_path = bld_save_fs[-1].path(bld_locs)

    print('Build Path for Partition Functions')
    print(pf_path)

    return pf_path, nasa_path
