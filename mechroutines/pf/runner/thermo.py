"""
  Run Programs for thermo param calculations
"""

import os
import sys
import subprocess
import shutil
import automol
import autofile
from mechanalyzer.inf import spc as sinfo
from mechlib.amech_io import printer as ioprinter


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
CUR_PATH = os.path.dirname(os.path.realpath(__file__))


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
    try:
        shutil.copyfile(pf_outfile, pf_datfile)
    except shutil.SameFileError:
        pass

    # Check for the existance of ThermP input and PF output
    assert os.path.exists(thermp_file), 'ThermP file does not exist'
    assert os.path.exists(pf_outfile), 'PF file does not exist'

    # Run thermp
    subprocess.check_call(['thermp', thermp_file])


def run_pac(formula, nasa_path):
    """
    Run pac99 for a given species name (formula)
    https://www.grc.nasa.gov/WWW/CEAWeb/readme_pac99.htm
    requires formula+'i97' and new.groups files
    """

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
        ioprinter.error_message(
            'PAC99 fit failed, maybe increase temperature ranges?')
        sys.exit()
    else:
        # Read the pac99 polynomial
        with open(os.path.join(nasa_path, formula+'.c97'), 'r') as pac99_file:
            pac99_str = pac99_file.read()
        if not pac99_str:
            ioprinter.error_message(
                'No polynomial produced from PAC99 fits, check for errors')
            sys.exit()


def thermo_paths(spc_dct_i, run_prefix, idx):
    """ Set up the path for saving the pf input and output.
        Placed in a MESSPF, NASA dirs high in run filesys.
    """

    # Get the formula and inchi key
    spc_info = sinfo.from_dct(spc_dct_i)
    spc_formula = automol.inchi.formula_string(spc_info[0])
    ich_key = automol.inchi.inchi_key(spc_info[0])

    # PF
    bld_locs = ['PF', idx]
    bld_save_fs = autofile.fs.build(run_prefix)
    bld_save_fs[-1].create(bld_locs)
    bld_path = bld_save_fs[-1].path(bld_locs)
    ioprinter.debug_message(
        'preparing thermo for:', spc_info[0],
        bld_path, spc_formula, ich_key)
    spc_pf_path = os.path.join(bld_path, spc_formula, ich_key)

    # NASA polynomials
    bld_locs = ['NASA', idx]
    bld_save_fs = autofile.fs.build(run_prefix)
    bld_save_fs[-1].create(bld_locs)
    bld_path = bld_save_fs[-1].path(bld_locs)
    spc_nasa_path = os.path.join(bld_path, spc_formula, ich_key)

    return spc_pf_path, spc_nasa_path
