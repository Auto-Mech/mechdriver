"""
  Run Programs for thermo param calculations
"""

import os
import automol
from routines.pf import thermo
import autofile

# New Libs
from lib.runner import script
from lib.filesystem import orb as fsorb


# MESSPF
def run_pf(pf_path, pf_script_str=script.MESSPF):
    """ run messpf
    """
    script.run_script(pf_script_str, pf_path)


# THERMP
def write_thermp_inp(spc_dct_i):
    """ write the thermp input file
    """
    ich = spc_dct_i['ich']
    h0form = spc_dct_i['Hfs'][0]
    formula = automol.inchi.formula(ich)

    # Write thermp input file
    enthalpyt = 0.
    breakt = 1000.
    thermo.runner.write_thermp_input(
        formula=formula,
        delta_h=h0form,
        enthalpy_temp=enthalpyt,
        break_temp=breakt,
        thermp_file_name='thermp.dat')


def run_thermp(pf_path, nasa_path):
    """ run thermp to convert partition functions to thermochemical data
    """
    # Run thermp
    thermo.runner.run_thermp(
        pf_path=pf_path,
        thermp_path=nasa_path,
        thermp_file_name='thermp.dat',
        pf_file_name='pf.dat'
        )
    with open('thermp.out', 'r') as thermfile:
        lines = thermfile.readlines()
    line = lines[-1]
    hf298k = line.split()[-1]
    return hf298k


# PAC99
def run_pac(spc_dct_i, nasa_path):
    """ run pac99 to convert thermochemical data to nasa polynomials
    """
    ich = spc_dct_i['ich']
    formula = automol.inchi.formula(ich)

    # Run pac99
    thermo.runner.run_pac99(nasa_path, formula)

    # Read the pac99 polynomial
    with open(os.path.join(nasa_path, formula+'.c97'), 'r') as pac99_file:
        pac99_str = pac99_file.read()
    pac99_poly_str = thermo.nasapoly.get_pac99_polynomial(pac99_str)

    return pac99_poly_str


# PATH CONTROL
def go_to_path(path):
    """ change directory to path and return the original working directory
    """
    starting_path = os.getcwd()
    os.chdir(path)
    return starting_path


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
    orb_restr = fsorb.orbital_restriction(
        spc_info, har_level)
    har_levelp = har_level[1:3]
    har_levelp.append(orb_restr)

    thy_save_fs = autofile.fs.theory(spc_save_path)
    thy_save_fs[-1].create(har_levelp)
    thy_save_path = thy_save_fs[-1].path(har_levelp)
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
