""" Generate NASA polynomials for a set of species and write them
    into a Chemkin THERMO block

    Input required:
        (1) species.csv file with contents
            name,smiles,hof
            <spc1_name>,<spc1_smiles>,<spc1_heat-of-formation>
            <spc2_name>,<spc2_smiles>,<spc2_heat-of-formation>
        where the heat-of-formation corresponds to 0 K and is in kcal/mol

        (2) directories containing pf.dat files for each species
            - Each directory name must correspond to name in species.csv file
            - all pf.dat files must be named as pf_<model>.dat
              where <model> corresponds to some thermo mode input by the
              user at the command line with option -m
"""

import os
import sys
import argparse
from io import StringIO
import pandas
import automol
from ioformat import pathtools
from phydat import phycon
import autorun


# Set the path to the current working directory
CWD = os.getcwd()

# General input information for thermo for all species
ENTHALPYT = 0.
BREAKT = 1000.


def main():
    """ Loop over the species
    """
    print('\nGENERATING THERMO')
    print('-----------------\n')

    # Parse the input
    opts = _command_line()
    print(f'Generating NASA polynomials for model {opts["model"]}...')

    spc_lst = _read_csv(opts['input'])

    # Loop over list and build the thermo block
    full_ckin_str = ''
    for spc in spc_lst:
        # Unpack species
        name, fml, hform0 = spc
        print('  - Generating NASA polynomial '
              f'for species {name}...')

        # Grab the correct pf string and calculate thermo
        pf_str = _read_pfdat_str(name, opts['model'])
        run_dir = _make_run_dir_path(name, opts['model'])
        hform298, nasa_poly = _run_autorun(name, fml, hform0, pf_str, run_dir)

        # Write the ckin strings
        full_ckin_str += _write_ckin_str(hform0, hform298, nasa_poly)

    # Write the full ckin string
    _write_full_ckin_str(full_ckin_str, opts['model'])
    print('\nThermo generation script completed successfully. Exiting...\n\n')


# UTIL FUNCTONS FOR MAIN #
# Parser
def _command_line():
    """ parse the command line
    """
    # Package command line into options
    par = argparse.ArgumentParser()
    par.add_argument('-i', '--input', default='species.csv',
                     help='name of input species mechanism file (species.csv)')
    par.add_argument('-m', '--model', default='har', type=str,
                     help='add heat-of-formation species (har)')
    opts = vars(par.parse_args())

    # Generate the output name and add to options dictionary
    output_name = f'{opts["model"]}.ckin'
    opts['output'] = output_name

    # Check if the output exists
    output_path = os.path.join(CWD, output_name)
    if os.path.exists(output_path):
        print(f'Output file {output_name} already exists. Exiting...\n\n')
        sys.exit()

    return opts


def _read_csv(input_name):
    """ Read the CSV string
    """

    def _get_data(input_name):
        """ read csv into dataframe """
        csv_str = pathtools.read_file(CWD, input_name)
        csv_file = StringIO(csv_str)
        data = pandas.read_csv(csv_file, comment='!', quotechar="'")
        data.columns = data.columns.str.strip()
        data.columns = map(str.lower, data.columns)
        return data

    data = _get_data(input_name)
    spc_lst = ()
    for name, smi, hof in zip(data.name, data.smiles, data.hof):
        fml = automol.chi.formula(automol.smiles.chi(smi))
        _hof = hof * phycon.KCAL2EH

        spc_lst += ((name, fml, _hof),)

    return spc_lst


def _read_pfdat_str(name, model):
    """ Read the pf.dat file for the specified species and model
    """
    dat_file_path = os.path.join(CWD, name)
    dat_file_name = f'pf_{model}.dat'
    pf_str = pathtools.read_file(dat_file_path, dat_file_name)

    dat_path = os.path.join(dat_file_path, dat_file_name)
    print(f'    - Reading MESS pf.dat string at {dat_path}')

    return pf_str


# Runner
def _make_run_dir_path(name, model):
    """ Make a path to the run directory
    """
    _name = name.replace('(', '_').replace(')', '_')
    return os.path.join(CWD, 'RUN', _name, model)


def _run_autorun(name, fml, hform0, pf_str, run_dir):
    """ Generate a NASA polynomial for a given species by calling
        ThermP+PAC99
    """

    print(f'    - Running ThermP+PAC99 at {run_dir}')

    fml_str = automol.form.string(fml)

    thermp_script_str = autorun.SCRIPT_DCT['thermp']
    pac99_script_str = autorun.SCRIPT_DCT['pac99'].format(fml_str)

    hform298, nasa_poly = autorun.thermo(
        thermp_script_str, pac99_script_str, run_dir,
        pf_str, name, fml, hform0,
        enthalpyt=ENTHALPYT, breakt=BREAKT, convert=True)

    return hform298, nasa_poly


# Writer
def _write_ckin_str(hform0, hform298, nasa_poly):
    """ Write the THERMO block for the species"
    """
    return (
        f'! Hf(0 K) = {(hform0*phycon.EH2KCAL):.2f}, '
        f'! Hf(298 K) = {(hform298*phycon.EH2KCAL):.2f} kcal/mol\n'
        f'{nasa_poly}\n'
    )


def _write_full_ckin_str(spc_ckin_str, model):
    """ Write the final full string
    """
    ckin_file_name = f'{model}.ckin'
    full_ckin_str = (
        'THERMO\n\n'
        f'{spc_ckin_str}'
        'END'
    )
    print(f'\nWriting thermo mechanism file {ckin_file_name}')
    pathtools.write_file(full_ckin_str, CWD, ckin_file_name)


if __name__ == '__main__':
    main()
