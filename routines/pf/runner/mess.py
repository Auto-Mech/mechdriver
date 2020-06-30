"""
    Run MESS calculations
"""

import os
from lib.submission import run_script
from lib.submission import DEFAULT_SCRIPT_DCT


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
CUR_PATH = os.path.dirname(os.path.realpath(__file__))


# Set paths to MESS jobs
def messrate_path(prefix, pes_formula, sub_pes_idx):
    """ Build a simple mess path using the run prefix
    """
    pes_str = '{}_{}'.format(pes_formula, sub_pes_idx)
    return os.path.join(prefix, 'MESSRATE', pes_str)


def messpf_path(prefix, spc_info):
    """ Build a simple mess path using the run prefix
    """
    spc_formula = automol.inchi.formula_string(spc_info[0])
    ich_key = automol.inchi.inchi_key(spc_info[0])
    path = os.path.join(prefix, 'MESSPF', spc_formula, ich_key)
        
    # Build the filesystem
    if not os.path.exists(os.path.join(run_prefix, 'MESSRATE')):
        os.mkdir(os.path.join(run_prefix, 'MESSRATE'))
    if not os.path.exists(mess_path):
        os.mkdir(mess_path)

    return path


# Write MESS files
def write_mess_file(mess_inp_str, dat_str_dct, mess_path,
                    fname='mess.inp', overwrite=True):
    """ Write MESS file
    """

    # Write the MESS file
    print('Writing MESS input file...')
    with open(os.path.join(mess_path, fname), 'w') as mess_file:
        mess_file.write(mess_inp_str)

    # Write all of the data files needed
    if dat_str_dct:
        print('Writing the MESS data files...')
    for fname, fstring in dat_str_dct.items():
        dat_path = os.path.join(mess_path, fname)
        print('Writing file: {}'.format(dat_path)
        if string:
            data_file_path = os.path.join(mess_path, name)
            with open(data_file_path, 'w') as data_file:
                data_file.write(string)
    # print(' - WARNING: File will be overwriten.')
    # print('No additional MESS input file will be written.')


def write_cwd_mess_file():
    """ Write a copy of the MESS file in the current working directory
    """
    starting_path = os.getcwd()
    jobdir_mess_path = ''.join([starting_path, '/mess'])
    if not os.path.exists(jobdir_mess_path):
        os.mkdir(jobdir_mess_path)
    jobdir_mess_file_path = ktp_runner.get_mess_path(
        jobdir_mess_path, pes_formula, sub_pes_idx)
    with open(jobdir_mess_file_path, 'a') as file_obj:
        file_obj.write('\n\n! NEW MESS FILE\n')


# Read MESS files
def read_mess_file(mess_path):
    """ read a mess file
    """
    mess_file = os.path.join(mess_path, 'mess.inp')
    print('Searching for MESS input file at {}'.format(mess_path))
    if os.path.exists(mess_file):
        print(' - Found mess.inp at path.')
        with open(mess_file, 'r') as mfile:
            mess_inp_str = mfile.read()
        # Need to read all of the data files
        dat_str_lst = []
    else:
        print(' - No mess.inp file found at path.')
        mess_inp_str = ''
        dat_str_lst = []

    return mess_inp_str, dat_str_lst


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


# Run MESS jobs
def run_rates(mess_path, script_str=DEFAULT_SCRIPT_DCT['messrate']):
    """ Run the mess file that was wriiten
    """
    run_script(script_str, mess_path)


def run_pf(mess_path, script_str=DEFAULT_SCRIPT_DCT['messpf']):
    """ Run the mess file that was wriiten
    """
    run_script(script_str, mess_path)
