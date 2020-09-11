"""
    Run MESS calculations
"""

import os
import numpy
import automol
import mess_io
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


def messpf_path(prefix, inchi):
    """ Build a simple mess path using the run prefix
    """
    spc_formula = automol.inchi.formula_string(inchi)
    ich_key = automol.inchi.inchi_key(inchi)
    return os.path.join(prefix, 'MESSPF', spc_formula, ich_key)


# Write MESS files
def write_mess_file(mess_inp_str, dat_str_dct, mess_path,
                    filename='mess.inp'):
    """ Write MESS file
        add overwrite keyword to handle overwrite MESS input
    """

    # Write the MESS file
    if not os.path.exists(mess_path):
        os.makedirs(mess_path)
    print('\n\nWriting MESS input file...')
    print(' - Path: {}'.format(mess_path))
    with open(os.path.join(mess_path, filename), 'w') as mess_file:
        mess_file.write(mess_inp_str)

    # Write all of the data files needed
    if dat_str_dct:
        print('Writing the MESS data files...')
    for fname, fstring in dat_str_dct.items():
        # dat_path = os.path.join(mess_path, fname)
        if fstring:
            data_file_path = os.path.join(mess_path, fname)
            print(' - Writing file: {}'.format(data_file_path))
            with open(data_file_path, 'w') as file_obj:
                file_obj.write(fstring)
    # print(' - WARNING: File will be overwriten.')
    # print('No additional MESS input file will be written.')


def write_cwd_pf_file(mess_str, inchi, fname='pf.inp'):
    """ Write a copy of the MESS file in the current working directory
    """

    # Set starting path
    starting_path = os.getcwd()

    # Set the MESS paths and build dirs if needed
    jobdir_mess_path = messpf_path(starting_path, inchi)
    if not os.path.exists(jobdir_mess_path):
        os.makedirs(jobdir_mess_path)

    # Write the files
    file_path = os.path.join(jobdir_mess_path, fname)
    if os.path.exists(file_path):
        for i in range(1, 51):
            if not os.path.exists(file_path+str(i+1)):
                fin_path = file_path+str(i+1)
                break
    else:
        fin_path = file_path
    with open(fin_path, 'w') as file_obj:
        file_obj.write(mess_str)

    print('Saving MESS input copy at {}'.format(fin_path))
    return fin_path

def write_cwd_rate_file(mess_str, pes_formula, sub_pes_idx, fname='mess.inp'):
    """ Write a copy of the MESS file in the current working directory
    """

    # Set starting path
    starting_path = os.getcwd()

    # Set the MESS paths and build dirs if needed
    jobdir_mess_path = messrate_path(
        starting_path, pes_formula, sub_pes_idx)
    if not os.path.exists(jobdir_mess_path):
        os.makedirs(jobdir_mess_path)

    # Write the files
    file_path = os.path.join(jobdir_mess_path, fname)
    if os.path.exists(file_path):
        for i in range(1, 51):
            if not os.path.exists(file_path+str(i+1)):
                fin_path = file_path+str(i+1)
                break
    else:
        fin_path = file_path
    with open(fin_path, 'w') as file_obj:
        file_obj.write(mess_str)

    print('Saving MESS input copy at {}'.format(fin_path))


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


def write_mess_output(formulastr, final_pf, mess_path, filename='pf.dat'):
    """ Write a mess output file for a pf file
    """
    mess_out_str = 'Natural log of the partition function and its derivatives:\n'
    mess_out_str += ' T, K            {}'.format(formulastr)
    # Write MESS output string
    temps, logq, dq_dt, d2q_dt2 = final_pf
    for idx in range(len(temps)):
        mess_out_str += '\n'
        mess_out_str += str(temps[idx]) + '\t'
        mess_out_str += str(logq[idx]) + '\t'
        mess_out_str += str(dq_dt[idx]) + '\t'
        mess_out_str += str(d2q_dt2[idx]) + '\t'
    # Write the MESS file
    if not os.path.exists(mess_path):
        os.makedirs(mess_path)
    print('\n\nWriting MESS Output file...')
    print(' - Path: {}'.format(mess_path))
    with open(os.path.join(mess_path, filename), 'w') as mess_file:
        mess_file.write(mess_out_str)


def read_messpf_temps(pf_path):
    """ Obtain the temperatures from the MESSPF file
    """

    # Obtain the temperatures, remove the 298.2 value
    temps, _, _, _ = read_messpf(pf_path)
    temps = [temp for temp in temps if not numpy.isclose(temp, 298.2)]

    return temps


def read_messpf(pf_path):
    """ Obtain the log partition functions from the MESSPF file
    """
    # Read MESSPF file
    messpf_file = os.path.join(pf_path, 'pf.dat')
    with open(messpf_file, 'r') as pffile:
        output_string = pffile.read()
    temps, logq, dq_dt, dq2_dt2 = mess_io.reader.pfs.partition_fxn(
        output_string)
    return temps, logq, dq_dt, dq2_dt2


def multiply_pfs(pfa, pfb, coeff):
    """ Obtain the pf information of the multiplication of pfa and pfb
    """
    tempsa, logqa, dq_dta, d2q_dt2a = pfa
    tempsb, logqb, dq_dtb, d2q_dt2b = pfb
    logq = [a+b+numpy.log(coeff) for a, b in zip(logqa, logqb)]
    dq_dt = [a+b+numpy.log(coeff) for a, b in zip(dq_dta, dq_dtb)]
    d2q_dt2 = [a+b+numpy.log(coeff) for a, b in zip(d2q_dt2a, d2q_dt2b)]
    return tempsa, logq, dq_dt, d2q_dt2


def divide_pfs(pfa, pfb, coeff):
    """ Obtain the pf information of the multiplication of pfa and pfb
    """
    tempsa, logqa, dq_dta, d2q_dt2a = pfa
    tempsb, logqb, dq_dtb, d2q_dt2b = pfb
    logq = [a-b-numpy.log(coeff) for a, b in zip(logqa, logqb)]
    dq_dt = [a-b-numpy.log(coeff) for a, b in zip(dq_dta, dq_dtb)]
    d2q_dt2 = [a-b-numpy.log(coeff) for a, b in zip(d2q_dt2a, d2q_dt2b)]
    return tempsa, logq, dq_dt, d2q_dt2


def run_rates(mess_path, script_str=DEFAULT_SCRIPT_DCT['messrate']):
    """ Run the mess file that was wriiten
    """
    run_script(script_str, mess_path)


def run_pf(mess_path, script_str=DEFAULT_SCRIPT_DCT['messpf']):
    """ Run the mess file that was wriiten
    """
    if os.path.exists(mess_path):
    #if os.path.exists(os.path.join(mess_path, 'pf.inp')):
        print('Running MESS input file...')
        print(' - Path: {}'.format(mess_path))
        run_script(script_str, mess_path)
    else:
        print('No MESS input file at path: {}'.format(mess_path))
