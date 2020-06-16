"""
  Run rates with MESS
"""

import os
from lib.submission import run_script
from lib.submission import DEFAULT_SCRIPT_DCT


def get_mess_path(prefix, pes_formula, sub_pes_idx):
    """ Build a simple mess path using the run prefix
    """
    pes_str = '{}_{}'.format(pes_formula, sub_pes_idx)
    return os.path.join(prefix, 'MESSRATE', pes_str)


def write_mess_file(mess_inp_str, dat_str_lst, mess_path, fname='mess.inp'):
    """ Write MESS file
    """

    # Write the MESS file
    with open(os.path.join(mess_path, fname), 'w') as mess_file:
        mess_file.write(mess_inp_str)

    # Write all of the data files needed
    # for dct in dat_str_lst:
    #     for data in dct.values():
    #         string, name = data
    #         # print('dat test', string, name)
    #         if string:
    #             data_file_path = os.path.join(mess_path, name)
    #             with open(data_file_path, 'w') as data_file:
    #                 data_file.write(string)


def read_mess_file(mess_path):
    """ read a mess file
    """
    dat_str_lst = []
    mess_file = os.path.join(mess_path, 'mess.inp')
    print('Searching for MESS input file at {}'.format(mess_path))
    if os.path.exists(mess_file):
        print(' - Found mess.inp at path.')
        print(' - WARNING: File will be overwriten.')
        # print('No additional MESS input file will be written...')
        with open(mess_file, 'r') as mfile:
            mess_inp_str = mfile.read()
    else:
        print(' - No mess.inp file found at path.')
        mess_inp_str = ''
        dat_str_lst = []

    return mess_inp_str, dat_str_lst


def run_rates(mess_path):
    """ Run the mess file that was wriiten
    """
    run_script(DEFAULT_SCRIPT_DCT['messrate'], mess_path)
