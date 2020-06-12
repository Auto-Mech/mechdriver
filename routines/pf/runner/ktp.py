"""
  Run rates with MESS
"""

import os
import autofile
from lib import filesys
from lib.amech_io import parser
from lib.submission import run_script
from lib.submission import DEFAULT_SCRIPT_DCT


def get_mess_path(prefix, pes_formula, sub_pes_idx):
    """ Build a simple mess path using the run prefix
    """
    pes_str = '{}_{}'.format(pes_formula, sub_pes_idx)
    return os.path.join(prefix, 'MESSRATE', pes_str)


def get_mess_path2(
        tsdct, rxn_save_path,
        model_dct, thy_dct, run_model):
    """ Set the path ot the MESS job
    """

    # Set information about the TS and theory methods
    ts_info = (tsdct['ich'], tsdct['chg'], tsdct['mul'])
    thy_info = parser.model.set_es_model_info(
        model_dct[run_model]['es'], thy_dct)[0]

    orb_restr = filesys.inf.orbital_restriction(ts_info, thy_info)
    ref_level = thy_info[1:3]
    ref_level.append(orb_restr)
    thy_save_fs = autofile.fs.theory(rxn_save_path)
    thy_save_fs[-1].create(ref_level)
    thy_save_path = thy_save_fs[-1].path(ref_level)

    # Set paths to the MESS file
    bld_locs = ['MESS', 0]
    bld_save_fs = autofile.fs.build(thy_save_path)
    bld_save_fs[-1].create(bld_locs)
    mess_path = bld_save_fs[-1].path(bld_locs)
    print('Build Path for MESS rate files:')
    print(mess_path)

    return mess_path


def write_mess_file(mess_inp_str, dat_str_lst, mess_path, fname='mess.inp'):
    """ Write MESS file
    """
    # Write total MESS input string
    print('Writing the MESS input file at {}'.format(mess_path))
    print(mess_inp_str)

    # Write the MESS file
    with open(os.path.join(mess_path, fname), 'w') as mess_file:
        mess_file.write(mess_inp_str)

    # Write all of the data files needed
    for dct in dat_str_lst:
        for data in dct.values():
            string, name = data
            # print('dat test', string, name)
            if string:
                data_file_path = os.path.join(mess_path, name)
                with open(data_file_path, 'w') as data_file:
                    data_file.write(string)


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
