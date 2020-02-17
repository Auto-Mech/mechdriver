"""
  Run rates with MESS
"""

import os
import autofile

# New Libs
from lib.load import model as loadmodel
from lib.runner import script
from lib.filesystem import orb as fsorb
from lib.filesystem import inf as finf


def run_rates(
        header_str, energy_trans_str,
        well_str, bim_str, ts_str, dat_str_lst,
        tsdct, rxn_save_path,
        model_dct, thy_dct, run_model):
    """ Generate k(T,P) by first compiling all the MESS strings
        and then running MESS
    """

    # Set information about the TS and theory methods
    ts_info = (tsdct['ich'], tsdct['chg'], tsdct['mul'])
    thy_info = loadmodel.set_es_model_info(
        model_dct[run_model]['es'], thy_dct)[0]

    orb_restr = fsorb.orbital_restriction(ts_info, thy_info)
    ref_level = thy_info[1:3]
    ref_level.append(orb_restr)
    thy_save_fs = autofile.fs.theory(rxn_save_path)
    thy_save_fs[-1].create(ref_level)
    thy_save_path = thy_save_fs[-1].path(ref_level)

    mess_inp_str = '\n'.join(
        [header_str, energy_trans_str, well_str, bim_str, ts_str])

    print('mess input file')
    print(mess_inp_str)

    # Set paths to the MESS file
    bld_locs = ['MESS', 0]
    bld_save_fs = autofile.fs.build(thy_save_path)
    bld_save_fs[-1].create(bld_locs)
    mess_path = bld_save_fs[-1].path(bld_locs)
    print('Build Path for MESS rate files:')
    print(mess_path)

    # Write the MESS file
    with open(os.path.join(mess_path, 'mess.inp'), 'w') as mess_file:
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

    # Run MESS on the written file    
    script.run_script(script.MESSRATE, mess_path)
    return mess_path
