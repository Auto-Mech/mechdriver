""" pf stuff
"""

import os
import mess_io


def get_pf_header(temp_step, ntemps):
    """ prepare partition function header string
    """
    global_pf_str = mess_io.writer.global_pf(
        [], temp_step, ntemps, rel_temp_inc=0.001,
        atom_dist_min=0.6)
    return global_pf_str


def get_pf_input(spc, spc_str, global_pf_str, zpe_str):
    """ prepare the full pf input string for running messpf
    """

    # create a messpf input file
    spc_head_str = 'Species ' + spc
    print('pf string test:', global_pf_str, spc_head_str, spc_str, zpe_str)
    pf_inp_str = '\n'.join(
        [global_pf_str, spc_head_str,
    #     spc_str, '\n'])
         spc_str, zpe_str, '\n'])
    return pf_inp_str


def write_pf_input(pf_inp_str, data_str_dct, pf_path):
    """ write the pf.inp file
    """
    # Write the data file
    for data in data_str_dct.values():
        string, name = data
        if string:
            data_file_path = os.path.join(pf_path, name)
            with open(data_file_path, 'w') as data_file:
                data_file.write(string)
            print('data string test')
            print(name)
            print(string)

    # Write the MESSPF input file
    with open(os.path.join(pf_path, 'pf.inp'), 'w') as pf_file:
        pf_file.write(pf_inp_str)
