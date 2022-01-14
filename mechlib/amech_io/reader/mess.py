"""
    Run MESS calculations
"""

import os
import ioformat.pathtools
import mess_io


def rate_strings(rate_paths_dct):
    """ Read the input and output for all of the PESs
    """

    rate_strs_dct, mess_paths_dct = {}, {}
    for pes_inf in rate_paths_dct:

        # Determine if MESS files read from standard/well-extended run
        # Defaults to well-extended if output found for that run
        base_mess_path = rate_paths_dct[pes_inf]['base']
        wext_mess_path = rate_paths_dct[pes_inf]['wext']
        if os.path.exists(os.path.join(wext_mess_path, 'rate.out')):
            mess_path = wext_mess_path
        else:
            mess_path = base_mess_path

        rate_strs_dct[pes_inf] = {
            'inp': ioformat.pathtools.read_file(mess_path, 'mess.inp'),
            'ktp_out': ioformat.pathtools.read_file(mess_path, 'rate.out'),
            'ke_out': ioformat.pathtools.read_file(mess_path, 'ke.out'),
            'ped': ioformat.pathtools.read_file(mess_path, 'ped.out'),
            'aux': ioformat.pathtools.read_file(mess_path, 'mess.aux'),
            'log': ioformat.pathtools.read_file(mess_path, 'mess.log')
        }
        mess_paths_dct[pes_inf] = mess_path

    return rate_strs_dct, mess_paths_dct


def messpf(pf_path):
    """ Obtain the log partition functions from the MESSPF file
    """
    # Read MESSPF file
    messpf_file = os.path.join(pf_path, 'pf.dat')
    with open(messpf_file, mode='r', encoding='utf-8') as pffile:
        output_string = pffile.read()
    temps, logq, dq_dt, dq2_dt2 = mess_io.reader.pfs.partition_function(
        output_string)
    return temps, logq, dq_dt, dq2_dt2
