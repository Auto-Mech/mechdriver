"""
    Run MESS calculations
"""

import os
from ioformat.pathtools import read_file
import mess_io


def rate_strings(rate_paths_dct):
    """ Read the input and output for all of the PESs

        returns dictionary for each pes group for a 'base' MESS run
        and a well-extended MESS run
    """

    rate_strs_dct, mess_paths_dct = {}, {}
    for pes_inf in rate_paths_dct:

        rate_strs_dct[pes_inf] = {}
        mess_paths_dct[pes_inf] = {}

        # Determine if MESS files read from standard/well-extended run
        # Defaults to well-extended if output found for that run
        for typ in ('base', 'wext'):

            mess_path = rate_paths_dct[pes_inf][typ]

            if os.path.exists(os.path.join(mess_path, 'rate.out')):
                rate_strs_dct[pes_inf][typ] = {
                    'inp': read_file(mess_path, 'mess.inp'),
                    'ktp_out': read_file(mess_path, 'rate.out'),
                    'ke_out': read_file(mess_path, 'ke.out'),
                    'ped': read_file(mess_path, 'ped.out'),
                    'aux': read_file(mess_path, 'mess.aux'),
                    'log': read_file(mess_path, 'mess.log')
                }
            else:
                rate_strs_dct[pes_inf][typ] = {}

            mess_paths_dct[pes_inf][typ] = mess_path

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
