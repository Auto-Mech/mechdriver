"""
    Run MESS calculations
"""

import os
import mess_io


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
