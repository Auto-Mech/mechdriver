"""
    Run MESS calculations
"""

import os
import numpy
import mess_io


def messpf_temps(pf_path):
    """ Obtain the temperatures from the MESSPF file
    """

    # Obtain the temperatures, remove the 298.2 value
    temps, _, _, _ = messpf(pf_path)
    temps = [temp for temp in temps if not numpy.isclose(temp, 298.2)]

    return temps


def messpf(pf_path):
    """ Obtain the log partition functions from the MESSPF file
    """
    # Read MESSPF file
    messpf_file = os.path.join(pf_path, 'pf.dat')
    with open(messpf_file, 'r') as pffile:
        output_string = pffile.read()
    temps, logq, dq_dt, dq2_dt2 = mess_io.reader.pfs.partition_function(
        output_string)
    return temps, logq, dq_dt, dq2_dt2
