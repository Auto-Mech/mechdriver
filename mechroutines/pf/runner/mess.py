"""
    Run MESS calculations
"""

import os
import numpy
import automol
import mess_io
from autorun import run_script
from mechlib.submission import DEFAULT_SCRIPT_DCT
from mechlib.amech_io import printer as ioprinter


# Read MESS files
def write_mess_output(formulastr, final_pf, mess_path, filename='pf.dat'):
    """ Write a mess output file for a pf file
    """
    
    mess_out_str = mess_io.writer.pf_output(FORMULA, *final_pf)

    ioprinter.messpf('write_output', mess_path)
    if not os.path.exists(mess_path):
        os.makedirs(mess_path)
    with open(os.path.join(mess_path, filename), 'w') as mess_file:
        mess_file.write(mess_out_str)


def read_messpf_temps(pf_path):
    """ Obtain the temperatures from the MESSPF file
    """
    
    messpf_file = os.path.join(pf_path, 'pf.dat')
    with open(messpf_file, 'r') as pffile:
        output_string = pffile.read()
    ret = mess_io.reader.pfs.partition_function(output_string)

    temps = [temp for temp in ret[0] if not numpy.isclose(temp, 298.2)]

    return temps


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
