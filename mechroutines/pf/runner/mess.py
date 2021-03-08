"""
    Run MESS calculations
"""

import os
import numpy
import mess_io
from mechlib.amech_io import printer as ioprinter
from mechlib.amech_io import messpf_path


# Read MESS files
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

    ioprinter.saving('MESS input copy', fin_path)

    return fin_path


def write_mess_output(formulastr, final_pf, mess_path, filename='pf.dat'):
    """ Write a mess output file for a pf file
    """

    mess_out_str = mess_io.writer.pf_output(formulastr, *final_pf)

    ioprinter.messpf('write_output', mess_path)
    if not os.path.exists(mess_path):
        os.makedirs(mess_path)
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
    temps, logq, dq_dt, dq2_dt2 = mess_io.reader.pfs.partition_function(
        output_string)
    return temps, logq, dq_dt, dq2_dt2


def multiply_pfs(pfa, pfb, coeff):
    """ Obtain the pf information of the multiplication of pfa and pfb
    """

    tempsa, logqa, dq_dta, d2q_dt2a = pfa
    _, logqb, dq_dtb, d2q_dt2b = pfb

    logq = [a+b+numpy.log(coeff) for a, b in zip(logqa, logqb)]
    dq_dt = [a+b+numpy.log(coeff) for a, b in zip(dq_dta, dq_dtb)]
    d2q_dt2 = [a+b+numpy.log(coeff) for a, b in zip(d2q_dt2a, d2q_dt2b)]

    return tempsa, logq, dq_dt, d2q_dt2


def divide_pfs(pfa, pfb, coeff):
    """ Obtain the pf information of the multiplication of pfa and pfb
    """

    tempsa, logqa, dq_dta, d2q_dt2a = pfa
    _, logqb, dq_dtb, d2q_dt2b = pfb

    logq = [a-b-numpy.log(coeff) for a, b in zip(logqa, logqb)]
    dq_dt = [a-b-numpy.log(coeff) for a, b in zip(dq_dta, dq_dtb)]
    d2q_dt2 = [a-b-numpy.log(coeff) for a, b in zip(d2q_dt2a, d2q_dt2b)]

    return tempsa, logq, dq_dt, d2q_dt2
