"""
    Run MESS calculations
"""

import os
import numpy
import mess_io
from mechlib.amech_io import printer as ioprinter


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


def combine_pfs(pfs, coeffs, operators):
    """ combine partition functions

        _pf = [temps, logq, dq_dt, d2q_dt2]
    """

    assert len(pfs) == len(coeffs)

    for midx, _pf in enumerate(pfs):
        if midx == 0:
            coeff = coeffs[midx]
            final_pf = _pf
        else:
            pf2 = _pf
            coeff = coeffs[midx]
            operator = operators[midx-1]
            if coeff < 0:
                coeff = abs(coeff)
                if operator == 'multiply':
                    operator = 'divide'
            final_pf = _combine_pfs(final_pf, pf2, coeff, operator)

    return final_pf


def _combine_pfs(pfa, pfb, coeff, operator):
    """ Obtain the pf information of the multiplication of pfa and pfb
    """

    tempsa, logqa, dq_dta, d2q_dt2a = pfa
    _, logqb, dq_dtb, d2q_dt2b = pfb

    if operator == 'multiply':
        logq = [a+b+numpy.log(coeff) for a, b in zip(logqa, logqb)]
        dq_dt = [a+b+numpy.log(coeff) for a, b in zip(dq_dta, dq_dtb)]
        d2q_dt2 = [a+b+numpy.log(coeff) for a, b in zip(d2q_dt2a, d2q_dt2b)]
    elif operator == 'divide':
        logq = [a-b-numpy.log(coeff) for a, b in zip(logqa, logqb)]
        dq_dt = [a-b-numpy.log(coeff) for a, b in zip(dq_dta, dq_dtb)]
        d2q_dt2 = [a-b-numpy.log(coeff) for a, b in zip(d2q_dt2a, d2q_dt2b)]

    return tempsa, logq, dq_dt, d2q_dt2
