"""
    Run MESS calculations
"""

import os
import mess_io
from mechlib.amech_io import printer as ioprinter


def output(formulastr, final_pf, mess_path, filename='pf.dat'):
    """ Write a mess output file for a pf file
    """

    mess_out_str = mess_io.writer.pf_output(formulastr, *final_pf)

    ioprinter.messpf('write_output', mess_path)
    if not os.path.exists(mess_path):
        os.makedirs(mess_path)
    with open(os.path.join(mess_path, filename), 'w') as mess_file:
        mess_file.write(mess_out_str)
