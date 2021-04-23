""" runtime messages
"""

from mechlib.amech_io.printer import message
from mechlib.amech_io.printer._pes import pes as print_pes
from mechlib.amech_io.printer._pes import channel as print_channel


def runlst(run_inf, run_lst):
    """ checks if run lst is a species lst
    """

    formula, pes_idx, sub_pes_idx = run_inf
    if formula != 'SPC':
        print_pes(pes_idx, formula, sub_pes_idx)
        for i, rxn in enumerate(run_lst):
           print_channel(i+1, rxn[0], rxn[1])
    else:
        for i, spc in enumerate(run_lst):
            message('Running SPC {}: {}'.format(i+1, spc))
