""" Driver for energy transfer evaluations including determination of
    Lennard-Jones parameters
"""

from mechroutines.trans import run_tsk
from mechlib.amech_io import parser


def run(spc_dct,
        thy_dct,
        rxn_lst,
        trans_tsk_lst,
        run_inp_dct):
    """ main driver for etransfer run
    """

    # Pull stuff from dcts for now
    save_prefix = run_inp_dct['save_prefix']
    run_prefix = run_inp_dct['run_prefix']

    # Build a list of the species to calculate thermochem for loops below
    spc_queue = parser.species.build_queue(rxn_lst)

    # Loop over Tasks
    for tsk_lst in trans_tsk_lst:
        [obj, tsk, etrans_keyword_dct] = tsk_lst
        run_tsk(tsk, spc_queue,
                spc_dct,
                thy_dct, etrans_keyword_dct,
                run_prefix, save_prefix)
