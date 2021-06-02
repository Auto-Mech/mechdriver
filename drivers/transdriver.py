""" Driver for determining and writing parameters that describe
    energy transfer processes.

    Main Loop of Driver:
        (1) PES or SPC list

    Main Workflow:
        (1) OneDMin:
            (1) Calculate Lennard-Jones sigma and epsilon parameters
        (2) Transport File:
            (1) Collate and process data from the SAVE filesystem
            (2) Format and write data into CHEMKIN tranport file
"""

from mechroutines.trans import run_tsk
from mechlib.amech_io import parser
from mechlib.amech_io import printer as ioprinter
from mechlib.reaction import split_unstable_full


def run(pes_rlst, spc_rlst,
        trans_tsk_lst,
        spc_mod_dct,
        spc_dct, thy_dct,
        run_prefix, save_prefix):
    """ main driver for etransfer run
    """

    # Print Header
    ioprinter.info_message('Calculating Transport:')
    ioprinter.runlst(('SPC', 0, 0), spc_rlst)

    # ---------------------------------------------------- #
    # PREPARE INFORMATION TO PASS TO TRANSPORTDRIVER TASKS #
    # ---------------------------------------------------- #

    spc_mods = list(spc_mod_dct.keys())  # hack
    spc_mod_dct_i = spc_mod_dct[spc_mods[0]]
    split_rlst = split_unstable_full(
        pes_rlst, spc_rlst, spc_dct, spc_mod_dct_i, save_prefix)
    spc_queue = parser.rlst.spc_queue(
        tuple(split_rlst.values())[0], 'SPC')

    # --------------------------------------- #
    # RUN THE REQUESTED TRANSPORTDRIVER TASKS #
    # --------------------------------------- #

    for tsk_lst in trans_tsk_lst:
        [_, tsk, etrans_keyword_dct] = tsk_lst
        run_tsk(tsk, spc_queue,
                spc_dct,
                thy_dct, etrans_keyword_dct,
                run_prefix, save_prefix)
