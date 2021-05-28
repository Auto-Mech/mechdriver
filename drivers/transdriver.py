""" Driver for energy transfer evaluations including determination of
    Lennard-Jones parameters
"""

from mechroutines.trans import run_tsk
from mechlib.amech_io import parser


def run(spc_rlst,
        trans_tsk_lst,
        spc_dct, thy_dct,
        run_prefix, save_prefix):
    """ main driver for etransfer run
    """

    # Print
    # for spc in RUN_SPC_LST_DCT:
    #     ioprinter.info_message(
    #         'Calculating Thermochem for species: {}'.format(spc),
    #         newline=1)

    # Build a list of the species to calculate thermochem for loops below
    spc_queue = parser.species.build_queue(spc_rlst, 'SPC')

    # Loop over Tasks
    for tsk_lst in trans_tsk_lst:
        [_, tsk, etrans_keyword_dct] = tsk_lst
        run_tsk(tsk, spc_queue,
                spc_dct,
                thy_dct, etrans_keyword_dct,
                run_prefix, save_prefix)
