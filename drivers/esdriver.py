""" electronic structure drivers
"""

from mechroutines.es import run_tsk
from mechlib.amech_io import parser
from mechlib.amech_io import printer as ioprinter


def run(pes_rlst, spc_rlst,
        es_tsk_lst,
        spc_dct, thy_dct,
        run_prefix, save_prefix):
    """ Central driver for all electronic structure tasks.

        :param pes_rlst: lst of PES-SUBPES-CHNLS ro tun
        :type pes_rlst: [dict[species, reacs, prods, model]]
        :param spc_rlst: lst of species to run
        :type spc_rlst: [dict[species, reacs, prods, model]]
        :param es_tsk_lst: list of the electronic structure tasks
        :type es_tsk_lst: list[[obj, tsk, keyword_dict]]
        :param spc_dct: species information
        :type spc_dct: dict[spc_name: spc_information]
        :param thy_dct: all of the theory information
        :type thy_dct: dict[]
    """

    # --------------------------------------------- #
    # PREPARE INFORMATION TO PASS TO ESDRIVER TASKS #
    # --------------------------------------------- #

    # Set the appropriate run lst; default to PES if any
    run_rlst = pes_rlst if pes_rlst is not None else spc_rlst

    # -------------------------------- #
    # RUN THE REQUESTED ESDRIVER TASKS #
    # -------------------------------- #

    for run_inf, run_lst in run_rlst.items():

        # Build a TS dictionary and add it to the spc dct if needed
        if any(tsk_lst[0] == 'ts' for tsk_lst in es_tsk_lst):
            pes_idx = run_inf[1]+1
            ts_dct, ts_queue = parser.spc.get_sadpt_dct(
                pes_idx, es_tsk_lst, run_lst,
                thy_dct, spc_dct,
                run_prefix, save_prefix)
            spc_dct = parser.spc.combine_sadpt_spc_dcts(
                ts_dct, spc_dct)

        # Print what is being run PESs that are being run
        ioprinter.runlst(run_inf, run_lst)

        for tsk_lst in es_tsk_lst:

            # Unpack the options
            [obj, tsk, es_keyword_dct] = tsk_lst

            # Build the queue of species based on user request
            if obj == 'all':
                obj_queue = parser.rlst.spc_queue(run_lst) + ts_queue
            if obj == 'spc':
                obj_queue = parser.rlst.spc_queue(run_lst)
            elif obj == 'ts':
                obj_queue = ts_queue
            elif obj == 'vdw':
                obj_queue = ()

            # Run the electronic structure task for all spc in queue
            for spc_name in obj_queue:
                run_tsk(tsk, spc_dct, spc_name,
                        thy_dct, es_keyword_dct,
                        run_prefix, save_prefix)
