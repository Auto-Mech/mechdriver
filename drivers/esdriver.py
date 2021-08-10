""" Driver which executes all of the electronic structure tasks
    for all species that comprise channels of input PES and and separate
    species.

    Two Main Loops of Driver:
        (1) species and transition states from input PES or SPC lists
        (2) electronic structure task list

    All tasks done following way:
        (1) Search for electronic structure data in save filesystem
        (2) Write, Run, Parse electronic structure job(s) in run fs
        (3) Write final data into save filesystem
"""

from mechroutines.es import run_tsk
from mechlib.amech_io import parser
from mechlib.amech_io import printer as ioprinter


def run(pes_rlst, spc_rlst,
        es_tsk_lst,
        spc_dct, glob_dct, thy_dct,
        run_prefix, save_prefix,
        print_debug=False):
    """ Executes all electronic structure tasks.

        :param pes_rlst: species from PESs to run
            [(PES formula, PES idx, SUP-PES idx)
            (CHANNEL idx, (REACS, PRODS))
        :type pes_rlst: tuple(dict[str: dict])
        :param spc_rlst: lst of species to run
        :type spc_rlst: tuple(dict[str: dict])
        :param es_tsk_lst: list of the electronic structure tasks
            tuple(tuple(obj, tsk, keyword_dict))
        :type es_tsk_lst: tuple(tuple(str, str, dict))
        :param spc_dct: species information
            dict[spc_name: spc_information]
        :type spc_dct: dict[str:dict]
        :param glob_dct: global information for all species
            dict[spc_name: spc_information]
        :type glob_dct: dict[str: dict]
        :param thy_dct: all of the theory information
            dict[thy name: inf]
        :type thy_dct: dict[str:dict]
        :param run_prefix: root-path to the run-filesystem
        :type run_prefix: str
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
    """

    # -------------------------------- #
    # RUN THE REQUESTED ESDRIVER TASKS #
    # -------------------------------- #

    # Set the appropriate run lst; default to PES if any
    # Runs through PESs, then SPC
    run_rlst = parser.rlst.combine(pes_rlst, spc_rlst)

    for (fml, pes_idx, subpes_idx), run_lst in run_rlst.items():

        # Print what is being run PESs that are being run
        ioprinter.runlst((fml, pes_idx, subpes_idx), run_lst)

        # Build a TS dictionary and add it to the spc dct if needed
        if (fml != 'SPC' and
           any(tsk_lst[0] == 'ts' for tsk_lst in es_tsk_lst)):
            ts_dct, ts_queue = parser.spc.ts_dct_from_estsks(
                pes_idx, es_tsk_lst, run_lst,
                thy_dct, spc_dct,
                run_prefix, save_prefix)
            spc_dct = parser.spc.combine_sadpt_spc_dcts(
                ts_dct, spc_dct, glob_dct)
        else:
            ts_queue = ()

        # Loop over the tasks
        for tsk_lst in es_tsk_lst:

            # Unpack the options
            [obj, tsk, es_keyword_dct] = tsk_lst

            # Build the queue of species based on user request
            if obj == 'all':
                obj_queue = parser.rlst.spc_queue(run_lst, fml) + ts_queue
            if obj == 'spc':
                obj_queue = parser.rlst.spc_queue(run_lst, fml)
            elif obj == 'ts':
                obj_queue = ts_queue
            elif obj == 'vdw':
                obj_queue = ()

            # Run the electronic structure task for all spc in queue
            for spc_name in obj_queue:
                run_tsk(tsk, spc_dct, spc_name,
                        thy_dct, es_keyword_dct,
                        run_prefix, save_prefix,
                        print_debug=print_debug)
