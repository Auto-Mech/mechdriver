""" Driver to read data from the SAVE filesystem and then
    processing and formatting it for the purposes for ease of
    reading and presentation.
"""

from mechroutines.proc import run_tsk
from mechlib.amech_io import parser
from mechlib.amech_io import printer as ioprinter


def run(pes_rlst, spc_rlst,
        proc_tsk_lst,
        spc_dct,
        pes_mod_dct, spc_mod_dct, thy_dct,
        run_prefix, save_prefix):
    """ Central driver for all output tasks.

        :param pes_idx: index for the PES where the channel/spc belong to
        :type pes_idx: int
        :param rxn_lst: species and models for all reactions being run
        :type rxn_lst: list[dict[species, reacs, prods, model]]
        :param spc_dct: species information
        :type spc_dct: dict[spc_name: spc_information]
        :param es_tsk_lst: list of the electronic structure tasks
        :type es_tsk_lst: list[[obj, tsk, keyword_dict]]
        :param thy_dct: all of the theory information
        :type thy_dct: dict[]
        :param run_inp_dct: information from input section of run.dat
        :type run_inp_dct: dict[]
    """

    # -------------------------------- #
    # RUN THE REQUESTED ESDRIVER TASKS #
    # -------------------------------- #

    run_rlst = parser.rlst.combine(pes_rlst, spc_rlst)

    for (fml, pes_idx, subpes_idx), run_lst in run_rlst.items():

        # Print what is being run PESs that are being run
        ioprinter.runlst((fml, pes_idx, subpes_idx), run_lst)

        # Loop over Tasks
        for tsk_lst in proc_tsk_lst:

            # Unpack the options
            [obj, tsk, prnt_keyword_dct] = tsk_lst

            # Get the model from the task
            spc_mods, pes_mod = parser.models.extract_models(tsk)
            spc_mod_dct_i = spc_mod_dct[spc_mods[0]]
            pes_mod_dct_i = pes_mod_dct[pes_mod]

            # Build the queue of species based on user request
            if obj == 'spc':
                obj_queue = parser.rlst.spc_queue(run_lst, fml)

            for spc_name in obj_queue:
                run_tsk(
                    tsk, spc_dct, spc_name,
                    thy_dct, prnt_keyword_dct,
                    pes_mod_dct_i, spc_mod_dct_i,
                    run_prefix, save_prefix)
