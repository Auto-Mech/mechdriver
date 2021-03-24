""" electronic structure drivers
"""

from mechroutines.es import run_tsk
from mechlib.amech_io import parser


def run(pes_idx,
        rxn_lst,
        spc_dct,
        es_tsk_lst,
        thy_dct,
        run_inp_dct):
    """ Central driver for all electronic structure tasks.

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

    # Pull stuff from dcts for now
    run_prefix = run_inp_dct['run_prefix']
    save_prefix = run_inp_dct['save_prefix']

    # Build a TS dictionary and add it to the spc dct if needed
    if any(tsk_lst[0] == 'ts' for tsk_lst in es_tsk_lst):
        ts_dct, ts_queue = parser.species.get_sadpt_dct(
            pes_idx, es_tsk_lst, rxn_lst,
            thy_dct, run_inp_dct, spc_dct,
            run_prefix, save_prefix,
            direction='forw')
        spc_dct = parser.species.combine_sadpt_spc_dcts(
            ts_dct, spc_dct)

    # Loop over Tasks
    for tsk_lst in es_tsk_lst:

        # Unpack the options
        [obj, tsk, es_keyword_dct] = tsk_lst

        # Build the queue of species based on user request
        if obj == 'all':
            obj_queue = parser.species.build_spc_queue(rxn_lst) + ts_queue
        if obj == 'spc':
            obj_queue = parser.species.build_spc_queue(rxn_lst)
        elif obj == 'ts':
            obj_queue = ts_queue
        elif obj == 'vdw':
            obj_queue = []

        # Run the electronic structure task for all spc in queue
        for spc_name, _ in obj_queue:
            run_tsk(tsk, spc_dct, spc_name,
                    thy_dct, es_keyword_dct,
                    run_prefix, save_prefix)
