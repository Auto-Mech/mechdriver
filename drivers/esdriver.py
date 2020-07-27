""" electronic structure drivers
"""

from routines.es import run_tsk
from lib import filesys
from lib.amech_io import parser
from lib.structure import instab


def run(pes_idx,
        rxn_lst,
        spc_dct,
        cla_dct,
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
        :param cla_dct: input to change class dict
        :type cla_dct: dict[]
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

    # Initialize variable for building a dct for ts
    built_dct = False

    # Loop over Tasks
    print('\nRunning electronic structure tasks given in the input...')
    for tsk_lst in es_tsk_lst:

        # Unpack the options
        [obj, tsk, es_keyword_dct] = tsk_lst

        # Set Task and theory information
        ini_thy_info = filesys.inf.get_es_info(
            es_keyword_dct['inplvl'], thy_dct)
        thy_info = filesys.inf.get_es_info(
            es_keyword_dct['runlvl'], thy_dct)
        var_scn_thy_info = None
        var_sp2_thy_info = None
        if es_keyword_dct['var_scnlvl'] is not None:
            var_scn_thy_info = filesys.inf.get_es_info(
                es_keyword_dct['var_scnlvl'], thy_dct)
        else:
            var_scn_thy_info = None
        var_sp1_thy_info = None
        if 'var_splvl1' in es_keyword_dct:
            if es_keyword_dct['var_splvl1'] is not None:
                var_sp1_thy_info = filesys.inf.get_es_info(
                    es_keyword_dct['var_splvl1'], thy_dct)
        # if es_keyword_dct['var_splvl2'] is not None:
        #     var_sp2_thy_info = filesys.inf.get_es_info(
        #         es_keyword_dct['var_splvl2'], thy_dct)
        # else:
        #     var_sp2_thy_info = None

        # Build the queue of species based on user request
        if obj == 'spc':
            spc_queue = parser.species.build_spc_queue(rxn_lst)
        elif obj == 'ts':
            if not built_dct:
                rxndirn = es_keyword_dct['rxndirn']
                ts_dct, ts_queue = parser.species.get_sadpt_dct(
                    pes_idx, es_tsk_lst, rxn_lst,
                    thy_dct, run_inp_dct, spc_dct, cla_dct,
                    direction=rxndirn)
                spc_dct = parser.species.combine_sadpt_spc_dcts(
                    ts_dct, spc_dct)
                built_dct = True
            spc_queue = ts_queue
        elif obj == 'vdw':
            spc_queue = []

        # Run the electronic structure task for all spc in queue
        for spc_name, _ in spc_queue:
            run_tsk(
                tsk, spc_dct, spc_name,
                thy_info, ini_thy_info,
                var_sp1_thy_info, var_sp2_thy_info, var_scn_thy_info,
                run_prefix, save_prefix,
                es_keyword_dct=es_keyword_dct)
