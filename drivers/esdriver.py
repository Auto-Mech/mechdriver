""" electronic structure drivers
"""

from routines.es import run_tsk
from lib import filesys
from lib.amech_io import parser


def run(pes_idx,
        rxn_lst,
        spc_dct,
        cla_dct,
        es_tsk_lst,
        thy_dct,
        run_inp_dct):
    """ driver for all electronic structure tasks
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
        mr_scn_thy_info = None
        mr_sp_thy_info = None
        # if es_keyword_dct['mr_scnlvl'] is not None:
        #     mr_scn_thy_info = filesys.inf.get_es_info(
        #         es_keyword_dct['mr_scnlvl'], thy_dct)
        # else:
        #     mr_scn_thy_info = None
        # if es_keyword_dct['mr_splvl'] is not None:
        #     mr_sp_thy_info = filesys.inf.get_es_info(
        #         es_keyword_dct['mr_splvl'], thy_dct)
        # else:
        #     mr_sp_thy_info = None

        # Build the queue of species based on user request
        if obj == 'spc':
            spc_queue = parser.mechanism.build_spc_queue(rxn_lst)
        elif obj == 'ts':
            if not built_dct:
                ts_dct, ts_queue = parser.species.get_sadpt_dct(
                    pes_idx, es_tsk_lst, rxn_lst,
                    thy_dct, run_inp_dct, spc_dct, cla_dct)
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
                mr_sp_thy_info, mr_scn_thy_info,
                run_prefix, save_prefix,
                es_keyword_dct=es_keyword_dct)
