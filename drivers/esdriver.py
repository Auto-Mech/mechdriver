""" electronic structure drivers
"""

import routines
from lib.load import run as loadrun
from lib.load import species as loadspc
from lib.load import mechanism as loadmech
from lib.filesystem import check as fcheck
from lib.filesystem import path as fpath
from lib.filesystem import inf as finf
from lib import printmsg


def run(rxn_lst,
        spc_dct,
        es_tsk_str,
        pes_model_dct, spc_model_dct,
        thy_dct,
        run_options_dct,
        run_inp_dct):
    """ driver for all electronic structure tasks
    """

    # Print the header message for the driver
    printmsg.program_header('es')

    # Pull stuff from dcts for now
    run_prefix = run_inp_dct['run_prefix']
    save_prefix = run_inp_dct['save_prefix']

    print('rxn_lst')
    print(rxn_lst)

    # Do some extra work to prepare the info to pass to the drivers
    es_tsk_lst = loadrun.build_run_es_tsks_lst(
        es_tsk_str, spc_model_dct, thy_dct)

    # Loop over Tasks
    for tsk_lst in es_tsk_lst:

        # Unpack the options
        [obj, tsk, es_run_key, es_ini_key, es_options] = tsk_lst

        # Set Task and theory information
        ini_thy_info = finf.get_es_info(es_ini_key, thy_dct)
        thy_info = finf.get_es_info(es_run_key, thy_dct)

        # Build the queue of species based on the requested obj
        if obj == 'spc':
            spc_queue = loadmech.build_spc_queue(rxn_lst)
        elif obj == 'ts':
            spc_queue = []
            ts_dct = loadspc.build_sadpt_dct(
                rxn_lst, thy_info, ini_thy_info,
                run_inp_dct, spc_dct, {})
            for sadpt in ts_dct:
                spc_queue.append((sadpt, ''))
                spc_dct.update({sadpt: ts_dct[sadpt]})
        elif obj == 'vdw':
            spc_queue = []

        # Loop over all requested species and run the task
        for spc in spc_queue:
            spc_name, _ = spc
            routines.es.run_tsk(
                tsk,
                spc_dct,
                spc_name,
                thy_info,
                ini_thy_info,
                run_prefix,
                save_prefix,
                es_options=es_options)
