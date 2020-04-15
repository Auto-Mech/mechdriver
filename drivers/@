""" electronic structure drivers
"""

import routines
from lib.load import species as loadspc
from lib.load import mechanism as loadmech
from lib.filesystem import inf as finf
from lib import printmsg


def run(pes_idx,
        rxn_lst,
        spc_dct,
        cla_dct,
        es_tsk_lst,
        thy_dct,
        run_inp_dct):
    """ driver for all electronic structure tasks
    """

    # Debug print statements
    print('rxn_lst')
    print(rxn_lst)

    # Print the header message for the driver
    printmsg.program_header('es')

    # Pull stuff from dcts for now
    run_prefix = run_inp_dct['run_prefix']
    save_prefix = run_inp_dct['save_prefix']

    # Loop over Tasks
    for tsk_lst in es_tsk_lst:

        # Unpack the options
        [obj, tsk, es_keyword_dct] = tsk_lst

        # Set Task and theory information
        es_ini_key = es_keyword_dct['inplvl']
        es_run_key = es_keyword_dct['runlvl']
        ini_thy_info = finf.get_es_info(es_ini_key, thy_dct)
        thy_info = finf.get_es_info(es_run_key, thy_dct)

        # Build the queue of species based on user request
        if obj == 'spc':
            spc_queue = loadmech.build_spc_queue(rxn_lst)
        elif obj == 'ts':
            spc_queue = []
            ts_dct = loadspc.build_sadpt_dct(
                pes_idx, rxn_lst, thy_info, ini_thy_info,
                run_inp_dct, spc_dct, cla_dct)
            for sadpt in ts_dct:
                if ts_dct[sadpt]['class'] != '':
                    spc_queue.append((sadpt, ''))
                    spc_dct = loadspc.combine_sadpt_spc_dcts(
                        ts_dct, spc_dct)
                else:
                    print('Skipping task for {} since no'.format(sadpt),
                          'class given/identified')
        elif obj == 'vdw':
            spc_queue = []

        # Run the electronic structure task for all spc in queue
        for spc_name, _ in spc_queue:
            routines.es.tsk.run(
                tsk,
                spc_dct,
                spc_name,
                thy_info,
                ini_thy_info,
                run_prefix,
                save_prefix,
                es_keyword_dct=es_keyword_dct)
