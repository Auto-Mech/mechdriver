""" Driver to read data from the SAVE filesystem and then
    processing and formatting it for the purposes for ease of
    reading and presentation.
"""

from mechroutines.proc import run_tsk
from mechlib.amech_io import parser


def run(pes_rlst, spc_rlst,
        proc_tsk_lst,
        spc_dct, thy_dct,
        pes_mod_dct, spc_mod_dct,
        run_prefix, save_prefix):
    """ Central driver for all output tasks.

        :param spc_dct: species information
        :type spc_dct: dict[spc_name: spc_information]
        :param es_tsk_lst: list of the electronic structure tasks
        :type es_tsk_lst: list[[obj, tsk, keyword_dict]]
        :param thy_dct: all of the theory information
        :type thy_dct: dict[]
        :param run_inp_dct: information from input section of run.dat
        :type run_inp_dct: dict[]
    """

    # Build the spc and ts queues
    spc_queue, ts_queue = (), ()
    run_rlst = parser.rlst.combine(pes_rlst, spc_rlst)

    ts_queue = ()
    for (fml, pes_idx, subpes_idx), run_lst_i in run_rlst.items():
        ioprinter.runlst((fml, pes_idx, subpes_idx), run_lst_i)
        if fml != 'SPC':
            if any(tsk_lst[0] == 'ts' for tsk_lst in proc_tsk_lst):
                ts_dct, _queue = parser.spc.ts_dct_from_proctsks(
                    pes_idx, proc_tsk_lst, run_lst_i,
                    thy_dct, spc_dct,
                    run_prefix, save_prefix)
                # doesnt allow for info from .dat file or internal defaults
                spc_dct.update(ts_dct)
                # Doesnt work it ts species from earlier part of for loop lost
                # spc_dct = parser.spc.combine_sadpt_spc_dcts(
                #     ts_dct, spc_dct, glob_dct)
                ts_queue += _queue
            spc_queue += parser.rlst.spc_queue(run_lst_i, fml)
        else:
            spc_queue += run_lst_i

    # Remove species from queue
    spc_queue = tuple(i for n, i in enumerate(spc_queue)
                      if i not in spc_queue[:n])

    # Loop over Tasks
    for tsk_lst in proc_tsk_lst:

        [obj, tsk, prnt_keyword_dct] = tsk_lst

        # Build the queue of species based on user request
        if obj == 'all':
            obj_queue = spc_queue + ts_queue
        if obj == 'spc':
            obj_queue = spc_queue
        elif obj == 'ts':
            obj_queue = ts_queue

        # Set up model dictionaries for the code to use
        spc_mod_dct_i = spc_mod_dct['global']
        pes_mod_dct_i = pes_mod_dct['global']

        run_tsk(
            tsk, obj_queue,
            prnt_keyword_dct,
            spc_dct, thy_dct,
            spc_mod_dct_i, pes_mod_dct_i,
            run_prefix, save_prefix)
