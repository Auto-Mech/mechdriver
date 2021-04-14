""" Formats all the strings in th run.dat file.

    Formats strings containing tasks for various drivers:
        ESDriver, ThermoDriver, kTPDriver, TransDriver, 

    into specifically formatted Python lists used to schedule and conduct
    all of the tasks requested by the user .
"""

import sys
import ioformat


# CHEM OBJECTS
def pes_idxs(string):
    """  Build a dictionary of the PESs.
            {pes_idx: [chn_idxs]}
            breaks if pes_idx is given on two lines
    """

    run_pes = {}
    for line in string.strip().splitlines():
        [pes_nums, chn_nums] = line.split(':')
        pes_idxs = ioformat.ptt.idx_lst_from_line(pes_nums)
        chn_idxs = ioformat.ptt.idx_lst_from_line(chn_nums)
        for idx in pes_idxs:
            run_pes.update({idx: chn_idxs})

    return run_pes


def spc_idxs(string):
    """  Build a dictionary of the PESs.
            {pes_idx: [chn_idxs]}
    """

    spc_idxs = ()
    for line in string.splitlines():
        spc_idxs += ioformat.ptt.idx_lst_from_line(line)

    return {1: spc_idxs}


# DRIVER TASK LISTS
def tsk_lst(tsk_str):
    """ Take the es tsk list string from input and set the tasks
        Right now, we presume the tasks given in the file are correct

        Needed for anything that takes a thy_dct
    """

    tsk_lst = _tsk_lst(es_tsk_str)
    mod_tsk_lst = _expand_tsks(tsk_lst)

    return mod_tsk_lst


def _tsk_lst(tsk_str):
    """ Set the sequence of electronic structure tasks for a given
        species or PESs
    """
    tsk_lst = []
    tsk_str = ioformat.remove_whitespace(tsk_str)
    for line in tsk_str.splitlines():
        try:
            tsk_line = line.split()
            obj, tsk, keyword_lst = tsk_line[0], tsk_line[1], tsk_line[2:]
            keyword_dct = ioformat.ptt.keyword_dct_from_block(
                '\n'.join(keyword_lst))
        except:
            print('Task line not formatted correctly:\n{}'.format(line))
            sys.exit()
        tsk_lst.append([obj, tsk, keyword_dct])

    return tsk_lst


EXPAND_DCT = {
    # 'find_ts': ('sadpt_scan', 'sadpt_opt', 'sadpt_hess')  # sadpt_check
}

def expand_tsks(tsks_lst):
    """ loop over tasks and expand tasks
    """

    # Expand the objects for running
    mod_tsks_lst = []
    for tsk_lst in tsks_lst:
        [obj, tsk, dct] = tsk_lst
        if obj == 'all':
            mod_tsks_lst.append(['spc', tsk, dct])
            mod_tsks_lst.append(['ts', tsk, dct])
        else:
            mod_tsks_lst.append([obj, tsk, dct])

    # Expand the tasks
    mod_tsks_lst2 = []
    for tsk_lst in mod_tsks_lst:
        [obj, tsk, dct] = tsk_lst
        expanded_tsks = EXPAND_DCT.get(tsk, None)
        if expanded_tsks is None:
            mod_tsks_lst2.append(tsk_lst)
        else:
            for tsk in expanded_tsks:
                mod_tsks_lst2.append([obj, tsk, dct])

    return mod_tsks_lst2
