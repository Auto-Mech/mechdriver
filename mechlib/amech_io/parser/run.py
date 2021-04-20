""" Formats all the strings in th run.dat file.

    Formats strings containing tasks for various drivers:
        ESDriver, ThermoDriver, kTPDriver, TransDriver,

    into specifically formatted Python lists used to schedule and conduct
    all of the tasks requested by the user .
"""

import sys
import ioformat


# Input block
def input_dictionary(run_str):
    """ dict
    """

    inp_block = ioformat.ptt.end_block(run_str, 'input', footer='input')

    inp_dct = ioformat.ptt.keyword_dct_from_block(inp_block)

    return inp_dct


# Chemistry objects
def pes_idxs(run_str):
    """  Build a dictionary of the PESs.
            {pes_idx: [chn_idxs]}
            breaks if pes_idx is given on two lines
    """

    pes_block = ioformat.ptt.end_block(run_str, 'pes', footer='pes')

    if pes_block is not None:
        run_pes = {}
        for line in pes_block.strip().splitlines():
            [pes_nums, chn_nums] = line.split(':')
            _pes_idxs = ioformat.ptt.idx_lst_from_line(pes_nums)
            _chn_idxs = ioformat.ptt.idx_lst_from_line(chn_nums)
            for idx in _pes_idxs:
                run_pes.update({idx: _chn_idxs})
    else:
        run_pes = None

    return run_pes


def spc_idxs(run_str):
    """  Build a dictionary of the PESs.
            {pes_idx: [chn_idxs]}
    """

    spc_block = ioformat.ptt.end_block(run_str, 'spc', footer='spc')

    if spc_block is not None:
        _idxs = ()
        for line in spc_block.splitlines():
            _idxs += ioformat.ptt.idx_lst_from_line(line)
        _spc_idxs = {1: _idxs}
    else:
        _spc_idxs = None

    return _spc_idxs


# DRIVER TASK LISTS
def tasks(run_str, thy_dct, kin_mod_dct, spc_mod_dct):
    """ runstr
    """

    es_block = ioformat.ptt.end_block(run_str, 'els', footer='els')
    trans_block = ioformat.ptt.end_block(run_str, 'trans', footer='trans')
    therm_block = ioformat.ptt.end_block(run_str, 'thermo', footer='thermo')
    ktp_block = ioformat.ptt.end_block(run_str, 'ktp', footer='ktp')
    proc_block = ioformat.ptt.end_block(run_str, 'proc', footer='proc')

    es_tsks = _tsk_lst(es_block, 3)
    therm_tsks = _tsk_lst(therm_block, 2)
    ktp_tsks = _tsk_lst(ktp_block, 2)
    trans_tsks = _tsk_lst(trans_block, 3)
    proc_tsks = _tsk_lst(proc_block, 3)

    # _new_check_dct(mod_tsk_lst, TSK_KEY_DCT, TSK_VAL_DCT, thy_dct)

    return {
        'es': es_tsks,
        'therm': therm_tsks,
        'ktp': ktp_tsks,
        'trans': trans_tsks,
        'proc': proc_tsks
    }


def _tsk_lst(tsk_str, num):
    """ Set the sequence of electronic structure tasks for a given
        species or PESs
    """

    def _split_line(line, num):
        """ Split a line
        """
        line = line.split()
        if num == 3:
            tsk, key_lst = line[:2], line[2:]
        elif num == 2:
            tsk, key_lst = line[:1], line[1:]
        key_dct = ioformat.ptt.keyword_dct_from_block('\n'.join(key_lst))

        return tsk + [key_dct]

    if tsk_str is not None:
        tsks = []
        tsk_str = ioformat.remove_whitespace(tsk_str)
        for line in tsk_str.splitlines():
            try:
                _tsk = _split_line(line, num)
            except:
                print('Task line not formatted correctly:\n{}'.format(line))
                sys.exit()
            tsks.append(_tsk)
        mod_tsks = _expand_tsks(tsks) if num == 3 else tsks
    else:
        mod_tsks = None

    return mod_tsks


EXPAND_DCT = {
    # 'find_ts': ('sadpt_scan', 'sadpt_opt', 'sadpt_hess')  # sadpt_check
}


def _expand_tsks(tsks_lst):
    """ loop over tasks and expand tasks
    """

    # Expand the objects for running
    mod_tsks_lst = []
    for _tsk_lst in tsks_lst:
        [obj, tsk, dct] = _tsk_lst
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


def extract_tsk(tsk, tsk_lst):
    """ Searches for a task in the tsk lst and pulls out the info
        if it is found. Only good if task in lst one time
    """

    tsk_inf = None
    for _tsk_inf in tsk_lst:
        if any(x == tsk for x in _tsk_inf):  # just looks in all pars
            tsk_inf = _tsk_inf
            break

    return tsk_inf
