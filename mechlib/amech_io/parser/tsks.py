""" Functions for setting the Electronic Structure tasks for
    an outlined procedure
"""

import sys
import ioformat
import automol
from mechlib.amech_io.parser._keywrd import build_tsk_default
from mechlib.amech_io.parser._keywrd import check_dictionary
from mechlib.amech_io.parser._keywrd import check_thy_lvls
from mechlib.amech_io.parser._keywrd import TSK_KEY_DCT
from mechlib.amech_io.parser._keywrd import TSK_VAL_DCT


def es_tsk_lst(es_tsk_str, thy_dct):
    """ Take the es tsk list string from input and set the tasks
        Right now, we presume the tasks given in the file are correct
    """

    tsk_lst = _tsk_lst(es_tsk_str)

    mod_tsk_lst = expand_tsks(tsk_lst)

    _new_check_dct(mod_tsk_lst, TSK_KEY_DCT, TSK_VAL_DCT, thy_dct)

    return mod_tsk_lst


def trans_tsk_lst(trans_tsk_str, thy_dct):
    """ Set the sequence of electronic structure tasks for a given
        species or PESs
    """

    tsk_lst = _tsk_lst(trans_tsk_str)
    # mod_tsk_lst = expand_tsks(tsk_lst)
    # _new_check_dct(mod_tsk_lst, TSK_KEY_DCT, TSK_VAL_DCT, thy_dct)
    _new_check_dct(tsk_lst, TSK_KEY_DCT, TSK_VAL_DCT, thy_dct)

    return tsk_lst


def therm_tsk_lst(therm_tsk_str):
    """ Set sequence of tasks for therm
    """
    raise NotImplementedError


def rate_tsk_lst(rate_tsk_str):
    """ Set sequence of tasks for rate
    """
    raise NotImplementedError


# def prnt_tsk_lst(prnt_tsk_str, thy_dct):
#     """ Set the sequence of electronic structure tasks for a given
#         species or PESs
#     """
#
#     # Set the task list using either the given models or supplied list
#     tsk_lst = prnt_tsks_from_lst(prnt_tsk_str)
#
#     # Ensure that all the tasks are in the supported tasks
#     check_prnt_tsks_supported(tsk_lst, thy_dct)
#
#     return mod_tsk_lst


def _tsk_lst(tsk_str):
    """ Set the sequence of electronic structure tasks for a given
        species or PESs
    """
    tsk_lst = []
    tsk_str = ioformat.remove_whitespace(tsk_str)
    for line in tsk_str.splitlines():
        tsk_line = line.split()
        if len(tsk_line) > 2:
            obj, tsk, keyword_lst = tsk_line[0], tsk_line[1], tsk_line[2:]
            keyword_dct = ioformat.ptt.keyword_dct_from_block(
                '\n'.join(keyword_lst))
        else:
            print('Task line not formatted correctly:\n{}'.format(line))
            sys.exit()
        tsk_lst.append([obj, tsk, keyword_dct])

    return tsk_lst


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


EXPAND_DCT = {
    # 'find_ts': ('sadpt_scan', 'sadpt_opt', 'sadpt_hess')  # sadpt_check
}


def _new_check_dct(tsk_lsts, tsk_key_dct, tsk_val_dct, thy_dct):
    """ Loop over all of the tasks, add default keywords and parameters
        and assesses if all the input is valid
    """

    for tsk_lst in tsk_lsts:

        # Unpack the task
        [obj, tsk, keyword_dct] = tsk_lst

        # Build the dictionary of default values for task
        default_dct = build_tsk_default(tsk, tsk_key_dct, tsk_val_dct)
        if obj not in tsk_key_dct[tsk][0]:  # correct
            print('tsk {}, not allowed for {}'.format(tsk, obj))
            print('')
            sys.exit()

        # Update the current task dct with the default
        new_key_dct = automol.util.dict_.right_update(default_dct, keyword_dct)

        # Check if the keyword values are allowed
        # need 2nd for anything that takes a string from the thy.dat file
        if check_dictionary(new_key_dct, tsk_key_dct, tsk_val_dct, 'ES_TSKS'):
            print('\n\nCHECK FAILED, QUITTING...')
            sys.exit()
        if check_thy_lvls(new_key_dct, thy_dct):
            print('\n\nCHECK FAILED, QUITTING...')
            sys.exit()
