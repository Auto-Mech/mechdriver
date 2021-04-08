""" Functions for setting the Electronic Structure tasks for
    an outlined procedure
"""

import sys
import ioformat
# from ioformat.ptt import format_tsk_keywords
from mechlib.amech_io.parser._keywrd import 


def es_tsk_lst(es_tsk_str, thy_dct):
    """ Take the es tsk list string from input and set the tasks
        Right now, we presume the tasks given in the file are correct
    """

    tsk_lst = _tsk_lst(es_tsk_str)

    mod_tsk_lst = expand_tsks(tsk_lst)

    # Ensure that all the tasks are in the supported tasks
    _new_check_dct(tsk_lst, thy_dct)

    # # Add defaults if they are missing
    # mod_tsk_lst = add_defaults_to_keyword_dct(
    #     tsk_lst,
    #     ES_TSK_KEYWORDS_DEFAULT_DCT,
    #     ES_TSK_KEYWORDS_SUPPORTED_DCT)

    return mod_tsk_lst


def trans_tsk_lst(trans_tsk_str):
    """ Set the sequence of electronic structure tasks for a given
        species or PESs
    """

    tsk_lst = _tsk_lst(trans_tsk_str)
    # mod_tsk_lst = add_defaults_to_keyword_dct(
    #     tsk_lst,
    #     TRANS_TSK_KEYWORDS_DEFAULT_DCT,
    #     TRANS_TSK_KEYWORDS_SUPPORTED_DCT)

    return mod_tsk_lst


def therm_tsk_lst(therm_tsk_str):
    """ Set sequence of tasks for therm
    """
    raise NotImplementedError


def rate_tsk_lst(rate_tsk_str):
    """ Set sequence of tasks for rate
    """
    raise NotImplementedError


def prnt_tsk_lst(prnt_tsk_str, thy_dct):
    """ Set the sequence of electronic structure tasks for a given
        species or PESs
    """

    # Set the task list using either the given models or supplied list
    tsk_lst = prnt_tsks_from_lst(prnt_tsk_str)

    # Ensure that all the tasks are in the supported tasks
    check_prnt_tsks_supported(tsk_lst, thy_dct)

    # Add defaults if they are missing
    # mod_tsk_lst = add_defaults_to_keyword_dct(
    #    tsk_lst,
    #    PRNT_TSK_KEYWORDS_DEFAULT_DCT,
    #    PRNT_TSK_KEYWORDS_SUPPORTED_DCT)

    return mod_tsk_lst


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
            # keyword_dct = format_tsk_keywords(keyword_lst)
            keyword_dct = ioformat.ptt.keyword_dct_from_block(inp_str)
        else:
            print('BAAD')
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
        expand_tsks = EXPAND_DCT.get(tsk, None)
        if expand_tsks is None:
            mod_tsks_lst2.append(tsk_lst)
        else:
            for tsk in expand_tsks:
                mod_tsks_lst2.append([obj, tsk, dct])

    return mod_tsks_lst2


EXPAND_DCT = {
    'find_ts': ('sadpt_scan', 'sadpt_opt', 'sadpt_hess')  # sadpt_check
}


def _new_check_dct(tsk_lst):
    """ Add defaults to tsks and check it
    """

    # First we add the defaults
    mod_tsk_lst = []
    for tsk_lst in tsk_lsts:

        # Unpack the task
        [obj, tsk, keyword_dct] = tsk_lst

        # Build the dictionary of default values for task
        supp_keywrds = ES_TSK_KEYWORDS_SUPPORTED_DCT[tsk]
        default_dct = dict(
            zip(keywrds, (NEW_ES_TSK_DCT[key][2] for key in keywrds)))

        # Update the current task dct with the default
        new_key_dct = automol.util.dict_.right_update(default_dct, keyword_dct)

        # Now check the defined in the dct are all supported
        assert set(new_key_dct.keys() <= supp_keywrds)

        # Check if the keyword values are allowed


# Check keywords
def check_es_tsks_supported(es_tsks, thy_dct):
    """ Check to see if the list of es tasks are supported by the code
    """
    tlvl = ('runlvl', 'inplvl',
            'var_splvl1', 'var_splvl2', 'var_scnlvl')
    if key in tlvl:
        if val not in thy_dct:
            print('*ERROR: tsk theory level',
                  '{} not given in theory.dat'.format(val))
            sys.exit()
            # elif key in ('mr_splvl', 'mr_scnlvl'):
            #     if val not in 'molpro2015':
            #         print('*ERROR: mr theory level only avail',
            #               'for molpro')
            #         sys.exit()
    elif key == 'cnf_range':
        if 'n' in val or 'e' in val:
            val2 = val[1:]
            try:
                val2 = int(val2)
            except ValueError:
                print('cnf_range should n<int> or e<float>')
        else:
            print('cnf_range should n<int> or e<float>')

    tlvl = ('geolvl', 'proplvl')
    if key in tlvl:
        if val not in thy_dct:
            print('*ERROR: tsk theory level',
                  '{} not given in theory.dat'.format(val))
            sys.exit()
