""" Functions for setting the Electronic Structure tasks for
    an outlined procedure
"""

import sys
import ioformat
from ioformat.ptt import set_value_type
from ioformat.ptt import format_tsk_keywords
from mechlib.amech_io.parser.keywords import ES_TSK_OBJ_SUPPORTED_LST
from mechlib.amech_io.parser.keywords import ES_TSK_SUPPORTED_DCT
from mechlib.amech_io.parser.keywords import ES_TSK_KEYWORDS_SUPPORTED_DCT
from mechlib.amech_io.parser.keywords import ES_TSK_KEYWORDS_VAL_SUPPORTED_DCT
from mechlib.amech_io.parser.keywords import ES_TSK_KEYWORDS_DEFAULT_DCT
from mechlib.amech_io.parser.keywords import TRANS_TSK_SUPPORTED_DCT
from mechlib.amech_io.parser.keywords import TRANS_TSK_KEYWORDS_SUPPORTED_DCT
from mechlib.amech_io.parser.keywords import TRANS_TSK_KEYWORDS_VAL_SUPPORTED_DCT
from mechlib.amech_io.parser.keywords import TRANS_TSK_KEYWORDS_VAL_DEFAULT_DCT
from mechlib.amech_io.parser.keywords import PRNT_TSK_OBJ_SUPPORTED_LST
from mechlib.amech_io.parser.keywords import PRNT_TSK_SUPPORTED_DCT
from mechlib.amech_io.parser.keywords import PRNT_TSK_KEYWORDS_SUPPORTED_DCT
from mechlib.amech_io.parser.keywords import PRNT_TSK_KEYWORDS_DEFAULT_DCT


def es_tsk_lst(es_tsk_str, thy_dct):
    """ Take the es tsk list string from input and set the tasks
        Right now, we presume the tasks given in the file are correct
    """
    
    tsk_lst = _tsk_lst(es_tsk_str)

    # mod_tsk_lst = expand_tsks(tsk_lst)
    
    # Ensure that all the tasks are in the supported tasks
    check_es_tsks_supported(tsk_lst, thy_dct)

    # Add defaults if they are missing
    mod_tsk_lst = add_defaults_to_keyword_dct(
        tsk_lst,
        ES_TSK_KEYWORDS_DEFAULT_DCT,
        ES_TSK_KEYWORDS_SUPPORTED_DCT)

    return mod_tsk_lst


def trans_tsk_lst(trans_tsk_str):
    """ Set the sequence of electronic structure tasks for a given
        species or PESs
    """

    tsk_lst = _tsk_lst(trans_tsk_str)
    mod_tsk_lst = add_defaults_to_keyword_dct(
        tsk_lst,
        TRANS_TSK_KEYWORDS_DEFAULT_DCT,
        TRANS_TSK_KEYWORDS_SUPPORTED_DCT)

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
    mod_tsk_lst = add_defaults_to_keyword_dct(
       tsk_lst,
       PRNT_TSK_KEYWORDS_DEFAULT_DCT,
       PRNT_TSK_KEYWORDS_SUPPORTED_DCT)

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
            keyword_dct = format_tsk_keywords(keyword_lst)
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
    'init_geom': ('test1', 'test2'),
    'find_ts': ('sadpt_scan', 'sadpt_opt', 'sadpt_hess')  # sadpt_check
}


# Add defaults
def add_defaults_to_keyword_dct(tsk_lsts, default_dct, supported_dct):
    """ Add default values to a checked keyword dct
    """
    mod_tsk_lst = []
    for tsk_lst in tsk_lsts:
        [obj, tsk, keyword_dct] = tsk_lst
        for dkey, dval in default_dct.items():
            if dkey not in keyword_dct:
                if dkey in supported_dct[tsk]:
                    keyword_dct[dkey] = dval
        mod_tsk_lst.append([obj, tsk, keyword_dct])

    return mod_tsk_lst


# Check keywords
def check_es_tsks_supported(es_tsks, thy_dct):
    """ Check to see if the list of es tasks are supported by the code
    """
    for tsk_lst in es_tsks:
        # try:
        # Unpack the list
        [obj, tsk, keyword_dct] = tsk_lst

        # Check the object
        if obj not in ES_TSK_OBJ_SUPPORTED_LST:
            print('*ERROR: object requested that is not allowed')
            print('Allowed objs')
            for supp_obj in ES_TSK_OBJ_SUPPORTED_LST:
                print(supp_obj)
            sys.exit()

        # Check the task
        if tsk not in ES_TSK_SUPPORTED_DCT[obj]:
            print('*ERROR: task requested not allowed for object')
            print('Allowed objs')
            for key, val in ES_TSK_SUPPORTED_DCT.items():
                print(key)
                print(val)
            sys.exit()

        # Check the requested es level
        for key, val in keyword_dct.items():
            if key not in ES_TSK_KEYWORDS_SUPPORTED_DCT[tsk]:
                print('*ERROR: option not allowed for given task')
                print('Allowed objs')
                for option in ES_TSK_KEYWORDS_SUPPORTED_DCT[tsk]:
                    print(tsk, option)
                sys.exit()
            else:
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
                elif key == 'hessmax':
                    if not isinstance(val, int):
                        print('{} must be set to an integer'.format(key))
                elif key == 'pot_thresh':
                    print(key, val, type(val))
                    if not isinstance(val, float):
                        print('{} must be set to an float'.format(key))
                elif key == 'cnf_range':
                    if 'n' in val or 'e' in val:
                        val2 = val[1:]
                        try:
                            val2 = int(val2)
                        except ValueError:
                            print('cnf_range should n<int> or e<float>')
                    else:
                        print('cnf_range should n<int> or e<float>')
                else:
                    if val not in ES_TSK_KEYWORDS_VAL_SUPPORTED_DCT[key]:
                        print('*ERROR: key {}'.format(key),
                              'not set to allowed value.')
                        print('Allowed values:')
                        for kval in ES_TSK_KEYWORDS_VAL_SUPPORTED_DCT[key]:
                            print(kval)
                        sys.exit()

        # Check  if inplvl and runlvl match for HR and IRC scans
        if tsk in ('hr_scan', 'irc_scan'):
            if keyword_dct['runlvl'] != keyword_dct['inplvl']:
                print('runlvl and inplvl MUST MATCH for',
                      'hr_scan and irc_scan tasks')
                sys.exit()
        # except:
        #     print('*ERROR: es_tsk not formatted correctly')
        #     print(tsk_lst)
        #     sys.exit()


def check_prnt_tsks_supported(es_tsks, thy_dct):
    """ Check to see if the list of es tasks are supported by the code
    """
    for tsk_lst in es_tsks:
        # try:
        # Unpack the list
        [obj, tsk, keyword_dct] = tsk_lst

        # Check the object
        if obj not in PRNT_TSK_OBJ_SUPPORTED_LST:
            print('*ERROR: object requested that is not allowed')
            print('Allowed objs')
            for supp_obj in PRNT_TSK_OBJ_SUPPORTED_LST:
                print(supp_obj)
            sys.exit()

        # Check the task
        if tsk not in PRNT_TSK_SUPPORTED_DCT[obj]:
            print('*ERROR: task requested not allowed for object')
            print('Allowed objs')
            for key, val in PRNT_TSK_SUPPORTED_DCT.items():
                print(key)
                print(val)
            sys.exit()

        # Check the requested es level
        for key, val in keyword_dct.items():
            if key not in PRNT_TSK_KEYWORDS_SUPPORTED_DCT[tsk]:
                print('*ERROR: option not allowed for given task')
                print('Allowed objs')
                for option in PRNT_TSK_KEYWORDS_SUPPORTED_DCT[tsk]:
                    print(tsk, option)
                sys.exit()
            else:
                tlvl = ('geolvl', 'proplvl')
                if key in tlvl:
                    if val not in thy_dct:
                        print('*ERROR: tsk theory level',
                              '{} not given in theory.dat'.format(val))
                        sys.exit()
