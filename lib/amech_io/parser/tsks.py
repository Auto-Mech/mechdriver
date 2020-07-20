""" Functions for setting the Electronic Structure tasks for
    an outlined procedure
"""

import sys
import ioformat
from lib.amech_io.parser.keywords import ES_TSK_OBJ_SUPPORTED_LST
from lib.amech_io.parser.keywords import ES_TSK_SUPPORTED_DCT
from lib.amech_io.parser.keywords import ES_TSK_KEYWORDS_SUPPORTED_DCT
from lib.amech_io.parser.keywords import ES_TSK_KEYWORDS_VAL_SUPPORTED_DCT
from lib.amech_io.parser.keywords import ES_TSK_KEYWORDS_DEFAULT_DCT


def es_tsk_lst(es_tsk_str, rxn_model_dct, thy_dct, saddle=False):
    """ Set the sequence of electronic structure tasks for a given
        species or PESs
    """

    # Set the task list using either the given models or supplied list
    if 'models' in es_tsk_str:
        tsk_lst = es_tsks_from_models(rxn_model_dct, saddle=saddle)
    else:
        tsk_lst = es_tsks_from_lst(es_tsk_str)

    # Ensure that all the tasks are in the supported tasks
    check_es_tsks_supported(tsk_lst, thy_dct)

    # Add defaults if they are missing
    mod_tsk_lst = add_defaults_to_tsk_keyword_dct(tsk_lst)

    return mod_tsk_lst


def es_tsks_from_lst(es_tsks_str):
    """ Take the es tsk list string from input and set the tasks
        Right now, we presume the tasks given in the file are correct
    """
    # Split the string into different strings of keywords
    tsk_lst = []
    es_tsks_str = ioformat.remove_whitespace(es_tsks_str)
    for line in es_tsks_str.splitlines():
        tsk_line = line.split()
        if len(tsk_line) > 2:
            obj, tsk, keyword_lst = tsk_line[0], tsk_line[1], tsk_line[2:]
            keyword_dct = format_tsk_keywords(keyword_lst)
        else:
            print('BAAD')
        tsk_lst.append([obj, tsk, keyword_dct])

    return tsk_lst


def format_tsk_keywords(keyword_lst):
    """ format keywords string
    """
    # Build initial keyword dct
    keyword_dct = {}
    for keyword in keyword_lst:
        [key, val] = keyword.split('=')
        keyword_dct[key] = format_val(val)

    return keyword_dct


def format_val(val):
    """ format val (probably switch to one in ptt later)
    """
    if val == 'True':
        formtd_val = True
    elif val == 'False':
        formtd_val = False
    elif val.isdigit():
        formtd_val = int(val)
    else:
        formtd_val = val

    return formtd_val


def add_defaults_to_tsk_keyword_dct(tsk_lsts):
    """ Add default values to a checked keyword dct
    """
    mod_tsk_lst = []
    for tsk_lst in tsk_lsts:
        [obj, tsk, keyword_dct] = tsk_lst
        for dkey, dval in ES_TSK_KEYWORDS_DEFAULT_DCT.items():
            if dkey not in keyword_dct:
                if dkey in ES_TSK_KEYWORDS_SUPPORTED_DCT[tsk]:
                    keyword_dct[dkey] = dval
        mod_tsk_lst.append([obj, tsk, keyword_dct])

    return mod_tsk_lst


def es_tsks_from_models(rxn_model_dct,
                        saddle=False, wells=False):
    """ Set a list of tasks using the model dictionary for
        setting up the electronic structure calculations
    """

    # Initialize task list
    tsk_lst = []

    # Tasks for getting the geoms for minima, ts, and wells
    tsk_lst.append('find_min')
    tsk_lst.append('min_samp')
    if saddle:
        tsk_lst.append('find_ts')
        tsk_lst.append('ts_samp')
    if wells:
        tsk_lst.append('find_wells')
        tsk_lst.append('well_samp')

    # Add hindered rotor calculations if requested
    if rxn_model_dct['tors'] == '1dhr':
        tsk_lst.append('hr_scan')

    # Calculations using the geometries
    tsk_lst.append('conf_hess')
    # if rxn_model_dct['vib'] == 'vpt2':
    #     tsk_lst.append('conf_vpt2')
    tsk_lst.append('conf_energy')

    # Add in the tasks for running ircs
    # if rxn_model_dct['ts_sadpt'] == 'vtst' and not barrierless:
    #    tsk_lst.append('run_irc')

    return tsk_lst


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
                if key in ('runlvl', 'inplvl', 'var_splvl1', 'var_splvl2, var_scn_lvl'):
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
                    print('check', key, val)
                    if not isinstance(val, int):
                        print('{} must be set to an integer'.format(key))
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
