""" Functions for setting the Electronic Structure tasks for
    an outlined procedure
"""

from lib.filesystem import inf as finf
from lib.load import ptt
from lib.load.keywords import ES_TSK_SUPPORTED_DCT
from lib.load.keywords import ES_TSK_OPTIONS_SUPPORTED_DCT


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
    assert check_es_tsks_supported(tsk_lst)

    return tsk_lst


def es_tsks_from_lst(es_tsks_str):
    """ Take the es tsk list string from input and set the tasks
        Right now, we presume the tasks given in the file are correct
    """
    # Split the string into different strings of keywords
    tsk_lst = []
    for line in es_tsks_str.splitlines():
        # Put a cleaner somehwere to get rid of blank lines
        tsk_line = line.strip().split()
        if len(tsk_line) == 4:
            [obj, tsk, run_lvl, in_lvl] = tsk_line
            formtd_options = []
        elif len(tsk_line) == 5:
            [obj, tsk, run_lvl, in_lvl, options] = tsk_line
            formtd_options = format_tsk_options(options)
        else:
            print('BAAD')
        tsk_lst.append([obj, tsk, run_lvl, in_lvl, formtd_options])

    return tsk_lst


def format_tsk_options(options_str):
    """ format options string
    """
    options_str = options_str.replace('[', '').replace(']', '')
    options_lst = [option.strip() for option in options_str.split(',')]
    return options_lst


def es_tsks_from_models(rxn_model_dct,
                        saddle=False, wells=False, barrierless=False):
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


def check_es_tsks_supported(es_tsks):
    """ Check to see if the list of es tasks are supported by the code
    """
    obj_good, tsk_good, opt_good = True, True, True
    for tsk_lst in es_tsks:
        [obj, tsk, _, _, options] = tsk_lst
        print(tsk_lst)
        if obj in ES_TSK_SUPPORTED_DCT:
            if tsk in ES_TSK_SUPPORTED_DCT[obj]:
                print(ES_TSK_OPTIONS_SUPPORTED_DCT[tsk])
                print(options)
                chk = all(option in ES_TSK_OPTIONS_SUPPORTED_DCT[tsk]
                          for option in options)
                if not chk:
                    print('opt not good')
                    opt_good = False
            else:
                print('tsk not good')
                tsk_good = False
        else:
            print('obj not good')
            obj_good = False

    return bool(obj_good and tsk_good and opt_good)
