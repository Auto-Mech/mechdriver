""" Functions for setting the Electronic Structure tasks for
    an outlined procedure
"""

from lib.filesystem import inf as finf
from lib.load import ptt
from lib.load.keywords import ES_TSK_SUPPORTED_LST


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
    # :assert check_es_tsks_supported(tsk_lst)
    mod_tsk_lst = []
    for tsk_info in tsk_lst:
        [tsk, es_run_key, es_ini_key, overwrite] = tsk_info
        ini_thy_info = finf.get_es_info(es_ini_key, thy_dct)
        thy_info = finf.get_es_info(es_run_key, thy_dct)
        mod_tsk_lst.append([tsk, thy_info, ini_thy_info, overwrite])

    return mod_tsk_lst


def es_tsks_from_lst(es_tsks_str):
    """ Take the es tsk list string from input and set the tasks
        Right now, we presume the tasks given in the file are correct
    """
    # Split the string into different strings of keywords
    tsk_lst = []
    for line in es_tsks_str.splitlines():
        # Put a cleaner somehwere to get rid of blank lines
        tmp = line.strip()
        if tmp:
            pvals = [string.strip() for string in tmp.split('=')]
            keyword, value_lst = ptt.format_param_vals(pvals)
            tsk_lst.append(format_tsk_lst(keyword, value_lst))

    # Hard to get the full tsk_lst for both function
    # tsk inp_geo out_info overwrite
    # Take the overwrite opt keyword plus el method from model

    return tsk_lst


def format_tsk_lst(keyword, value_lst):
    """ format keyword = [vals] string to the tsk_lst keywords
    """
    frmtd_keyword = keyword.split('tsk_')[1].lower()
    # setvals = [value.strip() for value in value_lst]
    tsk_lst = [frmtd_keyword, value_lst[0], value_lst[1], value_lst[2]]
    return tsk_lst


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
    return all(tsk in ES_TSK_SUPPORTED_LST for tsk in es_tsks)
