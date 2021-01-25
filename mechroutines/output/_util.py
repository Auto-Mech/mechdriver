""" util functions for output driver
"""

from mechlib.amech_io import parser
from mechlib import filesys

def freq_es_levels(print_keyword_dct):
    es_model = _default_es_levels(print_keyword_dct)
    es_model['harm'] = print_keyword_dct['proplvl']
    return es_model


def ene_es_levels(print_keyword_dct):
    es_model = _default_es_levels(print_keyword_dct)
    es_model['ene'] = print_keyword_dct['proplvl']
    return es_model


def _default_es_levels(print_keyword_dct):
    es_model = {'geo': print_keyword_dct['geolvl']}
    es_model['harm'] = print_keyword_dct['geolvl']
    es_model['ene'] = print_keyword_dct['geolvl']
    es_model['sym'] = print_keyword_dct['geolvl']
    es_model['tors'] = (
        print_keyword_dct['geolvl'],
        print_keyword_dct['geolvl'])
    es_model['vpt2'] = print_keyword_dct['geolvl']
    return es_model


def _set_conf_range(print_keyword_dct):
    cnf_range = print_keyword_dct['nconfs']
    if cnf_range == 'all':
        pass
    elif cnf_range != 'min':
        cnf_range = 'n{}'.format(cnf_range)
    else:
        cnf_range = print_keyword_dct['econfs']
        if cnf_range != 'min':
            cnf_range = 'e{}'.format(cnf_range)
    return cnf_range


def conformer_list(
        print_keyword_dct, save_prefix, run_prefix,
        spc_dct_i, thy_dct):
    # conformer range
    cnf_range = _set_conf_range(print_keyword_dct)

    # thy_info build
    thy_info = filesys.inf.get_es_info(
        print_keyword_dct['geolvl'], thy_dct)
    spc_info = filesys.inf.get_spc_info(spc_dct_i)
    thy_info = filesys.inf.modify_orb_restrict(spc_info, thy_info)

    # pf info build
    es_model = _default_es_levels(print_keyword_dct)
    pf_levels = parser.model.pf_level_info(
        es_model, thy_dct)
    pf_filesystems = filesys.models.pf_filesys(
        spc_dct_i, pf_levels, run_prefix, save_prefix, saddle=False)

    # confs
    [cnf_fs, _, _, _, _] = pf_filesystems['harm']
    cnf_locs_lst, cnf_locs_paths = filesys.mincnf.conformer_locators(
        cnf_fs, thy_info, cnf_range=cnf_range)
    return cnf_fs, cnf_locs_lst, cnf_locs_paths


def conformer_list_from_models(
        print_keyword_dct, save_prefix, run_prefix,
        spc_dct_i, thy_dct, pf_levels, pf_models):
    # conformer range
    cnf_range = _set_conf_range(print_keyword_dct)

    # thy_info build
    thy_info = pf_levels['geo'][1]
    spc_info = filesys.inf.get_spc_info(spc_dct_i)
    thy_info = filesys.inf.modify_orb_restrict(spc_info, thy_info)

    # pf info build
    pf_filesystems = filesys.models.pf_filesys(
        spc_dct_i, pf_levels, run_prefix, save_prefix, saddle=False)

    # confs
    [cnf_fs, _, _, _, _] = pf_filesystems['harm']
    cnf_locs_lst, cnf_locs_paths = filesys.mincnf.conformer_locators(
        cnf_fs, thy_info, cnf_range=cnf_range)
    return cnf_fs, cnf_locs_lst, cnf_locs_paths
