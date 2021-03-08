""" util functions for output driver
"""

import pandas

from mechlib.amech_io import parser
from mechlib import filesys
from autofile import io_ as io
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo


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
    """ Create a list of conformers based on the species name
        and run.dat geolvl/proplvl
    """
    # conformer range
    cnf_range = _set_conf_range(print_keyword_dct)

    # thy_info build
    thy_info = tinfo.from_dct(thy_dct.get(print_keyword_dct.get('geolvl')))
    spc_info = sinfo.from_dct(spc_dct_i)
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)

    # pf info build
    es_model = _default_es_levels(print_keyword_dct)
    pf_levels = parser.model.pf_level_info(
        es_model, thy_dct)
    pf_filesystems = filesys.models.pf_filesys(
        spc_dct_i, pf_levels, run_prefix, save_prefix, saddle=False)

    # confs
    [cnf_fs, _, _, _, _] = pf_filesystems['harm']
    if cnf_fs is not None:
        cnf_locs_lst, cnf_locs_paths = filesys.mincnf.conformer_locators(
            cnf_fs, mod_thy_info, cnf_range=cnf_range)
    else:
        cnf_locs_lst, cnf_locs_paths = [], []
    return cnf_fs, cnf_locs_lst, cnf_locs_paths


def conformer_list_from_models(
        print_keyword_dct, save_prefix, run_prefix,
        spc_dct_i, thy_dct, pf_levels, pf_models):
    """ Create a list of conformers based on the species name
        and model.dat info
    """
    # conformer range
    cnf_range = _set_conf_range(print_keyword_dct)

    # thy_info build
    thy_info = pf_levels['geo'][1]
    spc_info = sinfo.from_dct(spc_dct_i)
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)

    # pf info build
    pf_filesystems = filesys.models.pf_filesys(
        spc_dct_i, pf_levels, run_prefix, save_prefix, saddle=False)

    # confs
    [cnf_fs, _, _, _, _] = pf_filesystems['harm']
    cnf_locs_lst, cnf_locs_paths = filesys.mincnf.conformer_locators(
        cnf_fs, mod_thy_info, cnf_range=cnf_range)
    return cnf_fs, cnf_locs_lst, cnf_locs_paths


def set_csv_data(tsk):
    """ some tasks have nested dictionaries, prep for that
    """
    csv_data = {}
    if 'freq' in tsk:
        csv_data['freq'] = {}
        csv_data['tfreq'] = {}
        csv_data['allfreq'] = {}
        csv_data['scalefactor'] = {}

    return csv_data


def write_csv_data(tsk, csv_data, filelabel, spc_array):
    """ Write the csv data dictionary into the correct type of csv 
        or text file
    """
    if 'freq' in tsk:
        final_csv_data = {}
        for key in csv_data['freq']:
            final_csv_data[key] = csv_data['freq'][key]
        if csv_data['scalefactor']:
            final_csv_data['Scale Factor'] = []
            for key in csv_data['scalefactor']:
                final_csv_data[key + '_scalefactor'] = csv_data['scalefactor'][key]
            final_csv_data['Torsional Frequencies'] = []
            for key in csv_data['tfreq']:
                final_csv_data[key + '_tfreq'] = csv_data['tfreq'][key]
            final_csv_data['All RT Harmonic Frequencies'] = []
            for key in csv_data['allfreq']:
                final_csv_data[key + '_RTFreq'] = csv_data['allfreq'][key]
        print(final_csv_data)
        ncols = max([len(x) for x in final_csv_data.values()])
        df = pandas.DataFrame.from_dict(
            final_csv_data, orient='index',
            columns=['Path', 'ZPVE [A.U.]', *[''] * (ncols-2)])
        df.to_csv(filelabel, float_format='%.5f')
    if 'geo' in tsk:
        all_data = '\n'.join([spc_data for spc_data in csv_data.values()])
        io.write_file(filelabel, all_data)
    if 'zma' in tsk:
        all_data = '\n'.join([spc_data for spc_data in csv_data.values()])
        io.write_file(filelabel, all_data)
    if 'ene' in tsk:
        df = pandas.DataFrame.from_dict(
            csv_data, orient='index',
            columns=['Path', 'Energy [A.U.]'])
        df.to_csv(filelabel, float_format='%.8f')
    if 'enthalpy' in tsk:
        df = pandas.DataFrame.from_dict(
            csv_data, orient='index',
            columns=[
                'Path', 'ZPVE+Energy [A.U.]', 'Hf (0 K) [kcal/mol]',
                *spc_array])
        df.to_csv(filelabel, float_format='%.6f')
    if 'coeffs' in tsk:
        df = pandas.DataFrame.from_dict(
            csv_data, orient='index',
            columns=[
                *spc_array])
        df.to_csv(filelabel, float_format='%.2f')
