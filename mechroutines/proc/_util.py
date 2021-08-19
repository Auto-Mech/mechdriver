""" util functions for output driver
"""

import pandas
import ioformat
import automol
from autofile import io_ as io
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo
from mechlib import filesys
from mechlib.amech_io import printer as ioprinter


def freq_es_levels(print_keyword_dct):
    """ ?
    """
    es_levels = _default_es_levels(print_keyword_dct)
    es_levels['vib'] = print_keyword_dct['proplvl']
    return es_levels


def ene_es_levels(print_keyword_dct):
    """ ?
    """
    es_levels = _default_es_levels(print_keyword_dct)
    es_levels['ene'] = print_keyword_dct['proplvl']
    return es_levels


def generate_spc_model_dct(es_levels, thy_dct):
    """ Make a specie smodel dct to pass to pf_filesystems function
    """
    spc_model_dct_i = {}
    for prop in es_levels:
        spc_model_dct_i[prop] = {
            'geolvl': (
                es_levels['geo'],
                (1.00, tinfo.from_dct(thy_dct.get(es_levels['geo']))))}
        if prop == 'vib':
            spc_model_dct_i[prop]['mod'] = 'harm'
        elif prop == 'ene':
            spc_model_dct_i[prop]['lvl1'] = (
                es_levels['ene'],
                (1.00, tinfo.from_dct(thy_dct.get(es_levels['ene']))))
    return spc_model_dct_i


def _default_es_levels(print_keyword_dct):
    """ ?
    """
    es_model = {'geo': print_keyword_dct['geolvl']}
    es_model['vib'] = print_keyword_dct['geolvl']
    es_model['ene'] = print_keyword_dct['geolvl']
    es_model['sym'] = print_keyword_dct['geolvl']
    es_model['tors'] = (
        print_keyword_dct['geolvl'],
        print_keyword_dct['geolvl'])
    es_model['vpt2'] = print_keyword_dct['geolvl']
    return es_model

# ‘geolvl’: (‘lvl_wbt’, (1.00, thy_inf))
# es_model = {
#   ‘ene’: {
# ‘geolvl’: (‘geolvl’, parser.models.format_lvl(print_keyword_dct[‘geolvl’]))
#        }
#           ‘vib’: {
#               ‘geolvl’: print_keyword_dct[‘geolvl’]
#           }


def _set_sort_info_lst(sort_str, thy_dct, spc_info):
    """
    """
    sort_lvls = [None, None]
    sort_typ_lst = ['zpe', 'sp']
    if sort_str is not None:
        for sort_param in sort_str.split(','):
            idx = None
            for typ_idx, typ_str in enumerate(sort_typ_lst):
                if typ_str in sort_param:
                    lvl_key = sort_str.split(typ_str + '(')[1].split(')')[0]
                    idx = typ_idx
            if idx is not None:
                method_dct = thy_dct.get(lvl_key)
                if method_dct is None:
                    print(
                        'no {} in theory.dat, not using {} in sorting'.format(
                            lvl_key, sort_typ_lst[idx]))
                    continue
                thy_info = tinfo.from_dct(method_dct)
                mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
                sort_lvls[idx] = mod_thy_info
    return sort_lvls


def _set_conf_range(print_keyword_dct):
    """ ?
    """
    cnf_range = print_keyword_dct['cnf_range']
    #if cnf_range == 'all':
    #    pass
    #elif cnf_range != 'min':
    #    cnf_range = 'n{}'.format(cnf_range)
    #else:
    #    cnf_range = print_keyword_dct['econfs']
    #    if cnf_range != 'min':
    #        cnf_range = 'e{}'.format(cnf_range)
    return cnf_range


def conformer_list(
        spc_name, print_keyword_dct, save_prefix, run_prefix,
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
    sort_info_lst = _set_sort_info_lst(
        print_keyword_dct['sort'], thy_dct, spc_info)

    zrxn = spc_dct_i.get('zrxn', None)
    _root = filesys.root_locs(
        spc_dct_i,
        saddle=(zrxn is not None),
        name=spc_name)
    _, cnf_save_fs = filesys.build_fs(
        run_prefix, save_prefix, 'CONFORMER',
        thy_locs=mod_thy_info[1:],
        **_root)
    rng_cnf_locs_lst, rng_cnf_locs_path = filesys.mincnf.conformer_locators(
        cnf_save_fs, mod_thy_info, cnf_range=cnf_range, sort_info_lst=sort_info_lst)
    return cnf_save_fs, rng_cnf_locs_lst, rng_cnf_locs_path, mod_thy_info


def conformer_list_from_models(
        spc_name, print_keyword_dct, save_prefix, run_prefix,
        spc_dct_i, spc_mod_dct_i):
    """ Create a list of conformers based on the species name
        and model.dat info
    """
    # conformer range
    cnf_range = _set_conf_range(print_keyword_dct)

    # thy_info build
    thy_info = spc_mod_dct_i['vib']['geolvl'][1][1]
    spc_info = sinfo.from_dct(spc_dct_i)
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
    sort_info_lst = _set_sort_info_lst(
        print_keyword_dct['sort'], thy_dct, spc_info)

    zrxn = spc_dct_i.get('zrxn', None)
    _root = filesys.root_locs(
        spc_dct_i, name=spc_name, saddle=(zrxn is not None))
    _, cnf_save_fs = filesys.build_fs(
        run_prefix, save_prefix, 'CONFORMER',
        thy_locs=mod_thy_info[1:],
        **_root)
    rng_cnf_locs_lst, rng_cnf_locs_path = filesys.mincnf.conformer_locators(
        cnf_save_fs, mod_thy_info, cnf_range=cnf_range, sort_info_lst=sort_info_lst)
    return cnf_save_fs, rng_cnf_locs_lst, rng_cnf_locs_path, mod_thy_info


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
        fin_csv_data = {}
        for key in csv_data['freq']:
            fin_csv_data[key] = csv_data['freq'][key]
        if csv_data['scalefactor']:
            fin_csv_data['Scale Factor'] = []
            for key in csv_data['scalefactor']:
                fin_csv_data[key+'_scalefactor'] = csv_data['scalefactor'][key]
            fin_csv_data['Torsional Frequencies'] = []
            for key in csv_data['tfreq']:
                fin_csv_data[key+'_tfreq'] = csv_data['tfreq'][key]
            fin_csv_data['All RT Harmonic Frequencies'] = []
            for key in csv_data['allfreq']:
                fin_csv_data[key+'_RTFreq'] = csv_data['allfreq'][key]
        # print(fin_csv_data)
        ncols = max([len(x) for x in fin_csv_data.values()])
        dframe = pandas.DataFrame.from_dict(
            fin_csv_data, orient='index',
            columns=['Path', 'ZPVE [A.U.]', *[''] * (ncols-2)])
        dframe.to_csv(filelabel, float_format='%.5f')
    if 'geo' in tsk:
        all_data = '\n'.join(spc_data for spc_data in csv_data.values())
        io.write_file(filelabel, all_data)
    if 'molden' in tsk:
        all_data = '\n'.join(spc_data for spc_data in csv_data.values())
        io.write_file(filelabel, all_data)
    if 'zma' in tsk:
        for spc_data in csv_data.values():
            print(spc_data)
        all_data = '\n'.join(spc_data for spc_data in csv_data.values())
        io.write_file(filelabel, all_data)
    if 'ene' in tsk:
        dframe = pandas.DataFrame.from_dict(
            csv_data, orient='index',
            columns=['Path', 'Energy [A.U.]'])
        dframe.to_csv(filelabel, float_format='%.8f')
    if 'enthalpy' in tsk:
        dframe = pandas.DataFrame.from_dict(
            csv_data, orient='index',
            columns=[
                'Path', 'ZPVE+Energy [A.U.]', 'Hf (0 K) [kcal/mol]',
                *spc_array])
        dframe.to_csv(filelabel, float_format='%.6f')
    if 'coeffs' in tsk:
        dframe = pandas.DataFrame.from_dict(
            csv_data, orient='index',
            columns=[
                *spc_array])
        dframe.to_csv(filelabel, float_format='%.2f')


def write_missing_data_report(miss_data):
    """ Write all of the data collating the missing data
    """

    miss_data = automol.util.remove_duplicates_with_order(miss_data)

    print('\n\n\nMissing Data Requested by the User')
    print('{0:<20s}{1:<12s}{2:<12s}{3:<12s}'.format(
        'Name', 'Method', 'Basis', 'Data'))
    for data in miss_data:
        name, thy_info, dat = data
        method, basis = thy_info[1:3]
        print('{0:<20s}{1:<12s}{2:<12s}{3:<12s}'.format(
            name, method, basis, dat))


def write_data_dirs(data_dirs, prefix):
    """ write a series of data files for similar typ in a directory
    """

    for data_dir in data_dirs:
        name, string_dct = data_dir
        name = name.replace('(', '_').replace(')', '_')
        if string_dct:
            dir_name = ioformat.pathtools.prepare_path(
                [prefix, name], make=True)
            print('Writing files into {} directory'.format(dir_name))
            for key, string in string_dct.items():
                fname = '{}.dat'.format(key)
                print('Writing file {}'.format(fname))
                ioformat.pathtools.write_file(string, dir_name, fname)


def get_file_label(tsk, model_dct, proc_keyword_dct, spc_mod_dct_i):
    """ what is the name and extension for this processed file?
    """
    thylabel = ''  # used for data directory labels (e.g. displacements)
    if 'coeffs' in tsk:
        filelabel = 'coeffs'
        filelabel += '_{}'.format(model_dct['therm_fit']['ref_scheme'])
        filelabel += '.csv'
    elif 'freq' in tsk:
        filelabel = 'freq'
        geolvl = proc_keyword_dct.get('geolvl')
        if geolvl is not None:
            thylabel = proc_keyword_dct['geolvl']
            filelabel += '_{}'.format(proc_keyword_dct['geolvl'])
        else:
            thylabel = spc_mod_dct_i['vib']['geolvl'][0]
            filelabel += '_m{}'.format(spc_mod_dct_i['vib']['geolvl'][0])
        filelabel += '.csv'
    elif 'geo' in tsk:
        filelabel = 'molden'
        geolvl = proc_keyword_dct.get('geolvl')
        if geolvl is not None:
            filelabel += '_{}'.format(proc_keyword_dct['geolvl'])
        else:
            filelabel += '_m{}'.format(spc_mod_dct_i['vib']['geolvl'][0])
        filelabel += '.txt'
    elif 'molden' in tsk:
        filelabel = 'geo'
        geolvl = proc_keyword_dct.get('geolvl')
        if geolvl is not None:
            filelabel += '_{}'.format(proc_keyword_dct['geolvl'])
        else:
            filelabel += '_m{}'.format(spc_mod_dct_i['vib']['geolvl'][0])
        filelabel += '.txt'
    elif 'zma' in tsk:
        filelabel = 'zmat'
        geolvl = proc_keyword_dct.get('geolvl')
        if geolvl is not None:
            filelabel += '_{}'.format(proc_keyword_dct['geolvl'])
        if spc_mod_dct_i:
            filelabel += '_m{}'.format(spc_mod_dct_i['vib']['geolvl'][0])
        filelabel += '.txt'
    elif 'ene' in tsk:
        filelabel = 'ene'
        geolvl = proc_keyword_dct.get('geolvl')
        if geolvl is not None:
            filelabel += '_{}'.format(proc_keyword_dct['geolvl'])
            filelabel += '_{}'.format(proc_keyword_dct['proplvl'])
        else:
            filelabel += '_m{}'.format(spc_mod_dct_i['vib']['geolvl'][0])
            filelabel += '_{}'.format(spc_mod_dct_i['ene']['lvl1'][0])
        filelabel += '.csv'
    elif 'enthalpy' in tsk:
        filelabel = 'enthalpy'
        geolvl = proc_keyword_dct.get('geolvl')
        if geolvl is not None:
            filelabel += '_{}'.format(proc_keyword_dct['geolvl'])
            filelabel += '_{}'.format(proc_keyword_dct['proplvl'])
        else:
            filelabel += '_m{}'.format(spc_mod_dct_i['vib']['geolvl'][0])
            filelabel += '_{}'.format(spc_mod_dct_i['ene']['lvl1'][0])
        filelabel += '.csv'
    elif 'torsions' in tsk:
        filelabel = None
    return filelabel, thylabel


# def choose_theory(proc_keyword_dct, spc_mod_dct_i, thy_dct):
def choose_theory(proc_keyword_dct, spc_mod_dct_i):
    """ choose between theories set in models.dat and in run.dat
    """
    if proc_keyword_dct['geolvl']:
        # thy_info = tinfo.from_dct(thy_dct.get(
        #    proc_keyword_dct['geolvl']))
        spc_mod_dct_i = None
    # else:
    #    thy_info = spc_mod_dct_i['vib']['geolvl'][1][1]
    # return thy_info, spc_mod_dct_i
    return spc_mod_dct_i


def choose_conformers(
        spc_name, proc_keyword_dct, spc_mod_dct_i,
        save_prefix, run_prefix, spc_dct_i, thy_dct):
    """ get the locations (locs and paths) for the number of conformers
        set in the proc_keyword_dct and in the theory directory specified
        by either the same dct or by models.ddat
    """
    if spc_mod_dct_i is not None:
        ret = conformer_list_from_models(
            spc_name, proc_keyword_dct, save_prefix, run_prefix,
            spc_dct_i, spc_mod_dct_i)
        cnf_fs, rng_cnf_locs_lst, rng_cnf_locs_path, mod_thy_info = ret
    else:
        cnf_inf = conformer_list(
            spc_name, proc_keyword_dct, save_prefix, run_prefix,
            spc_dct_i, thy_dct)
        cnf_fs, rng_cnf_locs_lst, rng_cnf_locs_path, mod_thy_info = cnf_inf

    return cnf_fs, rng_cnf_locs_lst, rng_cnf_locs_path, mod_thy_info


# Queue manipulation functions
def remove_unstable(spc_queue, spc_dct, thy_dct, spc_mod_dct_i,
                    proc_key_dct, save_prefix):
    """ For each species in the queue see if there are files
        in the save filesystem denoting they are unstable. If so,
        that species is removed from the queue for collection tasks.
    """
    if spc_mod_dct_i is not None:
        thy_info = spc_mod_dct_i['vib']['geolvl'][1][1]
    else:
        thy_info = tinfo.from_dct(thy_dct.get(
            proc_key_dct['proplvl']))

    stable_queue = ()
    for spc_name in spc_queue:
        if 'ts_' in spc_name:
            stable_queue += (spc_name,)
        else:
            instab, path = filesys.read.instability_transformation(
                spc_dct, spc_name, thy_info, save_prefix)
            if instab is None:
                stable_queue += (spc_name,)
            else:
                ioprinter.info_message(
                    'Found instability file at path {}'.format(path),
                    newline=1)
                ioprinter.info_message(
                    'Removing {} from queue'.format(spc_name))

    return stable_queue


def remove_ts_missing(obj_queue, spc_dct):
    """ Remove TSs who have no information from the queue and
        include them in the missing data lists
    """

    new_queue = ()
    ts_miss_data = ()
    for obj in obj_queue:
        if 'ts_' in obj:
            ts_dct = spc_dct[obj]
            miss = ts_dct.get('missdata')
            if miss is not None:
                ts_miss_data += ((obj, ts_dct['missdata'], 'geometry'),)
            else:
                new_queue += (obj,)
        else:
            new_queue += (obj,)

    return new_queue, ts_miss_data
