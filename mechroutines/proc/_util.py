""" util functions for output driver
"""

import os
import pandas
import json
import ioformat
import automol
import autofile
from autofile import io_ as io
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo
from mechanalyzer.inf import rxn as rinfo
from mechlib import filesys
from mechlib.amech_io import printer as ioprinter


def freq_es_levels(proc_keyword_dct):
    """ ?
    """
    es_levels = _default_es_levels(proc_keyword_dct)
    es_levels['vib'] = proc_keyword_dct['proplvl']
    return es_levels


def ene_es_levels(proc_keyword_dct):
    """ ?
    """
    es_levels = _default_es_levels(proc_keyword_dct)
    es_levels['ene'] = proc_keyword_dct['proplvl']
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


def _default_es_levels(proc_keyword_dct):
    """ ?
    """
    es_model = {'geo': proc_keyword_dct['geolvl']}
    es_model['vib'] = proc_keyword_dct['geolvl']
    es_model['ene'] = proc_keyword_dct['geolvl']
    es_model['sym'] = proc_keyword_dct['geolvl']
    es_model['tors'] = (
        proc_keyword_dct['geolvl'],
        proc_keyword_dct['geolvl'])
    es_model['vpt2'] = proc_keyword_dct['geolvl']
    return es_model

# ‘geolvl’: (‘lvl_wbt’, (1.00, thy_inf))
# es_model = {
#   ‘ene’: {
# ‘geolvl’: (‘geolvl’, parser.models.format_lvl(proc_keyword_dct[‘geolvl’]))
#        }
#           ‘vib’: {
#               ‘geolvl’: proc_keyword_dct[‘geolvl’]
#           }


def _set_sort_info_lst(sort_str, thy_dct, spc_info):
    """ Return the levels to sort conformers by if zpve or sp
        levels were assigned in input

        if we ask for zpe(lvl_wbs),sp(lvl_b2t),gibbs(700)
        out sort_info_lst will be [('gaussian', 'wb97xd', '6-31*', 'RU'),
        ('gaussian', 'b2plypd3', 'cc-pvtz', 'RU'), None, None, 700.]
    """
    sort_lvls = [None, None, None, None, None]
    sort_typ_lst = ['freqs', 'sp', 'enthalpy', 'entropy', 'gibbs']
    if sort_str is not None:
        for sort_param in sort_str.split(','):
            idx = None
            for typ_idx, typ_str in enumerate(sort_typ_lst):
                if typ_str in sort_param:
                    lvl_key = sort_str.split(typ_str + '(')[1].split(')')[0]
                    idx = typ_idx
            if idx is not None:
                if idx < 2:
                    method_dct = thy_dct.get(lvl_key)
                    if method_dct is None:
                        ioprinter.warning_message(
                            f'no {lvl_key} in theory.dat, '
                            f'not using {sort_typ_lst[idx]} in sorting')
                        continue
                    thy_info = tinfo.from_dct(method_dct)
                    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
                    sort_lvls[idx] = mod_thy_info
                else:
                    sort_lvls[idx] = float(lvl_key)
    return sort_lvls


def choose_conformers(
        spc_name, proc_keyword_dct, spc_mod_dct_i,
        save_prefix, run_prefix, spc_dct, thy_dct):
    """ Create a list of conformers based on the species name
        and model.dat info
    """
    # conformer range
    spc_dct_i = spc_dct[spc_name]
    symm_dct = {}
    populate_symm = False
    cnf_range = proc_keyword_dct['cnf_range']
    if '_includesym' in cnf_range:
        cnf_range = cnf_range.split('_includesym')[0]
        populate_symm = True
    hbond_cutoffs = spc_dct_i['hbond_cutoffs']

    # thy_info build
    zrxn = spc_dct_i.get('zrxn', None)
    thy_info = spc_mod_dct_i['vib']['geolvl'][1][1]
    if zrxn is None:
        spc_info = sinfo.from_dct(spc_dct_i, canonical=True)
    else:
        spc_info = rinfo.ts_info(rinfo.from_dct(
            spc_dct_i['reacs'], spc_dct_i['prods'], spc_dct))
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
    sort_info_lst = _set_sort_info_lst(
        proc_keyword_dct['sort'], thy_dct, spc_info)

    _root = filesys.root_locs(
        spc_dct_i, name=spc_name, saddle=(zrxn is not None))
    _, cnf_save_fs = filesys.build_fs(
        run_prefix, save_prefix, 'CONFORMER',
        thy_locs=mod_thy_info[1:],
        **_root)
    rng_cnf_locs_lst, rng_cnf_locs_path = filesys.mincnf.conformer_locators(
        cnf_save_fs, mod_thy_info,
        cnf_range=cnf_range, sort_info_lst=sort_info_lst, hbond_cutoffs=hbond_cutoffs,
        print_enes=True, nprocs=proc_keyword_dct['nprocs'])
    if populate_symm:
        for cnf_locs, cnf_path in zip(rng_cnf_locs_lst, rng_cnf_locs_path):
            symm_fs = autofile.fs.symmetry(cnf_path)
            symm_locs_lst = symm_fs[-1].existing(ignore_bad_formats=True)
            if symm_locs_lst:
                symm_dct[':'.join(cnf_locs)] = symm_locs_lst
    cnf_info = (cnf_save_fs, rng_cnf_locs_lst, rng_cnf_locs_path,)
    return cnf_info, mod_thy_info, symm_dct


def set_csv_data(tsk):
    """ some tasks have nested dictionaries, prep for that
    """
    csv_data = {}
    if 'freq' in tsk:
        csv_data['freq'] = {}
        csv_data['tfreq'] = {}
        csv_data['allfreq'] = {}
        csv_data['scalefactor'] = {}
    elif 'weight' in tsk:
        csv_data['hf_array'] = []
        csv_data['pf_array'] = []
        csv_data['locs_lst'] = []
    elif 'gibbs' in tsk:
        csv_data['hf_array'] = []
        csv_data['locs_lst'] = []

    return csv_data


def write_csv_data(tsk, csv_data, filelabel, col_array, prefix):
    """ Write the csv data dictionary into the correct type of csv
        or text file
    """
    if filelabel is not None:
        filelabel = os.path.join(prefix, filelabel)
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
    if 'sidata' in tsk:
        all_data = '\n'.join(spc_data for spc_data in csv_data.values())
        io.write_file(filelabel, all_data)
    if 'json' in tsk:
        all_data = {}
        for label in csv_data:
            for ich in csv_data[label]:
                for cid in csv_data[label][ich]:
                    if ich not in all_data:
                        all_data[ich] = csv_data[label][ich]
                    else:
                        all_data[ich][cid] = csv_data[label][ich][cid]
        io.write_file(filelabel, json.dumps(all_data))
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
    if 'weight' in tsk:
        dframe = pandas.DataFrame.from_dict(
            csv_data, orient='index',
            columns=['weight'])
        dframe.to_csv(filelabel, float_format='%.8f')
    if 'gibbs' in tsk:
        dframe = pandas.DataFrame.from_dict(
            csv_data, orient='index',
            columns=['rel gibbs'])
        dframe.to_csv(filelabel, float_format='%.8f')
    if 'enthalpy' in tsk:
        dframe = pandas.DataFrame.from_dict(
            csv_data, orient='index',
            columns=[
                'Path', 'ZPVE+Energy [A.U.]', 'Hf (0 K) [kcal/mol]',
                *col_array])
        dframe.to_csv(filelabel, float_format='%.6f')
    if 'coeffs' in tsk:
        dframe = pandas.DataFrame.from_dict(
            csv_data, orient='index',
            columns=[
                *col_array])
        dframe.to_csv(filelabel, float_format='%.2f')
    if 'date' in tsk:
        dframe = pandas.DataFrame.from_dict(
            csv_data, orient='index',
            columns=[
                *col_array])
        dframe.to_csv(filelabel, float_format='%.2f')
    elif 'messpf_inp' in tsk:
        all_data = '\n'.join(spc_data for spc_data in csv_data.values())
        io.write_file(filelabel, all_data)

    elif 'pf' in tsk:
        dframe = pandas.DataFrame.from_dict(
            csv_data, orient='index',
            columns=[
                *col_array])
        dframe.to_csv(filelabel, float_format='%.2f')


def write_missing_data_report(miss_data, spc_dct):
    """ Write all of the data collating the missing data
    """

    # Get the indices of the species
    idx_dct = {name: str(idx+1) for idx, name in enumerate(spc_dct)}

    miss_data = automol.util.remove_duplicates_with_order(miss_data)

    print('\n\n\nMissing Data Requested by the User')
    print(f'{"Name":<20s}{"Method":<12s}{"Basis":<12s}{"Data":<12s}')
    for data in miss_data:
        name, thy_info, dat = data
        if 'ts_' not in name:
            _name = idx_dct[name] + ' ' + name
        else:
            _name = name
        method, basis = thy_info[1:3]
        print(f'{_name:<30s}{method:<12s}{basis:<12s}{dat:<12s}')


def write_data_dirs(data_dirs, prefix):
    """ write a series of data files for similar typ in a directory
    """

    for data_dir in data_dirs:
        name, string_dct = data_dir
        name = name.replace('(', '_').replace(')', '_')
        if string_dct:
            dir_name = ioformat.pathtools.prepare_path(
                [prefix, name], make=True)
            print(f'Writing files into {dir_name} directory')
            for key, string in string_dct.items():
                fname = f'{key}.dat'
                print(f'Writing file {fname}')
                ioformat.pathtools.write_file(string, dir_name, fname)


def get_file_label(tsk, model_dct, proc_keyword_dct, spc_mod_dct_i):
    """ what is the name and extension for this processed file?
    """
    thylabel = ''  # used for data directory labels (e.g. displacements)
    if 'coeffs' in tsk:
        filelabel = 'coeffs'
        filelabel += f'_{model_dct["therm_fit"]["ref_scheme"]}'
        filelabel += '.csv'
    elif 'date' in tsk:
        filelabel = 'date'
        geolvl = proc_keyword_dct.get('geolvl')
        if geolvl is not None:
            filelabel += f'_{proc_keyword_dct["geolvl"]}'
        else:
            filelabel += f'_m{spc_mod_dct_i["vib"]["geolvl"][0]}'
        filelabel += '.csv'
    elif 'freq' in tsk:
        filelabel = 'freq'
        geolvl = proc_keyword_dct.get('geolvl')
        if geolvl is not None:
            thylabel = proc_keyword_dct['geolvl']
            filelabel += f'_{proc_keyword_dct["geolvl"]}'
        else:
            thylabel = spc_mod_dct_i['vib']['geolvl'][0]
            filelabel += f'_m{spc_mod_dct_i["vib"]["geolvl"][0]}'
        filelabel += '.csv'
    elif 'geo' in tsk:
        filelabel = 'geo'
        geolvl = proc_keyword_dct.get('geolvl')
        if geolvl is not None:
            filelabel += f'_{proc_keyword_dct["geolvl"]}'
        else:
            filelabel += f'_m{spc_mod_dct_i["vib"]["geolvl"][0]}'
        filelabel += '.txt'
    elif 'sidata' in tsk:
        filelabel = 'si'
        geolvl = proc_keyword_dct.get('geolvl')
        if geolvl is not None:
            filelabel += f'_{proc_keyword_dct["geolvl"]}'
        else:
            filelabel += f'_m{spc_mod_dct_i["vib"]["geolvl"][0]}'
        filelabel += '.txt'
    elif 'molden' in tsk:
        filelabel = 'molden'
        geolvl = proc_keyword_dct.get('geolvl')
        if geolvl is not None:
            filelabel += f'_{proc_keyword_dct["geolvl"]}'
        else:
            filelabel += f'_m{spc_mod_dct_i["vib"]["geolvl"][0]}'
        filelabel += '.txt'
    elif 'zma' in tsk:
        filelabel = 'zmat'
        geolvl = proc_keyword_dct.get('geolvl')
        if geolvl is not None:
            filelabel += f'_{proc_keyword_dct["geolvl"]}'
        if spc_mod_dct_i:
            filelabel += f'_m{spc_mod_dct_i["vib"]["geolvl"][0]}'
        filelabel += '.txt'
    elif 'ene' in tsk:
        filelabel = 'ene'
        geolvl = proc_keyword_dct.get('geolvl')
        if geolvl is not None:
            filelabel += f'_{proc_keyword_dct["geolvl"]}'
            filelabel += f'_{proc_keyword_dct["proplvl"]}'
        else:
            filelabel += f'_m{spc_mod_dct_i["vib"]["geolvl"][0]}'
            filelabel += f'_m{spc_mod_dct_i["ene"]["lvl1"][0]}'
        filelabel += '.csv'
    elif 'weight' in tsk:
        filelabel = 'weight'
        geolvl = proc_keyword_dct.get('geolvl')
        if geolvl is not None:
            filelabel += f'_{proc_keyword_dct["geolvl"]}'
            filelabel += f'_{proc_keyword_dct["proplvl"]}'
        else:
            filelabel += f'_m{spc_mod_dct_i["vib"]["geolvl"][0]}'
            filelabel += f'_m{spc_mod_dct_i["ene"]["lvl1"][0]}'
        filelabel += '.csv'
    elif 'gibbs' in tsk:
        filelabel = 'gibbs'
        geolvl = proc_keyword_dct.get('geolvl')
        if geolvl is not None:
            filelabel += f'_{proc_keyword_dct["geolvl"]}'
            filelabel += f'_{proc_keyword_dct["proplvl"]}'
        else:
            filelabel += f'_m{spc_mod_dct_i["vib"]["geolvl"][0]}'
            filelabel += f'_m{spc_mod_dct_i["ene"]["lvl1"][0]}'
        filelabel += '.csv'
    elif 'enthalpy' in tsk:
        filelabel = 'enthalpy'
        geolvl = proc_keyword_dct.get('geolvl')
        if geolvl is not None:
            filelabel += f'_{proc_keyword_dct["geolvl"]}'
            filelabel += f'_{proc_keyword_dct["proplvl"]}'
        else:
            filelabel += f'_m{spc_mod_dct_i["vib"]["geolvl"][0]}'
            filelabel += f'_m{spc_mod_dct_i["ene"]["lvl1"][0]}'
        filelabel += '.csv'
    elif 'torsions' in tsk:
        filelabel = 'tors'
        geolvl = proc_keyword_dct.get('geolvl')
        if geolvl is not None:
            filelabel += f'_{proc_keyword_dct["geolvl"]}'
            filelabel += f'_{proc_keyword_dct["proplvl"]}'
        else:
            filelabel += f'_m{spc_mod_dct_i["tors"]["mod"][0]}'
            filelabel += f'_m{spc_mod_dct_i["tors"]["geolvl"][0]}'
    elif 'messpf_inp' in tsk:
        filelabel = 'messpf_input_global'
        filelabel += '.txt'
    elif 'hess_json' in tsk:
        filelabel = 'hess'
        filelabel += f'_m{spc_mod_dct_i["tors"]["mod"][0]}'
        filelabel += f'_m{spc_mod_dct_i["tors"]["geolvl"][0]}'
        filelabel += '.json'
    elif 'pf' in tsk:
        filelabel = 'pf_global'
        filelabel += '.csv'
    return filelabel, thylabel



# Queue manipulation functions
def remove_unstable(spc_queue, spc_dct, thy_dct, spc_mod_dct_i,
                    proc_key_dct, save_prefix):
    """ For each species in the queue see if there are files
        in the save filesystem denoting they are unstable. If so,
        that species is removed from the queue for collection tasks.
    """

    ioprinter.info_message(
        'Removing unstable species with instab.yaml file in SAVE filesystem',
        'from queue...\n')

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
                    f'Found instability file at path {path}', newline=1)
                ioprinter.info_message(
                    f'Removing {spc_name} from queue')

    return stable_queue


def remove_ts_missing(obj_queue, spc_dct):
    """ Remove TSs who have no information from the queue and
        include them in the missing data lists
    """
    ioprinter.info_message(
        'Removing TS with no zmat.r.yaml file in SAVE filesystem',
        'from queue...\n')

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


def remove_radrad_ts(obj_queue, spc_dct):
    """ Remove TSs who have no information from the queue and
        include them in the missing data lists
    """

    ioprinter.info_message(
        'Removing low-spin radical-radical TS from queue...\n')

    new_queue = ()
    for obj in obj_queue:
        if 'ts_' in obj:
            rxn_info = spc_dct[obj]['canon_rxn_info']
            ts_mul = rinfo.value(rxn_info, 'tsmult')
            high_ts_mul = rinfo.ts_mult(rxn_info, rxn_mul='high')
            if rinfo.radrad(rxn_info) and ts_mul != high_ts_mul:
                ioprinter.info_message(f'Removing {obj} from queue.')
            else:
                new_queue += (obj,)
        else:
            new_queue += (obj,)

    return new_queue


def reconcile_spc_mod_and_proc_keyword_dcts(
        tsk, spc_mod_dct, proc_keyword_dct, pes_mod_dct, thy_dct):
    spc_model = proc_keyword_dct.get('spc_model', 'global')
    kin_model = proc_keyword_dct.get('kin_model', 'global')
    if spc_model is None:
        print('No base species model specified, defaulting to "global" in models.dat')
        spc_model = 'global'
    if kin_model is None:
        print('No base kinetic model specified, defaulting to "global" in models.dat')
        kin_model = 'global'
    spc_mod_dct_i = spc_mod_dct[spc_model].copy()
    pes_mod_dct_i = pes_mod_dct[kin_model].copy()
    geo_lvl = proc_keyword_dct.get('geolvl', None)
    prop_lvl = proc_keyword_dct.get('proplvl', None)
    geo_tsks = [
        'date', 'geo', 'molden', 'zmatrix', 'torsions',
        'pf', 'messpf_inp', 'gibbs', 'weight', 'energy']
    if tsk in geo_tsks and geo_lvl:
        print('updating species model with specified geolvl', geo_lvl)
        geo_thy_info = (
            geo_lvl,
            (1.00, tinfo.from_dct(thy_dct.get(geo_lvl))))
        for key in ['symm', 'vib', 'tors']:
            if key in spc_mod_dct_i:
                spc_mod_dct_i[key]['geolvl'] = geo_thy_info
    if prop_lvl is not None:
        print('updating species model with specified proplvl', geo_lvl)
        if tsk in ['energy', 'enthalpy', 'gibbs', 'weight']:
            prop_thy_info = (
                prop_lvl,
                (1.00, tinfo.from_dct(thy_dct.get(prop_lvl))))
            spc_mod_dct_i['ene']['lvl1'] = prop_thy_info
        if tsk in ['freqs', 'hess_json']:
            prop_thy_info = (
                prop_lvl,
                (1.00, tinfo.from_dct(thy_dct.get(prop_lvl))))
            spc_mod_dct_i['vib']['geolvl'] = prop_thy_info
        if tsk in ['torsions']:
            prop_thy_info = (
                prop_lvl,
                (1.00, tinfo.from_dct(thy_dct.get(prop_lvl))))
            spc_mod_dct_i['tors']['enelvl'] = prop_thy_info
    return spc_mod_dct_i, pes_mod_dct_i

