""" read species
"""

import os
import chemkin_io
import automol
import autoparse.find as apf
import autofile
from lib.load import ptt
from lib.load import run as loadrun
from lib.phydat import symm, eleclvl
from lib.reaction import rxnid
from lib.filesystem import path as fpath
from lib.filesystem import read as fread
from lib.filesystem import inf as finf


CSV_INP = 'inp/species.csv'
DAT_INP = 'inp/species.dat'
CLA_INP = 'inp/class.csv'


def build_run_spc_dct(spc_dct, run_obj_dct):
    """ Get a dictionary of requested species matching the PES_DCT format
    """
    # Fix stuff, this function does not work right now

    spc_nums = [idx for idx in run_obj_dct['spc']]

    run_spc_lst = []
    for idx, spc in enumerate(spc_dct):
        if idx+1 in spc_nums:
            # model = spcmods[idx]
            run_spc_lst.append((spc, run_obj_dct['spc'][idx+1]))

    # Build the run dct
    run_dct = {}
    run_dct['all'] = {
        'species': run_spc_lst,
        'reacs': [],
        'prods': []
    }

    return run_dct


def build_spc_dct(job_path, spc_type):
    """ Get a dictionary of all the input species
        indexed by InChi string
    """
    spc_csv_str = ptt.read_inp_str(job_path, CSV_INP)
    if spc_type == 'csv':
        spc_dct = csv_dct(spc_csv_str, check_stereo=False)
    else:
        raise NotImplementedError

    # Modify spc dct with params from the AMech file
    mod_spc_dct = modify_spc_dct(job_path, spc_dct)

    return mod_spc_dct


def csv_dct(spc_str, check_stereo):
    """ read the species file in a .csv format
    """
    # Read in the initial CSV information
    smi_dct = chemkin_io.parser.mechanism.spc_name_dct(spc_str, 'smiles')
    ich_dct = chemkin_io.parser.mechanism.spc_name_dct(spc_str, 'inchi')
    mul_dct = chemkin_io.parser.mechanism.spc_name_dct(spc_str, 'mult')
    chg_dct = chemkin_io.parser.mechanism.spc_name_dct(spc_str, 'charge')
    sens_dct = chemkin_io.parser.mechanism.spc_name_dct(spc_str, 'sens')

    # Print message detailing if stereochemistry is used
    print_check_stero_msg(check_stereo)

    # Rebuild with stereochemistry if possible
    if check_stereo:
        spc_str = 'name,SMILES,InChI,mult,charge,sens \n'
        for name in ich_dct:
            ich = ich_dct[name]
            smi = smi_dct[name]
            mul = mul_dct[name]
            chg = chg_dct[name]
            sens = sens_dct[name]
            if not automol.inchi.is_complete(ich):
                print('adding stereochemistry for {0}, {1}, {2}'.format(
                    name, smi, ich))
                # this returns a list of ichs w/ different possible stereo vals
                # for now just taking the first of these
                ich = automol.inchi.add_stereo(ich)
                ich = ich[-1]
                ich_dct[name] = ich
            spc_str += '{0},\'{1}\',\'{2}\',{3},{4},{5} \n'.format(
                name, smi, ich, mul, chg, sens)

        stereo_path = 'species_stereo.csv'
        with open(stereo_path, 'w') as stereo_csv_file:
            stereo_csv_file.write(spc_str)

    # Build the final dictionary
    spc_names = []
    spc_dct = {}
    for name in mul_dct:
        spc_dct[name] = {}
        spc_dct[name]['smi'] = smi_dct[name]
        spc_dct[name]['ich'] = ich_dct[name]
        spc_dct[name]['chg'] = chg_dct[name]
        spc_dct[name]['mul'] = mul_dct[name]
        spc_names.append(name)

    return spc_dct


def read_spc_amech(job_path):
    """ Read an amech style input file for the species
    """

    spc_amech_str = ptt.read_inp_str(job_path, DAT_INP)

    spc_dct = {}
    if spc_amech_str:
        spc_sections = apf.all_captures(
            ptt.end_section_wname2('spc'), spc_amech_str)
        if spc_sections:
            for section in spc_sections:
                name = section[0]
                keyword_dct = ptt.build_keyword_dct(section[1])
                spc_dct[name] = keyword_dct

    return spc_dct


def modify_spc_dct(job_path, spc_dct):
    """ Modify the species dct using input from the additional AMech file
    """

    # Read in other dcts
    amech_dct = read_spc_amech(job_path)
    geom_dct = geometry_dictionary(job_path)

    mod_spc_dct = {}
    for spc in spc_dct:
        # Set the ich and mult
        ich = spc_dct[spc]['ich']
        mul = spc_dct[spc]['mul']
        chg = spc_dct[spc]['chg']
        mod_spc_dct[spc] = {}
        mod_spc_dct[spc]['ich'] = ich
        mod_spc_dct[spc]['mul'] = mul
        mod_spc_dct[spc]['chg'] = chg

        # Add the parameters from amech file
        if spc in amech_dct:
            if 'hind_inc' in amech_dct[spc]:
                mod_spc_dct[spc]['hind_inc'] = amech_dct[spc]['hind_inc']
            if 'hind_def' in amech_dct[spc]:
                mod_spc_dct[spc]['hind_def'] = amech_dct[spc]['hind_def']
            if 'elec_levs' in amech_dct[spc]:
                mod_spc_dct[spc]['elec_levs'] = amech_dct[spc]['elec_levs']
            if 'sym' in amech_dct[spc]:
                mod_spc_dct[spc]['sym'] = amech_dct[spc]['sym']
            if 'mc_nsamp' in amech_dct[spc]:
                mod_spc_dct[spc]['mc_nsamp'] = amech_dct[spc]['mc_nsamp']

        # Add the parameters from std lib if needed
        if 'elec_levs' not in mod_spc_dct[spc]:
            if (ich, mul) in eleclvl.DCT:
                mod_spc_dct[spc]['elec_levs'] = eleclvl.DCT[(ich, mul)]
            else:
                mod_spc_dct[spc]['elec_levs'] = [[0.0, mul]]
        if 'sym' not in mod_spc_dct[spc]:
            if (ich, mul) in symm.DCT:
                mod_spc_dct[spc]['sym'] = symm.DCT[(ich, mul)]

        # Add geoms from geo dct (prob switch to amech file)
        if ich in geom_dct:
            mod_spc_dct[spc]['geo_obj'] = geom_dct[ich]

    return mod_spc_dct


def read_class_dct(job_path):
    """ Read the class dictionary
    """
    cla_str = ptt.read_inp_str(job_path, CLA_INP)
    cla_str = cla_str if cla_str else 'REACTION,RCLASS\nEND'
    cla_dct = chemkin_io.parser.mechanism.reac_class_dct(cla_str, 'class')

    return cla_dct


def geometry_dictionary(job_path):
    """ read in dictionary of saved geometries
    """
    geom_path = os.path.join(job_path, 'data', 'geoms')
    geom_dct = {}
    for dir_path, _, file_names in os.walk(geom_path):
        for file_name in file_names:
            file_path = os.path.join(dir_path, file_name)
            if file_path.endswith('.xyz'):
                xyz_str = autofile.file.read_file(file_path)
                geo = automol.geom.from_xyz_string(xyz_str)
                ich = automol.geom.inchi(geo)
                if ich in geom_dct:
                    print('Warning: Dupilicate xyz geometry for ', ich)
                geom_dct[ich] = geo
    return geom_dct


def build_sadpt_dct(rxn_lst, model_dct, thy_dct, es_tsk_str,
                    run_inp_dct, run_options_dct, spc_dct, cla_dct):
    """ build dct
    """
    run_prefix = run_inp_dct['run_prefix']
    save_prefix = run_inp_dct['save_prefix']
    kickoff = run_options_dct['kickoff']

    print('\nBegin transition state prep')
    ts_dct = {}
    ts_idx = 0
    for _, rxn in enumerate(rxn_lst):
        es_tsk_lst = loadrun.build_run_es_tsks_lst(
            es_tsk_str, model_dct[rxn['model']]['es'], thy_dct, saddle=True)
        reacs = rxn['reacs']
        prods = rxn['prods']
        tsname = 'ts_{:g}'.format(ts_idx)
        ts_dct[tsname] = {}
        # Fix this rname build and class dct stuff
        # rname = rxn_name_lst[idx]
        # rname_eq = '='.join(rname.split('=')[::-1])
        # if rname in cla_dct:
        #     ts_dct[tsname]['given_class'] = cla_dct[rname]
        # elif rname_eq in cla_dct:
        #     ts_dct[tsname]['given_class'] = cla_dct[rname_eq]
        #     reacs = rxn['prods']
        #     prods = rxn['reacs']
        # else:
        #    ts_dct[tsname]['given_class'] = None
        ts_dct[tsname]['given_class'] = None
        if reacs and prods:
            ts_dct[tsname]['reacs'] = reacs
            ts_dct[tsname]['prods'] = prods
        ts_dct[tsname]['ich'] = ''
        ts_chg = 0
        for rct in reacs:  # rct_names_lst[idx]:
            ts_chg += spc_dct[rct]['chg']
        ts_dct[tsname]['chg'] = ts_chg
        mul_low, _, rad_rad = rxnid.ts_mul_from_reaction_muls(
            reacs, prods, spc_dct)
        ts_dct[tsname]['mul'] = mul_low
        ts_dct[tsname]['rad_rad'] = rad_rad
        ts_dct[tsname] = set_sadpt_info(
             es_tsk_lst, ts_dct, spc_dct, tsname,
             run_prefix, save_prefix, kickoff)
        ts_idx += 1

    return ts_dct


def set_sadpt_info(es_tsk_lst, ts_dct, spc_dct, sadpt,
                   run_prefix, save_prefix, kickoff):
    """ set the saddle point dct with info
    """
    for tsk in es_tsk_lst:
        if 'find_ts' in tsk:
            ini_thy_info = es_tsk_lst[0][2]
            thy_info = es_tsk_lst[0][1]
            break

    # Generate rxn data, reorder if necessary, and put in spc_dct for given ts
    rxn_ichs, rxn_chgs, rxn_muls, low_mul, high_mul = finf.rxn_info(
        save_prefix, sadpt, ts_dct, spc_dct, thy_info)
    ts_dct[sadpt]['low_mul'] = low_mul
    ts_dct[sadpt]['high_mul'] = high_mul

    # Generate rxn_fs from rxn_info stored in spc_dct
    [kickoff_size, kickoff_backward] = kickoff
    rxn_run_fs, rxn_save_fs, rxn_run_path, rxn_save_path = fpath.get_rxn_fs(
        run_prefix, save_prefix, rxn_ichs, rxn_chgs, rxn_muls, ts_dct[sadpt]['mul'])
    ts_dct[sadpt]['rxn_fs'] = [
        rxn_run_fs,
        rxn_save_fs,
        rxn_run_path,
        rxn_save_path]
    rct_zmas, prd_zmas, rct_cnf_save_fs, prd_cnf_save_fs = fread.get_zmas(
        ts_dct[sadpt]['reacs'], ts_dct[sadpt]['prods'], spc_dct,
        ini_thy_info, save_prefix, run_prefix, kickoff_size,
        kickoff_backward)
    ret = rxnid.ts_class(
        rct_zmas, prd_zmas, ts_dct[sadpt]['rad_rad'],
        ts_dct[sadpt]['mul'], low_mul, high_mul,
        rct_cnf_save_fs, prd_cnf_save_fs, ts_dct[sadpt]['given_class'])
    ret1, ret2 = ret
    if ret1:
        [rxn_class, spc_zma,
         dist_name, brk_name, grid,
         frm_bnd_key, brk_bnd_key,
         tors_names, update_guess] = ret1
        ts_dct[sadpt]['class'] = rxn_class
        ts_dct[sadpt]['grid'] = grid
        ts_dct[sadpt]['tors_names'] = tors_names
        ts_dct[sadpt]['zma'] = spc_zma
        # ts_dct[sadpt]['original_zma'] = spc_zma
        dist_info = [dist_name, 0., update_guess, brk_name]
        ts_dct[sadpt]['dist_info'] = dist_info
        ts_dct[sadpt]['frm_bnd_key'] = frm_bnd_key
        ts_dct[sadpt]['brk_bnd_key'] = brk_bnd_key
        # Adding in the rct and prd zmas for vrctst
        ts_dct[sadpt]['rct_zmas'] = rct_zmas
        ts_dct[sadpt]['prd_zmas'] = prd_zmas
        if ret2:
            ts_dct[sadpt]['bkp_data'] = ret2
        else:
            ts_dct[sadpt]['bkp_data'] = None
    else:
        ts_dct[sadpt]['class'] = None
        ts_dct[sadpt]['bkp_data'] = None

    return ts_dct[sadpt]


# Print messages
def print_check_stero_msg(check_stereo):
    """ Print a message detailing whether stereo is being checked
    """
    if check_stereo:
        print('Species will be treated with stereochemistry.',
              'Checking and adding stereo.')
    else:
        print('Stereochemistry will be ignored')
