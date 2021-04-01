""" read spcies
"""

import os
import automol
import autofile
import mechanalyzer
import ioformat
import autoparse.find as apf
from mechanalyzer.inf import thy as tinfo
from mechanalyzer.inf import rxn as rinfo
from phydat import symm
from phydat import eleclvl
from phydat import phycon
from mechlib.reaction import rxnid


CSV_INP = 'inp/species.csv'
DAT_INP = 'inp/species.dat'


# Build main dcts
def build_queue(rxn_lst):
    """ Build spc queue from the reaction lst for the drivers
        :return spc_queue: all the species and corresponding models in rxn
        :rtype: list[(species, model),...]
    """
    if 'all' in rxn_lst:
        # First check if rxn_lst is a bunch of species
        spc_queue = rxn_lst['all']['species']
    else:
        # Build the list from expanding the reacs and prods
        spc_queue = []
        for rxn in rxn_lst:
            model = rxn['model']
            spc_queue.extend(((reac, model) for reac in rxn['reacs']))
            spc_queue.extend(((prod, model) for prod in rxn['prods']))

    return spc_queue


def split_queue(spc_queue):
    new_queue = []
    op_dct = {'*': 'multiply', '+': 'add', '/': 'divide', '-': 'substract'}
    for (spc_name, (pes_model, spc_model)) in spc_queue:
        coeffs = []
        operators = []
        models = []
        coeff = ''
        model = ''
        for char in spc_model:
            if char == '.' or char.isdigit():
                coeff += char
            elif char.isalpha():
                model += char
            elif char in op_dct:
                operators.append(op_dct[char])
                if coeff:
                    coeffs.append(float(coeff))
                else:
                    coeffs.append(1)
                models.append(model)
                coeff = ''
                model = ''
        if coeff:
            coeffs.append(float(coeff))
        else:
            coeffs.append(1)
        models.append(model)
        new_queue.append((spc_name, (pes_model, models, coeffs, operators)))
    return new_queue


def build_spc_dct(job_path, spc_type):
    """ Build the species dct
    """

    spc_str = ioformat.ptt.read_inp_str(
        job_path, CSV_INP)
    spc_dct = mechanalyzer.parser.spc.build_spc_dct(spc_str, spc_type)

    # Modify spc dct with params from the AMech file
    mod_spc_dct = modify_spc_dct(job_path, spc_dct)

    return mod_spc_dct


def build_run_spc_dct(spc_dct, run_obj_dct):
    """ Get a dictionary of requested species matching the PES_DCT format
    """
    spc_nums = run_obj_dct['spc']
    run_spc_lst = []
    for idx, spc in enumerate(spc_dct):
        if spc != 'global':
            if idx+1 in spc_nums:
                run_spc_lst.append((spc, run_obj_dct['spc'][idx+1]))

    # Build the run dct
    run_dct = {}
    run_dct['all'] = {
        'species': run_spc_lst,
        'reacs': [],
        'prods': []
    }

    return run_dct


def build_spc_queue(rxn_lst):
    """ Build spc queue from the reaction lst for the drivers
        :return spc_queue: all the species and corresponding models in rxn
        :rtype: list[(species, model),...]
    """

    if 'all' in rxn_lst:
        # First check if rxn_lst is a bunch of species
        spc_queue = rxn_lst['all']['species']
    else:
        # Build the list from expanding the reacs and prods
        spc_queue = []
        for rxn in rxn_lst:
            model = rxn['model']
            spc_queue.extend(((reac, model) for reac in rxn['reacs']))
            spc_queue.extend(((prod, model) for prod in rxn['prods']))

    return spc_queue


def read_spc_amech(job_path):
    """ Read an amech style input file for the species
    """

    # Read the AMech species string
    if os.path.exists(os.path.join(job_path, DAT_INP)):
        spc_amech_str = ioformat.ptt.read_inp_str(
            job_path, DAT_INP, remove_comments='#')
        print('Found species.dat. Reading file...')
    else:
        spc_amech_str = ''
        print('No species.dat file...')

    # Build the keyword dcts
    amech_dct = {}
    if spc_amech_str:
        # Read each of the species sections and build the dcts
        spc_sections = apf.all_captures(
            ioformat.ptt.end_section_wname2('spc'), spc_amech_str)
        amech_dct = ptt.build_keyword_dct_msec(spc_sections)

    # Update the amech dct with the global
    new_amech_dct = {}
    glob = amech_dct.get('global', {})
    for spc in (x for x in amech_dct if x != 'global'):
        new_amech_dct[spc] = right_update(glob, amech_dct[spc])

    return amech_dct


def modify_spc_dct(job_path, spc_dct):
    """ Modify the species dct using input from the additional AMech file
    """

    # Read in other dcts
    amech_dct = read_spc_amech(job_path)
    geom_dct = geometry_dictionary(job_path)

    # Add in all of the species
    for spc in spc_dct:
        # Add stuff from the amech dct
        spc_dct[spc] = automol.util.dict_.right_update(
            spc_dct[spc], amech_dct.get(spc, {}))
        # Add the defaults
        spc_dct[spc] = automol.util.dict_.right_update(
            SPC_DEFAULT_DCT, mod_spc_dct[spc])

        # Add speciaized calls not in the default dct
        ich, mul = spc_dct[spc]['inchi'], spc_dct[spc]['mult']
        if 'elec_levels' not in spc_dct[spc]:
            spc_dct[spc]['elec_levels'] = eleclvl.DCT.get((ich, mul), (0.0, mul))
        if 'sym_factor' not in spc_dct[spc]:
            spc_dct[spc]['sym_factor'] = symm.DCT.get((ich, mul), 1.0)

    # Add the transitions states defined in species.dat that are not defined in spc_dct
    # they are not in the spc_dct currently, since we don't define TSs there; built later
    spc_dct.update(ts_dct)
    ts_dct = {}
    for spc in (x for x in amech_dct if 'ts' in x):
        ts_dct[spc] = {**amech_dct[spc]}
        
        # Add speciaized calls not in the default dct
        if 'active' not in mod_spc_dct[spc]:
            mod_spc_dct[spc]['active_space'] = None
        else:
            aspace = mod_spc_dct[spc].get('active')
            assert len(aspace) == 4, (
                'active must be length 4: {}'.format(aspace)
            )
            wfn_file = aspace[3]
            wfn_inp = os.path.join(os.path.join(job_path, 'inp/'+wfn_file))
            if os.path.exists(wfn_inp):
                wfn_str = ioformat.ptt.read_inp_str(job_path, wfn_inp)
                print('Found file: {}. Reading file...'.format(wfn_file))
            else:
                wfn_str = None
                print('No file: {}. Reading file...'.format(wfn_file))
            mod_spc_dct[spc]['active_space'] = (
                aspace[0], aspace[1], aspace[2], wfn_str)
    
    # add the TSs to the spc dct
    spc_dct.update(ts_dct)

        # Perform conversions as needed
        mod_spc_dct[spc]['hind_inc'] *= phycon.DEG2RAD

        mod_spc_dct[spc]['geo_inp'] = geom_dct.get(ich, None)

    return mod_spc_dct


def geometry_dictionary(job_path):
    """ read in dictionary of saved geometries
    """
    geom_path = os.path.join(job_path, 'data')
    geom_dct = {}
    for dir_path, _, file_names in os.walk(geom_path):
        for file_name in file_names:
            file_path = os.path.join(dir_path, file_name)
            if file_path.endswith('.xyz'):
                xyz_str = autofile.io_.read_file(file_path)
                geo = automol.geom.from_xyz_string(xyz_str)
                ich = automol.geom.inchi(geo)
                if ich in geom_dct:
                    print('Warning: Dupilicate xyz geometry for ', ich)
                geom_dct[ich] = geo
    return geom_dct


def get_sadpt_dct(pes_idx, es_tsk_lst, rxn_lst, thy_dct,
                  run_inp_dct, spc_dct, run_prefix, save_prefix,
                  direction='forw'):
    """ build a ts queue
    """

    print('\nTasks for transition states requested...')
    print('Identifying reaction classes for transition states...')

    # Build the ts_dct
    ts_dct = {}
    for tsk_lst in es_tsk_lst:
        [obj, _, es_keyword_dct] = tsk_lst
        if 'ts' in obj:
            method_dct = thy_dct.get(es_keyword_dct['runlvl'])
            ini_method_dct = thy_dct.get(es_keyword_dct['inplvl'])
            thy_info = tinfo.from_dct(method_dct)
            ini_thy_info = tinfo.from_dct(ini_method_dct)
            ts_dct = build_sadpt_dct(
                pes_idx, rxn_lst, thy_info, ini_thy_info,
                run_inp_dct, spc_dct, run_prefix, save_prefix,
                direction=direction)
            break

    # Build the queue
    if ts_dct:
        ts_queue = [(sadpt, '') for sadpt in ts_dct]
    else:
        ts_queue = []

    return ts_dct, ts_queue


def build_sadpt_dct(pes_idx, rxn_lst, thy_info, ini_thy_info,
                    run_inp_dct, spc_dct, run_prefix, save_prefix,
                    direction='forw'):
    """ build a dct for saddle points for all reactions in rxn_lst
    """

    ts_dct = {}
    for rxn in rxn_lst:
        ts_dct.update(
            build_sing_chn_sadpt_dct(
                pes_idx, rxn, thy_info, ini_thy_info,
                run_inp_dct, spc_dct, run_prefix, save_prefix,
                direction=direction)
        )

    return ts_dct


# def build_sing_chn_sadpt_dct(tsname, reaction, thy_info, ini_thy_info,
def build_sing_chn_sadpt_dct(pes_idx, reaction, thy_info, ini_thy_info,
                             run_inp_dct, spc_dct, run_prefix, save_prefix,
                             direction='forw'):
    """ build dct for single reaction
    """

    save_prefix = run_inp_dct['save_prefix']

    reacs = reaction['reacs']
    prods = reaction['prods']
    rxn_info = rinfo.from_dct(reacs, prods, spc_dct)
    print('  Preparing for reaction {} = {}'.format(
        '+'.join(reacs), '+'.join(prods)))

    # Set the reacs and prods for the desired direction
    reacs, prods = rxnid.set_reaction_direction(
        reacs, prods, rxn_info,
        thy_info, ini_thy_info, save_prefix, direction=direction)

    # Obtain the reaction object for the reaction
    zma_locs = (0,)
    zrxns, zmas = rxnid.build_reaction(
        rxn_info, ini_thy_info, zma_locs, save_prefix)

    # ts_dct = {}
    if zrxns is not None:
        ts_dct = {}
        for idx, (zrxn, zma) in enumerate(zip(zrxns, zmas)):
            rxn_run_fs = autofile.fs.reaction(run_prefix)
            rxn_save_fs = autofile.fs.reaction(save_prefix)
            rxn_run_path = rxn_run_fs[-1].path(rinfo.sort(rxn_info))
            rxn_save_path = rxn_save_fs[-1].path(rinfo.sort(rxn_info))
            rxn_fs = [rxn_run_fs, rxn_save_fs, rxn_run_path, rxn_save_path]
            tsname = 'ts_{:g}_{:g}_{:g}'.format(
                pes_idx, reaction['chn_idx'], idx)
            ts_dct[tsname] = {
                'zrxn': zrxn,
                'zma': zma,
                'reacs': reacs,
                'prods': prods,
                'rxn_info': rxn_info,
                'inchi': '',
                'charge': rinfo.value(rxn_info, 'charge'),
                'mult': rinfo.value(rxn_info, 'tsmult'),
                'elec_levels': [[0.0, rinfo.value(rxn_info, 'tsmult')]],
                'class': zrxn.class_,
                'rxn_fs': rxn_fs
            }
    else:
        ts_dct = {}
        print('Skipping reaction as class not given/identified')

    return ts_dct


def combine_sadpt_spc_dcts(sadpt_dct, spc_dct):
    """ Create a new dictionary that combines init spc_dct and sadpt dct
    """

    combined_dct = {}

    # Put all elements of spc_dct in combined dct that are NOT TSs
    for spc in spc_dct:
        if 'ts' not in spc:
            combined_dct[spc] = spc_dct[spc]

    # Now put in the TSs pulling info from everywhere
    for sadpt in sadpt_dct:

        # Put in stuff from the global dct
        combined_dct[sadpt] = {}
        if 'global' in spc_dct:
            for key, val in spc_dct['global'].items():
                combined_dct[sadpt][key] = val

        # Update any sadpt keywords if they are in the spc_dct from .dat file
        if sadpt in spc_dct:
            combined_dct[sadpt].update(spc_dct[sadpt])

        # Put in stuff from the sadpt_dct build
        for key, val in sadpt_dct[sadpt].items():
            if key not in combined_dct[sadpt]:
                combined_dct[sadpt][key] = val

        # Put in defaults if they were not defined
        combined_dct[sadpt] = automol.util._dict.right_update(
           keyword.TS_DEFAULT_DCT, combined_dct[sadpt])

    return combined_dct
