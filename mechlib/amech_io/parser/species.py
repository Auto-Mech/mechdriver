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
from mechlib import filesys
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
        if spc_sections:
            # Get the global species section
            for section in spc_sections:
                if section[0] == 'global':
                    # Build the global species section
                    keyword_dct = ioformat.ptt.build_keyword_dct(section[1])
                    amech_dct['global'] = keyword_dct
                else:
                    # Build each species dct to overwrite global dct
                    name = section[0]
                    keyword_dct = ioformat.ptt.build_keyword_dct(section[1])
                    amech_dct[name] = keyword_dct

    return amech_dct


def modify_spc_dct(job_path, spc_dct):
    """ Modify the species dct using input from the additional AMech file
    """

    # Read in other dcts
    amech_dct = read_spc_amech(job_path)
    geom_dct = geometry_dictionary(job_path)

    mod_spc_dct = {}
    for spc in spc_dct:
        # Set the ich and mult
        ich = spc_dct[spc]['inchi']
        mul = spc_dct[spc]['mult']
        chg = spc_dct[spc]['charge']
        mod_spc_dct[spc] = {}
        mod_spc_dct[spc]['inchi'] = ich
        mod_spc_dct[spc]['mult'] = mul
        mod_spc_dct[spc]['charge'] = chg

        # Add params from global dct in amech dct
        if 'global' in amech_dct:
            for key, val in amech_dct['global'].items():
                mod_spc_dct[spc][key] = val

        # Add/Reset the parameters from the species specific dct
        if spc in amech_dct:
            for key, val in amech_dct[spc].items():
                mod_spc_dct[spc][key] = val

    # Second loop for transtion states defined in species.dat
    for spc in amech_dct:
        if 'ts' in spc:
            mod_spc_dct[spc] = {}
            # Add params from global dct in amech dct
            if 'global' in amech_dct:
                for key, val in amech_dct['global'].items():
                    mod_spc_dct[spc][key] = val
            # Add the ts stuff
            mod_spc_dct[spc] = amech_dct[spc]

    # Add global to mod_spc_dct for other TS stuff later
    if 'global' in amech_dct:
        mod_spc_dct['global'] = amech_dct['global']

    # Final loop to add in things that are needed but could be missing
    for spc in mod_spc_dct:
        if spc != 'global' and 'ts_' not in spc:
            ich, mul = mod_spc_dct[spc]['inchi'], mod_spc_dct[spc]['mult']
            if 'elec_levels' not in mod_spc_dct[spc]:
                if (ich, mul) in eleclvl.DCT:
                    mod_spc_dct[spc]['elec_levels'] = eleclvl.DCT[(ich, mul)]
                else:
                    mod_spc_dct[spc]['elec_levels'] = [[0.0, mul]]
            if 'sym_factor' not in mod_spc_dct[spc]:
                if (ich, mul) in symm.DCT:
                    mod_spc_dct[spc]['sym_factor'] = symm.DCT[(ich, mul)]

        # Add defaults here for now
        if 'hind_inc' not in mod_spc_dct[spc]:
            mod_spc_dct[spc]['hind_inc'] = 30.0
        if 'kickoff' not in mod_spc_dct[spc]:
            mod_spc_dct[spc]['kickoff'] = [0.1, False]
        if 'tau_nsamp' not in mod_spc_dct[spc]:
            mod_spc_dct[spc]['tau_nsamp'] = [True, 12, 1, 3, 100, 25]
        if 'mc_nsamp' not in mod_spc_dct[spc]:
            mod_spc_dct[spc]['mc_nsamp'] = [True, 12, 1, 3, 100, 25]
        if 'pst_params' not in mod_spc_dct[spc]:
            mod_spc_dct[spc]['pst_params'] = [1.0, 6]

        # Set active space stuff
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

        # Perform conversions as needed
        mod_spc_dct[spc]['hind_inc'] *= phycon.DEG2RAD
        # if 'ts' in spc:
        #     print(mod_spc_dct[spc]['hind_inc'])

        # Add geoms from geo dct (prob switch to amech file)
        if ich in geom_dct:
            mod_spc_dct[spc]['geo_inp'] = geom_dct[ich]

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
                  run_inp_dct, spc_dct, cla_dct,
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
                run_inp_dct, spc_dct, cla_dct,
                direction=direction)
            break

    # Build the queue
    if ts_dct:
        ts_queue = [(sadpt, '') for sadpt in ts_dct]
    else:
        ts_queue = []

    return ts_dct, ts_queue


def build_sadpt_dct(pes_idx, rxn_lst, thy_info, ini_thy_info,
                    run_inp_dct, spc_dct, cla_dct, run_prefix, save_prefix,
                    direction='forw'):
    """ build a dct for saddle points for all reactions in rxn_lst
    """

    ts_dct = {}
    for rxn in rxn_lst:
        tsname = 'ts_{:g}_{:g}'.format(pes_idx, rxn['chn_idx'])
        ts_dct[tsname] = build_sing_chn_sadpt_dct(
            tsname, rxn, thy_info, ini_thy_info, run_prefix, save_prefix,
            run_inp_dct, spc_dct, cla_dct, direction=direction)

    return ts_dct


def build_sing_chn_sadpt_dct(tsname, reaction, thy_info, ini_thy_info,
                             run_inp_dct, spc_dct, cla_dct, run_prefix, save_prefix,
                             direction='forw'):
    """ build dct for single reaction
    """

    save_prefix = run_inp_dct['save_prefix']

    reacs = reaction['reacs']
    prods = reaction['prods']
    rxn_info = rinfo.from_dct(reacs, prods, spc_dct)
    print('  Preparing {} for reaction {} = {}'.format(
        tsname, '+'.join(reacs), '+'.join(prods)))

    # Set the reacs and prods for the desired direction
    reacs, prods, _ = rxnid.set_reaction_direction(
        reacs, prods, rxn_info, cla_dct,
        thy_info, ini_thy_info, save_prefix, direction=direction)

    # Obtain the reaction object for the reaction
    zma_locs = (0,)
    zrxn, zma = rxnid.build_reaction(
        rxn_info, ini_thy_info, zma_locs, save_prefix)

    if zrxn is not None:
        rxn_save_fs = autofile.fs.reaction(save_prefix)
        rxn_save_fs[-1].create(rinfo.sort(rxn_info))
        rxn_save_path = rxn_save_fs[-1].path(rinfo.sort(rxn_info))
        rxn_run_fs = autofile.fs.reaction(run_prefix)
        rxn_run_fs[-1].create(rinfo.sort(rxn_info))
        rxn_run_path = rxn_run_fs[-1].path(rinfo.sort(rxn_info))
        rxn_fs = [rxn_run_fs, rxn_save_fs, rxn_run_path, rxn_save_path]
        ts_dct = {
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
        if 'kickoff' not in combined_dct[sadpt]:
            combined_dct[sadpt]['kickoff'] = [0.1, False]
        if 'mc_nsamp' not in combined_dct[sadpt]:
            combined_dct[sadpt]['mc_nsamp'] = [True, 12, 1, 3, 100, 25]
        if 'tau_nsamp' not in combined_dct[sadpt]:
            combined_dct[sadpt]['tau_nsamp'] = [True, 12, 1, 3, 100, 25]
        if 'irc_idxs' not in combined_dct[sadpt]:
            combined_dct[sadpt]['irc_idxs'] = [
                -10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0,
                -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0,
                8.0, 9.0, 10.0]
        if 'pst_params' not in combined_dct[sadpt]:
            combined_dct[sadpt]['pst_params'] = [1.0, 6]

        # Perform conversions as needed
        # combined_dct[spc]['hind_inc'] *= phycon.DEG2RAD

    return combined_dct
