""" Build species dictionary for all spc and ts
"""

import os
import automol
import autofile
import ioformat
from mechanalyzer.inf import thy as tinfo
from mechanalyzer.inf import rxn as rinfo
from phydat import symm, eleclvl, phycon
from mechlib.reaction import rxnid
from mechlib.amech_io.parser._keywrd import defaults_from_val_dct
from mechlib.amech_io.parser._keywrd import SPC_VAL_DCT, TS_VAL_DCT


def modify_spc_dct(spc_dct, amech_dct, geo_dct):
    """ Modify the species dct using input from the additional AMech file
    """

    # Separate the global dct
    dat_dct, glob_dct = automol.util.dict_.separate_subdct(
        amech_dct, key='global')

    # Add in all of the species
    for spc in spc_dct:

        # Add stuff from the main amech dct and global dct
        spc_dct[spc] = automol.util.dict_.right_update(
            spc_dct[spc], glob_dct)
        spc_dct[spc] = automol.util.dict_.right_update(
            spc_dct[spc], amech_dct.get(spc, {}))

        # Add the defaults
        spc_dct[spc] = automol.util.dict_.right_update(
            defaults_from_val_dct(SPC_VAL_DCT), spc_dct[spc])

        # Add speciaized calls not in the default dct
        ich, mul = spc_dct[spc]['inchi'], spc_dct[spc]['mult']
        if 'elec_levels' not in spc_dct[spc]:
            spc_dct[spc]['elec_levels'] = eleclvl.DCT.get(
                (ich, mul), (0.0, mul))
        if 'sym_factor' not in spc_dct[spc]:
            spc_dct[spc]['sym_factor'] = symm.DCT.get(
                (ich, mul), 1.0)

    # Add transitions states defined in species.dat not defined in spc_dct
    ts_dct = {}
    for tsname in (x for x in dat_dct if 'ts' in x):
        ts_dct[tsname] = {**dat_dct[tsname]}

        # Need to add the TS defaults
        ts_dct[tsname] = automol.util.dict_.right_update(
            defaults_from_val_dct(TS_VAL_DCT), ts_dct[tsname])

        # Add speciaized calls not in the default dct
        # _set_active_key()

    # add the TSs to the spc dct
    spc_dct.update(ts_dct)

    # Final loop for conversions additions
    for spc in spc_dct:
        spc_dct[spc]['hind_inc'] *= phycon.DEG2RAD
        # spc_dct[spc]['geo_inp'] = geom_dct.get(ich, None)

    return spc_dct


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


# Functions to the spc_dct contributions for TS
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
