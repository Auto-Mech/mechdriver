""" Build species dictionary for all spc and ts
"""

import automol
import ioformat
import mechanalyzer
from mechanalyzer.inf import thy as tinfo
from mechanalyzer.inf import rxn as rinfo
from phydat import symm, eleclvl, phycon
from mechlib.reaction import rxnid
from mechlib.filesys import reaction_fs
from mechlib.amech_io.parser._keywrd import defaults_from_val_dct
from mechlib.amech_io.parser._keywrd import check_dct1


# DCTS
SPC_REQ = ('inchi', 'mult', 'charge', 'elec_levels', 'mc_nsamp')
TS_REQ = ('rxndirn', 'kt_pst', 'temp_pst', 'n_pst')
SPC_VAL_DCT = {
    'mult': ((int,), (), None),
    'charge': ((int,), (), None),
    'inchi': ((str,), (), None),
    'inchikey': ((str,), (), None),
    'smiles': ((str,), (), None),
    'sens': ((float,), (), None),  # auto from CSV reader, not used
    'fml': ((dict,), (), None),  # auto from CSV reader, not used
    'pst_params': ((tuple,), (), (1.0, 6)),
    # ^^ shouldn't be in spc, but auto dat glob prob for all TS keys)
    'tors_names': ((tuple,), (), None),
    'elec_levels': ((tuple,), (), None),
    'geo': ((tuple,), (), None),
    'sym_factor': ((float,), (), None),
    'kickoff': ((tuple,), (), (True, (0.1, False))),
    'hind_inc': ((float,), (), 30.0),
    'hbond_cutoffs': ((tuple,), (), (4.55, 1.92)),
    'mc_nsamp': ((tuple,), (), (True, 12, 1, 3, 100, 25)),
    'tau_nsamp': ((tuple,), (), (True, 12, 1, 3, 100, 25)),
    'smin': ((float,), (), None),
    'smax': ((float,), (), None),
    'etrans_nsamp': ((int,), (), None),
    'bath': ((str,), (), None),
    'lj': ((str,), (), None),
    'edown': ((str,), (), None),
    'active': ((tuple,), (), None),
    'zma_idx': ((int,), (), 0)
}
TS_VAL_DCT = {
    'rxndirn': ((str,), (), 'forw'),
    'kt_pst': ((float,), (), 4.0e-10),
    'temp_pst': ((float,), (), 300.0),
    'n_pst': ((float,), (), 6.0),
    'active': ((str,), (), None),
    'ts_search': ((str,), (), None),
    'ts_idx': ((int,), (), 0)
    # vrc_dct = autorun.varecof.VRC_DCT
    # machine_dct = {}
}
TS_VAL_DCT.update(SPC_VAL_DCT)


# Build spc
def species_dictionary(spc_str, dat_str, geo_dct, spc_type):
    """ Read each of the species input files:
            (1) species.csv: CSV file with basic info like names,inchis,mults
            (2) species.dat:
            (3) *.xyz: XYZ-files with geometries

        :param job_path: directory path where the input file(s) exist
        :type job_path: str
        :rtype dict[str: dict]
    """

    # Parse out the dcts from the strings
    spc_dct = mechanalyzer.parser.spc.build_spc_dct(spc_str, spc_type)

    dat_blocks = ioformat.ptt.named_end_blocks(dat_str, 'spc', footer='spc')
    dat_dct = ioformat.ptt.keyword_dcts_from_blocks(dat_blocks)

    # Merge all of the species inputs into a dictionary
    mod_spc_dct, glob_dct = modify_spc_dct(spc_dct, dat_dct, geo_dct)

    # Assess if the species.dat information is valid
    for name, dct in mod_spc_dct.items():
        # last comment breaks since TS only partially built at this stage
        # i.e. there is not a mult, that comes from the build later
        # prolly fine, since we add required stuff for TS ourselves
        # probably just check if stuff provided that is not supported
        # req_lst = SPC_REQ if 'ts' not in name else SPC_REQ+TS_REQ
        req_lst = SPC_REQ if 'ts' not in name else ()
        val_dct = SPC_VAL_DCT if 'ts' not in name else TS_VAL_DCT
        check_dct1(dct, val_dct, req_lst, f'Spc-{name}')

    return mod_spc_dct, glob_dct


# Format spc
def modify_spc_dct(spc_dct, amech_dct, geo_dct):
    """ Modify the species dct using input from the additional AMech file
    """

    # Build defaults
    spc_default = defaults_from_val_dct(SPC_VAL_DCT)
    ts_default = defaults_from_val_dct(TS_VAL_DCT)

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
            spc_default, spc_dct[spc])

        # Add speciaized calls not in the default dct
        ich, mul = spc_dct[spc]['inchi'], spc_dct[spc]['mult']
        if spc_dct[spc]['elec_levels'] is None:
            spc_dct[spc]['elec_levels'] = eleclvl.DCT.get(
                (ich, mul), ((0.0, mul),))
        if spc_dct[spc]['sym_factor'] is None:
            spc_dct[spc]['sym_factor'] = symm.DCT.get(
                (ich, mul), None)

    # Add transitions states defined in species.dat not defined in spc_dct
    ts_dct = {}
    for tsname in (x for x in dat_dct if 'ts' in x):
        ts_dct[tsname] = {**dat_dct[tsname]}

        # Need to add the TS defaults
        ts_dct[tsname] = automol.util.dict_.right_update(
            ts_default, ts_dct[tsname])

        # Add speciaized calls not in the default dct
        # _set_active_key()

    # add the TSs to the spc dct
    spc_dct.update(ts_dct)

    # Final loop for conversions and additions
    for spc in spc_dct:
        spc_dct[spc]['hind_inc'] *= phycon.DEG2RAD
        spc_dct[spc]['geo'] = geo_dct.get(spc, None)

    # Perform similar conversions where needed for glob dct
    if 'hind_inc' in glob_dct:
        glob_dct['hind_inc'] *= phycon.DEG2RAD

    return spc_dct, glob_dct


def combine_sadpt_spc_dcts(sadpt_dct, spc_dct, glob_dct):
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
        if len(list(glob_dct.keys())) > 0:
            for key, val in glob_dct.items():
                combined_dct[sadpt][key] = val

        # Update any sadpt keywords if they are in the spc_dct from .dat file
        if sadpt in spc_dct:
            combined_dct[sadpt].update(spc_dct[sadpt])

        # Put in stuff from the sadpt_dct build
        for key, val in sadpt_dct[sadpt].items():
            if key not in combined_dct[sadpt]:
                combined_dct[sadpt][key] = val

        # Put in defaults if they were not defined
        # hindered rotor being set incorrectly here
        combined_dct[sadpt] = automol.util.dict_.right_update(
            defaults_from_val_dct(TS_VAL_DCT), combined_dct[sadpt])

    return combined_dct


# Functions to the spc_dct contributions for TS
def ts_dct_from_estsks(pes_idx, es_tsk_lst, rxn_lst, thy_dct,
                       spc_dct, run_prefix, save_prefix):
    """ build a ts queue
    """

    print('\nTasks for transition states requested...')
    print('Identifying reaction classes for transition states...')

    # Build the ts_dct
    ts_dct = {}
    for tsk_lst in es_tsk_lst:
        obj, es_keyword_dct = tsk_lst[0], tsk_lst[-1]
        if obj in ('ts', 'all'):
            # want print for task list
            method_dct = thy_dct.get(es_keyword_dct['runlvl'])
            ini_method_dct = thy_dct.get(es_keyword_dct['inplvl'])
            thy_info = tinfo.from_dct(method_dct)
            ini_thy_info = tinfo.from_dct(ini_method_dct)
            break

    # Discern if TS should be reidentified
    re_id = False
    for tsk_lst in es_tsk_lst:
        obj, es_keyword_dct = tsk_lst[:-1], tsk_lst[-1]
        if 'find_ts' in obj:
            re_id = es_keyword_dct.get('re_id', False)

    ts_dct = {}
    for rxn in rxn_lst:
        ts_dct.update(
            ts_dct_sing_chnl(
                pes_idx, rxn,
                spc_dct, run_prefix, save_prefix,
                thy_info=thy_info, ini_thy_info=ini_thy_info, re_id=re_id)
        )

    # Build the queue
    # ts_queue = tuple(sadpt for sadpt in ts_dct) if ts_dct else ()

    return ts_dct
    # return ts_dct, ts_queue


def ts_dct_from_ktptsks(pes_idx, rxn_lst, ktp_tsk_lst,
                        spc_model_dct,
                        spc_dct, run_prefix, save_prefix):
    """ Build ts dct from ktp tsks
    """

    for tsk_lst in ktp_tsk_lst:
        [tsk, ktp_keyword_dct] = tsk_lst
        if 'mess' in tsk or 'fit' in tsk:
            spc_model = ktp_keyword_dct['spc_model']
            ini_thy_info = spc_model_dct[spc_model]['vib']['geolvl'][1][1]
            thy_info = spc_model_dct[spc_model]['ene']['lvl1'][1][1]
            break

    ts_dct = {}
    for rxn in rxn_lst:
        ts_dct.update(
            ts_dct_sing_chnl(
                pes_idx, rxn,
                spc_dct, run_prefix, save_prefix,
                thy_info=thy_info, ini_thy_info=ini_thy_info)
        )

    return ts_dct


def ts_dct_from_proctsks(pes_idx, proc_tsk_lst, rxn_lst, spc_mod_dct_i,
                         thy_dct, spc_dct, run_prefix, save_prefix):
    """ build a ts queue
    """

    print('\nTasks for transition states requested...')
    print('Identifying reaction classes for transition states...')

    # Build the ts_dct
    ts_dct = {}
    for tsk_lst in proc_tsk_lst:
        obj, proc_keyword_dct = tsk_lst[:-1], tsk_lst[-1]
        if 'ts' in obj or 'all' in obj:
            # want print for task list
            if spc_mod_dct_i is not None:
                ini_thy_info = spc_mod_dct_i['vib']['geolvl'][1][1]
                thy_info = spc_mod_dct_i['vib']['geolvl'][1][1]
            else:
                ini_thy_info = tinfo.from_dct(thy_dct.get(
                    proc_keyword_dct['proplvl']))
                thy_info = tinfo.from_dct(thy_dct.get(
                    proc_keyword_dct['proplvl']))
            break

    ts_dct = {}
    for rxn in rxn_lst:
        ts_dct.update(
            ts_dct_sing_chnl(
                pes_idx, rxn,
                spc_dct, run_prefix, save_prefix,
                thy_info=thy_info, ini_thy_info=ini_thy_info,
                id_missing=False)
        )

    # Build the queue
    ts_queue = tuple(sadpt for sadpt in ts_dct) if ts_dct else ()

    return ts_dct, ts_queue


def ts_dct_sing_chnl(pes_idx, reaction,
                     spc_dct, run_prefix, save_prefix,
                     thy_info=None, ini_thy_info=None,
                     id_missing=True, re_id=False):
    """ build dct for single reaction
    """

    # Unpack the reaction object
    chnl_idx, (reacs, prods) = reaction

    rxn_info = rinfo.from_dct(reacs, prods, spc_dct)
    rct_str, prd_str = '+'.join(reacs), '+'.join(prods)
    print(f'  Preparing for reaction {rct_str} = {prd_str}')

    # Set the reacs and prods for the desired direction
    reacs, prods = rxnid.set_reaction_direction(
        reacs, prods, rxn_info,
        thy_info, ini_thy_info, save_prefix, direction='forw')

    # Obtain the reaction object for the reaction
    zma_locs = (0,)
    # is there a better way to get this hbond param out of spc_dct and does
    # it matter in getting the mincofs to build the reaction if we bother
    # to include it?
    # hbond_cutoffs = spc_dct[reacs[0]]['hbond_cutoffs']
    zrxns, zmas, rclasses, status = rxnid.build_reaction(
        rxn_info, ini_thy_info, zma_locs, save_prefix,
        id_missing=id_missing, re_id=re_id)
    # , hbond_cutoffs=hbond_cutoffs)

    # Could reverse the spc dct
    if status not in ('MISSING-SKIP', 'MISSING-ADD'):
        ts_dct = {}
        for idx, (zrxn, zma, cls) in enumerate(zip(zrxns, zmas, rclasses)):
            tsname = f'ts_{pes_idx+1:d}_{chnl_idx+1:d}_{idx:d}'
            ts_dct[tsname] = {
                'zrxn': zrxn,
                'zma': zma,
                'reacs': reacs,
                'prods': prods,
                'rxn_info': rxn_info,
                'inchi': '',
                'charge': rinfo.value(rxn_info, 'charge'),
                'mult': rinfo.value(rxn_info, 'tsmult'),
                'elec_levels': ((0.0, rinfo.value(rxn_info, 'tsmult')),),
                'hind_inc': 30.0*phycon.DEG2RAD,
                'hbond_cutoffs': (4.55, 1.92),
                'class': cls,
                'rxn_fs': reaction_fs(run_prefix, save_prefix, rxn_info)
            }
    elif status == 'MISSING-ADD':
        tsname = f'ts_{pes_idx+1:d}_{chnl_idx+1:d}_0'
        ts_dct = {}
        ts_dct[tsname] = {'missdata': ini_thy_info}
    elif status == 'MISSING-SKIP':
        ts_dct = {}
        print('Skipping reaction as class not given/identified')

    return ts_dct


def base_tsname(pes_idx, chnl_idx):
    """ get tsname that precludes the confiuraton number
    """
    return f'ts_{pes_idx+1:d}_{chnl_idx+1:d}'


def tsnames_in_dct(pes_idx, chnl_idx, spc_dct, config_idxs=None):
    """ Get the names of all configuratons of a transition state
         for the channel of a PES.
    """
    _tsname = f'ts_{pes_idx+1:d}_{chnl_idx+1:d}'
    _tsname = _tsname + '_'
    if config_idxs is None:
        _tsnames = tuple(name for name in spc_dct.keys()
                         if _tsname in name)
    else:
        _tsnames = tuple(f'{_tsname}{idx}' for idx in config_idxs)
    return _tsnames
