"""
Write and Read MESS files for Rates
"""

import importlib
import copy
import automol
import mess_io
from mechlib.amech_io.parser.spc import tsnames_in_dct, base_tsname
from mechlib.amech_io import printer as ioprinter
from mechlib import filesys
from mechroutines.models import blocks
from mechroutines.models import build
from mechroutines.models import etrans
from mechroutines.models import tunnel
from mechroutines.models.inf import make_rxn_str
from mechroutines.models.typ import need_fake_wells
from mechroutines.models.typ import is_abstraction_pes
from mechroutines.ktp._ene import set_reference_ene
from mechroutines.ktp._ene import sum_channel_enes

from mechroutines.ktp._newstuff import energy_dist_params
from mechroutines.ktp._newstuff import set_prod_density_param
from mechroutines.ktp._newstuff import set_hot_enes


BLOCK_MODULE = importlib.import_module('mechroutines.models.blocks')


# Headers
def make_header_str(spc_dct, rxn_lst, pes_idx, pesgrp_num,
                    pes_param_dct, hot_enes_dct, label_dct,
                    temps, pressures, float_type):
    """ Built the head of the MESS input file that contains various global
        keywords used for running rate calculations.

        Function determines certain input parameters for the well-extension
        methodology based on the reaction type stored in spc_dct.

        :param spc_dct:
        :type spc_dct: dict[]
        :param temps: temperatures for the rate calculations (in K)
        :type temps: tuple(float)
        :param pressures: pressures for the rate calculations (in atm)
        :type pressures: tuple(float)
        :rtype: str
    """

    ioprinter.messpf('global_header')

    keystr1 = (
        'EnergyStepOverTemperature, ExcessEnergyOverTemperature, ' +
        'ModelEnergyLimit'
    )
    keystr2 = (
        'CalculationMethod, WellCutoff, ' +
        'ChemicalEigenvalueMax, ReductionMethod, AtomDistanceMin'
    )
    ioprinter.debug_message(f'     {keystr1}')
    ioprinter.debug_message(f'     {keystr2}')

    # Set the well extension energy thresh
    if is_abstraction_pes(spc_dct, rxn_lst, pes_idx):
        well_extend = None
    else:
        well_extend = 'auto'
        ioprinter.debug_message('Including WellExtend in MESS input')

    # Set other parameters
    # Need the PES number to pull the correct params out of lists
    ped_spc_lst, hot_enes_dct, micro_out_params = energy_dist_params(
        pesgrp_num, pes_param_dct, hot_enes_dct, label_dct)

    header_str = mess_io.writer.global_rates_input(
        temps, pressures,
        calculation_method='direct',
        well_extension=well_extend,
        ped_spc_lst=ped_spc_lst,
        hot_enes_dct=hot_enes_dct,
        excess_ene_temp=None,
        micro_out_params=micro_out_params,
        float_type=float_type)

    return header_str


def make_global_etrans_str(rxn_lst, spc_dct, etrans_dct):
    """ Writes a string with defining global energy transfer parameters used
        for all wells on the PES that do not have parameters defined in their
        respective sections.

        As a default, the function will obtain parameters for the first well
        that appears on the PES.
    """

    ioprinter.messpf('transfer_section')

    # Determine the species for which you
    ioprinter.messpf('well_section')
    well_info = etrans.set_etrans_well(rxn_lst, spc_dct)

    # Determine the bath
    ioprinter.messpf('bath_section')
    bath_info = etrans.set_bath(spc_dct, etrans_dct)

    # Write the MESS energy transfer strings
    edown_str, collid_str = etrans.make_energy_transfer_strs(
        well_info, bath_info, etrans_dct)
    energy_trans_str = mess_io.writer.global_energy_transfer_input(
        edown_str, collid_str)

    return energy_trans_str


# Reaction Channel Writers for the PES
def make_pes_mess_str(spc_dct, rxn_lst, pes_idx, pesgrp_num,
                      unstable_chnls,
                      run_prefix, save_prefix, label_dct,
                      tsk_key_dct, pes_param_dct,
                      thy_dct, pes_model_dct_i, spc_model_dct_i,
                      spc_model):
    """ Write all the MESS input file strings for the reaction channels
    """

    ioprinter.messpf('channel_section')

    # Initialize empty MESS strings
    full_well_str, full_bi_str, full_ts_str = '', '', ''
    full_dat_str_dct = {}

    # Set the energy and model for the first reference species
    ioprinter.info_message('\nCalculating reference energy for PES')
    ref_ene = set_reference_ene(
        rxn_lst, spc_dct,
        pes_model_dct_i, spc_model_dct_i,
        run_prefix, save_prefix, ref_idx=0)

    # Loop over all the channels and write the MESS strings
    written_labels = []
    basis_energy_dct = {}
    for rxn in rxn_lst:

        chnl_idx, (reacs, prods) = rxn

        ioprinter.obj('vspace')
        ioprinter.reading('PES electrion structure data')
        ioprinter.channel(chnl_idx, reacs, prods)

        # Get the names for all of the configurations of the TS
        tsname = base_tsname(pes_idx, chnl_idx)
        tsname_allconfigs = tsnames_in_dct(pes_idx, chnl_idx, spc_dct)

        # Obtain all of the species data
        if spc_model not in basis_energy_dct:
            basis_energy_dct[spc_model] = {}

        # Pass in full ts class
        chnl_infs, chn_basis_ene_dct = get_channel_data(
            reacs, prods, tsname_allconfigs,
            spc_dct, tsk_key_dct,
            basis_energy_dct[spc_model],
            thy_dct, pes_model_dct_i, spc_model_dct_i,
            run_prefix, save_prefix)

        basis_energy_dct[spc_model].update(chn_basis_ene_dct)

        # Calculate the relative energies of all spc on the channel
        chnl_enes = sum_channel_enes(chnl_infs, ref_ene)

        # Set the hot energies using the relative enes that will be
        # written into the global key section of MESS input later
        hot_enes_dct = set_hot_enes(pesgrp_num, reacs, prods,
                                    chnl_enes, pes_param_dct,
                                    ene_range=(10.,))

        # Write the mess strings for all spc on the channel
        mess_strs, dat_str_dct, written_labels = _make_channel_mess_strs(
            tsname, reacs, prods, pesgrp_num,
            spc_dct, label_dct, written_labels,
            pes_param_dct, chnl_infs, chnl_enes, spc_model_dct_i,
            unstable_chnl=(chnl_idx in unstable_chnls))

        # Append to full MESS strings
        [well_str, bi_str, ts_str] = mess_strs
        full_well_str += well_str
        full_bi_str += bi_str
        full_ts_str += ts_str
        full_dat_str_dct.update(dat_str_dct)

        ioprinter.debug_message('rxn', rxn)
        ioprinter.debug_message('enes', chnl_enes)
        ioprinter.debug_message('label dct', label_dct)
        ioprinter.debug_message('written labels', written_labels)

    # Combine all the reaction channel strings
    rxn_chan_str = '\n'.join([full_well_str, full_bi_str, full_ts_str])

    return rxn_chan_str, full_dat_str_dct, hot_enes_dct


def _make_channel_mess_strs(tsname, reacs, prods, pesgrp_num,
                            spc_dct, label_dct, written_labels,
                            pes_param_dct, chnl_infs, chnl_enes,
                            spc_model_dct_i,
                            unstable_chnl=False):
    """ For each reaction channel on the PES: take all of the pre-read and
        pre-processed information from the save filesys for the
        reactants, products, and transition state and write the appropriately
        formatted MESS input strings that will eventually be combined into the
        entire MESS input file.

        Also returns dictionary for all additional auxiliary data files,
        formatted as {file name: file string}, required by MESS.

        List of labels corresponding to MESS strings that have already been
        written and added to master string, meaning that species string does
        not need to be written again. Required since species appear on multiple
        channels.

        :param tsname: mechanism name of the transition state
        :param reacs: mechanisms name for the reactants of the reaction channel
        :type reacs: tuple(str)
        :param prods: mechanisms name for the products of the reaction channel
        :type prods: tuple(str)
        :param label_dct: mapping between mechanism name and MESS input label
        :type label_dct: dict[str: str]
        :param written_labels:
        :type written_labels:
        :param chnl_infs: collated molecular info obtained from save filesys
        :type chnl_infs: dict[str:__]
        :param chnl_enes: energies for channel, relative to PES reference
        :type chnl_enes: dict[str:float]
        :rtype: (str, str, str), str, dict[str:str]

    """

    # Initialize empty strings
    bi_str, well_str, ts_str = '', '', ''
    full_dat_dct = {}

    # Write the MESS string for the channel reactant(s) and product(s)
    for side in ('reacs', 'prods'):

        # Get information from relevant dictionaries
        rgt_names = reacs if side == 'reacs' else prods
        rgt_infs = chnl_infs[side]
        rgt_ene = chnl_enes[side]

        # Build the species string for reactant(s)/product(s)
        # Skip molec string building for termolecular species (may need agn)
        spc_strs = []
        if len(rgt_names) < 3:
            for inf in rgt_infs:
                spc_str, dat_dct = _make_spc_mess_str(inf)
                spc_strs.append(spc_str)
                full_dat_dct.update(dat_dct)

        # Set the labels to put into the file
        spc_label = [automol.inchi.smiles(spc_dct[name]['inchi'])
                     for name in rgt_names]
        _rxn_str = make_rxn_str(rgt_names)
        _rxn_str_rev = make_rxn_str(rgt_names[::-1])
        if _rxn_str in label_dct:
            chn_label = label_dct[_rxn_str]
        elif _rxn_str_rev in label_dct:
            chn_label = label_dct[_rxn_str_rev]
        else:
            ioprinter.warning_message(f'no {_rxn_str} in label dct')

        # Write the strings
        if chn_label not in written_labels:
            written_labels.append(chn_label)
            if len(rgt_names) == 3:
                bi_str += (
                    f'\n! {rgt_names[0]} + {rgt_names[1]} + {rgt_names[2]}\n')
                bi_str += mess_io.writer.dummy(chn_label, zero_ene=rgt_ene)
                # bi_str += '\n! DUMMY FOR UNSTABLE SPECIES\n'
                # bi_str += mess_io.writer.dummy(chn_label, zero_ene=None)
            elif len(rgt_names) == 2:
                # Determine if product densities should be calc'd
                if side == 'prods':
                    calc_dens = set_prod_density_param(
                        rgt_names, pesgrp_num, pes_param_dct)
                else:
                    calc_dens = (False, False)

                # bi_str += mess_io.writer.species_separation_str()
                bi_str += f'\n! {rgt_names[0]} + {rgt_names[1]}\n'
                bi_str += mess_io.writer.bimolecular(
                    chn_label, spc_label[0], spc_strs[0],
                    spc_label[1], spc_strs[1], rgt_ene,
                    calc_spc1_density=calc_dens[0],
                    calc_spc2_density=calc_dens[1])
            else:
                edown_str = rgt_infs[0].get('edown_str', None)
                collid_freq_str = rgt_infs[0].get('collid_freq_str', None)

                # well_str += mess_io.writer.species_separation_str()
                well_str += f'\n! {rgt_names[0]}\n'
                well_str += mess_io.writer.well(
                    chn_label, spc_strs[0],
                    zero_ene=rgt_ene,
                    edown_str=edown_str,
                    collid_freq_str=collid_freq_str)

        # Initialize the reactant and product MESS label
        if side == 'reacs':
            reac_label = chn_label
            inner_reac_label = chn_label
        else:
            prod_label = chn_label
            inner_prod_label = chn_label

    # For abstractions: Write MESS strings for fake reac and prod wells and TS
    if chnl_infs.get('fake_vdwr', None) is not None:

        # Write all the MESS Strings for Fake Wells and TSs
        fwell_str, fts_str, fake_lbl, fake_dct = _make_fake_mess_strs(
            (reacs, prods), 'reacs', chnl_infs['fake_vdwr'],
            chnl_enes, label_dct, reac_label)

        # Append the fake strings to overall strings
        well_str += fwell_str
        ts_str += fts_str

        # Re-set the reactant label for the inner transition state
        inner_reac_label = fake_lbl

        # Update the data string dct if necessary
        full_dat_dct.update(fake_dct)

    if chnl_infs.get('fake_vdwp', None) is not None:

        # Write all the MESS Strings for Fake Wells and TSs
        fwell_str, fts_str, fake_lbl, fake_dct = _make_fake_mess_strs(
            (reacs, prods), 'prods', chnl_infs['fake_vdwp'],
            chnl_enes, label_dct, prod_label)

        # Append the fake strings to overall strings
        well_str += fwell_str
        ts_str += fts_str

        # Reset the product labels for the inner transition state
        inner_prod_label = fake_lbl

        # Update the data string dct if necessary
        full_dat_dct.update(fake_dct)

    # Write MESS string for the inner transition state; append full
    # Label has to correspond only to base name (ignores configuration)
    ts_label = label_dct[tsname]
    rclass = spc_dct[tsname+'_0']['class']
    sts_str, ts_dat_dct = _make_ts_mess_str(
        chnl_infs, chnl_enes, spc_model_dct_i, rclass,
        ts_label, inner_reac_label, inner_prod_label,
        unstable_chnl=unstable_chnl)
    ts_str += sts_str
    full_dat_dct.update(ts_dat_dct)

    return [well_str, bi_str, ts_str], full_dat_dct, written_labels


def _make_spc_mess_str(inf_dct):
    """  Writes all processed save filesys data for a species and
         into an appropriately formatted MESS input string. Takes the
         pre-identified writer designation and calls the approprate
         MESS-block writer function in models/build module.

         :param inf_dct: save filesys data for species
         :type inf_dct: dict[]
         :rtype: str
    """
    mess_writer = getattr(BLOCK_MODULE, inf_dct['writer'])
    return mess_writer(inf_dct)


def _make_ts_mess_str(chnl_infs, chnl_enes, spc_model_dct_i, ts_class,
                      ts_label, inner_reac_label, inner_prod_label,
                      unstable_chnl=False):
    """  Writes all processed save filesys data for a transition state and
         into an appropriately formatted MESS input string. Takes the
         pre-identified writer designation and calls the approprate
         MESS-block writer function in models/build module.

        ^ slightly off, maybe add additional block function for variational,
        union, sadpt writing...

         Prior to writing, function does some additional data processing
         to write additional flux files and tunneling file strings for
         the input transition state.

         :param inf_dct: save filesys data for species
         :type inf_dct: dict[]
         :rtype: str
    """

    # Unpack info objects
    ts_mod = spc_model_dct_i['ts']

    # Write the initial data string and dat str dct with mdhr str
    mess_strs = []
    tunnel_strs = []
    ts_dat_dct = {}
    for idx, ts_inf_dct in enumerate(chnl_infs['ts']):

        # Build initial data block
        mstr, mdhr_dat, flux_dat = blocks.barrier_dat_block(
            ts_inf_dct, chnl_infs['reacs'], chnl_infs['prods'])

        # Write the appropriate string for the tunneling model
        tunnel_str, sct_dct = tunnel.write_mess_tunnel_str(
            ts_inf_dct, chnl_enes, ts_mod, ts_class, idx,
            unstable_chnl=unstable_chnl)

        # Update master TS list
        mess_strs.append(mstr)
        tunnel_strs.append(tunnel_str)
        if mdhr_dat:
            ts_dat_dct.update(mdhr_dat)
        if flux_dat:
            ts_dat_dct.update(flux_dat)
        if sct_dct:
            ts_dat_dct.update(sct_dct)

    # Place intermediate sadpt/rpath data into a MESS Barrier Block
    if len(mess_strs) == 1:
        mess_str = mess_strs[0]

        ts_sadpt, ts_nobar = ts_mod['sadpt'], ts_mod['nobar']
        radrad = bool('radical radical' in ts_class)
        write_ts_pt_str = bool(
            (not radrad and ts_sadpt != 'rpvtst') or
            (radrad and ts_nobar != 'rpvtst')
        )
        if write_ts_pt_str:
            ts_ene = chnl_enes['ts'][0]
            ts_str = '\n' + mess_io.writer.ts_sadpt(
                ts_label, inner_reac_label, inner_prod_label,
                mess_str, ts_ene, tunnel_str)
        else:
            ts_enes = chnl_enes['ts']
            ts_str = '\n' + mess_io.writer.ts_variational(
                ts_label, inner_reac_label, inner_prod_label,
                mess_str, ts_enes, tunnel_str)
    else:
        ts_enes = chnl_enes['ts']
        mess_str = mess_io.writer.configs_union(
            mess_strs, ts_enes, tunnel_strs=tunnel_strs)

        ts_str = '\n' + mess_io.writer.ts_sadpt(
            ts_label, inner_reac_label, inner_prod_label, mess_str)

    return ts_str, ts_dat_dct


def _make_fake_mess_strs(chnl, side, fake_inf_dcts,
                         chnl_enes, label_dct, side_label):
    """ write the MESS strings for the fake wells and TSs
    """

    # Set vars based on the reacs/prods
    reacs, prods = chnl
    if side == 'reacs':
        well_key = 'fake_vdwr'
        ts_key = 'fake_vdwr_ts'
        prepend_key = 'FRB'
        side_idx = 0
    elif side == 'prods':
        well_key = 'fake_vdwp'
        ts_key = 'fake_vdwp_ts'
        side_idx = 1
        if reacs in (prods, prods[::-1]):
            prepend_key = 'FRB'
        else:
            prepend_key = 'FPB'

    # Initialize well and ts strs and data dcts
    fake_dat_dct = {}
    well_str, ts_str = '', ''

    # Build a fake TS dct
    ts_inf_dct = {
        'n_pst': 6.0,
        'cn_pst': 10.0
    }

    # MESS string for the fake reactant side well
    well_dct_key = make_rxn_str(chnl[side_idx], prepend='F')
    well_dct_key_rev = make_rxn_str(chnl[side_idx][::-1], prepend='F')
    if well_dct_key in label_dct:
        fake_well_label = label_dct[well_dct_key]
    elif well_dct_key_rev in label_dct:
        fake_well_label = label_dct[well_dct_key_rev]
    else:
        ioprinter.warning_message(f'No label {well_dct_key} in label dict')
    # well_str += mess_io.writer.species_separation_str()
    _side_str = '+'.join(chnl[side_idx])
    well_str += f'\n! Fake Well for {_side_str}\n'
    fake_well, well_dat = blocks.fake_species_block(*fake_inf_dcts)
    well_str += mess_io.writer.well(
        fake_well_label, fake_well, chnl_enes[well_key])

    # MESS PST TS string for fake reactant side well -> reacs
    pst_dct_key = make_rxn_str(chnl[side_idx], prepend=prepend_key)
    pst_dct_key_rev = make_rxn_str(chnl[side_idx][::-1], prepend=prepend_key)
    if pst_dct_key in label_dct:
        pst_label = label_dct[pst_dct_key]
    elif pst_dct_key_rev in label_dct:
        pst_label = label_dct[pst_dct_key_rev]
    else:
        ioprinter.warning_message(f'No label {pst_dct_key} in label dict')
    pst_ts_str, pst_ts_dat = blocks.pst_block(ts_inf_dct, *fake_inf_dcts)
    ts_str += '\n' + mess_io.writer.ts_sadpt(
        pst_label, side_label, fake_well_label, pst_ts_str,
        chnl_enes[ts_key], tunnel='')

    # Build the data dct
    if well_dat:
        fake_dat_dct.update(well_dat)
    if pst_ts_dat:
        fake_dat_dct.update(pst_ts_dat)

    return well_str, ts_str, fake_well_label, fake_dat_dct


# Data Retriever Functions
def get_channel_data(reacs, prods, tsname_allconfigs,
                     spc_dct, tsk_key_dct,
                     model_basis_energy_dct,
                     thy_dct, pes_model_dct_i, spc_model_dct_i,
                     run_prefix, save_prefix):
    """ For all species and transition state for the channel and
        read all required data from the save filesys, then process and
        format it to be able to write it into a MESS filesystem.

        :param tsname: mechanism name of the transition state
        :param reacs: mechanisms name for the reactants of the reaction channel
        :type reacs: tuple(str)
        :param prods: mechanisms name for the products of the reaction channel
        :type prods: tuple(str)
    """

    # Initialize the dict
    chnl_infs = {}

    # Get the data for conformer sorting for reading the filesystem
    cnf_range = tsk_key_dct['cnf_range']
    sort_info_lst = filesys.mincnf.sort_info_lst(tsk_key_dct['sort'], thy_dct)

    # Determine the MESS data for the reactants and products
    # Gather data or set fake information for dummy reactants/products
    chnl_infs['reacs'], chnl_infs['prods'] = [], []
    for rgts, side in zip((reacs, prods), ('reacs', 'prods')):
        for rgt in rgts:
            spc_locs_lst = filesys.models.get_spc_locs_lst(
                spc_dct[rgt], spc_model_dct_i,
                run_prefix, save_prefix, saddle=False,
                cnf_range=cnf_range, sort_info_lst=sort_info_lst,
                name=rgt)
            chnl_infs_i, model_basis_energy_dct = build.read_spc_data(
                spc_dct, rgt,
                pes_model_dct_i, spc_model_dct_i,
                run_prefix, save_prefix, model_basis_energy_dct,
                spc_locs=spc_locs_lst[0])
            chnl_infs[side].append(chnl_infs_i)

    # Get data for all configurations for a TS
    chnl_infs['ts'] = []
    for name in tsname_allconfigs:
        spc_locs_lst = filesys.models.get_spc_locs_lst(
            spc_dct[name], spc_model_dct_i,
            run_prefix, save_prefix, saddle=True,
            cnf_range=cnf_range, sort_info_lst=sort_info_lst,
            name=name)
        inf_dct, model_basis_energy_dct = build.read_ts_data(
            spc_dct, name, reacs, prods,
            pes_model_dct_i, spc_model_dct_i,
            run_prefix, save_prefix, model_basis_energy_dct,
            spc_locs=spc_locs_lst[0])
        chnl_infs['ts'].append(inf_dct)

    # Set up the info for the wells
    rwell_model = spc_model_dct_i['ts']['rwells']
    pwell_model = spc_model_dct_i['ts']['pwells']
    rxn_class = spc_dct[tsname_allconfigs[0]]['class']
    if need_fake_wells(rxn_class, rwell_model):
        chnl_infs['fake_vdwr'] = copy.deepcopy(chnl_infs['reacs'])
    if need_fake_wells(rxn_class, pwell_model):
        chnl_infs['fake_vdwp'] = copy.deepcopy(chnl_infs['prods'])

    return chnl_infs, model_basis_energy_dct
