"""
Write and Read MESS files for Rates
"""

import importlib
import copy
import ioformat
import automol
import autorun
import mess_io
from mechlib.amech_io.parser.spc import tsnames_in_dct, base_tsname
from mechlib.amech_io import reader
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

from mechroutines.ktp._multipes import energy_dist_params
from mechroutines.ktp._multipes import set_prod_density_param
from mechroutines.ktp._multipes import set_hot_enes

BLOCK_MODULE = importlib.import_module('mechroutines.models.blocks')


# Create full string by writing the appropriate header, accounting for
# (1) MESS Version and (2) Use of Well-Extension
# And include the global_etrans and reaction channel strings
def make_full_str(energy_trans_str, rxn_chan_str, dats,
                  pesgrp_num, pes_param_dct, hot_enes_dct,
                  rate_paths_dct, pes_inf,
                  pes_mod_dct_i,
                  spc_dct, rxn_lst, pes_idx, tsk_key_dct):
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

    # Pull from PES model dct
    temps, pressures = pes_mod_dct_i['rate_temps'], pes_mod_dct_i['pressures']
    float_type = tsk_key_dct['float_precision']

    # Set other parameters
    # Need the PES number to pull the correct params out of lists
    ped_spc_lst, micro_out_params, pes_param_dct = energy_dist_params(
        pesgrp_num, pes_param_dct, hot_enes_dct, rxn_chan_str)
    
    ioprinter.messpf('global_header')

    # Write the header string
    if tsk_key_dct['mess_version'] == 'v1':
        _full_mess_v1(
            energy_trans_str, rxn_chan_str, dats,
            temps, pressures,
            ped_spc_lst, hot_enes_dct,
            micro_out_params,
            float_type,
            pes_mod_dct_i, spc_dct,
            rate_paths_dct, pes_inf,
            rxn_lst, pes_idx, tsk_key_dct)
    else:
        _full_mess_v2(
            energy_trans_str, rxn_chan_str, dats,
            temps, pressures,
            ped_spc_lst, hot_enes_dct,
            micro_out_params,
            float_type,
            pes_mod_dct_i, spc_dct,
            rate_paths_dct, pes_inf,
            rxn_lst, pes_idx)

    return pes_param_dct


def _full_mess_v1(energy_trans_str, rxn_chan_str, dats,
                  temps, pressures,
                  ped_spc_lst, hot_enes_dct,
                  micro_out_params,
                  float_type,
                  pes_mod_dct_i, spc_dct,
                  rate_paths_dct, pes_inf,
                  rxn_lst, pes_idx, tsk_key_dct):
    """ Make the global header string for MESS version 1

        last line of arguments only used to determine well-extension
    """

    ioprinter.debug_message(
        'EnergyStepOverTemperature, ExcessEnergyOverTemperature, ' +
        'ModelEnergyLimit')
    ioprinter.debug_message(
        'CalculationMethod, WellCutoff, ' +
        'ReductionMethod, AtomDistanceMin')

    if tsk_key_dct['well_lumping'] and not tsk_key_dct['well_extension']:
        print('*Warning: well lumping was selected, but well extension inactive:')
        print('Well Extension activated')
        tsk_key_dct['well_extension'] = True
        
    is_abstraction = is_abstraction_pes(spc_dct, rxn_lst, pes_idx)
    if is_abstraction and tsk_key_dct['well_extension']:
        well_extend = None # overwrite
    elif not is_abstraction and tsk_key_dct['well_extension']:
        well_extend = 0.001
    else:
        well_extend = None

    globkey_str = mess_io.writer.global_rates_input_v1(
        temps, pressures,
        calculation_method='direct',
        model_ene_limit=800.0,
        ene_stepover_temp=0.2, excess_ene_temp=None,
        well_extension=well_extend,
        well_reduction_thresh=10.0,
        chem_eig_max=0.2,
        ped_spc_lst=ped_spc_lst,
        hot_enes_dct=hot_enes_dct,
        micro_out_params=micro_out_params,
        float_type=float_type,
        ktp_outname='rate.out',
        ke_outname='ke.out',
        ped_outname='ped.out',
    )

    # Write base MESS input string into the RUN filesystem
    mess_inp_str = mess_io.writer.messrates_inp_str(
        globkey_str, rxn_chan_str,
        energy_trans_str=energy_trans_str,
        well_lump_str=None,
        use_short_names=True)
    # comment line of wellext
    mess_inp_str = mess_inp_str.replace('WellExtension', '!WellExtension')
    
    # print('rate_paths_dct test\n', rate_paths_dct)
    base_mess_path = rate_paths_dct[pes_inf]['base-v1']
    ioprinter.obj('line_plus')
    ioprinter.writing('MESS input file', base_mess_path)
    ioprinter.debug_message('MESS Input:\n\n'+mess_inp_str)
    autorun.write_input(
        base_mess_path, mess_inp_str,
        aux_dct=dats, input_name='mess.inp')


    # Write the second MESS string (well extended), if needed
    if not is_abstraction and tsk_key_dct['well_extension']:
        print('User requested well extension scheme for rates...')

        # Run the base MESSRATE
        print(f'  - Running MESS base job at path {base_mess_path}')
        autorun.run_script(autorun.SCRIPT_DCT['messrate-v1'], base_mess_path)

        # Write the well-extended MESSRATE file
        rate_strs_dct, mess_paths_dct = reader.mess.rate_strings(
            rate_paths_dct)
        read_mess_path = mess_paths_dct[pes_inf]['base-v1']
        print('  - Reading input and output from '
              f'base MESSRATE job at {read_mess_path}')

        wext_p = pes_mod_dct_i['well_extension_pressure']
        wext_t = pes_mod_dct_i['well_extension_temp']

        print('  - Setting up the well-extended MESSRATE input')
        wext_mess_inp_str_nolump = mess_io.well_lumped_input_file(
            rate_strs_dct[pes_inf]['base-v1']['inp'],
            rate_strs_dct[pes_inf]['base-v1']['ktp_out'],
            rate_strs_dct[pes_inf]['base-v1']['aux'],
            rate_strs_dct[pes_inf]['base-v1']['log'],
            wext_p,
            wext_t, lump = False)

        base_mess_path = rate_paths_dct[pes_inf]['base-v1']
        ioprinter.obj('line_plus')
        ioprinter.writing('New Well-Extended MESS input file '
                          f'at path {base_mess_path}')
        print('  - Warning, old base input overwritten.')
        ioprinter.debug_message('MESS Input:\n\n'+wext_mess_inp_str_nolump)
        autorun.write_input(
            base_mess_path, wext_mess_inp_str_nolump,
            aux_dct=dats, input_name='mess.inp')
        
        print(f'  - Running MESS base job at path {base_mess_path}')
        print('  - Warning, old base results overwritten.')
        autorun.run_script(autorun.SCRIPT_DCT['messrate-v1'], base_mess_path)
        
        if tsk_key_dct['well_lumping']:
            print('  - Setting up the well-Lumped MESSRATE input with')
            print(f'   lumping/extension Scheme for P={wext_p} atm, T={wext_t} K')
            wext_mess_inp_str = mess_io.well_lumped_input_file(
                rate_strs_dct[pes_inf]['base-v1']['inp'],
                rate_strs_dct[pes_inf]['base-v1']['ktp_out'],
                rate_strs_dct[pes_inf]['base-v1']['aux'],
                rate_strs_dct[pes_inf]['base-v1']['log'],
                wext_p,
                wext_t)

            wext_mess_path = rate_paths_dct[pes_inf]['wext-v1']
            ioprinter.obj('line_plus')
            ioprinter.writing('New Well-Lumped MESS input file '
                            f'at path {wext_mess_path}')
            ioprinter.debug_message('MESS Input:\n\n'+wext_mess_inp_str)
            autorun.write_input(
                wext_mess_path, wext_mess_inp_str,
                aux_dct=dats, input_name='mess.inp')

def _full_mess_v2(energy_trans_str, rxn_chan_str, dats,
                  temps, pressures,
                  ped_spc_lst, hot_enes_dct,
                  micro_out_params,
                  float_type,
                  pes_mod_dct_i, spc_dct,
                  rate_paths_dct, pes_inf,
                  rxn_lst, pes_idx):
    """ Make the global header string for MESS version 2
    """

    if is_abstraction_pes(spc_dct, rxn_lst, pes_idx):
        well_extend = None
    else:
        well_extend = 0.001
        ioprinter.debug_message('Including WellExtend in MESS input')

    globkey_str = mess_io.writer.global_rates_input_v2(
            temps, pressures,
            ref_temperature=pes_mod_dct_i['well_extension_temp'],
            ref_pressure=pes_mod_dct_i['well_extension_pressure'],
            model_ene_limit=800.0,
            ene_stepover_temp=0.2, ene_cutoff_temp=20.0, excess_ene_temp=10.0,
            chem_tol=1.0e-10, chem_thresh=0.1,
            well_pojection_thresh=0.1, well_reduction_thresh=10.0,
            time_propagation_limit=50.0, time_propagation_step=0.02,
            well_extension=well_extend,
            ped_spc_lst=ped_spc_lst, hot_enes_dct=hot_enes_dct,
            micro_out_params=micro_out_params,
            float_type=float_type,
            ktp_outname='rate.out',
            ke_outname='ke.out',
            ped_outname='ped.out',
    )

    # Write base MESS input string into the RUN filesystem
    mess_inp_str = mess_io.writer.messrates_inp_str(
        globkey_str, rxn_chan_str,
        energy_trans_str=energy_trans_str,
        well_lump_str=None,
        use_short_names=True)

    base_mess_path = rate_paths_dct[pes_inf]['base-v2']
    ioprinter.obj('line_plus')
    ioprinter.writing('MESS input file', base_mess_path)
    ioprinter.debug_message('MESS Input:\n\n'+mess_inp_str)
    autorun.write_input(
        base_mess_path, mess_inp_str,
        aux_dct=dats, input_name='mess.inp')


def make_global_etrans_str(rxn_lst, spc_dct, etrans_dct):
    """ Writes a string with defining global energy transfer parameters used
        for all wells on the PES that do not have parameters defined in their
        respective sections.

        As a default, the function will obtain parameters for the first well
        that appears on the PES.
    """

    ioprinter.messpf('global_transfer_section')

    # Determine the representative well and bath for global keys
    ioprinter.messpf('well_section')
    well_info = etrans.set_etrans_well(rxn_lst, spc_dct)
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
                      spc_model, nprocs=1):
    """ Write all the MESS input file strings for the reaction channels
    """

    ioprinter.messpf('channel_section')

    # Initialize data carrying objects and empty MESS strings
    basis_energy_dct = {}
    basis_energy_dct[spc_model] = {}

    full_well_str, full_bi_str, full_ts_str = '', '', ''
    full_dat_str_dct = {}

    # Set the energy and model for the first reference species
    ioprinter.info_message('\nCalculating reference energy for PES')
    ref_ene, model_basis_energy_dct = set_reference_ene(
        rxn_lst, spc_dct, tsk_key_dct,
        basis_energy_dct[spc_model],
        thy_dct, pes_model_dct_i, spc_model_dct_i,
        run_prefix, save_prefix, ref_idx=0, nprocs=nprocs)
    basis_energy_dct[spc_model].update(model_basis_energy_dct)

    # Loop over all the channels and write the MESS strings
    written_labels = []
    hot_enes_dct = {}
    for rxn in rxn_lst:

        chnl_idx, (reacs, prods) = rxn

        ioprinter.obj('vspace')
        ioprinter.reading('PES electronic structure data')
        ioprinter.channel(chnl_idx+1, reacs, prods)

        # Get the names for all of the configurations of the TS
        tsname = base_tsname(pes_idx, chnl_idx)
        tsname_allconfigs = tsnames_in_dct(pes_idx, chnl_idx, spc_dct)
        chnl_infs, chn_basis_ene_dct = get_channel_data(
            reacs, prods, tsname_allconfigs,
            spc_dct, tsk_key_dct,
            basis_energy_dct[spc_model],
            thy_dct, pes_model_dct_i, spc_model_dct_i,
            run_prefix, save_prefix, nprocs=nprocs)

        basis_energy_dct[spc_model].update(chn_basis_ene_dct)

        # Calculate the relative energies of all spc on the channel
        chnl_enes = sum_channel_enes(chnl_infs, ref_ene)

        # Set the hot energies using the relative enes that will be
        # written into the global key section of MESS input later
        hot_enes_dct = set_hot_enes(hot_enes_dct, pesgrp_num, reacs, prods,
                                    chnl_enes, pes_param_dct)

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

    # Combine all the reaction channel strings; remove empty lines
    rxn_chan_str = '\n'.join([full_well_str, full_bi_str, full_ts_str])
    rxn_chan_str = ioformat.remove_empty_lines(rxn_chan_str)

    if not hot_enes_dct:
        hot_enes_dct = None

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

        # Generate auxiliary labels corresponding to SMILES for quick IDs
        aux_labels = tuple(automol.chi.smiles(spc_dct[name]['inchi'])
                           for name in rgt_names)

        # old MESS channel labels system
        # Set the labels to put into the file
        # spc_labels = ()
        # for name in rgt_names:
        #     if name in label_dct:
        #         spc_labels += (label_dct[name],)
        #     else:
        #         spc_labels += (name,)
        #
        # _rxn_str = make_rxn_str(rgt_names)
        # _rxn_str_rev = make_rxn_str(rgt_names[::-1])
        # if _rxn_str in label_dct:
        #     chn_label = label_dct[_rxn_str]
        # elif _rxn_str_rev in label_dct:
        #     chn_label = label_dct[_rxn_str_rev]
        # else:
        #     ioprinter.warning_message(f'no {_rxn_str} in label dct')

        # new MESS channel labels system
        spc_labels = rgt_names+tuple()

        # always write as A+B
        _rxn_str = make_rxn_str(rgt_names)
        _rxn_str_rev = make_rxn_str(rgt_names[::-1])
        if _rxn_str_rev in written_labels:
            chn_label = _rxn_str_rev
        else:
            chn_label = _rxn_str

        if any(lbl in written_labels for lbl in (_rxn_str, _rxn_str_rev)):
            write_string = False
        else:
            write_string = True

        # Write the strings
        if write_string:
            # Append unwritten label to master list for future loops

            # Write appropriate string for Dummy, Bimol, Well
            written_labels.append(chn_label)
            if len(rgt_names) == 3:
                aux_str = (
                    f'[{aux_labels[0]} + {aux_labels[1]} + {aux_labels[2]}]'
                )
                bi_str += mess_io.writer.dummy(
                    chn_label,
                    aux_id_label=aux_str,
                    zero_ene=rgt_ene)
            elif len(rgt_names) == 2:
                # Determine if product densities should be calc'd
                if side == 'prods':
                    calc_dens = set_prod_density_param(
                        rgt_names, pesgrp_num, pes_param_dct)
                else:
                    calc_dens = (False, False)

                aux_str = (
                    f'[{aux_labels[0]} + {aux_labels[1]}]'
                )
                bi_str += mess_io.writer.bimolecular(
                    chn_label,
                    spc_labels[0], spc_strs[0],
                    spc_labels[1], spc_strs[1],
                    rgt_ene,
                    auxbimol_id_label=aux_str,
                    aux1_id_label=aux_labels[0],
                    aux2_id_label=aux_labels[1],
                    calc_spc1_density=calc_dens[0],
                    calc_spc2_density=calc_dens[1]) + '\n'
            else:
                edown_str = rgt_infs[0].get('edown_str', None)
                collid_freq_str = rgt_infs[0].get('collid_freq_str', None)

                aux_str = (
                    f'[{aux_labels[0]}]'
                )
                well_str += mess_io.writer.well(
                    chn_label, spc_strs[0],
                    aux_id_label=aux_str,
                    zero_ene=rgt_ene,
                    edown_str=edown_str,
                    collid_freq_str=collid_freq_str) + '\n'

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
            tsname, (reacs, prods), 'reacs', chnl_infs['fake_vdwr'],
            chnl_enes, reac_label)

        # Append the fake strings to overall strings
        if fake_lbl not in written_labels:
            well_str += fwell_str + '\n'
            ts_str += fts_str
            written_labels.append(fake_lbl)

        # Re-set the reactant label for the inner transition state
        inner_reac_label = fake_lbl

        # Update the data string dct if necessary
        full_dat_dct.update(fake_dct)

    if chnl_infs.get('fake_vdwp', None) is not None:

        # Write all the MESS Strings for Fake Wells and TSs
        fwell_str, fts_str, fake_lbl, fake_dct = _make_fake_mess_strs(
            tsname, (reacs, prods), 'prods', chnl_infs['fake_vdwp'],
            chnl_enes, prod_label)

        # Append the fake strings to overall strings
        if fake_lbl not in written_labels:
            well_str += fwell_str + '\n'
            ts_str += fts_str
            written_labels.append(fake_lbl)

        # Reset the product labels for the inner transition state
        inner_prod_label = fake_lbl

        # Update the data string dct if necessary
        full_dat_dct.update(fake_dct)

    # Write MESS string for the inner transition state; append full
    # Label has to correspond only to base name (ignores configuration)
    # ts_label = label_dct[tsname]
    ts_label = tsname   # change MESS labels
    rclass = spc_dct[tsname+'_0']['class']
    sts_str, ts_dat_dct = _make_ts_mess_str(
        chnl_infs, chnl_enes, spc_model_dct_i, rclass,
        ts_label, inner_reac_label, inner_prod_label,
        unstable_chnl=unstable_chnl)
    ts_str += sts_str
    full_dat_dct.update(ts_dat_dct)

    return (
        (well_str, bi_str, ts_str),
        full_dat_dct,
        written_labels
    )


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

        # If writer was previously overwritten with something else, use instead
        if chnl_infs['ts'][0]['writer'] != 'rpvtst_block':
            ts_ene = chnl_enes['ts'][0]
            ts_str = '\n' + mess_io.writer.ts_sadpt(
                ts_label, inner_reac_label, inner_prod_label, mess_str,
                aux_id_label=None,
                zero_ene=ts_ene,
                tunnel=tunnel_str)
        else:
            ts_enes = chnl_enes['ts']
            ts_str = '\n' + mess_io.writer.ts_variational(
                ts_label, inner_reac_label, inner_prod_label, mess_str,
                aux_id_label=None,
                zero_enes=ts_enes,
                tunnel=tunnel_str)
    else:
        ts_enes = chnl_enes['ts']
        mess_str = mess_io.writer.configs_union(
            mess_strs, ts_enes, tunnel_strs=tunnel_strs)

        ts_str = '\n' + mess_io.writer.ts_sadpt(
            ts_label, inner_reac_label, inner_prod_label, mess_str,
            aux_id_label=None,
            zero_ene=None, tunnel='')

    return ts_str, ts_dat_dct


def _make_fake_mess_strs(tsname, chnl, side, fake_inf_dcts,
                         chnl_enes, side_label):
    """ write the MESS strings for the fake wells and TSs
    """

    # Set vars based on the reacs/prods
    reacs, prods = chnl
    if side == 'reacs':
        well_key = 'fake_vdwr'
        ts_key = 'fake_vdwr_ts'
        #prepend_key = 'FakeRB'
        side_idx = 0
    elif side == 'prods':
        well_key = 'fake_vdwp'
        ts_key = 'fake_vdwp_ts'
        side_idx = 1
        #if reacs in (prods, prods[::-1]):
        #    prepend_key = 'FakeRB'
        #else:
        #    prepend_key = 'FakePB'

    # Initialize well and ts strs and data dcts
    fake_dat_dct = {}
    well_str, ts_str = '', ''

    # Build a fake TS dct
    ts_inf_dct = {
        'n_pst': 6.0,
        'cn_pst': 10.0
    }

    # MESS string for the fake reactant side well

    # Old MESS label code
    # well_dct_key = make_rxn_str(chnl[side_idx], prepend='F')
    # well_dct_key_rev = make_rxn_str(chnl[side_idx][::-1], prepend='F')
    # if well_dct_key in label_dct:
    #     fake_well_label = label_dct[well_dct_key]
    # elif well_dct_key_rev in label_dct:
    #     fake_well_label = label_dct[well_dct_key_rev]
    # else:
    #     ioprinter.warning_message(f'No label {well_dct_key} in label dict')

    # New MESS label code
    fake_well_label = make_rxn_str(chnl[side_idx], prepend='FakeW-')
    chn_idx = tsname.split('_')[2]  # ts_pesidx_chnidx_sadpt_idx
    _side_str = '+'.join(chnl[side_idx])
    aux_str = f'Fake Well for {_side_str}'
    fake_well, well_dat = blocks.fake_species_block(*fake_inf_dcts)
    well_str += mess_io.writer.well(
        fake_well_label, fake_well,
        aux_id_label=aux_str,
        zero_ene=chnl_enes[well_key])

    # MESS PST TS string for fake reactant side well -> reacs

    # Old MESS label code
    # pst_dct_key = make_rxn_str(chnl[side_idx], prepend=prepend_key)
    # pst_dct_key_rev = make_rxn_str(chnl[side_idx][::-1], prepend=prepend_key)
    # if pst_dct_key in label_dct:
    #     pst_label = label_dct[pst_dct_key]
    # elif pst_dct_key_rev in label_dct:
    #     pst_label = label_dct[pst_dct_key_rev]
    # else:
    #     ioprinter.warning_message(f'No label {pst_dct_key} in label dict')

    # New MESS label code (use channel index for PST barrier label)
    #pst_label = f'{prepend_key}{chn_idx}'
    pst_label = make_rxn_str(chnl[side_idx], prepend='FakeB-')
    pst_ts_str, pst_ts_dat = blocks.pst_block(ts_inf_dct, *fake_inf_dcts)
    ts_str += '\n' + mess_io.writer.ts_sadpt(
        pst_label, side_label, fake_well_label, pst_ts_str,
        aux_id_label=None,
        zero_ene=chnl_enes[ts_key],
        tunnel='')

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
                     run_prefix, save_prefix, nprocs=1):
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
        _need_ene_trans = bool(len(rgts) == 1)
        for rgt in rgts:
            spc_locs_lst = filesys.models.get_spc_locs_lst(
                spc_dct[rgt], spc_model_dct_i,
                run_prefix, save_prefix, saddle=False,
                cnf_range=cnf_range, sort_info_lst=sort_info_lst,
                name=rgt, nprocs=nprocs)
            chnl_infs_i, model_basis_energy_dct = build.read_spc_data(
                spc_dct, rgt,
                pes_model_dct_i, spc_model_dct_i,
                run_prefix, save_prefix, model_basis_energy_dct,
                calc_ene_trans=_need_ene_trans,
                spc_locs=spc_locs_lst[0])
            chnl_infs[side].append(chnl_infs_i)

    # Get data for all configurations for a TS
    chnl_infs['ts'] = []
    for name in tsname_allconfigs:
        spc_locs_lst = filesys.models.get_spc_locs_lst(
            spc_dct[name], spc_model_dct_i,
            run_prefix, save_prefix, saddle=True,
            cnf_range=cnf_range, sort_info_lst=sort_info_lst,
            name=name, nprocs=nprocs)
        spc_locs = spc_locs_lst[0] if spc_locs_lst else None
        inf_dct, model_basis_energy_dct = build.read_ts_data(
            spc_dct, name, reacs, prods,
            pes_model_dct_i, spc_model_dct_i,
            run_prefix, save_prefix, model_basis_energy_dct,
            spc_locs=spc_locs)
        chnl_infs['ts'].append(inf_dct)

    # Set up the info for the wells
    rwell_model = spc_model_dct_i['ts']['rwells']
    pwell_model = spc_model_dct_i['ts']['pwells']
    #rxn_class = spc_dct[tsname_allconfigs[0]]['class']  # no longer needed
    if need_fake_wells(reacs, rwell_model):
        chnl_infs['fake_vdwr'] = copy.deepcopy(chnl_infs['reacs'])
    if need_fake_wells(prods, pwell_model):
        chnl_infs['fake_vdwp'] = copy.deepcopy(chnl_infs['prods'])

    return chnl_infs, model_basis_energy_dct
