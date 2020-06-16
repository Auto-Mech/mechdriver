"""
Write and Read MESS files for Rates
"""

import importlib
import copy
import automol
import mess_io
from routines.pf.messf import blocks
from routines.pf.messf import models
from routines.pf.messf import set_reference_ene
from routines.pf.messf import calc_channel_enes
from routines.pf.messf import _tunnel as tunnel
from routines.pf.messf import is_atom
from routines.pf.ktp._util import set_pf_info
from routines.pf.ktp._util import set_ts_cls_info
from routines.pf.ktp._util import make_rxn_str
from routines.pf.ktp._util import treat_tunnel
from routines.pf.ktp._util import pst_ts
from routines.pf.ktp._util import need_fake_wells
from routines.pf.ktp._util import var_radrad
from routines.pf.ktp._util import print_pf_info


BLOCK_MODULE = importlib.import_module('routines.pf.messf.blocks')


# Writer
def make_header_str(temps, press):
    """ makes the standard header and energy transfer sections for MESS input file
    """

    print('\nPreparing global keywords section for MESS input...')

    print(' - Using temperatures and pressures defined by user')
    print(' - Using internal AutoMech defaults for other MESS keywords:')
    keystr1 = (
        'EnergyStepOverTemperature, ExcessEnergyOverTemperature, ' +
        'ModelEnergyLimit'
    )
    keystr2 = (
        'CalculationMethod, WellCutoff, ChemicalEigenvalueMax, ' +
        'ReductionMethod, AtomDistanceMin'
    )
    print('     {}'.format(keystr1))
    print('     {}'.format(keystr2))

    header_str = mess_io.writer.global_reaction(temps, press)

    return header_str


def make_etrans_str(spc_dct, rxn_lst,
                    exp_factor, exp_power, exp_cutoff,
                    eps1, eps2, sig1, sig2, mass1):
    """ makes the standard header and energy transfer sections for MESS input file
    """

    print('\nPreparing energy transfer section for MESS input...')

    # Get masses for energy transfer section
    print(' - Using masses of reactants of first channel')
    tot_mass = 0.
    for rct in rxn_lst[0]['reacs']:
        geo = automol.inchi.geometry(spc_dct[rct]['ich'])
        masses = automol.geom.masses(geo)
        for mass in masses:
            tot_mass += mass

    # Write MESS-format energy transfer section string
    print(' - Using Lennard-Jones sigma and epsilon parameters',
          'defined by the user.')
    print(' - Using exponential-down energy-transfer model parameters',
          'defined by the user.')
    # Energy transfer section
    energy_trans_str = mess_io.writer.energy_transfer(
        exp_factor, exp_power, exp_cutoff,
        eps1, eps2, sig1, sig2, mass1, tot_mass)

    return energy_trans_str


def make_pes_mess_str(spc_dct, rxn_lst, pes_idx,
                      run_prefix, save_prefix, label_dct,
                      model_dct, thy_dct):
    """ Write all the MESS input file strings for the reaction channels
    """

    print('\nPreparing reaction channel section for MESS input... ')

    # Initialize empty MESS strings
    full_well_str, full_bi_str, full_ts_str = '', '', ''
    full_dat_str_lst = []

    # Set the energy and model for the first reference species
    # print('\nCalculating reference energy for PES')
    ref_ene, ref_model = set_reference_ene(
        rxn_lst, spc_dct, thy_dct, model_dct,
        run_prefix, save_prefix, ref_idx=0)

    # Loop over all the channels and write the MESS strings
    written_labels = []
    for rxn in rxn_lst:

        print('Reading PES electronic structure data ' +
              'from save filesystem for')
        print('Channel {}: {} = {}...'.format(
            rxn['chn_idx'],
            '+'.join(rxn['reacs']),
            '+'.join(rxn['prods'])))

        # Set the TS name and channel model
        tsname = 'ts_{:g}_{:g}'.format(pes_idx, rxn['chn_idx'])
        chn_model = rxn['model'][1]

        # Obtain useful info objects
        pf_info = set_pf_info(model_dct, thy_dct, chn_model, ref_model)
        ts_cls_info = set_ts_cls_info(spc_dct, model_dct, tsname, chn_model)

        # Print models
        ref_ene_lvl = model_dct[ref_model]['es']['ene']
        print_pf_info(pf_info[1], pf_info[0], chn_model, ref_ene_lvl)

        # Obtain all of the species data
        chnl_infs = get_channel_data(rxn, tsname, spc_dct,
                                     pf_info, ts_cls_info,
                                     run_prefix, save_prefix)

        # Calculate the energies of all spc on the channel
        chnl_enes = calc_channel_enes(chnl_infs, ref_ene,
                                      chn_model, ref_model)

        # Write the mess strings for all spc on the channel
        mess_strs, dat_str_lst, written_labels = _make_channel_mess_strs(
            tsname, rxn, spc_dct, label_dct, written_labels,
            chnl_infs, chnl_enes, ts_cls_info)

        # Append to full MESS strings
        [well_str, bi_str, ts_str] = mess_strs
        full_well_str += well_str
        full_bi_str += bi_str
        full_ts_str += ts_str
        full_dat_str_lst.extend(dat_str_lst)

    # Add final End statement for the end of the MESS file
    full_ts_str += '\nEnd\n'

    return full_well_str, full_bi_str, full_ts_str, full_dat_str_lst


def _make_channel_mess_strs(tsname, rxn, spc_dct, label_dct, written_labels,
                            chnl_infs, chnl_enes, ts_cls_info):
    """ make the partition function strings for each of the channels
    includes strings for each of the unimolecular wells, bimolecular fragments,
    and transition states connecting them.
    It also includes a special treatment for abstraction to include phase space
    blocks and coupling bimolecular fragments to fake van der Waals wells
    """

    # Initialize empty strings
    bi_str, well_str, ts_str = '', '', ''

    # Write the MESS string for the channel reactant(s) and product(s)
    for side in ('reacs', 'prods'):

        # Get information from relevant dictionaries
        rgt_names = rxn[side]
        rgt_infs = chnl_infs[side]
        rgt_ene = chnl_enes[side]

        # Build the species string for reactant(s)/product(s)
        spc_str = [_make_spc_mess_str(inf) for inf in rgt_infs]

        # Set the labels to put into the file
        spc_label = [automol.inchi.smiles(spc_dct[name]['ich'])
                     for name in rgt_names]
        chn_label = label_dct[make_rxn_str(rgt_names)]

        # Write the strings
        if chn_label not in written_labels:
            written_labels.append(chn_label)
            if len(rgt_names) > 1:
                # bi_str += mess_io.writer.species_separation_str()
                bi_str += '\n! {} + {}\n'.format(rgt_names[0], rgt_names[1])
                bi_str += mess_io.writer.bimolecular(
                    chn_label, spc_label[0], spc_str[0],
                    spc_label[1], spc_str[1], rgt_ene)
            else:
                # well_str += mess_io.writer.species_separation_str()
                well_str += '\n! {}\n'.format(rgt_names[0])
                well_str += mess_io.writer.well(
                    chn_label, spc_str[0], rgt_ene)

        # Initialize the reactant and product MESS label
        if side == 'reacs':
            reac_label = chn_label
            inner_reac_label = chn_label
        else:
            prod_label = chn_label
            inner_prod_label = chn_label

    # For abstractions: Write MESS strings for fake reac and prod wells and TS
    print(list(chnl_infs.keys()))
    if 'fake_vdwr' in chnl_infs:
        print('after if')
        # Write all the MESS Strings for Fake Wells and TSs
        fwell_str, fts_str, fwellr_lbl, fwellp_lbl = _make_fake_mess_strs(
            rxn, chnl_infs['fake_vdwr'], chnl_infs['fake_vdwp'],
            chnl_enes, label_dct, reac_label, prod_label)

        # Append the fake strings to overall strings
        well_str += fwell_str
        ts_str += fts_str

        # Reset the reactant and product labels for the inner transition state
        inner_reac_label = fwellr_lbl
        inner_prod_label = fwellp_lbl

    # Write MESS string for the inner transition state; append full
    ts_label = label_dct[tsname]
    sts_str, dat_str_lst = _make_ts_mess_str(
        chnl_infs, chnl_enes, ts_cls_info,
        ts_label, inner_reac_label, inner_prod_label)
    ts_str += sts_str

    return [well_str, bi_str, ts_str], dat_str_lst, written_labels


def _make_spc_mess_str(inf_dct):
    """ makes the main part of the MESS species block for a given species
    """
    mess_writer = getattr(BLOCK_MODULE, inf_dct['writer'])
    return mess_writer(inf_dct)


def _make_ts_mess_str(chnl_infs, chnl_enes, ts_cls_info,
                      ts_label, inner_reac_label, inner_prod_label):
    """ makes the main part of the MESS species block for a given species
    """

    # Unpack info objects
    [_, ts_sadpt, ts_nobarrier, tunnel_model] = ts_cls_info

    # Initialize empty data string
    flux_str, mdhr_str, sct_str = '', '', ''

    # Write the initial data string
    mess_writer = getattr(BLOCK_MODULE, chnl_infs['ts']['writer'])
    mess_str = mess_writer(chnl_infs['ts'])
    # mess_str = mess_writer(*chnl_infs['ts'])

    # Write the appropriate string for the tunneling model
    tunnel_str, sct_str = '', ''
    if treat_tunnel(tunnel_model, ts_sadpt, ts_nobarrier):
        if tunnel_model == 'eckart':
            tunnel_str = tunnel.write_mess_eckart_str(
                chnl_enes['ts'], chnl_enes['reacs'], chnl_enes['prods'],
                chnl_infs['ts']['imag'])
        # elif tun_model == 'sct':
        #     tunnel_file = tsname + '_sct.dat'
        #     path = 'cat'
        #     tunnel_str, sct_str = tunnel.write_mess_sct_str(
        #         spc_dct[tsname], pf_levels, path,
        #         imag, tunnel_file,
        #         cutoff_energy=2500, coord_proj='cartesian')
    else:
        pass

    # Write the MESS string for the TS
    if ts_sadpt == 'vtst' or ts_nobarrier in ('vtst', 'vrctst'):
        ts_str = mess_io.writer.ts_variational(
            ts_label, inner_reac_label, inner_prod_label,
            mess_str, tunnel_str)
    else:
        ts_str = '\n' + mess_io.writer.ts_sadpt(
            ts_label, inner_reac_label, inner_prod_label,
            mess_str, chnl_enes['ts'], tunnel_str)

    # Combine dat strings together
    dat_str_lst = [flux_str, mdhr_str, sct_str]

    return ts_str, dat_str_lst


def _make_fake_mess_strs(rxn, fake_wellr_inf_dcts, fake_wellp_inf_dcts,
                         chnl_enes, label_dct, reac_label, prod_label):
    """ write the MESS strings for the fake wells and TSs
    """
    print('HERE')

    # Initialize well and ts strs
    well_str, ts_str = '', ''

    # MESS string for the fake reactant side well
    well_dct_key = make_rxn_str(rxn['reacs'], prepend='F')
    fake_wellr_label = label_dct[well_dct_key]
    # well_str += mess_io.writer.species_separation_str()
    well_str += '\n! Fake Well for {}\n'.format(
        '+'.join(rxn['reacs']))
    fake_wellr = blocks.fake_species_block(*fake_wellr_inf_dcts)
    well_str += mess_io.writer.well(
        fake_wellr_label, fake_wellr, chnl_enes['fake_vdwr'])

    # MESS PST TS string for fake reactant side well -> reacs
    well_dct_key = make_rxn_str(rxn['reacs'], prepend='FRB')
    pst_r_label = label_dct[well_dct_key]
    pst_r_ts_str = blocks.pst_block(*fake_wellr_inf_dcts)
    ts_str += '\n' + mess_io.writer.ts_sadpt(
        pst_r_label, reac_label, fake_wellr_label, pst_r_ts_str,
        chnl_enes['fake_vdwr_ts'], tunnel='')

    # MESS string for the fake product side well
    well_dct_key = make_rxn_str(rxn['prods'], prepend='F')
    fake_wellp_label = label_dct[well_dct_key]
    well_str += '\n! Fake Well for {}\n'.format(
        '+'.join(rxn['prods']))
    fake_wellp = blocks.fake_species_block(*fake_wellp_inf_dcts)
    well_str += mess_io.writer.well(
        fake_wellp_label, fake_wellp, chnl_enes['fake_vdwp'])

    # MESS PST TS string for fake product side well -> prods
    well_dct_key = make_rxn_str(rxn['prods'], prepend='FPB')
    pst_p_label = label_dct[well_dct_key]
    pst_p_ts_str = blocks.pst_block(*fake_wellp_inf_dcts)
    ts_str += '\n' + mess_io.writer.ts_sadpt(
        pst_p_label, prod_label, fake_wellp_label, pst_p_ts_str,
        chnl_enes['fake_vdwp_ts'], tunnel='')

    return well_str, ts_str, fake_wellr_label, fake_wellp_label


# Data Retriever Functions
def get_channel_data(rxn, tsname, spc_dct, pf_info, ts_cls_info,
                     run_prefix, save_prefix):
    """ generate dcts with the models
    """

    # Unpack info objects
    [chn_pf_levels, chn_pf_models, ref_pf_levels, ref_pf_models] = pf_info
    [ts_class, ts_sadpt, ts_nobarrier, _] = ts_cls_info

    # Determine the MESS data for the channel
    chnl_infs = {}
    chnl_infs['reacs'] = []
    chnl_infs['prods'] = []
    chnl_infs['ts'] = []
    for rct in rxn['reacs']:
        inf_dct = read_spc_data(spc_dct[rct], rct,
                                chn_pf_models, chn_pf_levels,
                                ref_pf_models, ref_pf_levels,
                                run_prefix, save_prefix)
        chnl_infs['reacs'].append(inf_dct)
    for prd in rxn['prods']:
        inf_dct = read_spc_data(spc_dct[prd], prd,
                                chn_pf_models, chn_pf_levels,
                                ref_pf_models, ref_pf_levels,
                                run_prefix, save_prefix)
        chnl_infs['prods'].append(inf_dct)

    # Set up data for TS
    if pst_ts(ts_class, ts_sadpt, ts_nobarrier):
        chnl_infs['ts'] = {'writer': 'blocks.pst_block'}
    else:
        inf_dct = read_ts_data(spc_dct[tsname], tsname,
                               chn_pf_models, chn_pf_levels,
                               ref_pf_models, ref_pf_levels,
                               run_prefix, save_prefix,
                               ts_class, ts_sadpt, ts_nobarrier)
        chnl_infs['ts'] = inf_dct

    # Set up the info for the wells
    if need_fake_wells(ts_class):
        chnl_infs['fake_vdwr'] = copy.deepcopy(chnl_infs['reacs'])
        chnl_infs['fake_vdwp'] = copy.deepcopy(chnl_infs['prods'])

    return chnl_infs


def read_spc_data(spc_dct_i, spc_name,
                  chn_pf_models, chn_pf_levels,
                  ref_pf_models, ref_pf_levels,
                  run_prefix, save_prefix):
    """ Determines which block writer to use tau
    """
    print(('\n++++++++++++++++++++++++++++++++++++++++++++++++' +
           '++++++++++++++++++++++++++++++++++++++'))
    print('\nReading filesystem info for {}'.format(spc_name))

    vib_model, tors_model = chn_pf_models['vib'], chn_pf_models['tors']
    if is_atom(spc_dct_i):
        inf_dct = models.atm_data(
            spc_dct_i,
            chn_pf_models, chn_pf_levels,
            ref_pf_models, ref_pf_levels,
            run_prefix, save_prefix)
        writer = 'atom_block'
    else:
        if vib_model == 'tau' or tors_model == 'tau':
            pass
            # inf_dct = models.tau_data(
            #     spc_dct_i,
            #     chn_pf_models, chn_pf_levels,
            #     ref_pf_models, ref_pf_levels,
            #     run_prefix, save_prefix, saddle=False)
            # writer = 'tau_block'
        else:
            inf_dct = models.mol_data(
                spc_dct_i,
                chn_pf_models, chn_pf_levels,
                ref_pf_models, ref_pf_levels,
                run_prefix, save_prefix, saddle=False, tors_wgeo=True)
            writer = 'species_block'

    # Add writer to inf dct
    inf_dct['writer'] = writer

    return inf_dct


def read_ts_data(spc_dct_i, tsname,
                 chn_pf_models, chn_pf_levels,
                 ref_pf_models, ref_pf_levels,
                 run_prefix, save_prefix,
                 ts_class, ts_sadpt, ts_nobarrier):
    """ Determine which block function to useset block functions
    """

    print(('\n++++++++++++++++++++++++++++++++++++++++++++++++' +
           '++++++++++++++++++++++++++++++++++++++'))
    print('\nReading filesystem info for {}'.format(tsname))

    # Get all of the information for the filesystem
    if not var_radrad(ts_class):
        # Build MESS string for TS at a saddle point
        if ts_sadpt == 'vtst':
            inf_dct = 'rpvtst_data'
            writer = 'vtst_saddle_block'
        else:
            inf_dct = models.mol_data(
                spc_dct_i,
                chn_pf_models, chn_pf_levels,
                ref_pf_models, ref_pf_levels,
                run_prefix, save_prefix, saddle=True, tors_wgeo=True)
            writer = 'species_block'
    else:
        # Build MESS string for TS with no saddle point
        if ts_nobarrier == 'vtst':
            inf_dct = 'rpvtst_data'
            writer = 'vtst_no_saddle_block'
        elif ts_nobarrier == 'vrctst':
            inf_dct = 'vrctst_data'
            writer = 'vrctst_block'

    # Add writer to inf dct
    inf_dct['writer'] = writer

    return inf_dct
