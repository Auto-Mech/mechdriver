"""
Write and Read MESS files for Rates
"""

import copy
import numpy
import automol
import mess_io
import ratefit
import autofile

# New libs
from lib.phydat import phycon
from lib.runner import script
from lib.load import model as loadmodel
from routines.pf.messf import blocks
from routines.pf.messf import get_fs_ene_zpe
from routines.pf.messf import calc_channel_enes
from routines.pf.messf import _tunnel as tunnel


# Writer
def make_header_str(spc_dct, rxn_lst,
                    temps, press,
                    exp_factor, exp_power, exp_cutoff,
                    eps1, eps2, sig1, sig2, mass1):
    """ makes the standard header and energy transfer sections for MESS input file
    """
    # Global Keys section
    globkey_str = mess_io.writer.global_reaction(temps, press)
    tot_mass = 0.
    for rct in rxn_lst[0]['reacs']:
        geo = automol.inchi.geometry(spc_dct[rct]['ich'])
        masses = automol.geom.masses(geo)
        for mass in masses:
            tot_mass += mass

    # Energy transfer section
    energy_trans_str = mess_io.writer.energy_transfer(
        exp_factor, exp_power, exp_cutoff,
        eps1, eps2, sig1, sig2, mass1, tot_mass)

    return globkey_str, energy_trans_str


def make_pes_mess_str(spc_dct, rxn_lst, pes_idx,
                      save_prefix, label_dct,
                      model_dct, thy_dct):
    """ Write all the MESS input file strings for the reaction channels
    """

    # Initialize empty MESS strings
    full_well_str, full_bim_str, full_ts_str = '', '', ''
    full_dat_str_lst = []

    # Get the model for the first reference species
    first_ground_model = rxn_lst[0]['model'][1]

    # Get the elec+zpe energy for the reference species
    first_ground_ene = 0.0
    first_spc = rxn_lst[0]['reacs']
    for rct in first_spc:
        first_ground_ene += get_fs_ene_zpe(
            spc_dct, rct,
            thy_dct, model_dct, first_ground_model,
            save_prefix, saddle=False)

    # Loop over all the channels and write the MESS strings
    written_labels = []
    for idx, rxn in enumerate(rxn_lst):

        # Set the TS name and channel model
        tsname = 'ts_{:g}_{:g}'.format(pes_idx, idx+1)
        chn_model = rxn['model'][1]

        # Calculate the energies of all spc on the channel
        channel_enes = calc_channel_enes(
            spc_dct, rxn, tsname,
            thy_dct, model_dct,
            chn_model, first_ground_model,
            save_prefix)

        # Write the mess strings for all spc on the channel
        mess_strs, dat_str_lst, written_labels = _make_channel_mess_strs(
            tsname, rxn, spc_dct, label_dct, written_labels,
            first_ground_ene, channel_enes,
            model_dct, thy_dct, save_prefix)

        # Append to full MESS strings
        [well_str, bim_str, ts_str] = mess_strs
        full_well_str += well_str
        full_bim_str += bim_str
        full_ts_str += ts_str
        full_dat_str_lst.extend(dat_str_lst)

    # Add final End statement for the end of the MESS file
    full_ts_str += '\nEnd\n'

    return full_well_str, full_bim_str, full_ts_str, full_dat_str_lst


def _make_channel_mess_strs(tsname, rxn, spc_dct, label_dct, written_labels,
                            first_ground_ene, channel_enes,
                            model_dct, thy_dct, save_prefix):
    """ make the partition function strings for each of the channels
    includes strings for each of the unimolecular wells, bimolecular fragments,
    and transition states connecting them.
    It also includes a special treatment for abstraction to include phase space
    blocks and coupling bimolecular fragments to fake van der Waals wells
    """

    # Initialize empty strings
    bim_str, well_str, ts_str = '', '', ''

    # Set the model and info for the reaction
    chn_model = rxn['model'][1]
    pf_levels = loadmodel.set_es_model_info(
        model_dct[chn_model]['es'], thy_dct)
    spc_model = loadmodel.set_pf_model_info(
        model_dct[chn_model]['pf'])
    # pst_params = (1.0, 6)
    ts_class = spc_dct[tsname]['class']

    # Unpack the energy dictionary and put energies in kcal
    first_ground_ene *= phycon.EH2KCAL
    reac_ene = channel_enes['reacs'] * phycon.EH2KCAL
    prod_ene = channel_enes['prods'] * phycon.EH2KCAL
    ts_ene = channel_enes['ts'] * phycon.EH2KCAL
    print('first_ground_ene', first_ground_ene)
    print('reac_ene', reac_ene)
    print('prod_ene', prod_ene)
    print('ts_ene', ts_ene)

    # Write the MESS string for the channel reactant(s) and product(s)
    rinfo = zip((rxn['reacs'], rxn['prods']), (reac_ene, prod_ene))
    for rct, rct_ene in rinfo:

        # Build the species data
        spc_data = []
        for spc in rct:
            spc_data.append(_make_spc_mess_str(
                spc_dct[spc], rxn, save_prefix, spc_model, pf_levels))

        # Write the MESS strings
        spc_label = [automol.inchi.smiles(spc_dct[spc]['ich']) for spc in rct]
        chn_label = label_dct[_make_rxn_str(rct)]
        if chn_label not in written_labels:
            written_labels.append(chn_label)
            if len(rct) > 1:
                ground_energy = rct_ene - first_ground_ene
                bim_str += '\n! {} + {}\n'.format(rct[0], rct[1])
                bim_str += mess_io.writer.bimolecular(
                    chn_label, spc_label[0], spc_data[0],
                    spc_label[1], spc_data[1], ground_energy)
            else:
                zero_energy = rct_ene - first_ground_ene
                well_str += '\n! {}\n'.format(rct)
                well_str += mess_io.writer.well(
                    chn_label, spc_data[0], zero_energy)

        # Initialize the reactant and product MESS label
        if rct == rxn['reacs']:
            reac_label = chn_label
            inner_reac_label = chn_label
        else:
            prod_label = chn_label
            inner_prod_label = chn_label

    # For abstractions: Write MESS strings for fake reac and prod wells and TS
    if _need_fake_wells(spc_dct[tsname]['class']):

        # Write all the MESS Strings for Fake Wells and TSs
        fwell_str, fts_str, fwellr_lbl, fwellp_lbl = _make_fake_mess_strs(
            spc_dct, rxn, label_dct, spc_model, pf_levels,
            save_prefix, first_ground_ene, reac_ene, prod_ene,
            reac_label, prod_label)

        # Append the fake strings to overall strings
        well_str += fwell_str
        ts_str += fts_str

        # Reset the reactant and product labels for the inner transition state
        inner_reac_label = fwellr_lbl
        inner_prod_label = fwellp_lbl

    # Write MESS string for the inner transition state; append full
    ts_label = label_dct[tsname]
    sts_str, dat_str_lst = _make_ts_mess_str(
        spc_dct, tsname, ts_class, rxn, save_prefix,
        model_dct, chn_model, spc_model, pf_levels,
        first_ground_ene, reac_ene, prod_ene, ts_ene,
        ts_label, inner_reac_label, inner_prod_label)
    ts_str += sts_str

    return [well_str, bim_str, ts_str], dat_str_lst, written_labels


def _make_spc_mess_str(spc_dct_i, rxn, save_prefix, spc_model, pf_levels):
    """ makes the main part of the MESS species block for a given species
    """
    tors_model, vib_model, _ = spc_model

    if vib_model == 'tau' or tors_model == 'tau':
        species_data = blocks.tau_block()
    else:
        species_data = blocks.species_block(
            spc_dct_i=spc_dct_i,
            rxn=rxn,
            spc_model=spc_model,
            pf_levels=pf_levels,
            save_prefix=save_prefix,
            saddle=False
        )

    return species_data


def _make_ts_mess_str(spc_dct, tsname, ts_class, rxn, save_prefix,
                      model_dct, chn_model, spc_model, pf_levels,
                      first_ground_ene, reac_ene, prod_ene, ts_ene,
                      ts_label, inner_reac_label, inner_prod_label,
                      pst_params=(1.0, 6)):
    """ makes the main part of the MESS species block for a given species
    """

    # Set the models for how we are treating the transition state
    ts_sadpt = model_dct[chn_model]['pf']['ts_sadpt']
    ts_nobarrier = model_dct[chn_model]['pf']['ts_barrierless']
    tun_model = model_dct[chn_model]['pf']['tunnel']

    # Unpack models
    # tors_model, vib_model, sym_model = spc_model
    multi_info = []

    # Initialize empty data string
    flux_str, mdhr_str, sct_str = '', '', ''

    # Get all of the information for the filesystem
    if not _var_radrad(ts_class):
        # Build MESS string for TS at a saddle point
        if ts_sadpt == 'pst':
            spc_dct_i = spc_dct[rxn['reacs'][0]]
            spc_dct_j = spc_dct[rxn['reacs'][1]]
            spc_str = blocks.pst_block(
                spc_dct_i, spc_dct_j, spc_model=spc_model,
                pf_levels=pf_levels,
                save_prefix=save_prefix,
                pst_params=pst_params)
        elif ts_sadpt == 'vtst':
            rpath_str_lst = blocks.vtst_saddle_block(
                spc_dct[tsname], pf_levels,
                ts_label, inner_reac_label, inner_prod_label, first_ground_ene)
        else:
            spc_str, mdhr_str, imag = blocks.species_block(
                spc_dct_i=spc_dct[tsname],
                spc_model=spc_model,
                pf_levels=pf_levels,
                save_prefix=save_prefix,
                saddle=True
            )
    else:
        # Build MESS string for TS with no saddle point
        if ts_nobarrier == 'vtst':
            if 'P' in inner_reac_label:
                spc_ene = reac_ene - first_ground_ene
            else:
                spc_ene = prod_ene - first_ground_ene
            rpath_str_lst = blocks.vtst_with_no_saddle_block(
                spc_dct[tsname], ts_label, inner_reac_label, inner_prod_label,
                spc_ene, multi_info)
        else:
            spc_str, flux_str = blocks.vrctst_block()

    # Write the appropriate string for the tunneling model
    tunnel_str, sct_str = '', ''
    if _treat_tunnel(tun_model, ts_sadpt, ts_nobarrier, _var_radrad(ts_class)):
        if tun_model == 'eckart':
            tunnel_str = tunnel.write_mess_eckart_str(
                ts_ene, reac_ene, prod_ene, imag)
        elif tun_model == 'sct':
            tunnel_file = tsname + '_sct.dat'
            path = 'cat'
            tunnel_str, sct_str = tunnel.write_mess_sct_str(
                spc_dct[tsname], pf_levels, path,
                imag, tunnel_file,
                cutoff_energy=2500, coord_proj='cartesian')
    else:
        pass

    # Write the MESS string for the TS
    # First if statement logic is bad
    if ts_sadpt == 'vtst' or ts_nobarrier in ('vtst', 'vrctst'):
        mess_str = mess_io.writer.ts_variational(
            ts_label, inner_reac_label, inner_prod_label, rpath_str_lst)
    else:
        zero_energy = ts_ene - first_ground_ene
        mess_str = '\n' + mess_io.writer.ts_sadpt(
            ts_label, inner_reac_label, inner_prod_label,
            spc_str, zero_energy, tunnel_str)

    # Combine dat strings together
    dat_str_lst = [flux_str, mdhr_str, sct_str]

    return mess_str, dat_str_lst


def _make_fake_mess_strs(spc_dct, rxn, label_dct, spc_model, pf_levels,
                         save_prefix,
                         first_ground_ene, reac_ene, prod_ene,
                         reac_label, prod_label,
                         pst_params=(1.0, 6)):
    """ write the MESS strings for the fake wells and TSs
    """
    # Initialize well and ts strs
    well_str, ts_str = '', ''

    # MESS string for the fake reactant side well
    spc_dct_i = spc_dct[rxn['reacs'][0]]
    spc_dct_j = spc_dct[rxn['reacs'][1]]
    well_dct_key = _make_rxn_str(rxn['reacs'], prepend='F')
    fake_wellr_label = label_dct[well_dct_key]
    vdwr_ene = reac_ene - 1.0
    zero_energy = vdwr_ene - first_ground_ene
    well_str += '\n! Fake Well for {}\n'.format(
        '+'.join(rxn['reacs']))
    fake_wellr = blocks.fake_species_block(
        spc_dct_i, spc_dct_j, spc_model, pf_levels)
    well_str += mess_io.writer.well(
        fake_wellr_label, fake_wellr, zero_energy)

    # MESS PST TS string for fake reactant side well -> reacs
    well_dct_key = _make_rxn_str(rxn['reacs'], prepend='FRB')
    pst_r_label = label_dct[well_dct_key]
    pst_r_ts_str = blocks.pst_block(
        spc_dct_i, spc_dct_j, spc_model=spc_model,
        pf_levels=pf_levels,
        save_prefix=save_prefix,
        pst_params=pst_params)
    zero_energy = reac_ene - first_ground_ene
    ts_str += '\n' + mess_io.writer.ts_sadpt(
        pst_r_label, reac_label, fake_wellr_label, pst_r_ts_str,
        zero_energy, tunnel='')

    # MESS string for the fake product side well
    spc_dct_i = spc_dct[rxn['prods'][0]]
    spc_dct_j = spc_dct[rxn['prods'][1]]
    well_dct_key = _make_rxn_str(rxn['prods'], prepend='F')
    fake_wellp_label = label_dct[well_dct_key]
    vdwp_ene = prod_ene - 1.0
    zero_energy = vdwp_ene - first_ground_ene
    well_str += '\n! Fake Well for {}\n'.format(
        '+'.join(rxn['prods']))
    fake_wellp = blocks.fake_species_block(
        spc_dct_i, spc_dct_j, spc_model, pf_levels)
    well_str += mess_io.writer.well(
        fake_wellp_label, fake_wellp, zero_energy)

    # MESS PST TS string for fake product side well -> prods
    well_dct_key = _make_rxn_str(rxn['prods'], prepend='FPB')
    pst_p_label = label_dct[well_dct_key]
    pst_p_ts_str = blocks.pst_block(
        spc_dct_i, spc_dct_j,
        spc_model=spc_model,
        pf_levels=pf_levels,
        save_prefix=save_prefix,
        pst_params=pst_params)
    zero_energy = prod_ene - first_ground_ene
    ts_str += '\n' + mess_io.writer.ts_sadpt(
        pst_p_label, prod_label, fake_wellp_label, pst_p_ts_str,
        zero_energy, tunnel='')

    return well_str, ts_str, fake_wellr_label, fake_wellp_label


# Build dictionary relating species name to MESS label
def make_pes_label_dct(rxn_lst, pes_idx, spc_dct):
    """ Builds a dictionary that matches the mechanism name to the labels used
        in the MESS input and output files for the whole PES
    """
    pes_label_dct = {}
    for idx, rxn in enumerate(rxn_lst):
        tsname = 'ts_{:g}_{:g}'.format(pes_idx, idx+1)
        pes_label_dct.update(
            _make_channel_label_dct(tsname, idx+1, pes_idx_dct, rxn, spc_dct))

    return pes_label_dct


def _make_channel_label_dct(tsname, chn_idx, label_dct, rxn, spc_dct):
    """ Builds a dictionary that matches the mechanism name to the labels used
        in the MESS input and output files
    """

    # Initialize idxs for bimol, well, and fake species
    pidx, widx, fidx = 1, 1, 1
    for val in label_dct.values():
        if 'P' in val:
            pidx += 1
        elif 'W' in val:
            widx += 1
        elif 'F' in val:
            fidx += 1

    # Determine the idxs for the channel reactants
    reac_label = ''
    bimol = bool(len(rxn['reacs']) > 1)
    well_dct_key1 = '+'.join(rxn['reacs'])
    well_dct_key2 = '+'.join(rxn['reacs'][::-1])
    if well_dct_key1 not in label_dct:
        if well_dct_key2 in label_dct:
            well_dct_key1 = well_dct_key2
        else:
            if bimol:
                reac_label = 'P' + str(pidx)
                pidx += 1
                label_dct[well_dct_key1] = reac_label
            else:
                reac_label = 'W' + str(widx)
                widx += 1
                label_dct[well_dct_key1] = reac_label
    if not reac_label:
        reac_label = label_dct[well_dct_key1]

    # Determine the idxs for the channel products
    prod_label = ''
    bimol = bool(len(rxn['prods']) > 1)
    well_dct_key1 = '+'.join(rxn['prods'])
    well_dct_key2 = '+'.join(rxn['prods'][::-1])
    if well_dct_key1 not in label_dct:
        if well_dct_key2 in label_dct:
            well_dct_key1 = well_dct_key2
        else:
            if bimol:
                prod_label = 'P' + str(pidx)
                label_dct[well_dct_key1] = prod_label
            else:
                prod_label = 'W' + str(widx)
                label_dct[well_dct_key1] = prod_label
    if not prod_label:
        prod_label = label_dct[well_dct_key1]

    # Set label for the inner transition state
    label_dct[tsname] = 'B' + str(chn_idx)

    # Determine idxs for any fake wells if they are needed
    fake_wellr_label = ''
    fake_wellp_label = ''
    if _need_fake_wells(spc_dct[tsname]['class']):
        well_dct_key1 = 'F' + '+'.join(rxn['reacs'])
        well_dct_key2 = 'F' + '+'.join(rxn['reacs'][::-1])
        if well_dct_key1 not in label_dct:
            if well_dct_key2 in label_dct:
                well_dct_key1 = well_dct_key2
            else:
                fake_wellr_label = 'F' + str(fidx)
                fidx += 1
                label_dct[well_dct_key1] = fake_wellr_label

                pst_r_label = 'FRB' + str(chn_idx)
                label_dct[well_dct_key1.replace('F', 'FRB')] = pst_r_label
            if not fake_wellr_label:
                fake_wellr_label = label_dct[well_dct_key1]
                pst_r_label = label_dct[well_dct_key1.replace('F', 'FRB')]
        else:
            fake_wellr_label = label_dct[well_dct_key1]
        well_dct_key1 = 'F' + '+'.join(rxn['prods'])
        well_dct_key2 = 'F' + '+'.join(rxn['prods'][::-1])
        if well_dct_key1 not in label_dct:
            if well_dct_key2 in label_dct:
                well_dct_key1 = well_dct_key2
            else:
                fake_wellp_label = 'F' + str(fidx)
                fidx += 1
                label_dct[well_dct_key1] = fake_wellp_label

                pst_p_label = 'FPB' + str(chn_idx)
                label_dct[well_dct_key1.replace('F', 'FPB')] = pst_p_label
            if not fake_wellp_label:
                fake_wellp_label = label_dct[well_dct_key1]
                pst_p_label = label_dct[well_dct_key1.replace('F', 'FPB')]
        else:
            fake_wellp_label = label_dct[well_dct_key1]

    return label_dct


# Various checkers to decide what kinds of MESS strings to write
def _need_fake_wells(tsclass):
    """ Return boolean to see if fake wells are needed
    """
    abst_rxn = bool('abstraction' in tsclass)
    # addn_rxn = bool('addition' in tsclass)
    subs_rxn = bool('substitution' in tsclass)
    return bool(abst_rxn or subs_rxn)
    # return bool(abst_rxn or addn_rxn or subs_rxn)


def _var_radrad(tsclass):
    """ Return boolean to see if fake wells are needed
    """
    rad_rad = 'radical radical' in tsclass
    low_spin = 'high spin' not in tsclass
    addn_rxn = 'addition' in tsclass
    return bool(rad_rad and low_spin and addn_rxn)


def _treat_tunnel(tunnel_model, ts_sadpt, ts_barrierless, var_radrad):
    """ Discern if tunneling will be treated
    """
    treat = False
    if tunnel_model == 'none':
        pass
    else:
        if var_radrad and ts_barrierless != 'pst':
            treat = True
        elif ts_sadpt != 'pst':
            treat = True

    return treat


# Formatting functions
def _make_rxn_str(rlst, prepend=''):
    """ convert list to string
    """
    return prepend + '+'.join(rlst)
