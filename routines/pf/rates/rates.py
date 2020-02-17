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


# Writer
def rate_headers(spc_dct, rxn_lst,
                 temps, press,
                 exp_factor, exp_power, exp_cutoff,
                 eps1, eps2, sig1, sig2, mass1):
    """ makes the standard header and energy transfer sections for MESS input file
    """
    # Header section
    header_str = mess_io.writer.global_reaction(temps, press)
    tot_mass = 0.
    for rct in rxn_lst[0]['reacs']:
        geo = automol.convert.inchi.geometry(spc_dct[rct]['ich'])
        masses = automol.geom.masses(geo)
        for mass in masses:
            tot_mass += mass

    # Energy transfer section
    energy_trans_str = mess_io.writer.energy_transfer(
        exp_factor, exp_power, exp_cutoff,
        eps1, eps2, sig1, sig2, mass1, tot_mass)

    return header_str, energy_trans_str


def write_channel_mess_strs(spc_dct, rxn_lst, pes_formula,
                            multi_info, pst_params,
                            save_prefix, idx_dct,
                            model_dct, thy_dct):
    """ Write all the MESS input file strings for the reaction channels
    """
    mess_strs = ['', '', '']

    # Get the model for the first reference species
    first_ground_model = rxn_lst[0]['model']

    # Get the elec+zpe energy for the reference species
    first_ground_ene = 0.0
    first_spc = rxn_lst[0]['reacs']
    for rct in first_spc:
        first_ground_ene += get_fs_ene_zpe(
            spc_dct, rct,
            thy_dct, model_dct, first_ground_model,
            save_prefix, saddle=False)

    # Write the MESS data strings for all the species; no ene
    species, dat_str_lst = make_all_species_data(
        rxn_lst, spc_dct, model_dct, thy_dct, save_prefix)

    # Loop over all the channels and write the MESS strings
    for idx, rxn in enumerate(rxn_lst):
        tsname = 'ts_{:g}'.format(idx)
        tsform = automol.formula.string(
            automol.geom.formula(
                automol.zmatrix.geometry(
                    spc_dct[tsname]['zma'])))
        # spc_dct[tsname]['original_zma'])))
        if tsform != pes_formula:
            print('Reaction ist contains reactions on different potential',
                  'energy surfaces: {} and {}'.format(tsform, pes_formula))
            print('Will proceed to construct only {}'.format(pes_formula))
            continue
        chn_model = rxn['model']
        channel_enes = calc_channel_enes(
            spc_dct, rxn, tsname,
            thy_dct, model_dct,
            chn_model, first_ground_model,
            save_prefix)
        # print('channel_enes')
        # print(channel_enes)
        mess_strs = make_channel_pfs(
            tsname, rxn, species, spc_dct, idx_dct, mess_strs,
            first_ground_ene, channel_enes,
            model_dct, thy_dct, multi_info, save_prefix,
            pst_params=pst_params)
    well_str, bim_str, ts_str = mess_strs
    ts_str += '\nEnd\n'

    return well_str, bim_str, ts_str, dat_str_lst


def make_all_species_data(rxn_lst, spc_dct, model_dct, thy_dct, save_prefix):
    """ generate the MESS species blocks for all the species
    """
    species = {}
    dat_str_lst = []

    spc_save_fs = autofile.fs.species(save_prefix)
    for idx, rxn in enumerate(rxn_lst):
        tsname = 'ts_{:g}'.format(idx)
        specieslist = rxn['reacs'] + rxn['prods']
        rxn_model = rxn['model']
        # Gather PF model and theory level info
        pf_levels = loadmodel.set_es_model_info(
            model_dct[rxn_model]['es'], thy_dct)
        pf_model = loadmodel.set_pf_model_info(
            model_dct[rxn_model]['pf'])

        # Get PF input header
        for name in specieslist:
            if name not in species:
                species_data, dat_str_dct, _ = make_species_data(
                    name, spc_dct[name], spc_save_fs,
                    pf_model, pf_levels)
                species[name] = species_data
                dat_str_lst.append(dat_str_dct)
        if 'radical radical addition' not in spc_dct[tsname]['class']:
            ret1, ret2, ret3 = make_species_data(
                tsname, spc_dct[tsname], save_prefix,
                pf_model, pf_levels)
            species[tsname] = ret1
            dat_str_lst.append(ret2)
            spc_dct[tsname]['imag_freq'] = ret3
    return species, dat_str_lst


def make_species_data(spc, spc_dct_i,
                      spc_save_fs, spc_model, pf_levels):
    """ makes the main part of the MESS species block for a given species
    """
    spc_info = (spc_dct_i['ich'], spc_dct_i['chg'], spc_dct_i['mul'])
    if 'ts_' in spc:
        save_path = spc_dct_i['rxn_fs'][3]
    else:
        spc_save_fs[-1].create(spc_info)
        save_path = spc_save_fs[-1].path(spc_info)
    species_data = blocks.species_block(
        spc=spc,
        spc_dct_i=spc_dct_i,
        spc_info=spc_info,
        spc_model=spc_model,
        pf_levels=pf_levels,
        save_prefix=save_path,
        )
    return species_data


def make_fake_species_data(spc_dct_i, spc_dct_j, spc_save_fs,
                           spc_model, pf_levels):
    """ make a fake MESS species block to represent the van der Waals well
    arising from the combination of two fragment species
    """
    spc_info_i = (spc_dct_i['ich'], spc_dct_i['chg'], spc_dct_i['mul'])
    spc_info_j = (spc_dct_j['ich'], spc_dct_j['chg'], spc_dct_j['mul'])
    spc_save_fs[-1].create(spc_info_i)
    spc_save_fs[-1].create(spc_info_j)
    save_path_i = spc_save_fs[-1].path(spc_info_i)
    save_path_j = spc_save_fs[-1].path(spc_info_j)
    species_data = blocks.fake_species_block(
        spc_dct_i=spc_dct_i,
        spc_dct_j=spc_dct_j,
        spc_info_i=spc_info_i,
        spc_info_j=spc_info_j,
        spc_model=spc_model,
        pf_levels=pf_levels,
        save_prefix_i=save_path_i,
        save_prefix_j=save_path_j
        )
    return species_data


def make_channel_pfs(
        tsname, rxn, species_data, spc_dct, idx_dct, strs,
        first_ground_ene, channel_enes,
        model_dct, thy_dct, multi_info, save_prefix,
        pst_params=(1.0, 6),
        rad_rad_ts='pst'):
    """ make the partition function strings for each of the channels
    includes strings for each of the unimolecular wells, bimolecular fragments,
    and transition states connecting them.
    It also includes a special treatment for abstraction to include phase space
    blocks and coupling bimolecular fragments to fake van der Waals wells
    """
    projrot_script_str = script.PROJROT
    bim_str, well_str, ts_str = strs

    # Set the model and info for the reaction
    chn_model = rxn['model']
    pf_levels = loadmodel.set_es_model_info(
        model_dct[chn_model]['es'], thy_dct)
    spc_model = loadmodel.set_pf_model_info(
        model_dct[chn_model]['pf'])

    # Unpack the energy dictionary and put energies in kcal
    first_ground_ene *= phycon.EH2KCAL
    reac_ene = channel_enes['reacs'] * phycon.EH2KCAL
    prod_ene = channel_enes['prods'] * phycon.EH2KCAL
    ts_ene = channel_enes['ts'] * phycon.EH2KCAL
    print('first_grouned_ene', first_ground_ene)
    print('reac_ene', reac_ene)
    print('prod_ene', prod_ene)
    print('ts_ene', ts_ene)

    # Set filesys object
    spc_save_fs = autofile.fs.species(save_prefix)

    # Find the number of uni and bimolecular wells already in the dictionary
    pidx = 1
    widx = 1
    fidx = 1
    for val in idx_dct.values():
        if 'P' in val:
            pidx += 1
        elif 'W' in val:
            widx += 1
        elif 'F' in val:
            fidx += 1

    # Set up new well for the reactants if that combo isn't already in the dct
    reac_label = ''
    bimol = False
    if len(rxn['reacs']) > 1:
        bimol = True
    well_data = []
    spc_label = []
    for reac in rxn['reacs']:
        spc_label.append(automol.inchi.smiles(spc_dct[reac]['ich']))
        well_data.append(species_data[reac])
    # Wells and bimol stuff
    well_dct_key1 = '+'.join(rxn['reacs'])
    well_dct_key2 = '+'.join(rxn['reacs'][::-1])
    if well_dct_key1 not in idx_dct:
        if well_dct_key2 in idx_dct:
            well_dct_key1 = well_dct_key2
        else:
            if bimol:
                reac_label = 'P' + str(pidx)
                pidx += 1
                ground_energy = reac_ene - first_ground_ene
                bim_str += '\n! {} + {} \n'.format(
                    rxn['reacs'][0], rxn['reacs'][1])
                bim_str += mess_io.writer.bimolecular(
                    reac_label, spc_label[0], well_data[0], spc_label[1],
                    well_data[1], ground_energy)
                idx_dct[well_dct_key1] = reac_label
            else:
                reac_label = 'W' + str(widx)
                widx += 1
                zero_energy = reac_ene - first_ground_ene
                well_str += '\n! {} \n'.format(rxn['reacs'][0])
                well_str += mess_io.writer.well(
                    reac_label, well_data[0], zero_energy)
                idx_dct[well_dct_key1] = reac_label
    if not reac_label:
        reac_label = idx_dct[well_dct_key1]

    # Set up a new well for the products if that combo isn't already in the dct
    prod_label = ''
    bimol = False
    if len(rxn['prods']) > 1:
        bimol = True
    well_data = []
    spc_label = []
    for prod in rxn['prods']:
        spc_label.append(automol.inchi.smiles(spc_dct[prod]['ich']))
        well_data.append(species_data[prod])
    zero_energy = prod_ene - reac_ene
    well_dct_key1 = '+'.join(rxn['prods'])
    well_dct_key2 = '+'.join(rxn['prods'][::-1])
    if well_dct_key1 not in idx_dct:
        if well_dct_key2 in idx_dct:
            well_dct_key1 = well_dct_key2
        else:
            if bimol:
                prod_label = 'P' + str(pidx)
                ground_energy = prod_ene - first_ground_ene
                bim_str += '\n! {} + {} \n'.format(
                    rxn['prods'][0], rxn['prods'][1])
                bim_str += mess_io.writer.bimolecular(
                    prod_label, spc_label[0], well_data[0],
                    spc_label[1], well_data[1],
                    ground_energy)
                idx_dct[well_dct_key1] = prod_label
            else:
                prod_label = 'W' + str(widx)
                zero_energy = prod_ene - first_ground_ene
                well_str += '\n! {} \n'.format(rxn['prods'][0])
                well_str += mess_io.writer.well(
                    prod_label, well_data[0], zero_energy)
                idx_dct[well_dct_key1] = prod_label
    if not prod_label:
        prod_label = idx_dct[well_dct_key1]

    # For abstraction first make fake wells and PST TSs
    # Then for tight TS's make ts_str
    zero_energy = ts_ene - first_ground_ene
    ts_label = 'B' + str(int(tsname.replace('ts_', ''))+1)
    imag_freq = 0
    if 'imag_freq' in spc_dct[tsname]:
        imag_freq = abs(spc_dct[tsname]['imag_freq'])
    if not imag_freq:
        print('No imaginary freq for ts: {}'.format(tsname))

    fake_wellr_label = ''
    fake_wellp_label = ''
    rad_rad = 'radical radical' in spc_dct[tsname]['class']
    low_spin = 'high spin' not in spc_dct[tsname]['class']
    abst_rxn = 'abstraction' in spc_dct[tsname]['class']
    addn_rxn = 'addition' in spc_dct[tsname]['class']
    subs_rxn = 'substitution' in spc_dct[tsname]['class']
    if abst_rxn or subs_rxn:
        # Make fake wells and PST TSs as needed
        well_dct_key1 = 'F' + '+'.join(rxn['reacs'])
        well_dct_key2 = 'F' + '+'.join(rxn['reacs'][::-1])
        pst_r_ts_str = ''
        if well_dct_key1 not in idx_dct:
            if well_dct_key2 in idx_dct:
                well_dct_key1 = well_dct_key2
            else:
                fake_wellr_label = 'F' + str(fidx)
                fidx += 1
                vdwr_ene = reac_ene - 1.0
                zero_energy = vdwr_ene - first_ground_ene
                well_str += '\n! Fake Well for {}\n'.format(
                    '+'.join(rxn['reacs']))
                fake_wellr = make_fake_species_data(
                    spc_dct[rxn['reacs'][0]], spc_dct[rxn['reacs'][1]],
                    spc_save_fs, spc_model, pf_levels)
                well_str += mess_io.writer.well(
                    fake_wellr_label, fake_wellr, zero_energy)
                idx_dct[well_dct_key1] = fake_wellr_label

                pst_r_label = 'FRB' + str(int(tsname.replace('ts_', ''))+1)
                idx_dct[well_dct_key1.replace('F', 'FRB')] = pst_r_label
                spc_dct_i = spc_dct[rxn['reacs'][0]]
                spc_dct_j = spc_dct[rxn['reacs'][1]]
                pst_r_ts_str = blocks.pst_block(
                    spc_dct_i, spc_dct_j, spc_model=spc_model,
                    pf_levels=pf_levels,
                    spc_save_fs=spc_save_fs,
                    pst_params=pst_params)
            if not fake_wellr_label:
                fake_wellr_label = idx_dct[well_dct_key1]
                pst_r_label = idx_dct[well_dct_key1.replace('F', 'FRB')]
            zero_energy = reac_ene - first_ground_ene
            tunnel_str = ''
            ts_str += '\n' + mess_io.writer.ts_sadpt(
                pst_r_label, reac_label, fake_wellr_label, pst_r_ts_str,
                zero_energy, tunnel_str)
        else:
            fake_wellr_label = idx_dct[well_dct_key1]
        well_dct_key1 = 'F' + '+'.join(rxn['prods'])
        well_dct_key2 = 'F' + '+'.join(rxn['prods'][::-1])
        if well_dct_key1 not in idx_dct:
            if well_dct_key2 in idx_dct:
                well_dct_key1 = well_dct_key2
            else:
                fake_wellp_label = 'F' + str(fidx)
                fidx += 1
                vdwp_ene = prod_ene - 1.0
                zero_energy = vdwp_ene - first_ground_ene
                well_str += '\n! Fake Well for {}\n'.format(
                    '+'.join(rxn['prods']))
                fake_wellp = make_fake_species_data(
                    spc_dct[rxn['prods'][0]], spc_dct[rxn['prods'][1]],
                    spc_save_fs, spc_model, pf_levels)
                well_str += mess_io.writer.well(
                    fake_wellp_label, fake_wellp, zero_energy)
                idx_dct[well_dct_key1] = fake_wellp_label

                pst_p_label = 'FPB' + str(int(tsname.replace('ts_', ''))+1)
                idx_dct[well_dct_key1.replace('F', 'FPB')] = pst_p_label
                spc_dct_i = spc_dct[rxn['prods'][0]]
                spc_dct_j = spc_dct[rxn['prods'][1]]
                pst_p_ts_str = blocks.pst_block(
                    spc_dct_i, spc_dct_j, spc_model=spc_model,
                    pf_levels=pf_levels,
                    spc_save_fs=spc_save_fs,
                    pst_params=pst_params)
            if not fake_wellp_label:
                fake_wellp_label = idx_dct[well_dct_key1]
                pst_p_label = idx_dct[well_dct_key1.replace('F', 'FPB')]
            zero_energy = prod_ene - first_ground_ene
            tunnel_str = ''
            ts_str += '\n' + mess_io.writer.ts_sadpt(
                pst_p_label, prod_label, fake_wellp_label, pst_p_ts_str,
                zero_energy, tunnel_str)
        else:
            fake_wellp_label = idx_dct[well_dct_key1]

        # Print inner TS data for radical radical call vtst or vrctst
        if rad_rad and low_spin and rad_rad_ts == 'pst':
            # zero_energy = SOMETHING
            pst_str = blocks.pst_block(
                spc_dct_i, spc_dct_j, spc_model=spc_model,
                pf_levels=pf_levels,
                spc_save_fs=spc_save_fs,
                pst_params=pst_params)
            ts_str += '\n' + mess_io.writer.ts_sadpt(
                ts_label, reac_label, prod_label, pst_str, zero_energy)
        elif rad_rad and low_spin:
            ts_label = 'B' + str(int(tsname.replace('ts_', ''))+1)
            if 'P' in reac_label:
                spc_ene = reac_ene - first_ground_ene
                spc_zpe = (
                    spc_dct[rxn['reacs'][0]]['zpe'] +
                    spc_dct[rxn['reacs'][1]]['zpe'])
            else:
                spc_ene = prod_ene - first_ground_ene
                spc_zpe = (
                    spc_dct[rxn['prods'][0]]['zpe'] +
                    spc_dct[rxn['prods'][1]]['zpe'])
            ts_str += '\n' + blocks.vtst_with_no_saddle_block(
                spc_dct[tsname], ts_label, fake_wellr_label, fake_wellp_label,
                spc_ene, spc_zpe, projrot_script_str, multi_info)
        else:
            vdwr_ene = reac_ene - 1.0
            vdwp_ene = prod_ene - 1.0
            zero_energy = ts_ene - first_ground_ene
            ts_reac_barr = ts_ene - vdwr_ene
            ts_prod_barr = ts_ene - vdwp_ene
            if ts_reac_barr < 0.:
                ts_reac_barr = 0.1
            if ts_prod_barr < 0.:
                ts_prod_barr = 0.1
            tunnel_str = mess_io.writer.tunnel_eckart(
                imag_freq, ts_reac_barr, ts_prod_barr)
            ts_str += '\n' + mess_io.writer.ts_sadpt(
                ts_label, fake_wellr_label, fake_wellp_label,
                species_data[tsname], zero_energy, tunnel_str)
    elif rad_rad and addn_rxn and low_spin:
        # For radical radical addition call vtst or vrctst
        ts_label = 'B' + str(int(tsname.replace('ts_', ''))+1)
        if 'P' in reac_label:
            spc_ene = reac_ene - first_ground_ene
            spc_zpe = (
                spc_dct[rxn['reacs'][0]]['zpe'] +
                spc_dct[rxn['reacs'][1]]['zpe'])
        else:
            spc_ene = prod_ene - first_ground_ene
            spc_zpe = (
                spc_dct[rxn['prods'][0]]['zpe'] +
                spc_dct[rxn['prods'][1]]['zpe'])
        ts_str += '\n' + blocks.vtst_with_no_saddle_block(
            spc_dct[tsname], ts_label, reac_label, prod_label,
            spc_ene, spc_zpe, projrot_script_str,
            multi_info)
    elif rad_rad and addn_rxn and low_spin and rad_rad_ts == 'pst':
        # zero_energy = SOMETHING
        pst_str = blocks.pst_block(
            spc_dct_i, spc_dct_j, spc_model=spc_model,
            pf_levels=pf_levels,
            spc_save_fs=spc_save_fs,
            pst_params=pst_params)
        ts_str += '\n' + mess_io.writer.ts_sadpt(
            ts_label, reac_label, prod_label, pst_str, zero_energy)
    else:
        ts_reac_barr = ts_ene - reac_ene
        ts_prod_barr = ts_ene - prod_ene
        if ts_reac_barr < 0.:
            ts_reac_barr = 0.1
        if ts_prod_barr < 0.:
            ts_prod_barr = 0.1
        tunnel_str = mess_io.writer.tunnel_eckart(
            imag_freq, ts_reac_barr, ts_prod_barr)
        ts_str += '\n' + mess_io.writer.ts_sadpt(
            ts_label, reac_label, prod_label,
            species_data[tsname], zero_energy, tunnel_str)

    return [well_str, bim_str, ts_str]


# Readers
def read_rates(rct_lab, prd_lab, mess_path, assess_pdep_temps,
               pdep_tolerance=20.0, no_pdep_pval=1.0,
               pdep_low=None, pdep_high=None, bimol=False):
    """ Read the rate constants from the MESS output and
        (1) filter out the invalid rates that are negative or undefined
        and obtain the pressure dependent values
    """

    # Dictionaries to store info; indexed by pressure (given in fit_ps)
    calc_k_dct = {}
    valid_calc_tk_dct = {}
    ktp_dct = {}

    # Read the MESS output file into a string
    with open(mess_path+'/rate.out', 'r') as mess_file:
        output_string = mess_file.read()
    # with open(mess_path+'/mess.inp', 'r') as mess_file:
    #     input_string = mess_file.read()

    # Read the temperatures and pressures out of the MESS output
    mess_temps, _ = mess_io.reader.rates.get_temperatures(
        output_string)
    mess_pressures, punit = mess_io.reader.rates.get_pressures(
        output_string)
    # mess_pressures, punit = mess_io.reader.rates.get_pressures_input(
    #     input_string)

    # Loop over the pressures obtained from the MESS output
    for pressure in mess_pressures:

        # Read the rate constants
        if pressure == 'high':
            rate_ks = mess_io.reader.highp_ks(
                output_string, rct_lab, prd_lab)
        else:
            rate_ks = mess_io.reader.pdep_ks(
                output_string, rct_lab, prd_lab, pressure, punit)

        # Store in a dictionary
        calc_k_dct[pressure] = rate_ks

    # Remove k(T) vals at each P where where k is negative or undefined
    # If ANY valid k(T,P) vals at given pressure, store in dct
    for pressure, calc_ks in calc_k_dct.items():
        filtered_temps, filtered_ks = ratefit.fit.get_valid_tk(
            mess_temps, calc_ks, bimol)
        if filtered_ks.size > 0:
            valid_calc_tk_dct[pressure] = numpy.concatenate(
                (filtered_temps, filtered_ks))

    # Filter the ktp dictionary by assessing the presure dependence
    if list(valid_calc_tk_dct.keys()) == ['high']:
        ktp_dct['high'] = valid_calc_tk_dct['high']
    else:
        rxn_is_pdependent = ratefit.calc.assess_pressure_dependence(
            valid_calc_tk_dct, assess_pdep_temps,
            tolerance=pdep_tolerance, plow=pdep_low, phigh=pdep_high)
        if rxn_is_pdependent:
            # Set dct to fit as copy of dct to do PLOG fits at all pressures
            ktp_dct = copy.deepcopy(valid_calc_tk_dct)
        else:
            # Set dct to have single set of k(T, P) vals: P is desired pressure
            if no_pdep_pval in valid_calc_tk_dct:
                ktp_dct['high'] = valid_calc_tk_dct[no_pdep_pval]
    if 'high' not in ktp_dct and 'high' in valid_calc_tk_dct.keys():
        ktp_dct['high'] = valid_calc_tk_dct['high']

    return ktp_dct
