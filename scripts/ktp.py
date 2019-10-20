""" reaction list test
"""

import os
import copy
import thermo
import automol
import autofile
import moldr
import mess_io
import mess_io.writer
import ratefit
import scripts
from datalibs import phycon
from submission import substr


def pf_headers(
        rct_ichs, temps, press, exp_factor, exp_power, exp_cutoff, eps1, eps2,
        sig1, sig2, mass1):
    """ makes the standard header and energy transfer sections for MESS input file
    """
    # header section
    header_str = mess_io.writer.global_reaction(temps, press)
    print(header_str)
    tot_mass = 0.
    for rct_ich in rct_ichs:
        geo = automol.convert.inchi.geometry(rct_ich)
        masses = automol.geom.masses(geo)
        for mass in masses:
            tot_mass += mass

    # energy transfer section
    energy_trans_str = mess_io.writer.energy_transfer(
        exp_factor, exp_power, exp_cutoff, eps1, eps2, sig1, sig2, mass1, tot_mass)

    return header_str, energy_trans_str


def make_all_species_data(rxn_lst, spc_dct, save_prefix, model_info, pf_info, ts_found, projrot_script_str):
    """ generate the MESS species blocks for all the species
    """
    species = {}
    spc_save_fs = autofile.fs.species(save_prefix)
    for idx, rxn in enumerate(rxn_lst):
        tsname = 'ts_{:g}'.format(idx)
        if tsname in ts_found:
            ts = spc_dct[tsname]
            specieslist = rxn['reacs'] + rxn['prods']
            for name in specieslist:
                if not name in species:
                    species[name], _ = make_species_data(
                        name, spc_dct[name], spc_save_fs, model_info, pf_info, projrot_script_str)
            if not 'radical radical addition' in spc_dct[tsname]['class']:
                species[tsname], spc_dct[tsname]['imag_freq'] = make_species_data(
                    tsname, ts, save_prefix, model_info, pf_info, projrot_script_str)
    return species


def make_species_data(spc, spc_dct_i, spc_save_fs, spc_model, pf_levels, projrot_script_str):
    """ makes the main part of the MESS species block for a given species
    """
    spc_info = (spc_dct_i['ich'], spc_dct_i['chg'], spc_dct_i['mul'])
    if 'ts_' in spc:
        save_path = spc_dct_i['rxn_fs'][3]
    else:
        spc_save_fs.leaf.create(spc_info)
        save_path = spc_save_fs.leaf.path(spc_info)
    species_data = moldr.pf.species_block(
        spc=spc,
        spc_dct_i=spc_dct_i,
        spc_info=spc_info,
        spc_model=spc_model,
        pf_levels=pf_levels,
        projrot_script_str=projrot_script_str,
        elec_levels=[[0., 1]], sym_factor=1.,
        save_prefix=save_path,
        )
    return species_data


def make_fake_species_data(spc_dct_i, spc_dct_j, spc_save_fs, pf_levels, projrot_script_str):
    """ make a fake MESS species block to represent the van der Waals well
    arising from the combination of two fragment species
    """
    spc_info_i = (spc_dct_i['ich'], spc_dct_i['chg'], spc_dct_i['mul'])
    spc_info_j = (spc_dct_j['ich'], spc_dct_j['chg'], spc_dct_j['mul'])
    spc_save_fs.leaf.create(spc_info_i)
    spc_save_fs.leaf.create(spc_info_j)
    save_path_i = spc_save_fs.leaf.path(spc_info_i)
    save_path_j = spc_save_fs.leaf.path(spc_info_j)
    if pf_levels[3]:
        spc_model = ['RIGID', 'HARM', 'SAMPLING']
    else:
        spc_model = ['RIGID', 'HARM', '']
    species_data = moldr.pf.fake_species_block(
        spc_dct_i=spc_dct_i,
        spc_dct_j=spc_dct_j,
        spc_info_i=spc_info_i,
        spc_info_j=spc_info_j,
        spc_model=spc_model,
        pf_levels=pf_levels,
        projrot_script_str=projrot_script_str,
        elec_levels=[[0., 1]], sym_factor=1.,
        save_prefix_i=save_path_i,
        save_prefix_j=save_path_j
        )
    return species_data


def make_channel_pfs(
        tsname, rxn, species_data, spc_dct, idx_dct, strs, first_ground_ene,
        spc_save_fs, pf_levels, multi_info, projrot_script_str):
    """ make the partition function strings for each of the channels
    includes strings for each of the unimolecular wells, bimolecular fragments, and
    transition states connecting them.
    It also includes a special treatment for abstraction to include phase space blocks
    coupling bimolecular fragments to fake van der Waals wells
    """
    bim_str, well_str, ts_str = strs
#Find the number of uni and bimolecular wells already in the dictionary
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
#Set up a new well for the reactants if that combo isn't already in the dct
    reac_label = ''
    reac_ene = 0.
    bimol = False
    if len(rxn['reacs']) > 1:
        bimol = True
    well_data = []
    spc_label = []
    for reac in rxn['reacs']:
        spc_label.append(automol.inchi.smiles(spc_dct[reac]['ich']))
        well_data.append(species_data[reac])
        reac_ene += scripts.thermo.spc_energy(spc_dct[reac]['ene'], spc_dct[reac]['zpe']) * phycon.EH2KCAL
    well_dct_key1 = '+'.join(rxn['reacs'])
    well_dct_key2 = '+'.join(rxn['reacs'][::-1])
    if not well_dct_key1 in idx_dct:
        if well_dct_key2 in idx_dct:
            well_dct_key1 = well_dct_key2
        else:
            if bimol:
                reac_label = 'P' + str(pidx)
                pidx += 1
                if not first_ground_ene:
                    first_ground_ene = reac_ene
                ground_energy = reac_ene - first_ground_ene
                bim_str += ' \t ! {} + {} \n'.format(rxn['reacs'][0], rxn['reacs'][1])
                bim_str += mess_io.writer.bimolecular(
                    reac_label, spc_label[0], well_data[0], spc_label[1],
                    well_data[1], ground_energy)
                idx_dct[well_dct_key1] = reac_label
            else:
                if not first_ground_ene:
                    first_ground_ene = reac_ene
                reac_label = 'W' + str(widx)
                widx += 1
                zero_energy = reac_ene - first_ground_ene
                well_str += ' \t ! {} \n'.format(rxn['reacs'][0])
                well_str += mess_io.writer.well(reac_label, well_data[0], zero_energy)
                idx_dct[well_dct_key1] = reac_label
    if not reac_label:
        reac_label = idx_dct[well_dct_key1]
#Set up a new well for the products if that combo isn't already in the dct
    prod_label = ''
    prod_ene = 0.
    bimol = False
    if len(rxn['prods']) > 1:
        bimol = True
    well_data = []
    spc_label = []
    for prod in rxn['prods']:
        spc_label.append(automol.inchi.smiles(spc_dct[prod]['ich']))
        well_data.append(species_data[prod])
        prod_ene += scripts.thermo.spc_energy(spc_dct[prod]['ene'], spc_dct[prod]['zpe']) * phycon.EH2KCAL
    zero_energy = prod_ene - reac_ene
    well_dct_key1 = '+'.join(rxn['prods'])
    well_dct_key2 = '+'.join(rxn['prods'][::-1])
    if not well_dct_key1 in idx_dct:
        if well_dct_key2 in idx_dct:
            well_dct_key1 = well_dct_key2
        else:
            if bimol:
                prod_label = 'P' + str(pidx)
                ground_energy = prod_ene - first_ground_ene
                bim_str += ' \t ! {} + {} \n'.format(
                    rxn['prods'][0], rxn['prods'][1])
                bim_str += mess_io.writer.bimolecular(
                    prod_label, spc_label[0], well_data[0], spc_label[1], well_data[1],
                    ground_energy)
                idx_dct[well_dct_key1] = prod_label
            else:
                prod_label = 'W' + str(widx)
                zero_energy = prod_ene - first_ground_ene
                well_str += ' \t ! {} \n'.format(rxn['prods'][0])
                well_str += mess_io.writer.well(prod_label, well_data[0], zero_energy)
                idx_dct[well_dct_key1] = prod_label
    if not prod_label:
        prod_label = idx_dct[well_dct_key1]

    #Set up a new well connected to ts
    if 'radical radical' in spc_dct[tsname]['class']:
        # for radical radical call vtst or vrctst
        ts_label = 'B' + str(int(tsname.replace('ts_', ''))+1)
        if 'P' in reac_label:
            spc_ene = reac_ene - first_ground_ene
            spc_zpe = spc_dct[rxn['reacs'][0]]['zpe'] + spc_dct[rxn['reacs'][1]]['zpe']
        else:
            spc_ene = prod_ene - first_ground_ene
            spc_zpe = spc_dct[rxn['prods'][0]]['zpe'] + spc_dct[rxn['prods'][1]]['zpe']
        ts_str += '\n' + moldr.pf.vtst_with_no_saddle_block(
            spc_dct[tsname], ts_label, reac_label, prod_label, spc_ene, spc_zpe, projrot_script_str,
            multi_info, elec_levels=[[0., 1]], sym_factor=1.
            )

    else:
        # for tight TS's make ts_str
        # for abstraction make fake wells and PST TSs
        ts_ene = scripts.thermo.spc_energy(spc_dct[tsname]['ene'], spc_dct[tsname]['zpe']) * phycon.EH2KCAL
        zero_energy = ts_ene - first_ground_ene
        ts_label = 'B' + str(int(tsname.replace('ts_', ''))+1)
        imag_freq = 0
        if 'imag_freq' in spc_dct[tsname]:
            imag_freq = abs(spc_dct[tsname]['imag_freq'])
        if not imag_freq:
            print('No imaginary freq for ts: {}'.format(tsname))


        fake_wellr_label = ''
        fake_wellp_label = ''
        if 'abstraction' in spc_dct[tsname]['class']:
        #Make fake wells and PST TSs as needed
            well_dct_key1 = 'F' + '+'.join(rxn['reacs'])
            well_dct_key2 = 'F' + '+'.join(rxn['reacs'][::-1])
            if not well_dct_key1 in idx_dct:
                if well_dct_key2 in idx_dct:
                    well_dct_key1 = well_dct_key2
                else:
                    fake_wellr_label = 'F' + str(fidx)
                    fidx += 1
                    vdwr_ene = reac_ene - 1.0
                    zero_energy = vdwr_ene - first_ground_ene
                    well_str += ' \t ! Fake Well for {}\n'.format('+'.join(rxn['reacs']))
                    fake_wellr = make_fake_species_data(
                        spc_dct[rxn['reacs'][0]], spc_dct[rxn['reacs'][1]],
                        spc_save_fs, pf_levels, projrot_script_str)
                    well_str += mess_io.writer.well(fake_wellr_label, fake_wellr, zero_energy)
                    idx_dct[well_dct_key1] = fake_wellr_label

                    pst_r_label = 'FRB' + str(int(tsname.replace('ts_', ''))+1)
                    idx_dct[well_dct_key1.replace('F', 'FRB')] = pst_r_label
                    spc_dct_i = spc_dct[rxn['reacs'][0]]
                    spc_dct_j = spc_dct[rxn['reacs'][1]]
                    if pf_levels[3]:
                        spc_model = ['RIGID', 'HARM', 'SAMPLING']
                    else:
                        spc_model = ['RIGID', 'HARM', '']
                    pst_r_ts_str = moldr.pf.pst_block(
                        spc_dct_i, spc_dct_j, spc_model=spc_model,
                        pf_levels=pf_levels, projrot_script_str=projrot_script_str,
                        spc_save_fs=spc_save_fs)
            print('fake_wellr_label test:', fake_wellr_label)
            if not fake_wellr_label:
                print('well_dct_key1 test:', well_dct_key1)
                fake_wellr_label = idx_dct[well_dct_key1]
                pst_r_label = idx_dct[well_dct_key1.replace('F', 'FRB')]
            zero_energy = reac_ene - first_ground_ene
            tunnel_str = ''
            ts_str += '\n' + mess_io.writer.ts_sadpt(
                pst_r_label, reac_label, fake_wellr_label, pst_r_ts_str,
                zero_energy, tunnel_str)

            well_dct_key1 = 'F' + '+'.join(rxn['prods'])
            well_dct_key2 = 'F' + '+'.join(rxn['prods'][::-1])
            if not well_dct_key1 in idx_dct:
                if well_dct_key2 in idx_dct:
                    well_dct_key1 = well_dct_key2
                else:
                    fake_wellp_label = 'F' + str(fidx)
                    fidx += 1
                    vdwp_ene = prod_ene - 1.0
                    zero_energy = vdwp_ene - first_ground_ene
                    well_str += ' \t ! Fake Well for {}\n'.format('+'.join(rxn['prods']))
                    fake_wellp = make_fake_species_data(
                        spc_dct[rxn['prods'][0]], spc_dct[rxn['prods'][1]],
                        spc_save_fs, pf_levels, projrot_script_str)
                    well_str += mess_io.writer.well(fake_wellp_label, fake_wellp, zero_energy)
                    idx_dct[well_dct_key1] = fake_wellp_label

                    pst_p_label = 'FPB' + str(int(tsname.replace('ts_', ''))+1)
                    idx_dct[well_dct_key1.replace('F', 'FPB')] = pst_p_label
                    spc_dct_i = spc_dct[rxn['prods'][0]]
                    spc_dct_j = spc_dct[rxn['prods'][1]]
                    if pf_levels[3]:
                        spc_model = ['RIGID', 'HARM', 'SAMPLING']
                    else:
                        spc_model = ['RIGID', 'HARM', '']
                    pst_p_ts_str = moldr.pf.pst_block(
                        spc_dct_i, spc_dct_j, spc_model=spc_model,
                        pf_levels=pf_levels, projrot_script_str=projrot_script_str,
                        spc_save_fs=spc_save_fs)
            if not fake_wellp_label:
                fake_wellp_label = idx_dct[well_dct_key1]
                pst_p_label = idx_dct[well_dct_key1.replace('F', 'FPB')]
            zero_energy = prod_ene - first_ground_ene
            tunnel_str = ''
            ts_str += '\n' + mess_io.writer.ts_sadpt(
                pst_p_label, prod_label, fake_wellp_label, pst_p_ts_str,
                zero_energy, tunnel_str)

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
                ts_label, fake_wellr_label, fake_wellp_label, species_data[tsname], zero_energy, tunnel_str)
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
                ts_label, reac_label, prod_label, species_data[tsname], zero_energy, tunnel_str)

    return [well_str, bim_str, ts_str], first_ground_ene


def run_rates(
        header_str, energy_trans_str, well_str, bim_str, ts_str, tsdct,
        thy_info, rxn_save_path):
    """ Generate k(T,P) by first compiling all the MESS strings and then running MESS
    """
    ts_info = (tsdct['ich'], tsdct['chg'], tsdct['mul'])
    orb_restr = moldr.util.orbital_restriction(ts_info, thy_info)
    ref_level = thy_info[1:3]
    ref_level.append(orb_restr)
    print('ref level test:', ref_level)
    thy_save_fs = autofile.fs.theory(rxn_save_path)
    thy_save_fs.leaf.create(ref_level)
    thy_save_path = thy_save_fs.leaf.path(ref_level)

    mess_inp_str = '\n'.join([header_str, energy_trans_str, well_str, bim_str, ts_str])
    print('mess input file')
    print(mess_inp_str)

    bld_locs = ['MESS', 0]
    bld_save_fs = autofile.fs.build(thy_save_path)
    bld_save_fs.leaf.create(bld_locs)
    mess_path = bld_save_fs.leaf.path(bld_locs)
    print('Build Path for MESS rate files:')
    print(mess_path)
    with open(os.path.join(mess_path, 'mess.inp'), 'w') as mess_file:
        mess_file.write(mess_inp_str)
    moldr.util.run_script(substr.MESSRATE, mess_path)
    return mess_path


def read_rates(rct_lab, prd_lab, mess_path, assess_pdep_temps,
               pdep_tolerance=20.0, no_pdep_pval=1.0,
               pdep_low=None, pdep_high=None):
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
        filtered_temps, filtered_ks = ratefit.fit.util.get_valid_tk(
            mess_temps, calc_ks)
        if filtered_ks.size > 0:
            valid_calc_tk_dct[pressure] = [filtered_temps, filtered_ks]

    # Filter the ktp dictionary by assessing the presure dependence
    rxn_is_pdependent = ratefit.err.assess_pressure_dependence(
        valid_calc_tk_dct, assess_pdep_temps,
        tolerance=pdep_tolerance, plow=pdep_low, phigh=pdep_high)

    if rxn_is_pdependent:
        # Set dct to fit as copy of dct to do PLOG fits at all pressures
        ktp_dct = copy.deepcopy(valid_calc_tk_dct)
    else:
        # Set dct to have single set of k(T, P) vals: P is desired pressure
        ktp_dct['high'] = valid_calc_tk_dct[no_pdep_pval]
        # premap = no_pdep_pval

    return ktp_dct


def mod_arr_fit(ktp_dct, mess_path, fit_type='single', fit_method='dsarrfit',
                t_ref=1.0, a_conv_factor=1.0):
    """
    Routine for a single reaction:
        (1) Grab high-pressure and pressure-dependent rate constants
            from a MESS output file
        (2) Fit rate constants to an Arrhenius expression
    """

    assert fit_type in ('single', 'double')
    assert fit_method == 'dsarrfit'

    # Dictionaries to store info; indexed by pressure (given in fit_ps)
    fit_param_dct = {}
    fit_k_dct = {}
    fit_err_dct = {}

    # Calculate the fitting parameters from the filtered T,k lists
    for pressure, tk_lsts in ktp_dct.items():

        # Set the temperatures and rate constants
        temps = numpy.array(tk_lsts[0])
        rate_constants = numpy.array(tk_lsts[1])

        # Fit rate constants using desired Arrhenius fit
        if fit_type == 'single':
            fit_params = ratefit.fit.arrhenius.single(
                temps, rate_constants, t_ref, fit_method,
                dsarrfit_path=mess_path, a_conv_factor=a_conv_factor)
        elif fit_type == 'double':
            fit_params = ratefit.fit.arrhenius.double(
                temps, rate_constants, t_ref, fit_method,
                dsarrfit_path=mess_path, a_conv_factor=a_conv_factor)

        # Store the fitting parameters in a dictionary
        fit_param_dct[pressure] = fit_params

    # Calculate fitted rate constants using the fitted parameters
    for pressure, params in fit_param_dct.items():

        # Set the temperatures
        # temps = valid_calc_tk_dct[pressure][0]
        temps = numpy.array(ktp_dct[pressure][0])

        # Calculate fitted rate constants, based on fit type
        if fit_type == 'single':
            fit_ks = ratefit.fxns.single_arrhenius(
                params[0], params[1], params[2],
                t_ref, temps)
        elif fit_type == 'double':
            fit_ks = ratefit.fxns.double_arrhenius(
                params[0], params[1], params[2],
                params[3], params[4], params[5],
                t_ref, temps)

        # Store the fitting parameters in a dictionary
        fit_k_dct[pressure] = fit_ks / a_conv_factor

    # Calculute the error between the calc and fit ks
    for pressure, fit_ks in fit_k_dct.items():

        # if pressure == 'high' and not rxn_is_pdependent:
        #     pressure = premap
        # calc_ks = valid_calc_tk_dct[pressure][1]
        calc_ks = ktp_dct[pressure][1]
        mean_avg_err, max_avg_err = ratefit.err.calc_sse_and_mae(
            calc_ks, fit_ks)

        # Store in a dictionary
        fit_err_dct[pressure] = [mean_avg_err, max_avg_err]

    # If no high-pressure is found in dct, add fake fitting parameters
    if 'high' not in fit_param_dct:
        if fit_type == 'single':
            fit_param_dct['high'] = [1.00, 0.00, 0.00]
        elif fit_type == 'double':
            fit_param_dct['high'] = [1.00, 0.00, 0.00, 1.00, 0.00, 0.00]

    return fit_param_dct, fit_err_dct


def species_thermo(
        spc_names,
        spc_info,
        spc_ref_names,
        elc_deg_dct,
        temp_par,
        ref_high_level,
        run_high_levels,
        high_level_coeff,
        spc_models,
        pf_levels,
        save_prefix,
        projrot_script_str,
        pf_script_str,
        ):
    """Evaluate the thermodynamics properties for all the species
    """

    # initializations
    temp_step = temp_par[0]
    ntemps = temp_par[1]
    chemkin_poly_strs = {}
    for hl_idx, _ in enumerate(run_high_levels):
        chemkin_poly_strs[hl_idx] = ''
    spc_ref_ich = []
    for ref_name in spc_ref_names:
        spc_ref_ich.append(spc_info[ref_name][0])

    # generate and store the energies so that the reference energies can be used as needed
    ene_hl = {}
    print('Evaluating thermo for following species')
    print(spc_names)
    for name in spc_names:
        # set up species information
        ich = spc_info[name][0]
        smi = automol.inchi.smiles(ich)
        print("smiles: {}".format(smi), "inchi: {}".format(ich))

        # get and store the high level energies
        print('name in high level energy routine')
        print(name)
        ene_cbs[name] = 0.
        for hl_idx, _ in enumerate(run_high_levels):
            min_ene = moldr.pf.get_high_level_energy(
                spc_info=spc_info[name],
                theory_low_level=ref_high_level,
                theory_high_level=run_high_levels[hl_idx],
                save_prefix=save_prefix,)
            ene_hl[(name, hl_idx)] = min_ene
            if high_level_coeff:
                ene_cbs[name] += min_ene*high_level_coeff(hl_idx)
        ene_hl[(name, hl_idx+1)] = ene_cbs[name]

    spc_save_fs = autofile.fs.species(save_prefix)
    for tors_model, vib_model in spc_models:
    # rot_model = RRHO or 1DHR or MDHR or TAU
        print('tors_model:', tors_model)
        print('vib_model:', vib_model)
        for har_level, tors_level, vpt2_level in pf_levels:
            # get the zero-point energy for each species
            print('harmonic_level:', har_level)
            print('tors_level:', tors_level)
            print('vpt2_level:', vpt2_level)
            print('Calculating zpe')
            spc_zpe = {}
            is_atom = {}
            zero_energy_str = {}
            for name in spc_names:
                ich = spc_info[name][0]
                smi = automol.inchi.smiles(ich)
                print("smiles: {}".format(smi), "inchi: {}".format(ich))
                spc_save_fs.leaf.create(spc_info[name])
                spc_save_path = spc_save_fs.leaf.path(spc_info[name])

                spc_zpe[name], is_atom[name] = moldr.pfdriver.get_zero_point_energy(
                    name,
                    spc_info[name],
                    tors_model, vib_model,
                    har_level, tors_level, vpt2_level,
                    script_str=pf_script_str,
                    elec_levels=[[0., 1]], sym_factor=1.,
                    save_prefix=spc_save_path)
                print(name, spc_zpe[name])
                zpe_str = '{0:<8.2f}\n'.format(spc_zpe[name])
                if is_atom[name]:
                   zero_energy_str[name] = 'End'
                else:
                   zero_energy_str[name] = ' ZeroEnergy[kcal/mol] ' + zpe_str
                   zero_energy_str[name] += 'End'
                print('zero_energy_str:', zero_energy_str)

            pf_inp_str = {}
            print('finished zpe')
            # get the partition function for each species
            
            chemkin_poly_strs = ['' for i in range(len(run_high_levels))]
            for name in spc_names:
                # set up species information
                ich = spc_info[name][0]
                smi = automol.inchi.smiles(ich)
                mul = spc_info[name][2]
                print("smiles: {}".format(smi), "inchi: {}".format(ich))
                spc_save_fs.leaf.create(spc_info[name])
                spc_save_path = spc_save_fs.leaf.path(spc_info[name])

                # to be generalized
                sym_factor = 1.
                elec_levels = [[0., mul]]
                if (ich, mul) in elc_deg_dct:
                    elec_levels = elc_deg_dct[(ich, mul)]

                # cycle through the low levels generating partition functions  for each
                spc_str = moldr.pfdriver.species_block(
                    spc_info=spc_info[name],
                    tors_model=tors_model,
                    vib_model=vib_model,
                    har_level=har_level,
                    tors_level=tors_level,
                    vpt2_level=vpt2_level,
                    script_str=projrot_script_str,
                    elec_levels=elec_levels,
                    sym_factor=sym_factor,
                    save_prefix=spc_save_path,
                    )

                # create a messpf input file
                global_pf_str = mess_io.writer.write_global_pf(
                    [], temp_step, ntemps, rel_temp_inc=0.001,
                    atom_dist_min=0.6)
                print(global_pf_str)
                spc_head_str = 'Species ' + name
                print(spc_head_str)
                pf_inp_str[name] = '\n'.join(
                    [global_pf_str, spc_head_str,
                     spc_str, zero_energy_str[name], '\n'])
                print(spc_str)
                print(pf_inp_str[name])

                orb_restr = moldr.util.orbital_restriction(
                        spc_info[name], tors_level)
                tors_levelp = tors_level[1:3]
                tors_levelp.append(orb_restr)

                thy_save_fs = autofile.fs.theory(spc_save_path)
                thy_save_fs.leaf.create(tors_levelp)
                thy_save_path = thy_save_fs.leaf.path(tors_levelp)
                bld_locs = ['PF', 0]
                bld_save_fs = autofile.fs.build(thy_save_path)
                bld_save_fs.leaf.create(bld_locs)
                pf_path = bld_save_fs.leaf.path(bld_locs)
                print('Build Path for Partition Functions')
                print(pf_path)

                # run messpf
                with open(os.path.join(pf_path, 'pf.inp'), 'w') as pf_file:
                    pf_file.write(pf_inp_str[name])
                moldr.util.run_script(pf_script_str, pf_path)
                print('finished partition function')

                formula = automol.inchi.formula(ich)
                print('\nformula:')
                print(formula)

                # Get atom count dictionary
                atom_dict = thermo.util.get_atom_counts_dict(formula)
                print('\natom dict:')
                print(atom_dict)

                # Get the list of the basis
                spc_bas = thermo.heatform.get_reduced_basis(spc_ref_ich, formula)
                print('\nbasis:')
                print(spc_bas)

                # Get the coefficients for the balanced heat-of-formation eqn
                coeff = thermo.heatform.calc_coefficients(spc_bas, atom_dict)
                print('\ncoeff:')
                print(coeff)

                # prepare NASA polynomials
                nasa_inp_str = ('nasa')
                bld_locs = ['NASA_POLY', 0]
                bld_save_fs.leaf.create(bld_locs)
                nasa_path = bld_save_fs.leaf.path(bld_locs)
                print('NASA build path')
                print(nasa_path)

                h_basis = []
                if high_level_coeff:
                    len_high_levels = len(run_high_levels)+1
                else:
                    len_high_levels = len(run_high_levels)
                for hl_idx in range(len_high_levels):
#                for hl_idx, _ in enumerate(run_high_levels):
                    for  name_ref in spc_ref_names:
                        ich = spc_info[name_ref][0]
                        if ich in spc_bas:
                            tmp = ene_hl[(name_ref, hl_idx)] + spc_zpe[name_ref]/phycon.EH2KCAL
                            h_basis.append(tmp)
                    print('\ne_basis:')
                    print(h_basis)

                    # Get the 0 K heat of formation
                    spc_ene = ene_hl[(name, hl_idx)] + spc_zpe[name]/phycon.EH2KCAL
                    spc_ene_cbs = ene_cbs + spc_zpe[name]/phycon.EH2KCAL
                    h0form = thermo.heatform.calc_hform_0k(spc_ene, h_basis, spc_bas, coeff, ref_set='ATcT')
                    print('h0form = ',h0form)

                    # need to change back to starting directory after running thermp and pac99 or rest of code is confused
                    starting_path = os.getcwd()
                    os.chdir(nasa_path)

                    # Write thermp input file
                    ENTHALPYT = 0.
                    BREAKT = 1000.
                    thermo.runner.write_thermp_input(
                        formula=formula,
                        deltaH=h0form,
                        enthalpyT=ENTHALPYT,
                        breakT=BREAKT,
                        thermp_file_name='thermp.dat')

                    # Run thermp
                    thermo.runner.run_thermp(
                        pf_path=pf_path,
                        thermp_path=nasa_path,
                        thermp_file_name='thermp.dat',
                        pf_file_name='pf.dat'
                        )

                    # Run pac99
                    print('formula test')
                    print(formula)
                    print(nasa_path)
                    thermo.runner.run_pac99(nasa_path, formula)

                    with open(os.path.join(nasa_path, 'thermp.out'), 'r') as thermp_outfile:
                        thermp_out_str = thermp_outfile.read()

                    with open(os.path.join(nasa_path, formula+'.c97'), 'r') as pac99_file:
                        pac99_str = pac99_file.read()

                    # Get the pac99 polynomial
                    pac99_poly_str = thermo.nasapoly.get_pac99_polynomial(pac99_str)
                    print('\nPAC99 Polynomial:')
                    print(pac99_poly_str)

                    # Convert the pac99 polynomial to chemkin polynomial
                    comment_str = '! tors model: {0}\n'.format(tors_model)
                    comment_str += '! vib model: {0}\n'.format(vib_model)
                    comment_str += '! har level: {0}{1}/{2}\n'.format(
                        har_level[3], har_level[1], har_level[2])
                    comment_str += '! tors level: {0}{1}/{2}\n'.format(
                        tors_level[3], tors_level[1], tors_level[2])
                    comment_str += '! vpt2 level: {0}{1}/{2}\n'.format(
                        vpt2_level[3], vpt2_level[1], vpt2_level[2])
                    comment_str += '! ref level for energy: {0}{1}/{2}\n'.format(
                        ref_high_level[3], ref_high_level[1], ref_high_level[2])
                    if hl_idx == len_high_levels-1:
                        comment_str += '! energy level from following sum of high levels:\n'
                        for hl_idx2, _ in enumerate(run_high_levels):
                            comment_str += '! energy level: {0}/{1}\n'.format(
                            run_high_levels[hl_idx2][1], run_high_levels[hl_idx2][2])
                    else:
                        comment_str += '! energy level: {0}/{1}\n'.format(
                        run_high_levels[hl_idx][1], run_high_levels[hl_idx][2])

                    chemkin_poly_str = thermo.nasapoly.convert_pac_to_chemkin(
                        name, atom_dict, comment_str, pac99_poly_str)
                    chemkin_set_str = thermo.nasapoly.convert_pac_to_chemkin(
                        name, atom_dict, '', pac99_poly_str)
                    print('\nCHEMKIN Polynomial:')
                    print(chemkin_poly_str)
                    print(hl_idx)
                    if chemkin_poly_strs[hl_idx] == '':
                        chemkin_poly_strs[hl_idx] += chemkin_poly_str
                    else:
                        chemkin_poly_strs[hl_idx] += chemkin_set_str
                    print('startig_path in thermo')
                    print(starting_path)
                    ckin_path = ''.join([starting_path, '/ckin'])
                    print(ckin_path)
                    os.chdir(starting_path)
                    with open(os.path.join(nasa_path, formula+'.ckin'), 'w') as nasa_file:
                        nasa_file.write(chemkin_poly_str)
                    with open(os.path.join(ckin_path, formula+'.ckin'), 'w') as nasa_file:
                        nasa_file.write(chemkin_poly_str)

            for hl_idx in range(len_high_levels):
#            for hl_idx, _ in enumerate(run_high_levels):
                hl_idx_str = str(hl_idx)
                with open(os.path.join(ckin_path, 'SPECIES'+hl_idx_str+'.ckin'), 'w') as nasa_file:
                    nasa_file.write(chemkin_poly_strs[hl_idx])


def reaction_rates(
        rxn_info_lst,
        smi_dct,
        chg_dct,
        mul_dct,
        elc_deg_dct,
        temperatures,
        pressures,
        etsfr_par,
        run_opt_levels,
        ref_high_level,
        run_high_levels,
        high_level_coeff,
        spc_models,
        pf_levels,
        save_prefix,
        projrot_script_str,
        rate_script_str,
    ):
    """Evaluate the reaction rates for given ab initio levels and theory models
    """

    idx = 0
    exp_factor = etsfr_par[idx]
    idx += 1

    exp_power = etsfr_par[idx]
    idx += 1

    exp_cutoff = etsfr_par[idx]
    idx += 1

    eps1 = etsfr_par[idx]
    idx += 1

    eps2 = etsfr_par[idx]
    idx += 1

    sig1 = etsfr_par[idx]
    idx += 1

    sig2 = etsfr_par[idx]
    idx += 1

    mass1 = etsfr_par[idx]
    idx += 1

    # set mass2 to mass of molecule
#    mass2 = etsfr_par[idx]
#    idx += 1

    # determine the number of pes

    formula_str_old = ''
    formula_pes = []
    pes_idx = -1
    nch_pes = []
    for formula_str, _, _, _  in rxn_info_lst:
        if formula_str != formula_str_old:
            pes_idx += 1
            formula_str_old = formula_str
            nch_pes.append(1)
            formula_pes.append(formula_str)
        else:
            nch_pes[pes_idx] += 1

    num_pes = pes_idx+1

    print('number of rxns:', sum(nch_pes))
    print('number of pes:', num_pes)
    for pes_idx in range(0, num_pes):
        print('There are ', nch_pes[pes_idx], 'channels on pes ', pes_idx+1)

    # spc_models, tors_model, and vib_model taken as outermost loops 
    # facilitates formation of full PES data and thus ME calculation or each model
    for tors_model, vib_model in spc_models:
        # tors_model = RRHO or 1DHR or MDHR or TAU
        print('tors_model:', tors_model)
        # vib_model = HARM or VPT2
        print('vib_model:', vib_model)
        for har_level, tors_level, vpt2_level in pf_levels:
            # get the zero-point energy for each species
            print('harmonic_level:', har_level)
            print('tors_level:', tors_level)
            print('vpt2_level:', vpt2_level)
            for opt_level_idx, _ in enumerate(run_opt_levels):
                # for prog, method, basis in run_opt_levels:
                prog = run_opt_levels[opt_level_idx][0]
                method = run_opt_levels[opt_level_idx][1]
                _, _, KWARGS, _ = (
                    moldr.util.run_qchem_par(prog, method))

                ts_orb_restr = moldr.util.orbital_restriction(
                    ts_info, run_opt_levels[opt_level_idx])
                thy_level = run_opt_levels[opt_level_idx][0:3]
                thy_level.append(ts_orb_restr)
                idx = -1
                for pes_idx in range(0, num_pes):
                    formula = formula_pes[pes_idx]
                    print()
                    print("Calcuating rate for {} Potential Energy Surface".format(formula))
                    print('nch test:', nch_pes[pes_idx])
                    for chn_idx in range(0,nch_pes[pes_idx]):
                        idx += 1
                        print('idx test:', idx)
                        _, rct_names, prd_names, rxn_name = rxn_info_lst[idx]
                        print('rct_names test:', rct_names)
                        print('prd_names test:', prd_names)
                        print('rxn_name test:', rxn_name)

                        # print the CHEMKIN reaction name for reference
                        print("Reaction: {}".format(rxn_name))

                        # determine inchis, charges, and multiplicities

                        rct_smis = list(map(smi_dct.__getitem__, rct_names))
                        prd_smis = list(map(smi_dct.__getitem__, prd_names))
                        rct_ichs = list(map(automol.smiles.inchi, rct_smis))
                        prd_ichs = list(map(automol.smiles.inchi, prd_smis))
                        rct_chgs = list(map(chg_dct.__getitem__, rct_names))
                        prd_chgs = list(map(chg_dct.__getitem__, prd_names))
                        rct_muls = list(map(mul_dct.__getitem__, rct_names))
                        prd_muls = list(map(mul_dct.__getitem__, prd_names))

                        ts_mul = automol.mult.ts.low(rct_muls, prd_muls)
                        ts_chg = sum(rct_chgs)
                        print('ts_chg test:',ts_chg)
                        ts_info = ('', ts_chg, ts_mul)

                        # print header and energy transfer sections only for the first channel
                        if chn_idx == 0:
                            # header section
                            header_str = mess_io.writer.global_reaction(temperatures, pressures)
                            print(header_str)

                            tot_mass = 0.
                            for rct_ich in rct_ichs:
                                geo = automol.convert.inchi.geometry(rct_ich)
                                masses = automol.geom.masses(geo)
                                for mass in masses:
                                    tot_mass += mass

                            # energy transfer section
                            energy_trans_str = mess_io.writer.energy_transfer(
                                exp_factor, exp_power, exp_cutoff, eps1, eps2, sig1, sig2, mass1, tot_mass)

                        # create channel sections for each channel
                        # check direction of reaction
                        rxn_ichs = [rct_ichs, prd_ichs]
                        rxn_chgs = [rct_chgs, prd_chgs]
                        rxn_muls = [rct_muls, prd_muls]
                        rxn_exo = moldr.util.reaction_energy(
                            save_prefix, rxn_ichs, rxn_chgs, rxn_muls, run_opt_levels[opt_level_idx])
                        print(rxn_exo)
                        if rxn_exo > 0:
                            rct_ichs, prd_ichs = prd_ichs, rct_ichs
                            rct_chgs, prd_chgs = prd_chgs, rct_chgs
                            rct_muls, prd_muls = prd_muls, rct_muls
                            print('ts search will be performed in reverse direction')

                        # construct the filesystem
                        rxn_ichs = [rct_ichs, prd_ichs]
                        rxn_chgs = [rct_chgs, prd_chgs]
                        rxn_muls = [rct_muls, prd_muls]

                        # set up the filesystem
                        is_rev = autofile.system.reaction_is_reversed(
                            rxn_ichs, rxn_chgs, rxn_muls)
                        rxn_ichs, rxn_chgs, rxn_muls = autofile.system.sort_together(
                            rxn_ichs, rxn_chgs, rxn_muls)
                        print(" - The reaction direction is {}"
                              .format('backward' if is_rev else 'forward'))

                        ts_mul = automol.mult.ts.low(rct_muls, prd_muls)

                        rxn_ichs = tuple(map(tuple, rxn_ichs))
                        rxn_chgs = tuple(map(tuple, rxn_chgs))
                        rxn_muls = tuple(map(tuple, rxn_muls))
     

                        print('rxn_save test', rxn_ichs, rxn_chgs, rxn_muls, ts_mul)
                        print(save_prefix)
                        rxn_save_fs = autofile.fs.reaction(save_prefix)
            #            rxn_save_fs.leaf.create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
                        rxn_save_path = rxn_save_fs.leaf.path(
                            [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

                        orb_restr = moldr.util.orbital_restriction(
                            ts_info, run_opt_levels[opt_level_idx])
                        ref_level = run_opt_levels[opt_level_idx][1:3]
                        ref_level.append(orb_restr)
                        print('ref level test:', ref_level)
                        thy_save_fs = autofile.fs.theory(rxn_save_path)
                        thy_save_fs.leaf.create(ref_level)
                        thy_save_path = thy_save_fs.leaf.path(ref_level)
                        spc_save_fs = autofile.fs.species(save_prefix)

                        # cycle over reactant and product species
                        # check if unimolecular or bimolecular species
                        print('Unimolecular or bimolecular test')
                        # initialize some strings if this is the first channel
                        if chn_idx == 0:
                            indxw = 0
                            indxp = 0
                            well_str = ''
                            bim_str = ''
                            ts_str = ''
                            well_ich = []
                            bim1_ich = []
                            bim2_ich = []
                        # first reactants
                        print(rct_ichs)
                        if len(rct_ichs) == 1:
                            # well
                            new_well = True
                            for idx in range(indxw):
                                if well_ich[idx] == rct_ichs[0]:
                                    new_well = False
                                    reac_label = 'W'+str(idx)
                            if new_well:
                                indxw += 1
                                well_ich.append = rct_ichs[0]
                                spc_info = (rct_ichs[0], rct_chgs[0], rct_muls[0])
                                spc_save_fs.leaf.create(spc_info)
                                spc_save_path = spc_save_fs.leaf.path(spc_info)
                                well_data = moldr.pfdriver.species_block(
                                    spc_info=spc_info,
                                    tors_model=tors_model,
                                    vib_model=vib_model,
                                    har_level=har_level,
                                    tors_level=tors_level,
                                    vpt2_level=vpt2_level,
                                    script_str=projrot_script_str,
                                    elec_levels=[[0., 1]], sym_factor=1.,
                                    save_prefix=spc_save_path,
                                    )
                                well_label = 'W'+str(indxw)
                                reac_label = well_label
                                zero_energy = 0.0
                                well_str += mess_io.writer.well(well_label, well_data, zero_energy)
                                print('well_str:', well_str)
                        else:
                            # bimolecular
                            new_bim = True
                            for idx in range(indxp):
                                if well_ich[idx] == rct_ichs[0]:
                                    new_bim = False
                                    reac_label = 'P'+str(idx)
                            if new_bim:
                                spc_data = ['', '']
                                spc_label = ['', '']
                                indxp += 1
                                bim_label = 'P'+str(indxp)
                                reac_label = bim_label
                                sym_factor = 1.
                                for idx, (rct_ich, rct_chg, rct_mul) in enumerate(zip(rct_ichs, rct_chgs, rct_muls)):
                                    spc_label[idx] = automol.inchi.smiles(rct_ich)
                                    spc_info = (rct_ich, rct_chg, rct_mul)
                                    spc_save_fs.leaf.create(spc_info)
                                    spc_save_path = spc_save_fs.leaf.path(spc_info)
                                    spc_data[idx] = moldr.pf.species_block(
                                        spc_info=spc_info,
                                        tors_model=tors_model,
                                        vib_model=vib_model,
                                        har_level=har_level,
                                        tors_level=tors_level,
                                        vpt2_level=vpt2_level,
                                        script_str=projrot_script_str,
                                        elec_levels=[[0., 1]], sym_factor=1.,
                                        save_prefix=spc_save_path,
                                        )
                                ground_energy = 0.0
                                bim_str += mess_io.writer.bimolecular(
                                    bim_label, spc_label[0], spc_data[0],
                                    spc_label[1], spc_data[1], ground_energy)
                                print('bim str:', bim_str)

                        # now products
                        print(prd_ichs)
                        if len(prd_ichs) == 1:
                            # well
                            new_well = True
                            for idx in range(indxw):
                                if well_ich[idx] == rct_ichs[0]:
                                    new_well = False
                                    reac_label = 'W'+str(idx)
                            if new_well:
                                indxw += 1
                                well_ich.append = rct_ichs[0]
                                for prd_ich, prd_chg, prd_mul in zip(prd_ichs, prd_chgs, prd_muls):
                                    spc_info = (prd_ich, prd_chg, prd_mul)
                                    spc_save_fs.leaf.create(spc_info)
                                    spc_save_path = spc_save_fs.leaf.path(spc_info)
                                    well_data = moldr.pf.species_block(
                                        spc_info=spc_info,
                                        tors_model=tors_model,
                                        vib_model=vib_model,
                                        har_level=har_level,
                                        tors_level=tors_level,
                                        vpt2_level=vpt2_level,
                                        script_str=projrot_script_str,
                                        elec_levels=[[0., 1]], sym_factor=1.,
                                        save_prefix=spc_save_path,
                                        )
                                well_label = 'W'+str(indxw)
                                prod_label = well_label
                                zero_energy = 0.0
                                well_str += mess_io.writer.well(well_label, well_data, zero_energy)
                                print('well_str:', well_str)
                # write W_indxw 
                        else:
                            new_bim = True
                            for idx in range(indxp):
                                if well_ich[idx] == rct_ichs[0]:
                                    new_bim = False
                                    prod_label = 'P'+str(idx)
                            if new_bim:
                            # bimolecular
                                indxp += 1
                                spc_data = ['', '']
                                spc_label = ['', '']
                                bimol_label = 'P'+str(indxp)
                                prod_label = bimol_label
                                for idx, (prd_ich, prd_chg, prd_mul) in enumerate(zip(prd_ichs, prd_chgs, prd_muls)):
                                    spc_info = (prd_ich, prd_chg, prd_mul)
                                    spc_save_fs.leaf.create(spc_info)
                                    spc_save_path = spc_save_fs.leaf.path(spc_info)
                                    spc_data[idx] = moldr.pf.species_block(
                                        spc_info=spc_info,
                                        tors_model=tors_model,
                                        vib_model=vib_model,
                                        har_level=har_level,
                                        tors_level=tors_level,
                                        vpt2_level=vpt2_level,
                                        script_str=projrot_script_str,
                                        elec_levels=[[0., 1]], sym_factor=1.,
                                        save_prefix=spc_save_path,
                                        )
                                    spc_label[idx] = automol.inchi.smiles(prd_ich)
                                ground_energy = 0.0
                                bim_str += mess_io.writer.bimolecular(
                                    bimol_label, spc_label[0], spc_data[0],
                                    spc_label[1], spc_data[1], ground_energy)
                                print('bim str:', bim_str)

                        elec_levels = [[0., ts_mul]]
                        symfactor = 1.
            #            prd_info=[prd_ichs[i], prd_chgs[i], prd_muls[i]]
                        # is it OK to use product info??
            #            spc_save_path = moldr.util.species_path(
            #                prd_ichs[i], prd_chgs[i], prd_muls[i], save_prefix)
            #            thy_save_path = moldr.util.theory_path(method, basis, ts_orb_restr, spc_save_path)
                        print('thy_save_path:', thy_save_path)
                        ts_data_str = moldr.pf.species_block(
                            spc_info=ts_info,
                            tors_model=tors_model,
                            vib_model=vib_model,
                            har_level=har_level,
                            tors_level=tors_level,
                            vpt2_level=vpt2_level,
                            script_str=projrot_script_str,
                            elec_levels=[[0., 1]], sym_factor=1.,
                            save_prefix=rxn_save_path,
                            )
                        ts_label = 'B1'
                        zero_energy = 50.

                        ts_str += mess_io.writer.ts_sadpt(ts_label, reac_label, prod_label, ts_data_str, zero_energy)
                        print(ts_str)

                        if chn_idx == nch_pes[pes_idx] - 1:
                            mess_inp_str = '\n'.join([header_str, energy_trans_str, well_str, bim_str, ts_str])
                            print('mess input file')
                            print(mess_inp_str)

                            bld_locs = ['MESS', 0]
                            bld_save_fs = autofile.fs.build(thy_save_path)
                            bld_save_fs.leaf.create(bld_locs)
                            mess_path = bld_save_fs.leaf.path(bld_locs)
                            print('Build Path for MESS rate files:')
                            print(mess_path)
                            with open(os.path.join(mess_path, 'mess.inp'), 'w') as mess_file:
                                mess_file.write(mess_inp_str)
                            moldr.util.run_script(rate_script_str, mess_path)
                

>>>>>>> clean ktp driver
