"""
  Fit the rate constants read from the MESS output to
  Arrhenius, Plog, or Troe expressions
"""

import os
import copy
import numpy
import ratefit
import chemkin_io
import mess_io
from lib.amech_io import writer
from lib.phydat import phycon


def fit_rates(inp_temps, inp_pressures, inp_tunit, inp_punit,
              pes_formula_str, idx_dct,
              es_info, pf_model,
              mess_path, fit_method, pdep_fit,
              arrfit_thresh):
    """ Parse the MESS output and fit the rates to
        Arrhenius expressions written as CHEMKIN strings
    """

    # Write header string containing thy information
    chemkin_header_str = writer.ckin.run_ckin_header(es_info, pf_model)
    chemkin_header_str += '\n'

    # Initialize full chemkin string and paths
    # chemkin_full_str = chemkin_header_str
    chemkin_full_str = ''
    starting_path = os.getcwd()
    ckin_path = ''.join([starting_path, '/ckin'])
    if not os.path.exists(ckin_path):
        os.mkdir(ckin_path)
    full_ckin_path = os.path.join(ckin_path, pes_formula_str+'.ckin')
    # if os.path.exists(full_ckin_path):
    #     os.remove(full_ckin_path)

    # Loop through reactions, fit rates, and write ckin strings
    labels = idx_dct.values()
    names = idx_dct.keys()
    for lab_i, name_i in zip(labels, names):
        a_conv_factor = phycon.NAVO if 'W' not in lab_i else 1.00
        if 'F' not in lab_i and 'B' not in lab_i:
            for lab_j, name_j in zip(labels, names):
                if 'F' not in lab_j and 'B' not in lab_j and lab_i != lab_j:

                    print(('\n--------------------------------------------' +
                           '------------------------------------------'))
                    # Set name
                    reaction = name_i + '=' + name_j
                    print('\nGetting Rates for {}'.format(
                        reaction))

                    # Initialize new chemkin str for reaction
                    chemkin_str = chemkin_header_str

                    # Read the rate constants out of the mess outputs
                    print('\nReading k(T,P)s from MESS output...')
                    ktp_dct = read_rates(
                        inp_temps, inp_pressures, inp_tunit, inp_punit,
                        lab_i, lab_j, mess_path, pdep_fit,
                        bimol=numpy.isclose(a_conv_factor, 6.0221e23))

                    # Move ahead in loop if no rates found
                    if not ktp_dct:
                        print('\nNo valid k(T,Ps)s from MESS output.')
                        print('\nSkipping to next reaction...')
                        continue

                    # Get the desired fits in the form of CHEMKIN strs
                    if fit_method == 'arrhenius':
                        print('\nFitting k(T,P)s to PLOG/Arrhenius Form....')
                        chemkin_str += perform_arrhenius_fits(
                            ktp_dct, reaction, mess_path,
                            a_conv_factor, arrfit_thresh)
                    elif fit_method == 'troe':
                        print('\nFitting k(T,P)s to Troe Form...')
                        # pass
                        # chemkin_str += perform_troe_fits(
                        #     ktp_dct, reaction, mess_path,
                        #     troe_param_fit_lst,
                        #     a_conv_factor, err_thresh)

                    # Write the CHEMKIN strings
                    chemkin_full_str += chemkin_str
                    chemkin_full_str += '\n\n'
                    # chemkin_str = chemkin_header_str + chemkin_str
                    print('\nFitting Parameters in CHEMKIN Format:')
                    print(chemkin_str)

                    # Print the results for each channel to a file
                    pes_chn_lab = pes_formula_str + '_' + name_i + '_' + name_j
                    ckin_name = os.path.join(ckin_path, pes_chn_lab+'.ckin')
                    with open(ckin_name, 'w') as cfile:
                        cfile.write(chemkin_str)

    # Print the results for the whole PES to a file
    with open(full_ckin_path, 'a') as cfile:
        cfile.write(chemkin_full_str)


# Functions to fit rates to Arrhenius/PLOG function
def perform_arrhenius_fits(ktp_dct, reaction, mess_path,
                           a_conv_factor, err_thresh):
    """ Read the rates for each channel and perform the fits
    """
    # Fit rate constants to single Arrhenius expressions
    sing_params_dct, sing_fit_temp_dct, sing_fit_success = mod_arr_fit(
        ktp_dct, mess_path, fit_type='single', fit_method='python',
        t_ref=1.0, a_conv_factor=a_conv_factor)
    if sing_fit_success:
        print('\nSuccessful fit to Single Arrhenius at all T, P')

    # Assess the errors of the single Arrhenius Fit
    sing_fit_err_dct = assess_arr_fit_err(
        sing_params_dct, ktp_dct,
        fit_type='single',
        t_ref=1.0, a_conv_factor=a_conv_factor)

    # Write a chemkin string for the single fit
    sing_chemkin_str = chemkin_io.writer.reaction.plog(
        reaction, sing_params_dct,
        err_dct=sing_fit_err_dct, temp_dct=sing_fit_temp_dct)

    # Assess single fitting errors:
    # are they within threshold at each pressure
    thresh, choice = arrfir_thresh
    if choice == 'max':
        test_val = max((
            vals[1] for vals in sing_fit_err_dct.values()))
    elif choice == 'mean':
        test_val = max((
            vals[0] for vals in sing_fit_err_dct.values()))
    sgl_fit_good = bool(test_val < thresh)

    # Put a double fit assess function
    # Skip first and last points in estimating maximum errors – or better yet, skip first or last data point, and try fit again.
    
    # Assess if a double Arrhenius fit is possible
    dbl_fit_poss = all(
        len(ktp_dct[p][0]) >= 6 for p in ktp_dct)

    # Write chemkin string for single/double fit, based on errors
    chemkin_str = ''
    if sgl_fit_good:
        print('\nSingle fit errors acceptable: Using single fits')
        chemkin_str += sing_chemkin_str
    elif not sgl_fit_good and not dbl_fit_poss:
        print('\nNot enough temperatures for a double fit:',
              ' Using single fits')
        chemkin_str += sing_chemkin_str
    elif not sgl_fit_good and dbl_fit_poss:
        print('\nSingle fit errs too large & double fit possible:',
              ' Trying double fit')

        # Fit rate constants to double Arrhenius expressions
        doub_params_dct, doub_fit_temp_dct, doub_fit_suc = mod_arr_fit(
            ktp_dct, mess_path, fit_type='double',
            fit_method='dsarrfit', t_ref=1.0,
            a_conv_factor=a_conv_factor,
            inp_params_dct=sing_params_dct)

        if doub_fit_suc:
            print('\nSuccessful fit to Double Arrhenius at all T, P')

            print('\nWriting fitting parameters and errors from',
                  'single arrhenius fit for comparison')
            print(sing_chemkin_str)

            # Assess the errors of the single Arrhenius Fit
            doub_fit_err_dct = assess_arr_fit_err(
                doub_params_dct, ktp_dct, fit_type='double',
                t_ref=1.0, a_conv_factor=a_conv_factor)
            chemkin_str += chemkin_io.writer.reaction.plog(
                reaction, doub_params_dct,
                err_dct=doub_fit_err_dct, temp_dct=doub_fit_temp_dct)
        else:
            print('\nDouble Arrhenius Fit failed for some reason:',
                  ' Using Single fits')
            chemkin_str += sing_chemkin_str

    return chemkin_str


def mod_arr_fit(ktp_dct, mess_path, fit_type='single', fit_method='dsarrfit',
                t_ref=1.0, a_conv_factor=1.0, inp_params_dct=None):
    """
    Routine for a single reaction:
        (1) Grab high-pressure and pressure-dependent rate constants
            from a MESS output file
        (2) Fit rate constants to an Arrhenius expression
    """

    assert fit_type in ('single', 'double'), 'Only single/double fits'
    if inp_param_dct is not None:
        assert set(list(ktp_dct.keys())) == set(list(inp_param_dct.keys())), (
                'Pressures for ktp and guess dcts must be the same')

    # Dictionaries to store info; indexed by pressure (given in fit_ps)
    fit_param_dct = {}
    fit_temp_dct = {}

    # Calculate the fitting parameters from the filtered T,k lists
    for pressure, tk_arr in ktp_dct.items():

        # Set the temperatures and rate constants
        temps = tk_arr[0]
        rate_constants = tk_arr[1]

        # Build guess params for a double fit
        if inp_param_dct is not None and fit_type == 'double':
            arr1_guess, arr2_guess = _generate_guess(
                inp_param_dct[pressure]) 

        # Fit rate constants using desired Arrhenius fit
        if fit_type == 'single':
            fit_params = ratefit.fit.arrhenius.single(
                temps, rate_constants, t_ref, fit_method,
                dsarrfit_path=mess_path, a_conv_factor=a_conv_factor)
        elif fit_type == 'double':
            fit_params = ratefit.fit.arrhenius.double(
                temps, rate_constants, t_ref, fit_method,
                arr1_guess=arr1_guess, arr2_guess=arr2_guess,
                dsarrfit_path=mess_path, a_conv_factor=a_conv_factor)

        # Store the fitting parameters in a dictionary
        fit_param_dct[pressure] = fit_params

        # Store the temperatures used to fit in a dictionary
        fit_temp_dct[pressure] = [min(temps), max(temps)]

    # Check if the desired fits were successful at each pressure
    fit_success = all(params for params in fit_param_dct.values())

    return fit_param_dct, fit_temp_dct, fit_success


def _generate_guess(params):
    """ Generate a set of double fit guess params based on input.
        Right now it is just a set of single params.
    """

    # Unpack input single params
    [a_inp, n_inp, ea_inp] = params

    # Generate new guesses
    a1_guess = a_inp * 0.5
    a2_guess = a_inp * 1.5
    n1_guess = n_inp * 0.9
    n2_guess = n_inp * 1.1
    ea1_guess = ea_inp
    ea2_guess = ea_inp

    arr1_guess = (a1_guess, n1_guess, ea1_guess)
    arr2_guess = (a2_guess, n2_guess, ea2_guess)

    return arr1_guess, arr2_guess


def _check_double_fit(sing_fit_dct, dbl_fit_dct,
                      temps=[2000.0], t_ref=1.0):
    """ Check if the double fit is bad
    """

    bad_dbl = False
    for sarr, darr in zip(sing_fit_dct.values(), dbl_fit_dct.values()):

        # Set the temperatures
        # temps = ktp_dct[pressure][0]
        sparams = sarr[1]
        dparams = darr[1]

        # Calculate fitted rate constants, based on fit type
        sgl_ks = ratefit.calc.single_arrhenius(
            sparams[0], sparams[1], sparams[2],
            t_ref, temps)
        dbl_ks = ratefit.calc.double_arrhenius(
            dparams[0], dparams[1], dparams[2],
            dparams[3], dparams[4], dparams[5],
            t_ref, temps)

        for sk, dk in zip(sgl_ks, dbl_ks):
            if dk >= sk * 10.0:
                bad_dbl = True

    return bad_dbl


def assess_arr_fit_err(fit_param_dct, ktp_dct, fit_type='single',
                       t_ref=1.0, a_conv_factor=1.0):
    """ Determine the errors in the rate constants that arise
        from the Arrhenius fitting procedure
    """

    fit_k_dct = {}
    fit_err_dct = {}

    assert fit_type in ('single', 'double')

    # Calculate fitted rate constants using the fitted parameters
    for pressure, params in fit_param_dct.items():

        # Set the temperatures
        temps = ktp_dct[pressure][0]

        # Calculate fitted rate constants, based on fit type
        if fit_type == 'single':
            fit_ks = ratefit.calc.single_arrhenius(
                params[0], params[1], params[2],
                t_ref, temps)
        elif fit_type == 'double':
            fit_ks = ratefit.calc.double_arrhenius(
                params[0], params[1], params[2],
                params[3], params[4], params[5],
                t_ref, temps)

        # Store the fitting parameters in a dictionary
        fit_k_dct[pressure] = fit_ks / a_conv_factor

    # Calculute the error between the calc and fit ks
    for pressure, fit_ks in fit_k_dct.items():

        calc_ks = ktp_dct[pressure][1]
        mean_avg_err, max_avg_err = ratefit.fit.fitting_errors(
            calc_ks, fit_ks)

        # Store in a dictionary
        fit_err_dct[pressure] = [mean_avg_err, max_avg_err]

    return fit_err_dct


# Functions to fit rates to Troe function
def perform_troe_fits(ktp_dct, reaction, mess_path,
                      troe_param_fit_lst, a_conv_factor, err_thresh):
    """ Fit rate constants to Troe parameters
    """

    # Dictionaries to store info; indexed by pressure (given in fit_ps)
    fit_param_dct = {}
    fit_temp_dct = {}

    # Calculate the fitting parameters from the filtered T,k lists
    new_dct = {}
    for key, val in ktp_dct.items():
        if key != 'high':
            new_dct[key] = val
    inv_ktp_dct = ratefit.fit.flip_ktp_dct(new_dct)
    fit_params = ratefit.fit.troe.std_form(
        inv_ktp_dct, troe_param_fit_lst, mess_path,
        highp_a=8.1e-11, highp_n=-0.01, highp_ea=1000.0,
        lowp_a=8.1e-11, lowp_n=-0.01, lowp_ea=1000.0,
        alpha=0.19, ts1=590, ts2=1.e6, ts3=6.e4,
        fit_tol1=1.0e-8, fit_tol2=1.0e-8,
        a_conv_factor=1.0)

    # Store the parameters in the fit dct
    for pressure, tk_arr in ktp_dct.items():

        # Store the fitting parameters in a dictionary
        fit_param_dct[pressure] = fit_params

        # Store the temperatures used to fit in a dictionary
        temps = tk_arr[0]
        fit_temp_dct[pressure] = [min(temps), max(temps)]

    # Check if the desired fits were successful at each pressure
    fit_success = all(params for params in fit_param_dct.values())

    # Calculate the errors from the Troe fits
    if fit_success:
        fit_err_dct = assess_troe_fit_err(
            fit_param_dct, ktp_dct, t_ref=1.0, a_conv_factor=a_conv_factor)
        fit_good = max((vals[1] for vals in fit_err_dct.values())) < err_thresh
        if fit_good:
            chemkin_str = chemkin_io.writer.reaction.troe(
                reaction,
                [fit_params[0], fit_params[1], fit_params[2]],
                [fit_params[3], fit_params[4], fit_params[5]],
                [fit_params[6], fit_params[8], fit_params[7], fit_params[9]],
                colliders=())

        # Store the fitting parameters in a dictionary
    else:
        chemkin_str = ''

    return chemkin_str


def assess_troe_fit_err(fit_param_dct, ktp_dct, t_ref=1.0, a_conv_factor=1.0):
    """ Determine the errors in the rate constants that arise
        from the Arrhenius fitting procedure
    """

    fit_k_dct = {}
    fit_err_dct = {}

    # Calculate fitted rate constants using the fitted parameters
    fit_pressures = fit_param_dct.keys()
    for pressure, params in fit_param_dct.items():

        # Set the temperatures
        temps = ktp_dct[pressure][0]

        # Calculate fitted rate constants, based on fit type
        highp_arrfit_ks = ratefit.calc.single_arrhenius(
            params[0], params[1], params[2],
            t_ref, temps)
        lowp_arrfit_ks = ratefit.calc.single_arrhenius(
            params[3], params[4], params[5],
            t_ref, temps)
        fit_ks = ratefit.calc.troe(
            highp_arrfit_ks, lowp_arrfit_ks, fit_pressures, temps,
            params[6], params[8], params[7], ts2=params[9])

        # Store the fitting parameters in a dictionary
        fit_k_dct[pressure] = fit_ks / a_conv_factor

    # Calculute the error between the calc and fit ks
    for pressure, fit_ks in fit_k_dct.items():

        calc_ks = ktp_dct[pressure][1]
        mean_avg_err, max_avg_err = ratefit.calc.fitting_errors(
            calc_ks, fit_ks)

        # Store in a dictionary
        fit_err_dct[pressure] = [mean_avg_err, max_avg_err]

    return fit_err_dct


# def check_single_fit_errors(sing_fit_err_dct, err_thresh):
#     """ Determine if fit errors
#     """
#     pass
#
#
# def print_fit_info(params_dct, err_dct):
#     """ print fitting infor
#     """
#     # inf_str = '{} {} {} {}'.format(
#     # for pressure, params in sing_params_dct.items():
#     #     inf_str += '{} {8.5} {} {}'.format(

# Readers
def read_rates(inp_temps, inp_pressures, inp_tunit, inp_punit,
               rct_lab, prd_lab, mess_path, pdep_fit, bimol=False):
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

    # Read the temperatures and pressures out of the MESS output
    mess_temps, tunit = mess_io.reader.rates.get_temperatures(
        output_string)
    mess_pressures, punit = mess_io.reader.rates.get_pressures(
        output_string)

    assert inp_temps <= mess_temps
    assert inp_pressures <= mess_pressures
    assert inp_tunit == tunit
    assert inp_punit == punit

    # Loop over the pressures obtained from the MESS output
    for pressure in mess_pressures:

        # Read the rate constants
        if pressure == 'high':
            rate_ks = mess_io.reader.highp_ks(
                output_string, rct_lab, prd_lab)
        else:
            rate_ks = mess_io.reader.pdep_ks(
                output_string, rct_lab, prd_lab, pressure)

        # Store in a dictionary
        calc_k_dct[pressure] = rate_ks

    # Remove k(T) vals at each P where where k is negative or undefined
    # If ANY valid k(T,P) vals at given pressure, store in dct
    print('\nRemoving invalid k(T,P)s from MESS output that are either:\n',
          '  (1) negative, (2) undefined [***], or (3) below 10**(-21) if',
          'reaction is bimolecular')
    for pressure, calc_ks in calc_k_dct.items():
        filtered_temps, filtered_ks = ratefit.fit.get_valid_tk(
            mess_temps, calc_ks, bimol)
        if filtered_ks.size > 0:
            valid_calc_tk_dct[pressure] = [filtered_temps, filtered_ks]

    # Filter the ktp dictionary by assessing the presure dependence
    if valid_calc_tk_dct:
        if list(valid_calc_tk_dct.keys()) == ['high']:
            print('\nValid k(T)s only found at High Pressure...')
            ktp_dct['high'] = valid_calc_tk_dct['high']
        else:
            if pdep_fit:
                print('\nUser requested to assess pressure dependence',
                      'of reaction.')
                assess_pdep_temps = pdep_fit['assess_pdep_temps']
                pdep_tolerance = pdep_fit['pdep_tolerance']
                no_pdep_pval = pdep_fit['no_pdep_pval']
                pdep_low = pdep_fit['pdep_low']
                pdep_high = pdep_fit['pdep_high']
                rxn_is_pdependent = ratefit.fit.assess_pressure_dependence(
                    valid_calc_tk_dct, assess_pdep_temps,
                    tolerance=pdep_tolerance, plow=pdep_low, phigh=pdep_high)
                if rxn_is_pdependent:
                    print('  Reaction found to be pressure dependent.',
                          'Fitting all k(T)s from all pressures',
                          'found in MESS.')
                    # Set dct as copy of dct to do PLOG fits at all pressures
                    ktp_dct = copy.deepcopy(valid_calc_tk_dct)
                else:
                    # Set dct w/ single set of k(T, P) at desired pressure
                    print('  No pressure dependence detected.',
                          'Grabbing k(T)s at {} {}'.format(
                              no_pdep_pval, punit))
                    if no_pdep_pval in valid_calc_tk_dct:
                        ktp_dct['high'] = valid_calc_tk_dct[no_pdep_pval]
            else:
                # Set dct fit as copy of dct to do PLOG fits at all pressures
                ktp_dct = copy.deepcopy(valid_calc_tk_dct)

        # Patchy way to get high-pressure rates in dct if needed
        if 'high' not in ktp_dct and 'high' in valid_calc_tk_dct.keys():
            ktp_dct['high'] = valid_calc_tk_dct['high']

    return ktp_dct
