"""
  Fit the rate constants read from the MESS output to
  Arrhenius expressions
"""

import ratefit
import chemkin_io
from mechroutines.pf.ktp.fit._util import pull_highp_from_dct
from mechlib.amech_io import printer as ioprinter


# Functions to fit rates to Arrhenius/PLOG function
def perform_fits(ktp_dct, reaction, mess_path,
                 a_conv_factor, arrfit_thresh):
    """ Read the rates for each channel and perform the fits
    """

    # Fit rate constants to single Arrhenius expressions
    sing_params_dct, sing_fit_temp_dct, sing_fit_success = mod_arr_fit(
        ktp_dct, mess_path, fit_type='single', fit_method='python',
        t_ref=1.0, a_conv_factor=a_conv_factor)
    if sing_fit_success:
        ioprinter.info_message(
            'Successful fit to Single Arrhenius at all T, P', newline=1)

    # Assess the errors of the single Arrhenius Fit
    sing_fit_err_dct = assess_arr_fit_err(
        sing_params_dct, ktp_dct,
        fit_type='single',
        t_ref=1.0, a_conv_factor=a_conv_factor)

    # Write a chemkin string for the single fit
    sing_highp, sing_plog_dct, pressures = pull_highp_from_dct(sing_params_dct)
    if sing_plog_dct:  # if PLOG
        sing_chemkin_str = chemkin_io.writer.reaction.plog(
            reaction, sing_plog_dct)
    else:  # if Arrhenius
        sing_chemkin_str = chemkin_io.writer.reaction.arrhenius(
            reaction, sing_highp)
    sing_chemkin_str += '\n'
    sing_chemkin_str += chemkin_io.writer.reaction.fit_info(
        pressures, sing_fit_temp_dct, sing_fit_err_dct)

    # Assess single fitting errors:
    # are they within threshold at each pressure
    thresh, choice = arrfit_thresh
    if choice == 'max':
        test_val = max((
            vals[1] for vals in sing_fit_err_dct.values()))
    elif choice == 'mean':
        test_val = max((
            vals[0] for vals in sing_fit_err_dct.values()))
    sgl_fit_good = bool(test_val < thresh)

    # Put a double fit assess function
    # Skip first and last points in estimating maximum errors –
    # or better yet, skip first or last data point, and try fit again.

    # Assess if a double Arrhenius fit is possible
    dbl_fit_poss = all(
        len(ktp_dct[p][0]) >= 6 for p in ktp_dct)

    # Write chemkin string for single/double fit, based on errors
    chemkin_str = ''
    if sgl_fit_good:
        ioprinter.info_message(
            'Single fit errors acceptable: Using single fits')
        chemkin_str += sing_chemkin_str
    elif not sgl_fit_good and not dbl_fit_poss:
        ioprinter.info_message(
            'Not enough temperatures for a double fit:',
            ' Using single fits')
        chemkin_str += sing_chemkin_str
    elif not sgl_fit_good and dbl_fit_poss:
        ioprinter.info_message(
            'Single fit errs too large & double fit possible:',
            ' Trying double fit')

        # Generate guess parameters
        ioprinter.info_message(
            'Generating Double Fit Guess from Single Fit Parameters')
        guess_params_dct = make_dbl_fit_guess(sing_params_dct)

        # Fit rate constants to double Arrhenius expressions
        doub_params_dct, doub_fit_temp_dct, doub_fit_suc = mod_arr_fit(
            ktp_dct, mess_path, fit_type='double',
            fit_method='dsarrfit', t_ref=1.0,
            a_conv_factor=a_conv_factor,
            inp_param_dct=guess_params_dct)

        if doub_fit_suc:
            ioprinter.info_message(
                'Successful fit to double Arrhenius at all T, P', newline=1)

            ioprinter.info_message(
                'Writing fitting parameters and errors from ',
                'single arrhenius fit for comparison', newline=1)
            ioprinter.info_message(sing_chemkin_str)

            # Assess the errors of the single Arrhenius Fit
            doub_fit_err_dct = assess_arr_fit_err(
                doub_params_dct, ktp_dct, fit_type='double',
                t_ref=1.0, a_conv_factor=a_conv_factor)
            doub_highp, doub_plog_dct, pressures = pull_highp_from_dct(
                doub_params_dct)
            if doub_plog_dct:  # if PLOG
                chemkin_str += chemkin_io.writer.reaction.plog(
                    reaction, doub_plog_dct)
            else:  # if Arrhenius
                chemkin_str = chemkin_io.writer.reaction.arrhenius(
                    reaction, doub_highp)

            chemkin_str += '\n'
            chemkin_str += chemkin_io.writer.reaction.fit_info(
                pressures, doub_fit_temp_dct, doub_fit_err_dct)
        else:
            ioprinter.warning_message(
                'Double Arrhenius fit failed for some reason:',
                ' Using single fits')
            chemkin_str += sing_chemkin_str

    return chemkin_str


def mod_arr_fit(ktp_dct, mess_path,
                fit_type='single', fit_method='dsarrfit',
                t_ref=1.0, a_conv_factor=1.0, inp_param_dct=None):
    """
    Routine for a single reaction:
        (1) Grab high-pressure and pressure-dependent rate constants
            from a MESS output file
        (2) Fit rate constants to an Arrhenius expression
    """

    assert fit_type in ('single', 'double'), 'Only single/double fits'
    if inp_param_dct is not None:
        assert set(list(ktp_dct.keys())) <= set(list(inp_param_dct.keys())), (
            'Pressure in ktp not in guess: \n {} \n{}'.format(
                list(ktp_dct.keys()), list(inp_param_dct.keys()))
            )

    # Dictionaries to store info; indexed by pressure (given in fit_ps)
    fit_param_dct = {}
    fit_temp_dct = {}

    # Calculate the fitting parameters from the filtered T,k lists
    for pressure, tk_arr in ktp_dct.items():

        # Set the temperatures and rate constants
        temps = tk_arr[0]
        rate_constants = tk_arr[1]

        # Fit rate constants using desired Arrhenius fit
        if fit_type == 'single':

            fit_params = ratefit.fit.arrhenius.single(
                temps, rate_constants, t_ref, fit_method,
                dsarrfit_path=mess_path, a_conv_factor=a_conv_factor)

        elif fit_type == 'double':

            # Generate guess parameters
            if inp_param_dct is not None:
                arr_guess = inp_param_dct[pressure]
                arr1_guess, arr2_guess = arr_guess[:3], arr_guess[3:]
            else:
                arr1_guess, arr2_guess = (8.1e-11, -0.01, 2000.0), ()

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


def make_dbl_fit_guess(params_dct):
    """ Make dbl fit term
    """

    dbl_guess_dct = {}
    for pressure, params in params_dct.items():
        dbl_guess_dct[pressure] = _generate_guess(params)

    return dbl_guess_dct


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

    return arr1_guess + arr2_guess


# Fit Assessments
def _check_double_fit(sing_fit_dct, dbl_fit_dct,
                      temps=(2000.0), t_ref=1.0):
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

        for sk_val, dk_val in zip(sgl_ks, dbl_ks):
            if dk_val >= sk_val * 10.0:
                bad_dbl = True

    return bad_dbl


def assess_arr_fit_err(fit_param_dct, ktp_dct, fit_type='single',
                       t_ref=1.0, a_conv_factor=1.0,
                       err_set='all'):
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

        # Assess the errors using some subset of the rate constants
        test_calc_ks, test_fit_ks = _gen_err_set(
            calc_ks, fit_ks, err_set=err_set)
        mean_avg_err, max_avg_err = ratefit.fit.fitting_errors(
            test_calc_ks, test_fit_ks)

        # Store in a dictionary
        fit_err_dct[pressure] = [mean_avg_err, max_avg_err]

    return fit_err_dct


def _gen_err_set(calc_ks, fit_ks, err_set='all'):
    """ look at err ranges
    """

    if err_set == 'skip':
        test_calc_ks = calc_ks[1:-2]
        test_fit_ks = fit_ks[1:-2]
    else:
        test_calc_ks = calc_ks
        test_fit_ks = fit_ks

    return test_calc_ks, test_fit_ks
