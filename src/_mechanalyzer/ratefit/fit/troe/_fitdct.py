"""
  Fit the rate constants read from the MESS output to
  Troe expressions
"""

import chemkin_io
import ratefit
# from ratefit.fit import fitting_temperatures_dct
from ratefit.fit import fitting_error_dct
from ratefit.fit import set_a_conversion_factor


def pes(ktp_dct, reaction, mess_path, params, tol):
    """ Fit rate constants to Troe parameters
    """

    # Initialize the a conversion factor
    a_conv_factor = set_a_conversion_factor(reaction)

    # Dictionaries to store info; indexed by pressure (given in fit_ps)
    fit_param_dct = {}
    # fit_temp_dct = {}

    # Calculate the fitting parameters from the filtered T,k lists
    fit_params = ratefit.fit.troe.reaction(
        ktp_dct, mess_path,
        troe_param_fit_lst=params,
        highp_guess=(8.1e-11, -0.01, 1000.0),
        lowp_guess=(8.1e-11, -0.01, 1000.0),
        alpha=0.19, ts1=590, ts2=1.e6, ts3=6.e4,
        fit_tol1=1.0e-8, fit_tol2=1.0e-8,
        a_conv_factor=a_conv_factor)

    # Store the parameters in the fit dct
    for pressure in ktp_dct:
        fit_param_dct[pressure] = fit_params

    # Build fit temp dct
    # fit_temp_dct = fitting_temperatures_dct(fit_param_dct)

    # Check if the desired fits were successful at each pressure
    fit_success = all(params for params in fit_param_dct.values())

    # Calculate the errors from the Troe fits
    if fit_success:
        fit_err_dct = assess_troe_fit_err(
            fit_param_dct, ktp_dct, t_ref=1.0, a_conv_factor=a_conv_factor)
        fit_good = max((vals[1] for vals in fit_err_dct.values())) < tol
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
        fit_k_dct[pressure] = (temps, fit_ks / a_conv_factor)

    # Calculute the error between the calc and fit ks
    fit_err_dct = fitting_error_dct(ktp_dct, fit_k_dct)

    return fit_err_dct
