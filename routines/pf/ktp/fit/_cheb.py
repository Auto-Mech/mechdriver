"""
  Fit the rate constants read from the MESS output to
  Arrhenius expressions
"""


import numpy
import ratefit
import chemkin_io


def perform_fits(ktp_dct, inp_temps, reaction, mess_path,
                 a_conv_factor=1.0, t_ref=1.0, arrfit_method='python',
                 tdeg=6, pdeg=4):
    """ Read the rates for each channel and perform the fits
    """

    # Obtain the fit paramts for the high-pressure rate constants, if needed
    if 'high' in ktp_dct:
        [temps, rate_constants] = ktp_dct['high']
        highp_params = ratefit.fit.arrhenius.single(
            temps, rate_constants, t_ref, arrfit_method,
            dsarrfit_path=mess_path, a_conv_factor=a_conv_factor)
    else:
        highp_params = [1.0, 0.0, 0.0]

    # If there are rates for more than one pressure
    pressures = tuple(pressure for pressure in ktp_dct.keys()
                      if pressure != 'high')

    # Fit rate constants to Chebyshev polynomial
    alpha, trange, prange = ratefit.fit.chebyshev.kfit(
        inp_temps, ktp_dct, tdeg=tdeg, pdeg=pdeg)
    tmin, tmax = trange
    pmin, pmax = prange

    # Calculate the fitted rate constants
    fit_ktps = ratefit.calc.chebyshev(
        alpha, tmin, tmax, pmin, pmax, inp_temps, pressures)

    # Calculate errors
    err_dct, temp_dct = {}, {}
    num_kts = []
    # print('reaction in rate fit test:', reaction)
    for pressure in pressures:
        rate_kts = ktp_dct[pressure][1]
        fit_kts = numpy.array(fit_ktps[pressure])
        mean_avg_err, max_avg_err = ratefit.fit.fitting_errors(
            rate_kts, fit_kts)
        # print('rate fit test:', pressure, rate_kts, fit_kts)

        # Add to the dct and lst
        err_dct[pressure] = [mean_avg_err, max_avg_err]
        num_kts.append(len(rate_kts))

    # Assess if fit is viable by looking at err and same num k(T) at each P
    fit_viable = True
    if len(set(num_kts)) == 1:
        print('Different number of k(T) values at different pressures...')
        fit_viable = False
    elif max((vals[1] for vals in err_dct.values())) < 20.0:
        print('Errors from Chebyshve fit too large (see string)...')
        fit_viable = False

    # Write the Chemkin strings
    chemkin_str = chemkin_io.writer.reaction.chebyshev(
        reaction, highp_params, alpha, tmin, tmax, pmin, pmax)
    chemkin_str += '\n'
    chemkin_str += chemkin_io.writer.reaction.fit_info(
        pressures, temp_dct, err_dct)

    if not fit_viable:
        # Print message and reset string to empty to trigger Arrhenius
        print('Chemkin string from Chebyshev fit')
        print(chemkin_str)
        chemkin_str = ''

    return chemkin_str
