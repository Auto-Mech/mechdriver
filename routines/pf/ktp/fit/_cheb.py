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

    # Fit the high-pressure rate constants, if needed
    if 'high' in ktp_dct:
        [temps, rate_constants] = ktp_dct['high']
        high_params = ratefit.fit.arrhenius.single(
            temps, rate_constants, t_ref, arrfit_method,
            dsarrfit_path=mess_path, a_conv_factor=a_conv_factor)
    else:
        high_params = [1.0, 0.0, 0.0]
    
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
    for pressure in pressures:
        rate_kts = ktp_dct[pressure][1]
        fit_kts = numpy.array(fit_ktps[pressure])
        mean_avg_err, max_avg_err = ratefit.fit.fitting_errors(
            rate_kts, fit_kts)

        # Add to the dct
        err_dct[pressure] = [mean_avg_err, max_avg_err]

    # Write the Chemkin strings
    chemkin_str = chemkin_io.writer.reaction.chebyshev(
        reaction, high_params, alpha, tmin, tmax, pmin, pmax)
    chemkin_str += '\n'
    chemkin_str += chemkin_io.writer.reaction.fit_info(
        pressures, temp_dct, err_dct)

    return chemkin_str


def _cheb_err():
    """ chebyshev error
    """

    fit_ktps = ratefit.calc.chebyshev(
        alpha, TMIN, TMAX, PMIN, PMAX, TEMPS, PRESSURES)



