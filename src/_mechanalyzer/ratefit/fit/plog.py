""" Fits rate constants to a PLOG form
"""

from autoreact.params import RxnParams
from ratefit.fit import arr


def get_params(ktp_dct, dbltol=15, dbl_iter=1, tref=1.0):
    """ Gets the fitting parameters for a PLOG fit to rate constant data.
        Also gets the errors of that fit. Performs either a single or double
        Arrhenius fit at each pressure

        :param ktp_dct: rate constants to be fitted; should be P-dependent
        :type ktp_dct: dict {pressure: (temps, kts)}
        :param dbltol: the tolerance at which to reject a single fit
        :type dbltol: float
        :param dbl_iter: max number of iterations for double fitting
        :type dbl_iter: int
        :param tref: reference temp for the modified Arrhenius form
        :type tref: float
        :return params: fitted Arrhenius parameters
        :rtype: autoreact.RxnParams object
        :return err_dct: fitting errors
        :rtype: dict {pressure: (temps, errs)}
    """

    # Get the pressures
    pressures = [pressure for pressure in ktp_dct
                 if pressure != 'high']

    # Run the Arrhenius fitter for each pressure
    plog_dct = {}
    err_dct = {}
    for pressure in pressures:
        # Create a ktp_dct with only one pressure for use with Arrhenius fitter
        temp_ktp_dct = {pressure: ktp_dct[pressure]}
        temp_params, temp_err_dct = arr.get_params(
            temp_ktp_dct, dbltol=dbltol, dbl_iter=dbl_iter, tref=tref)
        arr_params = temp_params.arr  # get the Arrhenius parameters
        plog_dct[pressure] = arr_params  # update the plog_dct
        err_dct[pressure] = temp_err_dct[pressure]  # update the err_dct

    params = RxnParams(plog_dct=plog_dct)  # instantiate RxnParams

    return params, err_dct
