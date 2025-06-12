""" Fits rate constants to a Chebyshev form
"""

import numpy
from scipy.special import eval_chebyt
from autoreact.params import RxnParams
from ratefit.fit import arr
from ratefit.fit import err

RCOND = -1  # a parameter for the Numpy least squares fitting


def get_params(ktp_dct, tdeg=4, pdeg=6, tol=20.0):
    """ Gets the fitting parameters for a Chebyshev fit to rate constant data.
        Also gets the errors of that fit.

        :param ktp_dct: rate constants to be fitted; should be P-dependent
        :type ktp_dct: dict {pressure: (temps, kts)}
        :param tdeg: number of temperature coefficients
        :type tdeg: int
        :param pdeg: number of pressure coefficients
        :type pdeg: int
        :param tol: percent error tolerance, above which to warn the user of a
            bad Chebyshev fit
        :type tol: float
        :return params: fitted Arrhenius parameters
        :rtype: autoreact.RxnParams object
        :return err_dct: fitting errors
        :rtype: dict {pressure: (temps, errs)}
    """

    # Check the viability of the ktp_dct
    is_viable = check_viability(ktp_dct)
    assert is_viable, '# of temps is not the same at all pressures'

    # Obtain fit params for the 1-atm rate constants (used for printing only)
    one_atm_arr = get_one_atm_arr(ktp_dct)  # tuple of tuples

    alpha, tlim, plim = get_alpha(ktp_dct, tdeg=tdeg, pdeg=pdeg)

    # Store items in a cheb_dct and instantiate RxnParams
    cheb_dct = {
        'alpha': alpha,
        'tlim': tlim,
        'plim': plim,
        'one_atm_arr': one_atm_arr}
    params = RxnParams(cheb_dct=cheb_dct)
    err_dct = err.get_err_dct(ktp_dct, params)

    # Print a warning if the max_err exceeds the input tol
    max_err = err.get_max_err(err_dct)
    if max_err > tol:
        print(f'The maximum Chebyshev fitting error, {max_err:.1f}, is greater'
              f' than the desired max error, {tol:.1f}')

    return params, err_dct


def get_alpha(ktp_dct, tdeg=4, pdeg=6):
    """ Performs the Chebyshev fit: gets the alpha matrix

        :param ktp_dct: rate constants to be fitted; should be P-dependent
        :type ktp_dct: dict {pressure: (temps, kts)}
        :param tdeg: number of temperature coefficients
        :type tdeg: int
        :param pdeg: number of pressure coefficients
        :type pdeg: int
        :return alpha: array of Chebyshev polynomial coefficients
        :rtype: numpy.ndarray of shape (tdeg, pdeg)
        :return tlim: minimum and maximum temperatures of fit
        :rtype: tuple (tmin, tmax)
        :return plim: minimum and maximum pressures of fit
        :rtype: tuple (pmin, pmax)
    """

    def conv_dct_to_array(ktp_dct):
        """ Converts the contents of a ktp_dct to a Numpy array

            :return tp_array: array of kts
            :rtype: Numpy array of shape (num_temps, num_pressures)
        """

        # Check the first pressure to get the number of temps (same for all)
        num_temps = len(ktp_dct[list(ktp_dct.keys())[0]][0])
        num_pressures = len(ktp_dct)
        pt_array = numpy.ndarray([num_pressures, num_temps])
        pressures = tuple(pressure for pressure in ktp_dct.keys()
                          if pressure != 'high')
        for pidx, pressure in enumerate(pressures):
            _, kts = ktp_dct[pressure]
            pt_array[pidx, :] = kts

        tp_array = numpy.transpose(pt_array)

        return tp_array

    pressures = tuple(pressure for pressure in ktp_dct.keys()
                      if pressure != 'high')
    temps = ktp_dct[pressures[0]][0]  # all temp vectors should be the same

    # Get info on temps and pressures
    tnum, pnum = len(temps), len(pressures)
    tmin, tmax = min(temps), max(temps)
    pmin, pmax = min(pressures), max(pressures)

    # Get reduced temperatures and pressures
    tred = (2*temps**(-1) - tmin**(-1) - tmax**(-1))\
        / (tmax**(-1) - tmin**(-1))
    pred = (2*numpy.log10(pressures) - numpy.log10(pmin) - numpy.log10(pmax))\
        / (numpy.log10(pmax) - numpy.log10(pmin))

    # Convert the ktp_dct to a Numpy array of shape (num_temps, num_pressures)
    tp_array = conv_dct_to_array(ktp_dct)

    # Create matrices for fits
    amat = numpy.zeros((tnum * pnum, tdeg * pdeg), numpy.float64)
    bvec = numpy.zeros((tnum * pnum), numpy.float64)
    nzero = 0
    for tidx1, temp in enumerate(tred):
        for pidx1, press in enumerate(pred):
            for tidx2 in range(tdeg):
                for pidx2 in range(pdeg):
                    idx1 = (pidx1 * tnum + tidx1)
                    idx2 = (pidx2 * tdeg + tidx2)
                    amat[idx1, idx2] = (
                        eval_chebyt(tidx2, temp) * eval_chebyt(pidx2, press))
            if tp_array[tidx1, pidx1] is not None:
                bvec[idx1] = numpy.log10(tp_array[tidx1, pidx1])
            else:
                bvec[idx1] = None
                nzero += 1

    nnonzero = tnum * pnum - nzero
    idxp = -1
    amatp = numpy.zeros((nnonzero, tdeg * pdeg), numpy.float64)
    bvecp = numpy.zeros((nnonzero), numpy.float64)
    for idx in range(tnum*pnum):
        if not numpy.isnan(bvec[idx]):
            idxp += 1
            bvecp[idxp] = bvec[idx]
            for idx2 in range(tdeg*pdeg):
                amatp[idxp, idx2] = amat[idx, idx2]

    # Perform least-squares fit to get alpha coefficients
    theta = numpy.linalg.lstsq(amatp, bvecp, rcond=RCOND)[0]

    alpha = numpy.zeros((tdeg, pdeg), numpy.float64)
    for tidx2 in range(tdeg):
        for pidx2 in range(pdeg):
            alpha[tidx2, pidx2] = theta[pidx2 * tdeg + tidx2]
    tlim = (tmin, tmax)
    plim = (pmin, pmax)

    return alpha, tlim, plim


def check_viability(ktp_dct):
    """ Checks if a Chebyshev fit is viable based on whether all the temp
        values are the same. Note: ignores 'high' if present.

        :param ktp_dct: rate constants to be fitted; should be P-dependent
        :type ktp_dct: dict {pressure: (temps, kts)}
        :return cheb_viable: whether or not a Chebyshev fit is viable
        :rtype: Bool
    """

    cheb_viable = True
    pressures = tuple(pressure for pressure in ktp_dct.keys()
                      if pressure != 'high')

    # Loop over each pressure
    for pidx, pressure in enumerate(pressures):
        if pidx == 0:
            baseline_temps = ktp_dct[pressure][0]  # store first temps

        temps = ktp_dct[pressure][0]
        if len(temps) != len(baseline_temps):  # check temps are same length
            cheb_viable = False
        elif not numpy.allclose(temps, baseline_temps):
            cheb_viable = False

    return cheb_viable


def get_one_atm_arr(ktp_dct):
    """ Gets Arrhenius parameters at 1 atm, if available

        :param ktp_dct: rate constants to be fitted; should be P-dependent
        :type ktp_dct: dict {pressure: (temps, kts)}
        :return one_atm_arr: single Arrhenius fitting parameters at 1 atm
        :rtype: tuple of tuples
    """

    # If there is rate constant data at 1 atm
    if 1 in ktp_dct.keys():
        high_ktp_dct = {'high': ktp_dct[1]}  # set to 'high' to use Arrhenius
        arr_params, _ = arr.get_params(high_ktp_dct, dbltol=500)  # only single
        one_atm_arr = arr_params.arr
    # Otherwise, just use fake parameters
    else:
        one_atm_arr = ((1.0, 0.0, 0.0),)

    return one_atm_arr
