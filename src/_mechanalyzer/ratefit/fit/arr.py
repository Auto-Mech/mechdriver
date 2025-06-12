""" Fits rate constants to a modified Arrhenius (single or double) form
"""

import numpy
from scipy.optimize import leastsq
from scipy.optimize import least_squares
from autoreact.params import RxnParams
from phydat import phycon
from ratefit.fit import err

RC = phycon.RC_CAL  # universal gas constant in cal/mol-K
GUESS_BNDS = ((1, 1e4), (0.1, 20), (1, 100))  # guess bounds for double fitting


def get_params(ktp_dct, dbltol=15, dbl_iter=1, tref=1.0):
    """ Gets the fitting parameters for an Arrhenius fit to rate constant data.
        Also gets the errors of that fit. Performs either a single or double
        Arrhenius fit.

        :param ktp_dct: rate constants to be fitted; should be P-independent
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

    assert len(ktp_dct) == 1, (
        f'ktp_dct for Arrhenius fit needs only 1 pressure, not {len(ktp_dct)}')
    pressure = list(ktp_dct.keys())[0]  # get first (and only) pressure

    # Perform single fit and assess its errors
    (temps, kts) = ktp_dct[pressure]  # read data
    sing_params = single_arr(temps, kts)
    sing_err_dct = err.get_err_dct(ktp_dct, sing_params)
    sing_max_err = err.get_max_err(sing_err_dct)

    # Put in checks to see if a single fit should be used
    # Sees if enough rate constants to do double fit, and error necessitates it
    use_single_fit = False
    if len(kts) < 6:
        print('There are less than six rate constants to fit. Using single fit.')
        use_single_fit = True
    else:
        if sing_max_err < dbltol:
            print(f'Single fit error is {sing_max_err:.1f}%, which is less than'
                  f' the input limit of {dbltol}%. Using single fit.')
            use_single_fit = True
            
    # Perform double fit if needed or just take single fit params
    if use_single_fit:
        params = sing_params
        err_dct = sing_err_dct
    else:
        print('Attempting double fit...')

        # Perfom a double fit
        doub_params, guess_idx = double_arr(temps, kts, sing_params, tref=tref,
                                            dbltol=dbltol, dbl_iter=dbl_iter)

        # Assess errors
        doub_err_dct = err.get_err_dct(ktp_dct, doub_params)
        doub_max_err = err.get_max_err(doub_err_dct)
        print(f'Double fit obtained with max error of {doub_max_err:.1f}% '
              f'after {guess_idx + 1} iteration(s).')

        # Use single fit if double fit is worse than single fit
        if doub_max_err > sing_max_err:
            print('Double fit error is worse than single; using single fit.')
            params = sing_params
            err_dct = sing_err_dct
        # Use single fit if double fit contains infinite values
        elif check_for_inf(doub_params):
            print('Double fit error contains inf values; using single fit.')
            params = sing_params
            err_dct = sing_err_dct
        # Otherwise, use the double fit
        else:
            params = doub_params
            err_dct = doub_err_dct

    return params, err_dct


def single_arr(temps, kts, tref=1.0):
    """ Fits pressure-independent rate constants to a single Arrhenius form

        :param temps: temperatures at which rate constants are defined (K)
        :type temps: numpy.ndarray of shape (num_temps,)
        :param kts: rate constants
        :type kts: numpy.ndarray of shape (num_temps,)
        :param tref: reference temp for the modified Arrhenius form
        :type tref: float
        :return params: fitted single Arrhenius parameters
        :rtype: autoreact.RxnParams object
    """

    # Consider several cases depending on the number of valid rate constants
    # no k is positive, so return all zeros
    if kts.size == 0:
        a_fit, n_fit, ea_fit = 0.0, 0.0, 0.0

    # if num(k) is 1: set A = k
    elif kts.size == 1:
        a_fit, n_fit, ea_fit = kts[0], 0.0, 0.0

    # if num(k) > 0 is 2,3: fit A and Ea
    elif kts.size in (2, 3):
        # Build vectors and matrices used for the fitting
        a_vec = numpy.ones(len(temps))
        ea_vec = (-1.0 / RC) * (1.0 / temps)
        coeff_mat = numpy.array([a_vec, ea_vec], dtype=numpy.float64)
        coeff_mat = coeff_mat.transpose()
        k_vec = numpy.log(kts)
        # Perform the least-squares fit
        theta = numpy.linalg.lstsq(coeff_mat, k_vec, rcond=None)[0]
        # Set the fitting parameters
        a_fit, n_fit, ea_fit = numpy.exp(theta[0]), 0.0, theta[1]

    # if num(k) > 0 is more than 3: fit A, n, and Ea
    elif kts.size > 3:
        # Build vectors and matrices used for the fitting
        a_vec = numpy.ones(len(temps))
        n_vec = numpy.log(temps / tref)
        ea_vec = (-1.0 / RC) * (1.0 / temps)
        coeff_mat = numpy.array([a_vec, n_vec, ea_vec], dtype=numpy.float64)
        coeff_mat = coeff_mat.transpose()
        k_vec = numpy.log(kts)
        # Perform the least-squares fit
        theta = numpy.linalg.lstsq(coeff_mat, k_vec, rcond=None)[0]
        # Set the fitting parameters
        a_fit, n_fit, ea_fit = numpy.exp(theta[0]), theta[1], theta[2]
    # check on nans
    if any(numpy.isnan([a_fit, n_fit, ea_fit])):
        print('*Warning: NaN params found- set all params to 0.')
        a_fit, n_fit, ea_fit = 0.0, 0.0, 0.0
    # Pack the parameters into an arr_dct and instantiate RxnParams
    arr_dct = {'arr_tuples': [[a_fit, n_fit, ea_fit]]}
    params = RxnParams(arr_dct=arr_dct)

    return params


def double_arr(temps, kts, sing_params, tref=1.0, dbltol=15, dbl_iter=1):
    """ Fit a set of T-dependent rate constants to a double Arrhenius form; can
        do so by iterating across a range of guesses by modifying the provided
        single Arrhenius fitting parameters

        :param temps: temperatures at which rate constants are defined (K)
        :type temps: numpy.ndarray of shape (num_temps,)
        :param kts: rate constants
        :type kts: numpy.ndarray of shape (num_temps,)
        :param sing_params: the best-fit single Arrhenius parameters
        :type sing_params: autoreact.RxnParams object
        :param tref: reference temp for the modified Arrhenius form
        :type tref: float
        :param dbltol: the tolerance at which to reject a single fit
        :type dbltol: float
        :param dbl_iter: max number of iterations for double fitting
        :type dbl_iter: int
        :return params: fitted double Arrhenius parameters
        :rtype: autoreact.RxnParams object
        :return guess_idx: number of double fits performed
        :rtype: int
    """

    def fit_doub_arr(temps, kts, sing_params, a_change, n_change, tref=1.0,
                     allow_neg=False):
        """ Performs one double Arrhenius fit by generating initial guesses
            using changes to the single Arrhenius fit

            :param a_change: amount by which to vary the A factor
            :type a_change: float
            :param n_change: amount by which to vary the temperature exponent
            :type n_change: float
            :return params: fitted double Arrhenius parameters
            :rtype: autoreact.RxnParams object
            (all other inputs same as parent function)
        """

        sing_a, sing_n, sing_ea = sing_params.arr[0]  # get first (only) entry

        # Get a new tref for the double fit: the logarithmic midpoint temp
        doub_tref = numpy.sqrt(max(temps) / min(temps)) * min(temps)

        # Generate guesses by varying single parameters
        sing_a = sing_a * (doub_tref / tref) ** sing_n  # convert to new basis
        init_guess = [(sing_a * a_change), (sing_n + n_change), sing_ea,
                      (sing_a * (1 - a_change)), (sing_n - n_change), sing_ea]

        # Set bounds: np.inf works better than setting e.g., 1e+300, but slower
        if allow_neg:  # no bounds
            bounds = ([-numpy.inf, -numpy.inf, -numpy.inf, -numpy.inf, 
                       -numpy.inf, -numpy.inf], [numpy.inf, numpy.inf,
                       numpy.inf, numpy.inf, numpy.inf, numpy.inf]) 
        else:  # constrain A factors to be positive
            bounds = ([0, -numpy.inf, -numpy.inf, 0, 
                       -numpy.inf, -numpy.inf], [numpy.inf, numpy.inf,
                       numpy.inf, numpy.inf, numpy.inf, numpy.inf]) 

        # Perform a least-squares fit
        # note: previous version scipy.optimize.leastsq (unbounded): used method='lm'
        # same or better results obtained with x_scale='jac' for new cases tested
        try:
            plsq = least_squares(_resid_func, init_guess, bounds=bounds,
                                args=(temps, kts, doub_tref), x_scale = 'jac', #method = 'lm',
                                ftol=1.0E-8, xtol=1.0E-8, max_nfev=100000)
        except ValueError: #test older method without bounds
            try:
                plsq = least_squares(_resid_func, init_guess,
                                    args=(temps, kts, doub_tref), method='lm',
                                    ftol=1.0E-8, xtol=1.0E-8, max_nfev=100000)
            except ValueError:
                try:
                    plsq = least_squares(_resid_func, init_guess,
                                        loss = 'arctan',
                                        args=(temps, kts, doub_tref),
                                        ftol=1.0E-8, xtol=1.0E-8, max_nfev=100000)
                except ValueError:
                    plsq = None
                    raw_params = [numpy.inf, numpy.inf,
                        numpy.inf, numpy.inf, numpy.inf, numpy.inf] 
        # Retrieve the fit params and convert A back to the input tref
        if plsq:
            raw_params = list(plsq.x)  # list of length 6
            raw_params[0] = raw_params[0] * (tref / doub_tref) ** raw_params[1]
            raw_params[3] = raw_params[3] * (tref / doub_tref) ** raw_params[4]

        # Instantiate RxnParams
        arr_dct = {'arr_tuples': [raw_params[:3], raw_params[3:]]}
        params = RxnParams(arr_dct=arr_dct)

        return params

    # Create an array of predefined values to guess if the initial guess fails
    a_changes = [0.1, 0.3, 0.5, 0.7, 0.9]
    n_changes = [1.2, 1.5, 1.9, 2.5, 3]
    predef_iter = len(a_changes) * len(n_changes) + 1  # +1 for the first guess

    # Make a maximum of dbl_iter attempts at a double fit
    max_errs = []
    prev_params = []
    # Note: using 'high' here in place of the actual pressure; doesn't matter
    # since this ref_ktp_dct never gets returned
    ref_ktp_dct = {'high': (temps, kts)}  # used for err_dct later
    for guess_idx in range(dbl_iter):
        # Use SJK's guesses as first try
        if guess_idx == 0:
            a_change = 0.5
            n_change = 2
        # Otherwise, loop over the predefined A and n changes
        else:
            a_change = a_changes[(guess_idx - 1) % 5]
            n_change = n_changes[int((guess_idx - 1) / 5)]

        # Perform a double fit
        params = fit_doub_arr(temps, kts, sing_params, a_change, n_change,
                              tref=tref)

        # Get an err_dct and the max error
        err_dct = err.get_err_dct(ref_ktp_dct, params)
        max_err = err.get_max_err(err_dct)
        max_errs.append(max_err)
        prev_params.append(params)

        # Exit the loop if tolerance is satisfied
        if max_err < dbltol:
            break
        # If tolerance not satisfied and guesses have been exhausted, grab
        # the result with the lowest error
        if guess_idx in (dbl_iter - 1, predef_iter - 1):
            min_idx = numpy.argmin(max_errs)
            params = prev_params[min_idx]
            break

    return params, guess_idx


def _resid_func(curr_guess, temps, kts, tref):
    """ Computes the residual between fit and data for double fitter

        :param curr_guess: current guess for double Arrhenius params
        :type curr_guess: list [A1, n1, Ea1, A2, n2, Ea2]
        :param temps: temperatures at which rate constants are defined (K)
        :type temps: numpy.ndarray of shape (num_temps,)
        :param kts: rate constants
        :type kts: Numpy.ndarray of shape (num_temps,)
        :param tref: reference temp for the modified Arrhenius form
        :type tref: float
        :return resid: residual between fit and data
        :rtype: Numpy.ndarray of shape (num_temps,)
    """

    # Compute the fitted rate constants
    k_fit1 = curr_guess[0] * numpy.exp(
        curr_guess[1] * numpy.log(temps / tref) - curr_guess[2] / (RC * temps))
    k_fit2 = curr_guess[3] * numpy.exp(
        curr_guess[4] * numpy.log(temps / tref) - curr_guess[5] / (RC * temps))
    k_fit = k_fit1 + k_fit2
    resid = numpy.log10(kts) - numpy.log10(k_fit)

    return resid


def check_for_inf(params):
    """ Checks for infinite values in fitted Arrhenius parameters

        :param params: fitted Arrhenius parameters
        :type params: autoreact.RxnParams object
        :return contains_inf: whether Arrhenius parameters contain infinity
        :rtype: Bool
    """

    arr_tuples = params.arr
    contains_inf = False
    for arr_tuple in arr_tuples:
        for element in arr_tuple:
            if numpy.isinf(element):
                contains_inf = True
                break

    return contains_inf
