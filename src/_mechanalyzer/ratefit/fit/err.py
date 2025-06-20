""" Calculates errors between reference rate constants and fitted parameters
"""

from mechanalyzer.calculator import rates


def get_err_dct(ref_ktp_dct, params):
    """ Calculates an err_dct by comparing a reference ktp_dct with a ktp_dct
        calculated from provided fitting parameters

        :param ref_ktp_dct: reference ktp_dct for comparison
        :type ref_ktp_dct: dict
        :param params: fitting parameters
        :type params: autoreact.RxnParams object
        :return err_dct: fitting errors for a single reaction
        :rtype: dict {pressure: (temps, errs)}
    """

    # Calculate a ktp_dct using the fitting params
    temps_lst, pressures = get_temps_pressures(ref_ktp_dct)
    fit_ktp_dct = rates.eval_params(params, temps_lst, pressures)

    # Calculate the err_dct
    err_dct = {}
    for pressure, (temps, fit_kts) in fit_ktp_dct.items():
        ref_kts = ref_ktp_dct[pressure][1]
        errs = 100 * (fit_kts - ref_kts) / ref_kts
        err_dct[pressure] = (temps, errs)

    return err_dct


def get_max_err(err_dct):
    """ Gets the singular max (absolute) error from an err_dct

        :param err_dct: fitting errors for a single reaction
        :type err_dct: dict
        :return max_err: maximum absolute error in an err_dct
        :rtype: float
    """

    max_err = 0
    for (_, errs) in err_dct.values():
        if max(abs(errs)) > max_err:
            max_err = max(abs(errs))

    return max_err


def get_temps_pressures(ktp_dct):
    """ Reads a ktp_dct and gets the list of pressure and corresponding list of
        temperature arrays

        :param ktp_dct: k(T,P) values
        :type ktp_dct: dict {pressure: (temps, kts)}
        :return temps_lst: list of temperature arrays at each pressure (K)
        :rtype: list [numpy.ndarray1, numpy.ndarray2, ...]
        :return pressures: list of pressures (atm)
        :rtype: list
    """

    temps_lst = []
    pressures = []
    for pressure, (temps, _) in ktp_dct.items():
        temps_lst.append(temps)
        pressures.append(pressure)

    return temps_lst, pressures
