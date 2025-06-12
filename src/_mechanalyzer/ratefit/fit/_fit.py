""" This is Clayton's new version
"""

import copy
import numpy
from ratefit.fit import arr
from ratefit.fit import plog
from ratefit.fit import cheb

DEFAULT_PDEP = {
    'temps': (500.0, 1000, 2000.0),
    'tol': 20.0,
    'plow': None,
    'phigh': None,
    'pval': 1.0}
DEFAULT_ARR = {  # also used for PLOG fitting
    'dbltol': 50.0,
    'dbl_iter': 30}
DEFAULT_TROE = {
    'params': ('ts1', 'ts2', 'ts3', 'alpha'),
    'tol': 20.0}
DEFAULT_CHEB = {
    'tdeg': 6,
    'pdeg': 4,
    'tol': 20.0}
NICKNAMES = {  # used for printing messages
    'arr': 'Arrhenius',
    'plog': 'PLOG',
    'cheb': 'Chebyshev',
    'troe': 'Troe'}
ALLOWED_FIT_METHODS = (
    'arr',
    'plog',
    'cheb',
    'troe')


def fit_rxn_ktp_dct(rxn_ktp_dct, fit_method, pdep_dct=None, arrfit_dct=None,
                    chebfit_dct=None, troefit_dct=None):
    """ Fits all reactions in a rxn_ktp_dct to some desired form

        :param rxn_ktp_dct: rate constants to be fitted, for multiple reactions
        :type rxn_ktp_dct: dict {rxn: ktp_dct}
        :param fit_method: desired fit form; 'arr', 'plog', 'cheb', or 'troe'
        :type fit_method: str
        :param pdep_dct: instructions for checking P dependence
        :type pdep_dct: dict
        :param arrfit_dct: instructions for Arrhenius fitting (also for PLOG)
        :type arrfit_dct: dict
        :param chebfit_dct: instructions for Chebyshev fitting
        :type chebfit_dct: dict
        :param troefit_dct: instructions for Troe fitting
        :type troefit_dct: dict
        :return rxn_param_dct: fitted parameters for each reaction
        :rtype: dict {rxn: params}
        :return rxn_err_dct: fitting errors for each reaction
        :rtype: dict {rxn: err_dct}
    """
    def rxn_name_str(rxn, newline=False):
        """ get a reaction name string
        """
        return ' = '.join((' + '.join(rxn[0]), ' + '.join(rxn[1])))

    rxn_param_dct = {}
    rxn_err_dct = {}
    for rxn, ktp_dct in rxn_ktp_dct.items():
        print(f'\nFitting Reaction: {rxn_name_str(rxn)}')
        params, err_dct = fit_ktp_dct(ktp_dct, fit_method, pdep_dct=pdep_dct,
                                      arrfit_dct=arrfit_dct,
                                      chebfit_dct=chebfit_dct,
                                      troefit_dct=troefit_dct)
        if all(x is not None for x in (params, err_dct)):
            rxn_param_dct[rxn] = params
            rxn_err_dct[rxn] = err_dct
        print('--------------------------------\n')

    return rxn_param_dct, rxn_err_dct


def fit_ktp_dct(ktp_dct, fit_method, pdep_dct=None, arrfit_dct=None,
                chebfit_dct=None, troefit_dct=None):
    """ Fits a single ktp_dct to some desired form

        :param ktp_dct: rate constants to be fitted
        :type ktp_dct: dict {pressure: (temps, kts)}
        :param fit_method: desired fit form; 'arr', 'plog', 'cheb', or 'troe'
        :type fit_method: str
        :param pdep_dct: instructions for checking P dependence
        :type pdep_dct: dict
        :param arrfit_dct: instructions for Arrhenius fitting (also for PLOG)
        :type arrfit_dct: dict
        :param chebfit_dct: instructions for Chebyshev fitting
        :type chebfit_dct: dict
        :param troefit_dct: instructions for Troe fitting
        :type troefit_dct: dict
        :return params: fitted parameters
        :rtype: autoreact.RxnParams object
        :return err_dct: fitting errors
        :rtype: dict {pressure: (temps, errs)}
    """

    # Load instructions, etc.
    assert fit_method in ALLOWED_FIT_METHODS, (
        f"Fit method should be in {ALLOWED_FIT_METHODS}, not '{fit_method}'")
    pdep_dct = pdep_dct or DEFAULT_PDEP
    arrfit_dct = arrfit_dct or DEFAULT_ARR
    chebfit_dct = chebfit_dct or DEFAULT_CHEB
    troefit_dct = troefit_dct or DEFAULT_TROE

    # If rates exist, fit them to desired functional form
    params, err_dct = None, None
    if ktp_dct:
        # Get the pdep_version of the ktp_dct
        pdep_ktp_dct = get_pdep_ktp_dct(
            ktp_dct, assess_temps=pdep_dct['temps'],
            tol=pdep_dct['tol'],
            plow=pdep_dct['plow'], phigh=pdep_dct['phigh'],
            pval=pdep_dct['pval'])

        # Check the ktp_dct and the fit_method to see how to fit rates
        actual_fit_method = assess_fit_method(pdep_ktp_dct, fit_method)

        # Get desired fit as instance of the RxnParams class, and the err_dct
        if actual_fit_method == 'troe':
            print('Troe not implemented. Not doing a fit')
        elif actual_fit_method == 'arr':
            # dbl_iter = arrfit_dct.get('dbl_iter')  # unused for now
            params, err_dct = arr.get_params(
                pdep_ktp_dct, dbltol=arrfit_dct['dbltol'])
        elif actual_fit_method == 'plog':
            # dbl_iter = arrfit_dct.get('dbl_iter')  # unused for now
            params, err_dct = plog.get_params(
                pdep_ktp_dct, dbltol=arrfit_dct['dbltol'])
        elif actual_fit_method == 'cheb':
            params, err_dct = cheb.get_params(
                pdep_ktp_dct,
                tdeg=chebfit_dct['tdeg'], pdeg=chebfit_dct['pdeg'],
                tol=chebfit_dct['tol'])
    else:
        print('No rate constants to fit.')

    return params, err_dct


def assess_fit_method(pdep_ktp_dct, fit_method):
    """ Checks a pdep_ktp_dct (i.e., already filtered for pressure dependence)
        to see if the selected fit method is acceptable. If something is amiss,
        corrects the fit method to the proper form.

        :param pdep_ktp_dct: pressure-dependent rate constants; either 'high'
            only (if P-independent) or with all pressures (if P-dependent)
        :type pdep_ktp_dct: dict {pressure: (temps, kts)}
        :param fit_method: desired fit form; 'arr', 'plog', 'cheb', or 'troe'
        :type fit_method: str
        :return actual_fit_method: fit method after checking P dependence; may
            be unchanged from original selection
        :rtype: str
    """

    pressures = list(pdep_ktp_dct.keys())
    npressures = len(pressures)
    if npressures == 1 or (npressures == 2 and 'high' in pressures):
        actual_fit_method = 'arr'
    else:
        # If the fit_method is Arrhenius but there are multiple pressures
        if fit_method == 'arr':
            actual_fit_method = 'plog'
        else:
            actual_fit_method = fit_method

    # If Chebyshev was requested and granted, check that it is viable to do
    # Chebyshev based on the number of k(T) values at each pressure
    if actual_fit_method == 'cheb':
        # Note: this check ignores the 'high' entry
        cheb_viable = cheb.check_viability(pdep_ktp_dct)
        if not cheb_viable:
            actual_fit_method = 'plog'

    # Print message to say what fitting will be done
    if actual_fit_method != fit_method:  # if the fit method was changed
        if actual_fit_method == 'plog':
            # If changed to PLOG from Arrhenius
            if fit_method == 'arr':
                print('\nPressure dependence found. Fitting to PLOG form...')
            elif fit_method == 'cheb':
                print('\nRates invalid for Chebyshev. Fitting to PLOG form...')
        else:  # if changed to Arrhenius from anything else
            print(f'\nNot enough pressures for {NICKNAMES[fit_method]} form.')
            print('\nFitting k(T,P)s to Arrhenius form...')
    elif fit_method is None:  # if the ktp_dct was empty
        print('\nNo valid k(T,Ps)s to fit. Skipping to next reaction...')
    else:  # if the fit method was unchanged
        print(f'\nFitting k(T,P)s to {NICKNAMES[fit_method]} form...')

    return actual_fit_method


def get_pdep_ktp_dct(ktp_dct, assess_temps=DEFAULT_PDEP['temps'],
                     tol=DEFAULT_PDEP['tol'], plow=DEFAULT_PDEP['plow'],
                     phigh=DEFAULT_PDEP['phigh'], pval=DEFAULT_PDEP['pval']):
    """ Returns a ktp_dct with either (i) P-dependent or (ii) P-independent
        rate constants, depending on the P dependence of the rate constants

        :param ktp_dct: rate constants to be fitted
        :type ktp_dct: dict {pressure: (temps, kts)}
        :param assess_temps: temperature(s) at which to assess P dependence
        :type assess_temps: tuple
        :param tol: percent difference threshold for determining P dependence
        :type tol: float
        :param plow: low pressure for assessing P dependence; if None, the
            lowest pressure will be used
        :type plow: float
        :param phigh: high pressure for assessing P dependence; if None, the
            highest pressure will be used
        :type phigh: float
        :param pval: pressure at which to get kts if ktp_dct is P-independent;
            if pval isn't in the ktp_dct, gets kts at the highest pressure
        :return pdep_ktp_dct: pressure-dependent rate constants; either 'high'
            only (if P-independent) or with all pressures (if P-dependent)
        :rtype: dict {pressure: (temps, kts)}
    """

    # Assess the pressure dependence of the rate constants
    rxn_is_pdependent = assess_pdep(ktp_dct, assess_temps, tol=tol, plow=plow,
                                    phigh=phigh)

    if rxn_is_pdependent:
        print('Reaction found to be pressure dependent.',
              'Fitting all k(T)s from all pressures.')
        pdep_ktp_dct = copy.deepcopy(ktp_dct)
    else:
        # If pval is in the dct, return kts at that pressure
        if pval in ktp_dct:
            pdep_ktp_dct = {pval: ktp_dct[pval]}
            print(f'No pressure dependence. Grabbing k(T)s at {pval} atm.')
        # Otherwise, get kts at some other pressure
        else:
            print(f'No pressure dependence, but no k(T)s at {pval} atm.')
            pressures = [pressure for pressure in ktp_dct
                         if pressure != 'high']
            if pressures:
                pressure = pressures[-1]
                pdep_ktp_dct = {pressure: ktp_dct[pressure]}
                print(f'Grabbing kts at {pressure} atm.')
            else:
                pdep_ktp_dct = {'high': ktp_dct['high']}
                print('No numerical pressures. Grabbing "high" kts.')

    return pdep_ktp_dct


def assess_pdep(ktp_dct, assess_temps=DEFAULT_PDEP['temps'],
                tol=DEFAULT_PDEP['tol'], plow=DEFAULT_PDEP['plow'],
                phigh=DEFAULT_PDEP['phigh']):
    """ Assesses if there are significant differences between k(T,P) values
        at low and high pressure, signaling pressure dependence

        :param ktp_dct: rate constants to be fitted
        :type ktp_dct: dict {pressure: (temps, kts)}
        :param assess_temps: temperature(s) at which to assess P dependence
        :type assess_temps: tuple
        :param tol: percent difference threshold for determining P dependence
        :type tol: float
        :param plow: low pressure for assessing P dependence; if None, the
            lowest pressure will be used
        :type plow: float
        :param phigh: high pressure for assessing P dependence; if None, the
            highest pressure will be used
        :type phigh: float
        :param pval: pressure at which to get kts if ktp_dct is P-independent;
            if pval isn't in the ktp_dct, gets kts at the highest pressure
        :return pdep_ktp_dct: pressure-dependent rate constants; either 'high'
            only (if P-independent) or with all pressures (if P-dependent)
        :rtype: dict {pressure: (temps, kts)}
    """

    def is_pdep_atT(ktp_dct, plow, phigh, assess_temp):
        is_pdep = False
        temps_low = ktp_dct[plow][0]
        temps_high = ktp_dct[phigh][0]
        temp_low_match = numpy.where(
            numpy.isclose(temps_low, assess_temp))[0]
        temp_high_match = numpy.where(
            numpy.isclose(temps_high, assess_temp))[0]
        if temp_low_match.size > 0 and temp_high_match.size > 0:
            temp_low_idx = temp_low_match[0]
            temp_high_idx = temp_high_match[0]
            # Grab k value for the appropriate temp and pressure
            kt_low = ktp_dct[plow][1][temp_low_idx]
            kt_high = ktp_dct[phigh][1][temp_high_idx]
            # Calculate the % difference and see if above threshold
            kt_dif = (abs(kt_low - kt_high) / kt_low) * 100.0
            if kt_dif > tol:
                is_pdep = True
        return is_pdep

    # Get list of pressures, ignoring the high-pressure limit rates
    pressures = [pressure for pressure in ktp_dct
                 if pressure != 'high']

    # Set the low- and high-pressure if not specified by user
    if plow is None and pressures != []:
        plow = min(pressures)
    if phigh is None and pressures != []:
        phigh = max(pressures)

    # Check % difference for k(T, P) vals
    is_pdep = False
    if plow in ktp_dct and phigh in ktp_dct:  # won't get run if only 'high'
        # Loop over temps to examine for large % dif in k(T) at low- and high-P
        for assess_temp in assess_temps:
            # For low and high P, find the idx for assess_temp
            is_pdep = is_pdep_atT(ktp_dct, plow, phigh, assess_temp)
            if is_pdep:
                break
        # additional check: highest T at each P compared to the highest P
        # unless you specifically indicated a temperature range
        if assess_temps == DEFAULT_PDEP['temps']:
            pressures.sort()
            for plow in pressures[:-1]:
                temps_p_max = max(ktp_dct[plow][0])
                is_pdep = is_pdep_atT(ktp_dct, plow, phigh, temps_p_max)
                if is_pdep:
                    break
    return is_pdep
