"""
Calculate rates with various fitting functions
"""

import sys
import copy
import numpy
import pandas
from phydat import phycon
from scipy.special import eval_chebyt
from mechanalyzer.calculator import thermo

RC = phycon.RC_CAL  # gas constant in cal/(mol.K)
RC2 = phycon.RC_ATM  # gas constant in cm^3.atm/(mol.K)


def eval_rxn_param_dct(rxn_param_dct, temps_lst, pressures, tref=1.0):
    """ Loops through all rxns in a rxn_param_dct and gets a ktp_dct for
        each one

        :param rxn_param_dct: rate parameters for all rxns in a mech
        :type rxn_param_dct: dict {rxn: params}
        :param temps_lst: list of temperature arrays used to get k(T,P)s (K)
        :type temps_lst: list [numpy.ndarray1, numpy.ndarray2, ...]
        :param pressures: pressures used to get k(T,P)s (atm)
        :type pressure: list
        :return rxn_ktp_dct: k(T,Ps) at all temps and pressures for each rxn
        :rtype: dict {rxn: ktp_dct}
    """

    temps_lst = check_p_t(temps_lst, pressures)  # enforce formatting rules

    rxn_ktp_dct = {}
    for rxn, params in rxn_param_dct.items():
        ktp_dct = eval_params(params, temps_lst, pressures, tref=tref)
        rxn_ktp_dct[rxn] = ktp_dct

    return rxn_ktp_dct


def eval_params(params, temps_lst, pressures, tref=1.0):
    """ Look through a params and evaluate k(T,P) based on the contents.
        Return a ktp_dct.

        :param params: object describing the functional fits
        :type params: autochem/autoreact RxnParams object
        :param temps_lst: list of temperature arrays used to get k(T,P)s (K)
        :type temps_lst: list [numpy.ndarray1, numpy.ndarray2, ...]
        :param pressures: pressures used to get k(T,P)s (atm)
        :type pressure: list
        :return ktp_dct: k(T,Ps) at all temps and pressures
        :rtype: dict {pressure: (temps, kts)}
    """

    temps_lst = check_p_t(temps_lst, pressures)  # enforce formatting rules

    # Get the list of existing forms
    forms = params.get_existing_forms()
    assert forms != (), 'The params object is empty'

    # Loop over each rxn type (usually only one, but trying to be general)
    ktp_dct = {}
    for form in forms:
        new_ktp_dct = {}
        if form == 'arr':
            # Pick the temperatures to use in the Arrhenius assessment
            if 'high' in pressures:
                temps = temps_lst[pressures.index('high')]
                pressure = 'high'
            else:
                temps = temps_lst[-1]  # use the last pressure
                pressure = pressures[-1]

            arr_tuples = params.arr
            kts = arr(arr_tuples, temps, tref)
            new_ktp_dct[pressure] = (temps, kts)

        elif form == 'plog':
            plog_dct = params.plog
            new_ktp_dct = plog(plog_dct, temps_lst, pressures, tref=tref)

        elif form == 'cheb':
            cheb_dct = params.cheb
            alpha = cheb_dct['alpha']
            tlim = cheb_dct['tlim']
            plim = cheb_dct['plim']
            new_ktp_dct = cheb(alpha, tlim, plim, temps_lst, pressures)

        elif form == 'troe':
            troe_dct = params.troe
            highp_arr = troe_dct['highp_arr']
            lowp_arr = troe_dct['lowp_arr']
            troe_params = troe_dct['troe_params']
            new_ktp_dct = troe(highp_arr, lowp_arr, troe_params, temps_lst,
                               pressures, tref=tref)

        elif form == 'lind':
            lind_dct = params.lind
            highp_arr = lind_dct['highp_arr']
            lowp_arr = lind_dct['lowp_arr']
            new_ktp_dct = lind(highp_arr, lowp_arr, temps_lst, pressures,
                               tref=tref)

        ktp_dct = add_ktp_dcts(new_ktp_dct, ktp_dct)

    # Deal with the unusual (i.e., wrong) case of duplicate fits of one form
    dup_ktp_dct = handle_duplicates(params, temps_lst, pressures, tref=tref)
    ktp_dct = add_ktp_dcts(dup_ktp_dct, ktp_dct)

    return ktp_dct


def handle_duplicates(params, temps_lst, pressures, tref=1.0):
    """ Calculates any unusual duplicate cases. These only occur in strange
        cases when duplicates of a functional form (e.g., two PLOGs)
        are described for the same reaction. These should not really occur,
        but are nonetheless handled here.

        :param params: object describing the functional fits
        :type params: autochem/autoreact RxnParams object
        :param temps_lst: list of temperature arrays used to get k(T,P)s (K)
        :type temps_lst: list [numpy.ndarray1, numpy.ndarray2, ...]
        :param pressures: pressures used to get k(T,P)s (atm)
        :type pressure: list
        :return dup_ktp_dct: k(T,Ps) at all temps and pressures for duplicates
        :rtype: dict {pressure: (temps, kts)}
    """

    _, dup_counts = params.check_for_dups()  # dups is a Boolean

    dup_ktp_dct = {}
    for form, dup_count in dup_counts.items():
        # Note that Arrhenius is absent because this is allowed to have dups
        # and is thus not an unusual case
        if form == 'plog':
            for dup_idx in range(dup_count):
                plog_dct = params.plog_dups[dup_idx]
                new_ktp_dct = plog(plog_dct, temps_lst, pressures, tref=tref)
                dup_ktp_dct = add_ktp_dcts(new_ktp_dct, dup_ktp_dct)

        elif form == 'cheb':
            for dup_idx in range(dup_count):
                cheb_dct = params.cheb_dups[dup_idx]
                alpha = cheb_dct['alpha']
                tlim = cheb_dct['tlim']
                plim = cheb_dct['plim']
                new_ktp_dct = cheb(alpha, tlim, plim, temps_lst, pressures)
                dup_ktp_dct = add_ktp_dcts(new_ktp_dct, dup_ktp_dct)

        elif form == 'troe':
            for dup_idx in range(dup_count):
                troe_dct = params.troe_dups[dup_idx]
                highp_arr = troe_dct['highp_arr']
                lowp_arr = troe_dct['lowp_arr']
                troe_params = troe_dct['troe_params']
                new_ktp_dct = troe(highp_arr, lowp_arr, troe_params, temps_lst,
                                   pressures, tref=tref)
                dup_ktp_dct = add_ktp_dcts(new_ktp_dct, dup_ktp_dct)

        elif form == 'lind':
            for dup_idx in range(dup_count):
                lind_dct = params.lind_dups[dup_idx]
                highp_arr = lind_dct['highp_arr']
                lowp_arr = lind_dct['lowp_arr']
                new_ktp_dct = lind(highp_arr, lowp_arr, temps_lst, pressures,
                                   tref=tref)
                dup_ktp_dct = add_ktp_dcts(new_ktp_dct, dup_ktp_dct)

    return dup_ktp_dct


def arr(arr_tuples, temps, tref, rval=RC):
    """ Calculates T-dependent rate constants [k(T)]s using a list of any number
            of Arrhenius parameters

        :param arr_tuples: Arrhenius fit parameters
        :type arr_tuples: tuple ((A1, n1, Ea1), (A2, n2, Ea2), ...)
        :param temps: temperature array used to get k(T)s (K)
        :type temps: numpy.ndarray
        :param tref: reference temperature used for modified Arrhenius (K)
        :type tref: float
        :return kts: k(T)s for the given temperatures
        :rtype: numpy.ndarray
    """

    kts = numpy.zeros(len(temps))
    for arr_tuple in arr_tuples:
        assert len(arr_tuple) == 3, (
            f'Length of Arrhenius tuple should be 3, not {len(arr_tuple)}')
        a_par, n_par, ea_par = arr_tuple
        kts += a_par * ((temps/tref)**n_par) * numpy.exp(-ea_par/(rval*temps))

    return kts


def plog(plog_dct, temps_lst, pressures, tref=1.0):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a PLOG functional expression.

        :param plog_dct: Arrhenius fitting parameters at several pressures
        :type plog_dct: dict {pressure: arr_tuples}
        :param temps_lst: list of temperature arrays used to get k(T,P)s (K)
        :type temps_lst: list [numpy.ndarray1, numpy.ndarray2, ...]
        :param pressures: pressures used to get k(T,P)s (atm)
        :type pressure: list
        :param tref: reference temperature used for modified Arrhenius (K)
        :type tref: float
        :return ktp_dct: k(T,Ps) at all temps and pressures
        :rtype: dict {pressure: (temps, kts)}
    """

    def plog_one_p(plog_dct, temps, pressure, tref=1.0):
        """ Calculates T,P-dependent rate constants [k(T,P)]s using
            a PLOG functional expression, at a single pressure,
            across several temperatures.
        """

        # Check if pressure in plog dct; use plog pressure for numerical stab
        pressure_defined = False
        plog_pressures = list(plog_dct.keys())
        for plog_pressure in plog_pressures:
            if numpy.isclose(pressure, plog_pressure, rtol=1.0e-2):
                pressure_defined = True
                plog_params = plog_dct[plog_pressure]

        # If pressure equals value use, arrhenius expression
        if pressure_defined:
            kts = arr(plog_params, temps, tref=tref)
        else:
            # Calculate pressure term for PLOG expression
            # Use two PLOG pressures our pressure of interest sits between
            for i, _ in enumerate(plog_pressures):
                if i != len(plog_pressures)-1:
                    if plog_pressures[i] < pressure < plog_pressures[i+1]:
                        plow = plog_pressures[i]
                        phigh = plog_pressures[i+1]
                        plow_params = plog_dct[plow]
                        phigh_params = plog_dct[phigh]
                        break
            # Note: log10 instead of ln; no difference
            pres_term = ((numpy.log10(pressure) - numpy.log10(plow)) /
                         (numpy.log10(phigh) - numpy.log10(plow)))

            # Calculate k(T)s at high-P and low-P with Arrhenius expressions
            kts_low = arr(plow_params, temps, tref=tref)
            kts_high = arr(phigh_params, temps, tref=tref)

            # Calculate K(T,P)s with PLOG expression
            log_kts = (
                numpy.log10(kts_low) +
                ((numpy.log10(kts_high)-numpy.log10(kts_low)) * pres_term)
            )
            kts = 10**(log_kts)

        return kts

    # Remove 'high' from pressures and the corresponding temperature array
    temps_lst, pressures = remove_high(temps_lst, pressures)

    # Set the plog pressures to see if pressure is in range
    max_plog = max(list(plog_dct.keys()))
    min_plog = min(list(plog_dct.keys()))

    # Build the ktp_dct
    ktp_dct = {}
    for pidx, pressure in enumerate(pressures):
        temps = temps_lst[pidx]
        # If input pressure falls between min and max pressures
        if min_plog <= pressure <= max_plog:
            kts = plog_one_p(plog_dct, temps, pressure, tref=tref)
            ktp_dct[pressure] = (temps, kts)
        # If input pressure falls below min pressure, use min pressure
        elif pressure < min_plog:
            kts = plog_one_p(plog_dct, temps, min_plog, tref=tref)
        # If input pressure falls above max pressure, use max pressure
        elif pressure > max_plog:
            kts = plog_one_p(plog_dct, temps, max_plog, tref=tref)
        ktp_dct[pressure] = (temps, kts)

    return ktp_dct


def cheb(alpha, tlim, plim, temps_lst, pressures):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a Chebyshev functional expression.

        :param alpha: Chebyshev coefficient matrix
        :type alpha: numpy.ndarray
        :param tlim: minimum temperature Chebyshev model is defined
        :type tlim: float
        :param plim: minimum pressure Chebyshev model is defined
        :type plim: float
        :param temps_lst: list of temperature arrays used to get k(T,P)s (K)
        :type temps_lst: list [numpy.ndarray1, numpy.ndarray2, ...]
        :param pressures: pressures used to get k(T,P)s (atm)
        :type pressure: list
        :return ktp_dct: k(T,Ps) at all temps and pressures
        :rtype: dict {pressure: (temps, kts)}
    """

    def cheb_one_p(alpha, tlim, plim, temps, pressure):
        """ Calculates T,P-dependent rate constants [k(T,P)]s using
            a Chebyshev functional expression, at a single pressure,
            across several temperatures.

            :param temps: temperature array used to get k(T)s (K)
            :type temps: numpy.ndarray
            :return kts: k(T)s at the given pressure
            :rtype: numpy.ndarray
            (other inputs unchanged from parent function)
        """

        tmin, tmax = tlim
        pmin, pmax = plim
        logp = numpy.log10(pressure)
        logpmin = numpy.log10(pmin)
        logpmax = numpy.log10(pmax)
        alpha_nrows, alpha_ncols = alpha.shape

        kts = numpy.zeros(len(temps))
        for i, temp in enumerate(temps):
            ctemp = (
                (2.0 * temp**(-1) - tmin**(-1) - tmax**(-1)) /
                (tmax**(-1) - tmin**(-1)))
            cpress = (
                (2.0 * logp - logpmin - logpmax) /
                (logpmax - logpmin))

            log_kt = 0.0
            for j in range(alpha_nrows):
                for k in range(alpha_ncols):
                    log_kt += (alpha[j][k] * eval_chebyt(j, ctemp) *
                               eval_chebyt(k, cpress))

            kts[i] = 10**(log_kt)

        return kts

    # Remove 'high' from pressures and the corresponding temperature array
    temps_lst, pressures = remove_high(temps_lst, pressures)

    ktp_dct = {}
    for pidx, pressure in enumerate(pressures):
        temps = temps_lst[pidx]
        kts = cheb_one_p(alpha, tlim, plim, temps, pressure)
        ktp_dct[pressure] = (temps, kts)

    return ktp_dct


def troe(highp_arr, lowp_arr, troe_params, temps_lst, pressures,
         collid_factor=1.0, tref=1.0):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a Troe functional expression.

        :param highp_arr: high-P limit Arrhenius parameters
        :type highp_arr: list [[A1, n1, Ea1], [A2, n2, Ea2,], ...]
        :param lowp_arr: low-P limit Arrhenius parameters
        :type lowp_arr: list [[A1, n1, Ea1], [A2, n2, Ea2,], ...]
        :param troe_params: 3 or 4 Troe fitting coefficients: alpha, T***, T*,
            and T** (T** is optional)
        :type troe_params: list
        :param temps_lst: list of temperature arrays used to get k(T,P)s (K)
        :type temps_lst: list [numpy.ndarray1, numpy.ndarray2, ...]
        :param pressures: pressures used to get k(T,P)s (atm)
        :type pressure: list
        :param collid_factor: collider efficiency factor
        :type collid_factor: float
        :param tref: reference temperature used for modified Arrhenius (K)
        :type tref: float
        :return ktp_dct: k(T,Ps) at all temps and pressures
        :rtype: dict {pressure: (temps, kts)}
    """

    def troe_one_p(highp_arr, lowp_arr, troe_params, temps, pressure,
                   collid_factor=1.0, tref=1.0):
        """ Calculates T,P-dependent rate constants [k(T,P)]s using
            a Troe functional expression, at a single pressure,
            across several temperatures.

            :param temps: temperature array used to get k(T)s (K)
            :type temps: numpy.ndarray
            :return kts: k(T)s at the given pressure
            :rtype: numpy.ndarray
            (other inputs unchanged from parent function)
        """

        def _f_broadening_term(pr_term, troe_params, temps):
            """ Calculates the F broadening factor term

                :param pr_term: reduced pressure terms at all temps
                :type pr_term: numpy.ndarray
                :return f_cent: F broadening term
                :rtype: float
                (other inputs unchanged from parent function)
            """

            if len(troe_params) == 3:
                alpha, ts3, ts1 = troe_params
                ts2 = numpy.nan  # prevents annoying Pylint warning
            else:
                alpha, ts3, ts1, ts2 = troe_params

            # Calculate Fcent term
            f_cent = ((1.0 - alpha) * numpy.exp(-temps / ts3) +
                      alpha * numpy.exp(-temps / ts1))
            if ts2 is not None and ts2 is not numpy.nan:  # catches nan & None
                f_cent += numpy.exp(-ts2 / temps)

            # Calculate the Log F term
            c_val = -0.4 - 0.67 * numpy.log10(f_cent)
            n_val = 0.75 - 1.27 * numpy.log10(f_cent)
            d_val = 0.14
            val = ((numpy.log10(pr_term) + c_val) /
                   (n_val - d_val * (numpy.log10(pr_term) + c_val)))**2
            logf = (1.0 + val)**(-1) * numpy.log10(f_cent)

            # Calculate F broadening term
            f_term = 10**(logf)

            return f_term

        # Calculate the high- and low-P rate constants
        highp_kts = arr(highp_arr, temps, tref=tref)
        lowp_kts = arr(lowp_arr, temps, tref=tref)

        # If the pressure is 'high', just use the high-P kts
        if pressure == 'high':
            kts = highp_kts
        # Otherwise, calculate the pressure-dependent kts
        else:
            # Calculate the pr term and broadening factor
            pr_term = _pr_term(highp_kts, lowp_kts, temps, pressure, collid_factor)
            f_term = _f_broadening_term(pr_term, troe_params, temps)

            # Calculate Troe rate constants
            kts = highp_kts * (pr_term / (1.0 + pr_term)) * f_term

        return kts

    ktp_dct = {}
    for pidx, pressure in enumerate(pressures):
        temps = temps_lst[pidx]
        kts = troe_one_p(highp_arr, lowp_arr, troe_params, temps, pressure,
                         collid_factor=collid_factor, tref=tref)
        ktp_dct[pressure] = (temps, kts)

    return ktp_dct


def lind(highp_arr, lowp_arr, temps_lst, pressures, collid_factor=1.0,
         tref=1.0):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a Lindemann functional expression.

        :param highp_arr: high-P limit Arrhenius parameters
        :type highp_arr: list [[A1, n1, Ea1], [A2, n2, Ea2,], ...]
        :param lowp_arr: low-P limit Arrhenius parameters
        :type lowp_arr: list [[A1, n1, Ea1], [A2, n2, Ea2,], ...]
        :param temps_lst: list of temperature arrays used to get k(T,P)s (K)
        :type temps_lst: list [numpy.ndarray1, numpy.ndarray2, ...]
        :param pressures: pressures used to get k(T,P)s (atm)
        :type pressure: list
        :param collid_factor: collider efficiency factor
        :type collid_factor: float
        :param tref: reference temperature used for modified Arrhenius (K)
        :type tref: float
        :return ktp_dct: k(T,Ps) at all temps and pressures
        :rtype: dict {pressure: (temps, kts)}
    """

    def lind_one_p(highp_arr, lowp_arr, temps, pressure, collid_factor=1.0,
                   tref=1.0):
        """ Calculates T,P-dependent rate constants [k(T,P)]s using
            a Lindemann functional expression, at a single pressure,
            across several temperatures.

            :param temps: temperature array used to get k(T)s (K)
            :type temps: numpy.ndarray
            :return kts: k(T)s at the given pressure
            :rtype: numpy.ndarray
            (other inputs unchanged from parent function)
        """

        # Calculate the high- and low-P rate constants
        highp_kts = arr(highp_arr, temps, tref=tref)
        lowp_kts = arr(lowp_arr, temps, tref=tref)

        # If the pressure is 'high', just use the high-P kts
        if pressure == 'high':
            kts = highp_kts
        # Otherwise, calculate the pressure-dependent kts
        else:
            # Calculate the pr term
            pr_term = _pr_term(highp_kts, lowp_kts, temps, pressure, collid_factor)

            # Calculate Lindemann rate constants
            kts = highp_kts * (pr_term / (1.0 + pr_term))

        return kts

    ktp_dct = {}
    for pidx, pressure in enumerate(pressures):
        temps = temps_lst[pidx]
        kts = lind_one_p(highp_arr, lowp_arr, temps, pressure,
                         collid_factor=collid_factor, tref=tref)
        ktp_dct[pressure] = (temps, kts)

    return ktp_dct


def merge_rxn_ktp_dcts(full_rxn_ktp_dct, rxn_ktp_dct):
    """ Merge a reaction ktp dictionary into an existing one.
        If a reaction currently exists in both, then we add the kts
        together for the reactions.
        If a reaction is a self-reaction A+B=A+B (and same bath gas):
        the reaction is not added
    """
    # Add dictionary to another ktp dictionary
    # We will add the rates together if rxn is prev. found
    for rxn, _ktp1 in rxn_ktp_dct.items():
        # if reactants and products are the same: do not add
        rcts = list(rxn[0])
        prds = list(rxn[1])
        rcts.sort(), prds.sort()
        if rcts == prds:
            print('self reaction {} skipped'.format(rxn))
            continue
        
        if rxn not in full_rxn_ktp_dct:
            # Simply place rates into the full dct
            full_rxn_ktp_dct[rxn] = _ktp1
        else:
            # Add the existing rates from the full and small dct together
            _ktp2 = full_rxn_ktp_dct[rxn]
            full_rxn_ktp_dct[rxn] = add_ktp_dcts(_ktp1, _ktp2)

    return full_rxn_ktp_dct


def add_ktp_dcts(ktp_dct1, ktp_dct2):
    """ Add the rates in two ktp_dcts
    """

    # If either starting dct is empty, simply copy the other
    if not ktp_dct1:
        added_dct = copy.deepcopy(ktp_dct2)
    elif not ktp_dct2:
        added_dct = copy.deepcopy(ktp_dct1)
    # Otherwise, add the dcts
    else:
        added_dct = {}

        # Determine pressures common to both dcts and only in p1/p2
        pr1 = set(ktp_dct1.keys())
        pr2 = set(ktp_dct2.keys())

        pr_common = sorted(list(pr1 & pr2))
        pr1_only = sorted(list(pr1 - pr2))
        pr2_only = sorted(list(pr2 - pr1))

        # Loop over common pressures and add the rate constants togther
        for pressure in pr_common:
            temps1, kts1 = ktp_dct1[pressure]
            temps2, kts2 = ktp_dct2[pressure]

            # Check if two arrays are the same
            if numpy.array_equal(temps1, temps2):
                added_kts = kts1 + kts2
                added_dct[pressure] = (temps1, added_kts)  # store added values
            else:
                added_temps = sorted(
                    list(set(temps1.tolist() + temps2.tolist())))
                added_kts = []
                for temp in added_temps:
                    if temp in temps1:
                        _kt1_idx = next(i for i, t in enumerate(temps1)
                                        if numpy.isclose(temp, t))
                        _kt1 = kts1[_kt1_idx]
                    else:
                        _kt1 = 0
                    if temp in temps2:
                        _kt2_idx = next(i for i, t in enumerate(temps2)
                                        if numpy.isclose(temp, t))
                        _kt2 = kts2[_kt2_idx]
                    else:
                        _kt2 = 0
                    added_kts.append(_kt1+_kt2)
                added_dct[pressure] = (numpy.array(added_temps),
                                       numpy.array(added_kts))

        # Add rate constants for pressures only in ktp_dct1
        for pressure in pr1_only:
            added_dct[pressure] = ktp_dct1[pressure]

        # Add rate constants for pressures only in ktp_dct2
        for pressure in pr2_only:
            added_dct[pressure] = ktp_dct2[pressure]

    return added_dct


def check_p_t(temps_lst, pressures):
    """ Enforces rules on the temps_lst and pressures. In the case where the
        user has only provided one temp array in temps_lst, duplicates temp
        to occur the same number of times as the number of pressures

        :param temps_lst: list of temperature arrays used to get k(T,P)s (K)
        :type temps_lst: list [numpy.ndarray1, numpy.ndarray2, ...]
        :param pressures: pressures used to get k(T,P)s (atm)
        :type pressure: list
    """

    # Accepted types for temps_lst and pressures
    good = (list, tuple)

    # Check that temps_lst is a list or tuple of Numpy arrays
    assert isinstance(temps_lst, good) and isinstance(temps_lst[0], numpy.ndarray), (
        'The temps_lst input should be a list or tuple of Numpy.ndarrays.')

    # Check that pressures is a list or tuple
    assert isinstance(pressures, good), (
        f'The pressures input should be a list or tuple, not {type(pressures)}')

    # In the case where the user only provided one temp array, duplicate it to
    # make the number of temp arrays match the number of pressures
    if len(temps_lst) == 1:
        new_temps_lst = []
        for _ in range(len(pressures)):
            new_temps_lst.append(temps_lst[0])
    # Otherwise, make sure the number of temp_arrays is the number of pressures
    # and that the temp arrays are all Numpy arrays
    else:
        assert len(temps_lst) == len(pressures), (
            f'There are {len(temps_lst)} temp_arrays and {len(pressures)}'
            ' pressures. These should be the same.')
        for temps in temps_lst:
            assert isinstance(temps, numpy.ndarray), (
                'The temps in temps_lst should be Numpy arrays')
        new_temps_lst = copy.copy(temps_lst)  # temps_lst is unchanged

    return new_temps_lst


def checks_temp_pressure_and_extend(ped_df, hotbf_df):
    """ compare T,P """
    temp_ped, pressure_ped = [ped_df.index, ped_df.columns]
    temp_hot, pressure_hot = [hotbf_df.index, hotbf_df.columns]
    # check that T of temp_hot are at least as many as those of temp_ped
    # if they're not, it doesn't make sense to continue
    if not all(temp_ped_i in temp_hot for temp_ped_i in temp_ped):
        print('*Error: temperature range in HOTenergies '
              'does not cover the full range')
        sys.exit()

    # if in pressure_ped not all pressures of hoten are available:
    # extend the range of pressures
    # ex.: for H abstractions, they will be pressure independent
    # but probably the successive decomposition is not
    if not all(pressure_hot_i in pressure_ped
               for pressure_hot_i in pressure_hot):
        print('*Warning: P range of PedOutput smaller than HOTenergies:\n')
        print('Energy distribution at other pressure '
              'approximated from available values \n')
        ped_df = extend_dct_with_pressure(
            pressure_hot, pressure_ped, ped_df)

    # check that pressures of pressure_ped are contained in hot energies
    # if they are not: extend pressure range assuming behavior is the same
    # ELIF no good here - might be that pressure_ped and pressure_hot
    #   only intersect in a limited range
    if not all(pressure_ped_i in pressure_hot
               for pressure_ped_i in pressure_ped):
        print('Warning: P range of HOTenergies smaller than PEDoutput:')
        print('Energy distribution at other pressures '
              'approximated from available values \n')
        hotbf_df = extend_dct_with_pressure(
            pressure_ped, pressure_hot, hotbf_df)

    return ped_df, hotbf_df


def extend_dct_with_pressure(pressure_all, pressure_red, dct_toextend):
    """ Takes a dictionary(dataframe) with pressure keys(columns) and extends it:
        approximate the values at missing pressures with available values

        :param pressure_all: all pressures
        :type pressure_all: list
        :param dct_toextend: dictionary (dataframe) w/ pressure_reduced as keys
        :type dct_toextend: dictionary (dataframe)
        :return dct_extended: dictionary with pressure_all as keys (columns)
        :rtype: dictionary (dataframe)
    """

    for pressure in pressure_all:
        if pressure not in pressure_red:
            # approximate pressure:
            # provides the minimum difference with pressure_ped_i
            pressure_approx = pressure_red[
                numpy.argmin(
                    [numpy.log(abs(pressure_red_i-pressure))
                     for pressure_red_i in pressure_red])]
            # extend original dataframe
            dct_toextend[pressure] = dct_toextend[pressure_approx]

    return dct_toextend


def p_to_m(pressure, temps, rval=RC2):
    """ Convert the pressure to the concentration of a gas, [M], assuming an
        ideal gas form where [M] ~ P/RT.

        :param pressure: pressure of gas (atm)
        :type pressure: float
        :return mconc: concentration of gas (mol/cm^3)
        :rtype: float
    """
    return pressure / (rval * temps)


def _pr_term(highp_kts, lowp_kts, temps, pressure, collid_factor=1.0, rval=RC2):
    """ Calculates the reduced pressure term for a single pressure
        used for Lindemann and Troe P-dependent functional expressions.

        :param highp_kts: k(T)s determined at the high-pressure limit
        :type highp_kts: numpy.ndarray
        :param lowp_kts: k(T)s determined at the low-pressure limit
        :type lowp_kts: numpy.ndarray
        :param temps: temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressure: pressure used to calculate reduced pressure
        :type pressure: float
        :rtype: numpy.ndarray
    """

    pr_term = (
        (lowp_kts / highp_kts) *
        p_to_m(pressure, temps, rval=rval) *
        collid_factor)

    return pr_term


def read_rxn_ktp_dct(rxn_ktp_dct, rxn, pressure, val):
    """ Reads the entries of a rxn_ktp_dct for a single rxn and single pressure

        :param rxn_ktp_dct: ktp_dcts indexed by rxn name
        :type rxn_ktp_dct: dict {rxn: ktp_dct}
        :param rxn: rxn key
        :type rxn: tuple
        :param pressure: pressure at which to get values (atm); can be 'high'
        :type pressure: float or str
        :param val: the quantity to return; either 'temps' or 'rates'
        :type val: str
        :return: the desired quantity
        :rtype: Numpy.ndarray
    """

    assert pressure == 'high' or isinstance(pressure, float)
    assert val in ('temps', 'rates'), (
        f"The val input should be 'temps' or 'rates', not '{val}'")

    # Deal with pressure floats
    if pressure == 'high':
        pval = 'high'
    else:
        for _pressure in (p for p in rxn_ktp_dct[rxn].keys() if p != 'high'):
            if numpy.isclose(pressure, _pressure, atol=1.0e-4):
                pval = _pressure

    # Set final val index for reading rates versus temp array
    val_idx = 1 if val == 'rates' else 0

    return rxn_ktp_dct[rxn][pval][val_idx]


def remove_high(temps_lst, pressures):
    """ Removes the 'high' entry from a list of pressures and removes the
        corresponding temps array from the temps_lst. If 'high' is not
        present, will return temps_lst and pressures unchanged.

        :param temps_lst: list of temperature arrays used to get k(T,P)s (K)
        :type temps_lst: list [numpy.ndarray1, numpy.ndarray2, ...]
        :param pressures: pressures used to get k(T,P)s (atm)
        :type pressure: list
        :return fixed_temps_lst: list of temperature arrays, with the array
            corresponding to 'high' removed (if it was present)
        :rtype: list [numpy.ndarray1, numpy.ndarray2, ...]
        :return fixed_pressures: list of pressures, with 'high' removed (if it
            was present)
        :rtype: list
    """

    fixed_temps_lst = copy.deepcopy(temps_lst)
    fixed_pressures = copy.deepcopy(pressures)

    if 'high' in fixed_pressures:
        high_idx = fixed_pressures.index('high')
        fixed_pressures.remove('high')
        fixed_temps_lst.pop(high_idx)

    return fixed_temps_lst, fixed_pressures


def mult_by_factor(ktp_dct, factor):
    """ Multiplies all k(T,P) values in a ktp_dct by a fixed constant
    """

    new_ktp_dct = {}
    for pressure, (temps, kts) in ktp_dct.items():
        new_ktp_dct[pressure] = (temps, kts * factor)

    return new_ktp_dct


def get_bw_ktp_dct(ktp_dct_entry, rcts, prds, therm_df):
    """get backward rate constants for entries in a ktp dct

    Args:
        ktp_dct_entry: usual ktp dct entry
        rcts (tuple): reactants
        prds (tuple): products
        therm_df (dataframe): dataframe with thermochemistry
    returns: ktp_dct_entry_bw
    """
    ktp_dct_entry_bw = {}
    for P, vals in ktp_dct_entry.items():
        Tvect, k = vals
        kt_series = pandas.Series(k, index=Tvect)
        dg_rxn = thermo.extract_deltaX_therm(therm_df, rcts, prds, 'G')
        kt_series_bw = get_bw_rate(kt_series, rcts, prds, dg_rxn)
        ktp_dct_entry_bw[P] = (Tvect, kt_series_bw.values)

    return ktp_dct_entry_bw


def get_bw_rate(kt_series, rcts, prds, dg_rxn):
    """get backward rate constat

    Args:
        kt_series (pd series): k values. index is temperature
        rcts (tuple): reactants
        prds (tuple): products
        dg_rxn (pd series): DG(T)
    """
    deltanu = len(prds) - len(rcts)
    # 1e+5 Pa / 1e+6 to get cm3
    kt_series *= numpy.power((1e+5/8.314/kt_series.index/numpy.power(10, 6)),-deltanu)
    kt_series *= numpy.exp(dg_rxn[kt_series.index]/1.987/kt_series.index)
    
    return kt_series
           
    