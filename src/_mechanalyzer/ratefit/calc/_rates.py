"""
Calculate rates with various fitting functions
"""

import numpy as np
from scipy.special import eval_chebyt
from phydat import phycon


RC = phycon.RC_CAL  # gas constant in cal/(mol.K)
RC2 = phycon.RC_ATM  # gas constant in cm^3.atm/(mol.K)


def single_arrhenius(a_par, n_par, ea_par,
                     t_ref, temps, rval=RC):
    """ Calculates T-dependent rate constants [k(T)]s using
        a single Arrhenius functional expression.

        :param a_par: pre-exponential A parameter
        :type a_par: float
        :param n_par: temperature exponent n parameter
        :type n_par: float
        :param ea_par: activation energy Ea parameter (kcal.mol-1)
        :type ea_par: float
        :param t_ref: Reference temperature (K)
         :type t_ref: float
        :param temps: List of Temperatures (K)
        :type temps: numpy.ndarray
        :return kts: T-dependent rate constants
        :rtype: numpy.ndarray
    """
    kts = a_par * ((temps / t_ref)**n_par) * np.exp(-ea_par/(rval*temps))
    return kts


def double_arrhenius(a_par1, n_par1, ea_par1,
                     a_par2, n_par2, ea_par2,
                     t_ref, temps, rval=RC):
    """ Calculates T-dependent rate constants [k(T)]s using
        a double Arrhenius functional expression.

        :param a_par1: 1st pre-exponential A parameter
        :type a_par1: float
        :param n_par1: 1st temperature exponent n parameter
        :type n_par1: float
        :param ea_par1: 1st activation energy Ea parameter (kcal.mol-1)
        :type ea_par1: float
        :param a_par2: 2nd pre-exponential A parameter
        :type a_par2: float
        :param n_par2: 2nd temperature exponent n parameter
        :type n_par2: float
        :param ea_par2: 2nd activation energy Ea parameter (kcal.mol-1)
        :type ea_par2: float
        :param t_ref: Reference temperature (K)
        :type t_ref: float
        :param temps: List of Temperatures (K)
        :type temps: numpy.ndarray
        :return kts: T-dependent rate constants
        :rtype: numpy.ndarray
    """
    kts = (
        a_par1 * ((temps / t_ref)**n_par1) * np.exp(-ea_par1/(rval*temps)) +
        a_par2 * ((temps / t_ref)**n_par2) * np.exp(-ea_par2/(rval*temps))
    )
    return kts


def arrhenius(arr_tuples, temps, t_ref, rval=RC):
    """ Calculates T-dependent rate constants [k(T)]s using a list of any number
            of Arrhenius parameters

        :param arr_tuples: Arrhenius fit parameters
        :type arr_tuples: tuple ((A1, n1, Ea1), (A2, n2, Ea2), ...)
        :param temps: temperatures at which to do calculations (K)
        :type temps: numpy.ndarray
        :param t_ref: reference temperature (K)
        :type t_ref: float
        :return kts: T-dependent rate constants
        :rtype: numpy.ndarray
    """

    kts = np.zeros_like(temps)
    for arr_tuple in arr_tuples:
        assert len(arr_tuple) == 3, (
            f'Length of Arrhenius tuples should be 3, not {len(arr_tuple)}')
        a_par, n_par, ea_par = arr_tuple
        new_kts = a_par * ((temps/t_ref)**n_par) * np.exp(-ea_par/(rval*temps))
        kts = np.add(kts, new_kts, out=kts, casting='unsafe')  # stops Np errs

    return kts


def lowp_limit(highp_kts, temps, pressures, collid_factor=1.0, rval=RC2):
    """ Calculates T,P-dependent rate constants [k(T,P)]s assuming
        the reaction occurs in the low-pressure regime where the
        rates are linear with pressure.

        :param highp_kts: k(T)s determined at high-pressure
        :type highp_kts: numpy.ndarray
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressures: Pressures used to calculate k(T,P)s
        :type pressures: list(float)
        :param collid_factor: Buffer enhancement collision factor
        :type collid_factor: float
        :return ktp_dct: k(T,Ps) at all temps and pressures
        :rtype: dict[pressure: temps]
    """
    kp_dct = {}
    for pressure in pressures:
        kp_dct[pressure] = lowp_limit_one_pressure(
            highp_kts, temps, pressure,
            collid_factor=collid_factor, rval=rval)

    ktp_dct = _ktp_dct(kp_dct, temps)

    return ktp_dct


def lowp_limit_one_pressure(highp_kts, temps, pressure,
                            collid_factor=1.0, rval=RC2):
    """ Calculates the reduced pressure term for a single pressure
        used for Lindemann and Troe P-dependent functional expressions.

        :param list highp_kts: k(T)s determined at high-pressure
        :type highp_kts: numpy.ndarray
        :param temps: Temps used to calculate high- and low-k(T)s
        :temps: numpy.ndarray
        :param pressure: Pressure used to calculate reduced pressure
        :type pressure: float
        :param collid_factor: Buffer enhancement collision factor
        :type collid_factor: float
        :rtype: numpy.ndarray
    """
    return highp_kts * p_to_m(pressure, temps, rval=rval) * collid_factor


def lindemann(highp_kts, lowp_kts, temps, pressures, collid_factor=1.0):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a Lindemann functional expression.

        :param highp_kts: k(T)s determined at high-pressure
        :type highp_kts: numpy.ndarray
        :param lowp_kts: k(T)s determined at low-pressure
        :type lowp_kts: numpy.ndarray
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressures: Pressures used to calculate k(T,P)s
        :type pressures: list(float)
        :param collid_factor: Buffer enhancement collision factor
        :type collid_factor: float
        :return ktp_dct: k(T,Ps) at all temps and pressures
        :rtype: dict[pressure: temps]
    """
    kp_dct = {}
    for pressure in pressures:
        kp_dct[pressure] = lindemann_one_pressure(
            highp_kts, lowp_kts, temps, pressure, collid_factor=collid_factor)

    ktp_dct = _ktp_dct(kp_dct, temps, highp_kts=highp_kts)

    return ktp_dct


def lindemann_one_pressure(highp_kts, lowp_kts, temps, pressure,
                           collid_factor=1.0):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a Lindemann functional expression, at a given pressure,
        across several temperatures.

        :param highp_kts: k(T)s determined at high-pressure
        :type highp_kts: numpy.ndarray
        :param lowp_kts: k(T)s determined at low-pressure
        :type lowp_kts: numpy.ndarray
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressure: Pressure used to calculate k(T,P)s
        :type pressure: float
        :param collid_factor: Buffer enhancement collision factor
        :type collid_factor: float
        :return ktps: Set of k(T,P)s at given pressure
        :rtype numpy.ndarray
    """
    # Calculate the pr term
    pr_terms = _pr_term(highp_kts, lowp_kts, temps, pressure,
                        collid_factor=collid_factor)

    # Calculate Lindemann rate constants
    ktps = highp_kts * (pr_terms / (1.0 + pr_terms))

    return ktps


def troe(highp_kts, lowp_kts, temps, pressures,
         alpha, ts3, ts1, ts2=None, collid_factor=1.0):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a Troe functional expression.

        :param highp_kts: k(T)s determined at high-pressure
        :type highp_kts: numpy.ndarray
        :param lowp_kts: k(T)s determined at low-pressure
        :type lowp_kts: numpy.ndarray
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressures: Pressures used to calculate k(T,P)s
        :type pressure: list(float)
        :param alpha: Troe alpha parameter
        :type alpha: float
        :param ts3: Troe T3 parameter
        :type ts3: float
        :param ts1: Troe T1 parameter
        :type ts1: float
        :param ts2: Troe T2 parameter
        :type ts2: float
        :param collid_factor: Buffer enhancement collision factor
        :type collid_factor: float
        :return ktp_dct: k(T,Ps) at all temps and pressures
        :rtype: dict[pressure: temps]
    """
    kp_dct = {}
    for pressure in pressures:
        kp_dct[pressure] = troe_one_pressure(
            highp_kts, lowp_kts, temps, pressure,
            alpha, ts3, ts1, ts2=ts2, collid_factor=collid_factor)

    ktp_dct = _ktp_dct(kp_dct, temps, highp_kts=highp_kts)

    return ktp_dct


def troe_one_pressure(highp_ks, lowp_ks, temps, pressure,
                      alpha, ts3, ts1, ts2=None, collid_factor=1.0):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a Troe functional expression, at a given pressure,
        across several temperatures.

        :param highp_ks: k(T)s determined at high-pressure
        :type highp_ks: numpy.ndarray
        :param lowp_ks: k(T)s determined at low-pressure
        :type lowp_ks: numpy.ndarray
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressures: Pressure used to calculate k(T,P)s
        :type pressure: float
        :param alpha: Troe alpha parameter
        :type alpha: float
        :param ts3: Troe T3 parameter
        :type ts3: float
        :param ts1: Troe T1 parameter
        :type ts1: float
        :param ts2: Troe T2 parameter
        :type ts2: float
        :param collid_factor: Buffer enhancement collision factor
        :type collid_factor: float
        :return ktps: Set of k(T,P)s at given pressure
        :rtype numpy.ndarray
        :return ktps: Set of k(T,P)s at given pressure
        :rtype: dict[pressure: temps]
    """

    # Calculate the pr term and broadening factor
    pr_term = _pr_term(highp_ks, lowp_ks, temps, pressure, collid_factor)
    f_term = _f_broadening_term(pr_term, alpha, ts3, ts1, ts2, temps)

    # Calculate Troe rate constants (collision factor could be wrong)
    ktps = highp_ks * (pr_term / (1.0 + pr_term)) * f_term

    return ktps


def plog(plog_dct, t_ref, temps, pressures):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a PLOG functional expression.

        :param plog_dct: Arrhenius fitting parameters at several pressures
        :type plog_dct: dict[pressure: [fit_params]]
        :param t_ref: Reference temperature (K)
        :type t_ref: float
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressures: Pressures used to calculate k(T,P)s
        :type pressures: list(float)
        :return ktp_dct: k(T,Ps) at all temps and pressures
        :rtype: dict[pressure: temps]
    """

    # Set the plog pressures to see if pressure is in range
    plog_pressures = list(plog_dct.keys())

    kp_dct = {}
    for pressure in pressures:
        if min(plog_pressures) <= pressure <= max(plog_pressures):
            kp_dct[pressure] = plog_one_pressure(
                plog_dct, t_ref, temps, pressure)
            # From mechanalyzer, might need fix
            # if np.ndim(temps) == 1:
            #     kp_dct[pressure] = plog_one_pressure(
            #         plog_dct, temps, pressure, t_ref)
            # else:
            #     kp_dct[pressure] = plog_one_pressure(
            #         plog_dct, temps[index], pressure, t_ref)

    ktp_dct = _ktp_dct(kp_dct, temps)

    return ktp_dct


def plog_one_pressure(plog_dct, t_ref, temps, pressure):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a PLOG functional expression, at a given pressure,
        across several temperatures.

        :param plog_dct: Arrhenius fitting parameters at several pressures
        :type plog_dct: dict[pressure: [fit_params]]
        :param t_ref: Reference temperature (K)
        :type t_ref: float
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressure: Pressure used to calculate k(T,P)s
        :type pressure: float
        :return ktps: Set of k(T,P)s at given pressure
        :rtype numpy.ndarray
    """

    plog_pressures = list(plog_dct.keys())

    # Check if pressure is in plog dct; use plog pressure for numerical stab
    pressure_defined = False
    for plog_pressure in plog_pressures:
        if np.isclose(pressure, plog_pressure, atol=1.0e-3):
            pressure_defined = True
            plog_params = plog_dct[plog_pressure]

    # If pressure equals value use, arrhenius expression
    if pressure_defined:
        ktps = arrhenius(plog_params, t_ref, temps)
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
        pres_term = (
            (np.log10(pressure) - np.log10(plow)) /
            (np.log10(phigh) - np.log10(plow))
        )

        # Calculate k(T)s at high-P and low-P with Arrhenius expressions
        kt_low = arrhenius(plow_params, t_ref, temps)
        kt_high = arrhenius(phigh_params, t_ref, temps)

        # Calculate K(T,P)s with PLOG expression
        logkt = (
            np.log10(kt_low) +
            ((np.log10(kt_high) - np.log10(kt_low)) * pres_term)
        )
        ktps = 10**(logkt)

    return ktps


def cheb(alpha, tlim, plim, temps, pressures):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a Chebyshev functional expression.

        :param alpha: Chebyshev coefficient matrix
        :type alpha: numpy.ndarray
        :param tmin: minimum temperature Chebyshev model is defined
        :type tmin: float
        :param pmin: minimum pressure Chebyshev model is defined
        :type pmin: float
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressures: Pressures used to calculate k(T,P)s
        :type pressures: list(float)
        :return ktp_dct: k(T,Ps) at all temps and pressures
        :rtype: dict[pressure: temps]
    """

    def cheb_one_pressure(alpha, tlim, plim, temps, pressure):
        """ Calculates T,P-dependent rate constants [k(T,P)]s using
            a Chebyshev functional expression, at a given pressure,
            across several temperatures.
        """

        tmin, tmax = tlim
        pmin, pmax = plim
        alpha_nrows, alpha_ncols = alpha.shape

        ktps = np.zeros(len(temps))
        for i, temp in enumerate(temps):
            ctemp = (
                (2.0 * temp**(-1) - tmin**(-1) - tmax**(-1)) /
                (tmax**(-1) - tmin**(-1)))
            cpress = (
                (2.0 * np.log10(pressure) - np.log10(pmin) - np.log10(pmax)) /
                (np.log10(pmax) - np.log10(pmin)))

            logktp = 0.0
            for j in range(alpha_nrows):
                for k in range(alpha_ncols):
                    logktp += (alpha[j][k] * eval_chebyt(j, ctemp) *
                               eval_chebyt(k, cpress))

            ktps[i] = 10**(logktp)

        return ktps

    kp_dct = {}
    for pressure in pressures:
        kp_dct[pressure] = cheb_one_pressure(
            alpha, tlim, plim, temps, pressure)

    ktp_dct = _ktp_dct(kp_dct, temps)

    return ktp_dct


# Functions for calculating terms in certain P-dependent expressions
def _pr_term(highp_rateks, lowp_rateks, temps, pressure,
             collid_factor=1.0, rval=RC2):
    """ Calculates the reduced pressure term for a single pressure
        used for Lindemann and Troe P-dependent functional expressions.

        :param list highp_ks: k(T)s determined at high-pressure
        :type highp_ks: numpy.ndarray
        :param list lowp_ks: k(T)s determined at low-pressure
        :type lowp_ks: numpy.ndarray
        :param temps: Temps used to calculate high- and low-k(T)s
        :temps: numpy.ndarray
        :param pressure: Pressure used to calculate reduced pressure
        :type pressure: float
        :rtype: numpy.ndarray
    """

    pr_term = (
        (lowp_rateks / highp_rateks) *
        p_to_m(pressure, temps, rval=rval) *
        collid_factor
    )

    return pr_term


def _f_broadening_term(pr_term, alpha, ts3, ts1, ts2, temp):
    """ Calculates the F broadening factor term used for Troe expressions.

        :param pr_term: reduced pressure
        :type pr_term: float
        :param alpha: Troe alpha parameter
        :type alpha: float
        :param ts3: Troe T3 parameter
        :type ts3: float
        :param ts1: Troe T1 parameter
        :type ts1: float
        :param ts2: Troe T2 parameter
        :type ts2: float
        :param temp: temperature
        :type temp: float
        :return f_cent: F broadening term
        :rtype: float
    """

    # Calculate Fcent term
    f_cent = ((1.0 - alpha) * np.exp(-temp / ts3) +
              alpha * np.exp(-temp / ts1))
    if ts2 is not None:
        f_cent += np.exp(-ts2 / temp)

    # Calculate the Log F term
    c_val = -0.4 - 0.67 * np.log10(f_cent)
    n_val = 0.75 - 1.27 * np.log10(f_cent)
    d_val = 0.14
    val = ((np.log10(pr_term) + c_val) /
           (n_val - d_val * (np.log10(pr_term) + c_val)))**2
    logf = (1.0 + val)**(-1) * np.log10(f_cent)

    # Calculate F broadening term
    f_term = 10**(logf)

    return f_term


def p_to_m(pressure, temps, rval=RC2):
    """ Convert the pressure to the concentration of a gas [M]
        assuming an ideal gas form where [M] ~ P/RT.

        :param pressure: Pressure of gas (in atm)
        :type pressure: float
        :return mconc: Conncentration of Gas (mol/cm^3)
        :rtype: float
    """
    return pressure / (rval * temps)


# Helper functions
def _ktp_dct(kp_dct, temps, highp_kts=None):
    """ Creates a kTP dictionary from a kP dictionary and an
        array of temperatures

        :param kp_dct: dct of the form {P:[k@T1, k@T2]}
        :type kp_dct: {float: numpy.array}
        :param temps: array of temperatures, either 1- or 2-D
        :type temps: numpy.ndarray
        :return ktp_dct:
        :rtype:

    """

    ktp_dct = {}
    for index, pressure in enumerate(kp_dct):
        kts = kp_dct[pressure]
        if np.ndim(temps) == 1:
            ktp_dct[pressure] = (temps, kts)
        else:
            ktp_dct[pressure] = (temps[index], kts)

    # Add the high-P k(T)s to the kTP dictionary if needed
    if highp_kts is not None:
        if np.ndim(temps) == 1:
            ktp_dct['high'] = (temps, highp_kts)
        else:
            ktp_dct['high'] = (temps[-1], highp_kts)

    return ktp_dct
