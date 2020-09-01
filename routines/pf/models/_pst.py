"""
  Calculate the pre-exponential parameter required to calculate some desired
  Phase-Space Theory rate constant, given an exponential parameter and
  temperature. Using rate expression constant gives functional
  dependence  of k(n, mu, C0, T).

  For (n=N, mu=MU, T=T), the target parameter (C0_TGT) needed to obtain
  the target rate constant (k_TGT) can be found via

    C0_TGT = [ k_TGT / k(n=N, mu=MU, C0=1.0, T=300.0) ] * 1.0
"""

import scipy
import numpy
import automol


def kt_pst(n_par, mred, cn_par, temp):

    """ Calculate a rate constant according to Phase-Space Theory.

        :param n_par: exponential parameter
        :param mred: reduced mass ()
        :param cn_par: pre-exponential potential coefficient [in Bohr]
        :param temp: temperature (K)
        :return: k(T)
        :rtype: float
    """

    kt_val = (
        (8.0 * numpy.pi)**(1.0/2.0) *
        ((n_par - 2) / 2)**(2.0/n_par) *
        scipy.special.gamma(1.0 - 2.0/n_par) *
        mred**(-1.0/2.0) *
        cn_par**(2.0/n_par) *
        temp**(1.0/2.0 - 2.0/n_par)
    )

    return kt_val


if __name__ == '__main__':

    # Input variables read from filesys
    GEO1 = 'geo'
    GEO2 = 'geo'

    # input variables from theuser
    KT_PST = 1e-12
    N_PST = 6.0
    T_PST = 300.0

    # Calculate reduced mass
    MRED = automol.geom.reduced_mass(GEO1, GEO2)

    # Calculate the pre-exponential needed to get specific PST k(T)
    C0 = KT_PST / kt_pst(N_PST, MRED, 1.0, T_PST)
