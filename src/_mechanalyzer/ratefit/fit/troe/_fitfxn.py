""" fit rate constants to Arrhenius expressions
"""

import os
from ratefit.fit.troe import troefit_io
from ratefit.fit import invert_ktp_dct
from ratefit.fit import pull_highp_from_dct


RC = 1.98720425864083e-3  # Gas Constant in kcal/mol.K


def reaction(ktp_dct, troefit_path,
             troe_param_fit_lst=('ts1', 'ts2', 'ts3', 'alpha'),
             highp_guess=(8.1e-11, -0.01, 1000.0),
             lowp_guess=(8.1e-11, -0.01, 1000.0),
             alpha=0.19, ts1=590.0, ts2=1.0e6, ts3=6.0e4,
             fit_tol1=1.0e-8, fit_tol2=1.0e-8,
             a_conv_factor=1.0):
    """ Fits a set of T,P-dependent rate constants [k(T,P)]s to a
        Troe functional expression using troefit fitting code.

        :param ktp_dct: k(T,P) values
        :type ktp_dct: dict[pressure:temps]
        :param troe_param_fit_lst: Troe parameters to include in fitting
        :param troe_param_fit_lst: list(str)
        :param highp_guess: A, n, Ea Arrhenius params for high-P limit
        :type higp_guess: list(float)
        :param lowp_guess: A, n, Ea Arrhenius params for low-P limit
        :type higp_guess: list(float)
        :param alpha: Troe alpha parameter
        :type alpha: float
        :param ts3: Troe T3 parameter
        :type ts3: float
        :param ts1: Troe T1 parameter
        :type ts1: float
        :param ts2: Troe T2 parameter
        :type ts2: float
        :param fit_tol1: tolerance for fitting ???
        :type fit_tol1: float
        :param fit_tol2: tolerance for fitting ???
        :type fit_tol2: float
        :param troe_path: path to run troe
        :type troe_path: str
        :param a_conv_factor: Conversion factor for A parameter
        :type a_conv_factor: float
        :return fit_params: fitting parameters for function
        :rtype: list
    """

    # Pull high pressures out
    _, num_ktp_dct, _ = pull_highp_from_dct(ktp_dct)

    # Invert the k(T,P) dct to k(P,T) dct
    kpt_dct = invert_ktp_dct(num_ktp_dct)

    # Write the input file for the ratefit code
    ratefit_inp_str = troefit_io.write_input(
        kpt_dct,
        troe_param_fit_lst=troe_param_fit_lst,
        highp_guess=highp_guess,
        lowp_guess=lowp_guess,
        alpha=alpha, ts1=ts1, ts2=ts2, ts3=ts3,
        fit_tol1=fit_tol1, fit_tol2=fit_tol2)
    troefit_inp_file = os.path.join(troefit_path, 'troefit.dat')
    with open(troefit_inp_file, 'w', encoding='utf-8') as troefit_infile:
        troefit_infile.write(ratefit_inp_str)

    # Run the ratefit program
    troefit_io.run_troefit(troefit_path)

    # Parse the ratefit files for the Arrhenius fit parameters
    troename = 'troefit.out'
    troefit_out_file = os.path.join(troefit_path, troename)
    with open(troefit_out_file, 'r', encoding='utf-8') as troefit_outfile:
        troefit_out_str = troefit_outfile.read()

    fit_params = troefit_io.read_params(troefit_out_str, a_conv_factor)

    return fit_params
