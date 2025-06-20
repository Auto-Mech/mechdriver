""" interface to troefit code by sjk
"""

import os
import subprocess
from mako.template import Template


RC = 1.98720425864083e-3  # Gas Constant in kcal/mol.K
RC2 = 0.0820573660809596 * 1000.0  # Gas Constant in cm^3.atm/mol.K

# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))


def write_input(kpt_dct,
                troe_param_fit_lst=('ts1', 'ts2', 'ts3', 'alpha'),
                highp_guess=(8.1e-11, -0.01, 1000.0),
                lowp_guess=(8.1e-11, -0.01, 1000.0),
                alpha=0.19, ts1=590.0, ts2=1.0e6, ts3=6.0e4, rval=RC2,
                fit_tol1=1.0e-8, fit_tol2=1.0e-8):
    """ Write the troefit input file.

        :param kpt_dct: k(T,P) values
        :type kpt_dct: dict[pressure:temps]
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
        :return troefit_str: string for the troefit input file
        :rtype: string
    """

    # Write the fitting parameters at the top of the file
    fit_tol_str = f'{fit_tol1:<12.10f} {fit_tol2:<12.10f}'

    # Set the indices for what is fitted
    fit_idxs = ['1', '1', '1', '1', '1', '1', None, None, None, None]
    fit_idxs[6] = '1' if 'alpha' in troe_param_fit_lst else '0'
    fit_idxs[7] = '1' if 'ts1' in troe_param_fit_lst else '0'
    fit_idxs[8] = '1' if 'ts2' in troe_param_fit_lst else '0'
    fit_idxs[9] = '1' if 'ts3' in troe_param_fit_lst else '0'
    fit_idx_str = ' '.join(fit_idxs)
    fit_option_idx = 1

    # Write the parameters strings
    highp_params_str = (
        f'{highp_guess[0]:>8.5E} '
        f'{highp_guess[1]:>8.5E} '
        f'{highp_guess[2]:>8.5E}'
    )
    lowp_params_str = (
        f'{lowp_guess[0]:>8.5E} '
        f'{lowp_guess[1]:>8.5E} '
        f'{lowp_guess[2]:>8.5E}'
    )
    troe_params_str = f'{alpha:>8.5E} {ts1:>8.5E} {ts2:>8.5E} {ts3:>8.5E}'

    # Write the rate constants
    num_tk = len(kpt_dct)
    kpt_str = ''
    for temp, pk_arr in kpt_dct.items():
        [pressures, rate_constants] = pk_arr
        kpt_str += f'{temp:<8.2f}{len(pressures):<4d}\n'
        for pressure, rate in zip(pressures, rate_constants):
            # use p_to_m function from calc/_rates.py  at some point
            density = pressure / (rval * temp)
            kpt_str += f'{density:<14.5E}{rate:<14.8E}\n'

    # Build the fill value dictionary
    fit_keys = {
        'fit_tol_str': fit_tol_str,
        'fit_idx_str': fit_idx_str,
        'fit_option_idx': fit_option_idx,
        'highp_params_str': highp_params_str,
        'lowp_params_str': lowp_params_str,
        'troe_params_str': troe_params_str,
        'num_tk': num_tk,
        'kpt_str': kpt_str
        }

    # Set template name and path for the global keywords section
    template_file_name = 'troefit.mako'
    template_file_path = os.path.join(SRC_PATH, template_file_name)

    # Build global section string
    troefit_str = Template(filename=template_file_path).render(**fit_keys)

    return troefit_str


def run_troefit(path):
    """ Run the troefit executable.

        :param path: Path where the troefit input file exists
        :type path: str
    """

    # Go to path
    start_path = os.getcwd()
    os.chdir(path)

    # Run the executable
    exe_cmd = 'troefit.x'
    subprocess.check_call([exe_cmd])

    # Return to starting dir
    os.chdir(start_path)


def read_params(output_string, conv_factor=1.000):
    """ Parse the output of the troefit code for the fitting
        parameters for a Troe fit.
        :param output_str: output of troefit code
        :type output_str: str
        :param fit_type:
        :type fit_type: str
        :param a_conv_factor: Conversion factor for A parameter
        :type a_conv_factor: float
    """

    # Loop over the lines and find the resulting fit params line
    lines = output_string.splitlines()
    lines.reverse()
    kinf_str, k0_str, ptroe_str = '', '', ''
    for line in lines:
        if line.startswith(' results for iteration'):
            kinf_str = lines[lines.index(line)-3]
            k0_str = lines[lines.index(line)-4]
            ptroe_str = lines[lines.index(line)-5]
            break

    # Assess status of fits (single always assumed True for now)
    fit_success = bool(kinf_str and k0_str and ptroe_str)

    # Grab the fitting parameters
    # Multiple A by given conversion factor and Ea/R term by R to get Ea
    _ = conv_factor
    if fit_success:
        high_params = [float(param) for param in kinf_str.split()]
        low_params = [float(param) for param in k0_str.split()]
        troe_params = [float(param) for param in ptroe_str.split()]
        temp_params = troe_params[:3]
        fcent = troe_params[3]
        fit_params = [high_params, low_params, temp_params, fcent]
    else:
        fit_params = []

    return fit_params
