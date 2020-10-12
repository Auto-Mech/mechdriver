"""
  Fit the rate constants read from the MESS output to
  Arrhenius, Plog, Troe, and Chebyshev expressions
"""

import copy
import numpy
import ratefit
import mess_io
from lib.amech_io import writer
from lib.phydat import phycon
from routines.pf.ktp.fit import _arr as arr
from routines.pf.ktp.fit import _cheb as cheb


def fit_rates(inp_temps, inp_pressures, inp_tunit, inp_punit,
              pes_formula, label_dct, es_info, pf_model,
              mess_path, inp_fit_method, pdep_fit,
              arrfit_thresh):
    """ Parse the MESS output and fit the rates to
        Arrhenius expressions written as CHEMKIN strings
    """

    # Initialize chemkin dct with header
    chemkin_str_dct = {
        'header': writer.ckin.model_header(es_info, pf_model)
    }

    # Loop through reactions, fit rates, and write ckin strings
    rxn_pairs = gen_reaction_pairs(label_dct)
    for (name_i, lab_i), (name_j, lab_j) in rxn_pairs:

        # Set the name and A conversion factor
        reaction = name_i + '=' + name_j
        a_conv_factor = phycon.NAVO if 'W' not in lab_i else 1.00

        print(('\n--------------------------------------------' +
               '------------------------------------------'))
        print('\nReading and Fitting Rates for {}'.format(
            reaction))

        # Read the rate constants out of the mess outputs
        print('\nReading k(T,P)s from MESS output...')
        ktp_dct = read_rates(
            inp_temps, inp_pressures, inp_tunit, inp_punit,
            lab_i, lab_j, mess_path, pdep_fit,
            bimol=numpy.isclose(a_conv_factor, 6.0221e23))

        # Check the ktp dct and fit_method to see how to fit rates
        fit_method = _assess_fit_method(ktp_dct, inp_fit_method)
        
        # Get the desired fits in the form of CHEMKIN strs
        if fit_method is None:
            continue
        elif fit_method == 'arrhenius':
            chemkin_str = arr.perform_fits(
                ktp_dct, reaction, mess_path,
                a_conv_factor, arrfit_thresh)
        elif fit_method == 'chebyshev':
            chemkin_str = cheb.perform_fits(
                ktp_dct, inp_temps, reaction, mess_path,
                a_conv_factor)
        # elif fit_method == 'troe':
        #     # chemkin_str += troe.perform_fits(
        #     #     ktp_dct, reaction, mess_path,
        #     #     troe_param_fit_lst,
        #     #     a_conv_factor, err_thresh)

        # Update the chemkin string dct
        print('\nFitting Parameters in CHEMKIN Format:')
        print(chemkin_str)
        ridx = pes_formula + '_' + reaction.replace('=', '_')
        chemkin_str_dct.update({ridx: chemkin_str})

    return chemkin_str_dct


def gen_reaction_pairs(label_dct):
    """ Generate pairs of reactions
    """
    rxn_pairs = ()
    for name_i, lab_i in label_dct.items():
        if 'F' not in lab_i and 'B' not in lab_i:
            for name_j, lab_j in label_dct.items():
                if 'F' not in lab_j and 'B' not in lab_j and lab_i != lab_j:
                    rxn_pairs += (((name_i, lab_i), (name_j, lab_j)),)

    sorted_rxn_pairs = ()
    for pair in rxn_pairs:
        rct, prd = pair
        if (rct, prd) in sorted_rxn_pairs or (prd, rct) in sorted_rxn_pairs:
            continue
        else:
            sorted_rxn_pairs += ((rct, prd),)

    return sorted_rxn_pairs


# Readers
def read_rates(inp_temps, inp_pressures, inp_tunit, inp_punit,
               rct_lab, prd_lab, mess_path, pdep_fit, bimol=False):
    """ Read the rate constants from the MESS output and
        (1) filter out the invalid rates that are negative or undefined
        and obtain the pressure dependent values
    """

    # Dictionaries to store info; indexed by pressure (given in fit_ps)
    calc_k_dct = {}
    valid_calc_tk_dct = {}
    ktp_dct = {}

    # Read the MESS output file into a string
    mess_file = mess_path+'/rate.out'
    print('mess file', mess_file)
    with open(mess_file, 'r') as mess_file:
        output_string = mess_file.read()

    # Read the temperatures and pressures out of the MESS output
    mess_temps, tunit = mess_io.reader.rates.get_temperatures(
        output_string)
    mess_pressures, punit = mess_io.reader.rates.get_pressures(
        output_string)

    assert inp_temps <= mess_temps
    assert inp_pressures <= mess_pressures
    assert inp_tunit == tunit
    assert inp_punit == punit

    # Loop over the pressures obtained from the MESS output
    for pressure in mess_pressures:

        # Read the rate constants
        if pressure == 'high':
            rate_ks = mess_io.reader.highp_ks(
                output_string, rct_lab, prd_lab)
        else:
            rate_ks = mess_io.reader.pdep_ks(
                output_string, rct_lab, prd_lab, pressure)

        # Store in a dictionary
        calc_k_dct[pressure] = rate_ks

    # Remove k(T) vals at each P where where k is negative or undefined
    # If ANY valid k(T,P) vals at given pressure, store in dct
    print('\nRemoving invalid k(T,P)s from MESS output that are either:\n',
          '  (1) negative, (2) undefined [***], or (3) below 10**(-21) if',
          'reaction is bimolecular')
    for pressure, calc_ks in calc_k_dct.items():
        filtered_temps, filtered_ks = ratefit.fit.get_valid_tk(
            mess_temps, calc_ks, bimol)
        if filtered_ks.size > 0:
            valid_calc_tk_dct[pressure] = [filtered_temps, filtered_ks]

    # Filter the ktp dictionary by assessing the presure dependence
    if valid_calc_tk_dct:
        if list(valid_calc_tk_dct.keys()) == ['high']:
            print('\nValid k(T)s only found at High Pressure...')
            ktp_dct['high'] = valid_calc_tk_dct['high']
        else:
            if pdep_fit:
                print('\nUser requested to assess pressure dependence',
                      'of reaction.')
                pkeys = ('assess_pdep_temps', 'tolerance', 'plow', 'phigh')
                dct = {k: pdep_fit[k] for k in pdep_fit if k in pkeys}
                rxn_is_pdependent = ratefit.fit.assess_pressure_dependence(
                    valid_calc_tk_dct, **dct)
                if rxn_is_pdependent:
                    print('  Reaction found to be pressure dependent.',
                          'Fitting all k(T)s from all pressures',
                          'found in MESS.')
                    ktp_dct = copy.deepcopy(valid_calc_tk_dct)
                else:
                    no_pdep_pval = pdep_fit['no_pdep_pval']
                    print('  No pressure dependence detected.',
                          'Grabbing k(T)s at {} {}'.format(
                              no_pdep_pval, punit))
                    if no_pdep_pval in valid_calc_tk_dct:
                        ktp_dct['high'] = valid_calc_tk_dct[no_pdep_pval]
            else:
                ktp_dct = copy.deepcopy(valid_calc_tk_dct)

        # Patchy way to get high-pressure rates in dct if needed
        if 'high' not in ktp_dct and 'high' in valid_calc_tk_dct.keys():
            ktp_dct['high'] = valid_calc_tk_dct['high']

    return ktp_dct


def _assess_fit_method(ktp_dct, inp_fit_method):
    """ Assess if there are any rates to fit and if so, check if
        the input fit method should be used, or just simple Arrhenius
        fits will suffice because there is only one pressure for which
        rates exist to be fit.
    """

    if ktp_dct:
        pressures = list(ktp_dct.keys())
        npressures = len(pressures)
        # If only one pressure (outside HighP limit), just run Arrhenius
        if npressures == 1:
            fit_method = 'arrhenius'
        elif npressures == 2 and 'high' in pressures:
            fit_method = 'arrhenius'
        else:
            fit_method = inp_fit_method
    else:
        fit_method = None

    # Print message to say what fitting will be done
    if fit_method == 'arrhenius':
        if inp_fit_method != 'arrhenius':
            print('\nNot enough pressure-dependent rates for Troe/Chebyshev.')
        print('\nFitting k(T,P)s to PLOG/Arrhenius Form....')
    elif fit_method == 'chebyshev':
        print('\nFitting k(T,P)s to Chebyshev Form...')
    elif fit_method == 'troe':
        print('\nFitting k(T,P)s to Troe Form...')
    elif fit_method is None:
        print('\nNo valid k(T,Ps)s from MESS output to fit.')
        print('\nSkipping to next reaction...')

    return fit_method
