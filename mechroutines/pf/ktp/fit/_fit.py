"""
  Fit the rate constants read from the MESS output to
  Arrhenius, Plog, Troe, and Chebyshev expressions
"""

import copy
import numpy
import ratefit
import mess_io
from phydat import phycon
from mechlib.amech_io import writer
from mechlib.amech_io import printer as ioprinter
from mechroutines.pf.ktp.fit import _arr as arr
from mechroutines.pf.ktp.fit import _cheb as cheb


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

        ioprinter.obj('vspace')
        ioprinter.obj('line_dash')
        ioprinter.reading(
            'Reading and Fitting Rates for {}'.format(
                reaction))

        # Read the rate constants out of the mess outputs
        ioprinter.reading(
            'Reading k(T,P)s from MESS output...', newline=1)
        ktp_dct = read_rates(
            inp_temps, inp_pressures, inp_tunit, inp_punit,
            lab_i, lab_j, mess_path, pdep_fit,
            bimol=numpy.isclose(a_conv_factor, 6.0221e23))

        # Check the ktp dct and fit_method to see how to fit rates
        fit_method = _assess_fit_method(ktp_dct, inp_fit_method)

        # Get the desired fits in the form of CHEMKIN strs
        if fit_method is None:
            continue
        if fit_method == 'arrhenius':
            chemkin_str = arr.perform_fits(
                ktp_dct, reaction, mess_path,
                a_conv_factor, arrfit_thresh)
        elif fit_method == 'chebyshev':
            chemkin_str = cheb.perform_fits(
                ktp_dct, inp_temps, reaction, mess_path,
                a_conv_factor)
            if not chemkin_str:
                chemkin_str = arr.perform_fits(
                    ktp_dct, reaction, mess_path,
                    a_conv_factor, arrfit_thresh)
        # elif fit_method == 'troe':
        #     # chemkin_str += troe.perform_fits(
        #     #     ktp_dct, reaction, mess_path,
        #     #     troe_param_fit_lst,
        #     #     a_conv_factor, err_thresh)

        # Update the chemkin string dct
        ioprinter.obj('vspace')
        ioprinter.info_message(
            'Final Fitting Parameters in CHEMKIN Format:\n', chemkin_str)
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
    filt_ktp_dct = {}
    ktp_dct = {}

    # Read the MESS output file into a string
    mess_file = mess_path+'/rate.out'
    ioprinter.debug_message('mess file', mess_file)
    with open(mess_file, 'r') as mess_file:
        output_string = mess_file.read()

    # Read the temperatures and pressures out of the MESS output
    # mess_temps, tunit = mess_io.reader.rates.get_temperatures(
    #     output_string)
    _, punit = mess_io.reader.rates.get_pressures(
        output_string)

    mess_temps = inp_temps
    mess_pressures = inp_pressures

    mess_temps = list(set(list(mess_temps)))
    mess_temps.sort()

    assert inp_temps <= mess_temps
    assert inp_pressures <= mess_pressures
    # assert inp_tunit == tunit
    # assert inp_punit == punit

    # Loop over the pressures obtained from the MESS output
    calc_ktp_dct = mess_io.reader.rates.ktp_dct(
        output_string, rct_lab, prd_lab)

    # Remove k(T) vals at each P where where k is negative or undefined
    # If ANY valid k(T,P) vals at given pressure, store in dct
    ioprinter.info_message(
        'Removing invalid k(T,P)s from MESS output that are either:\n',
        '  (1) negative, (2) undefined [***], or (3) below 10**(-21) if',
        'reaction is bimolecular', newline=1)
    filt_ktp_dct = ratefit.fit.filter_ktp_dct(calc_ktp_dct, bimol)

    # Filter the ktp dictionary by assessing the presure dependence
    if filt_ktp_dct:
        if list(filt_ktp_dct.keys()) == ['high']:
            ioprinter.info_message(
                'Valid k(T)s only found at High Pressure...', newline=1)
            ktp_dct['high'] = filt_ktp_dct['high']
        else:
            if pdep_fit:
                ioprinter.info_message(
                    'User requested to assess pressure dependence',
                    'of reaction.', newline=1)
                pkeys = ('assess_pdep_temps', 'tolerance', 'plow', 'phigh')
                dct = {k: pdep_fit[k] for k in pdep_fit if k in pkeys}
                rxn_is_pdependent = ratefit.fit.assess_pressure_dependence(
                    filt_ktp_dct, **dct)
                if rxn_is_pdependent:
                    ioprinter.info_message(
                        'Reaction found to be pressure dependent.',
                        'Fitting all k(T)s from all pressures',
                        'found in MESS.', indent=1/2.)
                    ktp_dct = copy.deepcopy(filt_ktp_dct)
                else:
                    no_pdep_pval = pdep_fit['no_pdep_pval']
                    ioprinter.info_message(
                        'No pressure dependence detected.',
                        'Grabbing k(T)s at {} {}'.format(
                            no_pdep_pval, punit), newline=1)
                    if no_pdep_pval in filt_ktp_dct:
                        ktp_dct['high'] = filt_ktp_dct[no_pdep_pval]
            else:
                ktp_dct = copy.deepcopy(filt_ktp_dct)

        # Patchy way to get high-pressure rates in dct if needed
        if 'high' not in ktp_dct and 'high' in filt_ktp_dct.keys():
            ktp_dct['high'] = filt_ktp_dct['high']

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
            ioprinter.error_message(
                'Rates at not enough pressures for Troe/Chebyshev.', newline=1)
        ioprinter.info_message(
            'Fitting k(T,P)s to PLOG/Arrhenius Form....', newline=1)
    elif fit_method == 'chebyshev':
        ioprinter.info_message(
            'Fitting k(T,P)s to Chebyshev Form...', newline=1)
    elif fit_method == 'troe':
        ioprinter.info_message(
            'Fitting k(T,P)s to Tree Form...', newline=1)
    elif fit_method is None:
        ioprinter.error_message(
            'No valid k(T,Ps)s from MESS output to fit.',
            'Skipping to next reaction...', newline=1)

    return fit_method
