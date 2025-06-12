"""
functions called in the sorter
to be sorted . . .
"""
import numpy
import copy
import pandas as pd

# functions for reaction keys of rxn/ktp dictionaries
def get_rxns_for_species(spc, all_rxns):
    """ filters keys of all_rxns and returns only those in which spc is involved in

    Args:
        spc (str): species of interest
        all_rxns (list): rxns tuples [((A,),(B,),(thrdby,)), ...]

    Returns:
        rxns: filtered reaction list
    """
    rxns = []
    for rxn in all_rxns:
        if spc in rxn[0] or spc in rxn[1]:
            rxns.append(rxn)
    return rxns
            
# functions for ktp dictionaries

def get_aligned_rxn_ratio_dct(aligned_rxn_dct_entry):
    """ converts the entry of the aligned_rxn_ktp_dictionary to the ratios
        between the reference rate and the given rate

    :param aligned_rxn_dct_entry: entry of aligned_rxn_ktp/ratio_dct
    :type aligned_rxn_dct_entry:
        list[dct{pressure: numpy.array(temps), numpy.array(values)}]
    :return aligned_ratio_dct_entry: aligned dictionary entry
    :rtype: list(dct)
    """

    ref_ktp_dct = aligned_rxn_dct_entry[0]
    ratio_dct_entry = []
    for mech_idx, ktp_dct in enumerate(aligned_rxn_dct_entry):
        # If (1) you are on the first ktp_dct, (2) the ref_ktp_dct is None,
        # or (3) the current_ktp_dct is None, set the ratio_dct to None
        if mech_idx == 0 or ref_ktp_dct is None or ktp_dct is None:
            ratio_dct = None
        # Otherwise, calculate the ratio_dct
        else:
            ratio_dct = {}
            for pressure, (temps, kts) in ktp_dct.items():
                # If pressure defined in ref ktp_dct: calculate and store ratio
                if pressure in ref_ktp_dct.keys():
                    _, ref_kts = ref_ktp_dct[pressure]
                    ratios = kts / ref_kts
                    ratio_dct[pressure] = (temps, ratios)
            if ratio_dct == {}:  # account for when no pressures contain ratios
                ratio_dct = None

        # Append the current ratio_dct
        ratio_dct_entry.append(ratio_dct)

    return ratio_dct_entry


def get_max_aligned_values(aligned_rxn_dct_entry):
    """ Gets the maximum values for each reaction from an entry (value) of
        either an aligned_rxn_ktp_dct or an aligned_rxn_ratio_dct

    :param aligned_rxn_dct_entry: entry of aligned_rxn_ktp/ratio_dct
                                  or of ktp dct
    :type aligned_rxn_dct_entry:
        list[dct{pressure: numpy.array(temps), numpy.array(values)}]
        or directly dct{pressure: numpy.array(temps), numpy.array(values)}
    :return max_val: max value
    :rtype: float
    """

    max_val = 0

    for single_dct in aligned_rxn_dct_entry:
        if single_dct is not None:
            for _, (_, values) in single_dct.items():
                if max(values) > max_val:
                    max_val = max(values)

    return max_val
