"""Driver for mechanism sorting
"""

import os
import mechanalyzer


def run(path, spc_name, mech_name, sortmech_name,
        sort_list, isolate_species=None):
    """ main driver for mechanism sorting
    path: path where the input files are
    *_name: filenames
    sort_list: list of sorting criteria
    isolate_species: list of species to be isolated
    """

    spc_dct_full, rxn_param_dct, elem_tuple = mechanalyzer.parser.mech.readfiles(
        os.path.join(path, spc_name), os.path.join(path, mech_name))

    # BUILD  MECH INFORMATION
    mech_info = mechanalyzer.parser.mech.build_dct(spc_dct_full, rxn_param_dct)

    # SORTING: sort the mech and build the sorted rxn param dct
    sorted_idx, cmts_dct, spc_dct = mechanalyzer.parser.mech.sort_mechanism(
        mech_info, spc_dct_full, sort_list, isolate_species)
    rxn_param_dct_sorted = mechanalyzer.parser.mech.reordered_mech(
        rxn_param_dct, sorted_idx)

    # WRITE THE NEW MECHANISM
    mechanalyzer.parser.mech.write_mech(
        elem_tuple, spc_dct, rxn_param_dct_sorted, sortmech_name, comments=cmts_dct)

    if isolate_species:
        rxn_param_dct_rest = mechanalyzer.parser.util.filter_keys(
            rxn_param_dct, rxn_param_dct_sorted)
        mechanalyzer.parser.mech.write_mech(
            elem_tuple, spc_dct_full, rxn_param_dct_rest,
            os.path.join(path, 'mechanism_rest.txt'))
