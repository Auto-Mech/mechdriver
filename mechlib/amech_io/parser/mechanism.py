"""
Read the mechanism file
"""

from mechanalyzer.parser import pes
from mechanalyzer.parser import mech
from ioformat import ptt
import autoparse.find as apf

MECH_INP = 'inp/mechanism.dat'
SORT_INP = 'inp/sort.dat'


def build_pes_dct(job_path, mech_type,
                  spc_dct, run_obj_dct):
    """Build the PES dct
    """

    # Read the string
    mech_str = ptt.read_inp_str(job_path, MECH_INP, remove_comments='!')

    # Build the total PES dct
    _, mech_info, _ = mech.parse_mechanism(
        mech_str, mech_type, spc_dct)
    pes_dct = pes.build_pes_dct(*mech_info[1:])

    # Build an index dct relating idx to formula
    idx_dct, form_dct = pes.build_pes_idx_dct(pes_dct)

    # Reduce the PES dct to only what the user requests
    if run_obj_dct:
        pesnums = [idx_pair[0] for idx_pair in run_obj_dct]
    else:
        pesnums = [idx_pair[0] for idx_pair in run_obj_dct]
    reduced_pes_dct = pes.reduce_pes_dct_to_user_inp(pes_dct, pesnums)

    # Get a dct for all of the connected channels with the PESs to run
    conn_chnls_dct = pes.connected_channels_dct(
        reduced_pes_dct)

    # Form the pes dct that has info formatted to run
    # Get the models in here
    run_pes_dct = pes.pes_dct_w_rxn_lsts(
        reduced_pes_dct, idx_dct, form_dct, conn_chnls_dct, run_obj_dct)

    # Print the channels for the whole mechanism file
    pes.print_pes_channels(pes_dct)

    return run_pes_dct


def read_sort_section(job_path):
    """ reads the options for sorting
        coverts them to two list:
            isolate_species
            sort_list
        for proper input to sortdriver
    """
    # read the sorting input
    sort_str = ptt.read_inp_str(job_path, SORT_INP, remove_comments='#')

    submech_section = apf.all_captures(
        ptt.end_section('isolate_submech'), sort_str)

    if submech_section is None:
        # empty section
        isolate_species = []
    else:
        # format the section
        species = apf.first_capture(
            ptt.paren_section('species'), submech_section[0])
        isolate_species = ptt.build_keyword_lst(species)

    sortmech_section = apf.all_captures(
        ptt.end_section('sort_mech'), sort_str)
    # this section is mandatory
    if sortmech_section is None:
        print('*ERROR: sort_mech section is not defined')
    else:
        criteria = apf.first_capture(
            ptt.paren_section('criteria'), sortmech_section[0])
        n_criteria = apf.first_capture(
            ptt.keyword_pattern('n_criteria_headers'), sortmech_section[0])
        sort_list = ptt.build_keyword_lst(criteria+n_criteria)

    return isolate_species, sort_list
