""" Runs the sorter
"""
from mechanalyzer.parser import mech as mparser
from mechanalyzer.parser import spc as sparser
from mechanalyzer.parser import new_spc as new_sparser
from mechanalyzer.builder import sort_fct


SPC_TYPE = 'csv'
MECH_TYPE = 'chemkin'
DFG = {
    'Tref': 1500,
    'keepfiltered': 0,
    'DH': 0,
    'H5H3ratio': 0,
    'kratio': 1e50, #depends too much on temperature to make it a default that makes sense
    'kabs': 1e50,
    'lookforpromptchains': 0 
}

# Functions to take the mechanism strings (may want to further simplify)

def sorted_pes_dct(spc_str, mech_str, isolate_spc, sort_lst):
    """ Function that extracts sorted subpes for a mech
    """

    srt_mch, _ = _sort_objs(spc_str, mech_str, sort_lst, isolate_spc)

    return srt_mch.return_pes_dct()


def sorted_mech(spc_str, mech_str, isolate_spc, sort_lst, spc_therm_dct=None, dct_flt_grps={}, stereo_optns=False):
    """ Function that conducts the sorting process for all of the above tests
    """

    # Build mech information
    srt_mch, rxn_param_dct = _sort_objs(
        spc_str, mech_str, sort_lst, isolate_spc, stereo_optns=stereo_optns)

    pes_groups = ''
    rxns_filter = ''
    # if prompt groups detected: retrieve grps info
    if 'submech_prompt' in sort_lst and spc_therm_dct and dct_flt_grps:
        DFG.update(dct_flt_grps)
        # check that Tref is contained in spc_therm_dct, otherwise raise exception
        # assume T is the same for any species
        T_thermdct = list(spc_therm_dct.values())[0][0]

        if DFG['Tref'] not in T_thermdct:
            raise ValueError('Tref {} K for calculations not among T of therm_dct {} \
                plese update ref value in dct_flt_grps[Tref] or add Tref to therm dct'.format(DFG['Tref'], T_thermdct))
        elif len(T_thermdct) < 4:
            raise ValueError('therm dictionary should contain at least 3 temperatures for interpolation purposes')
        
        srt_mch.filter_groups_prompt(
            spc_therm_dct, DFG)
        pes_groups = srt_mch.grps
        rxns_filter = srt_mch.rxns_dh

    elif 'submech_prompt' in sort_lst:
        pes_groups = srt_mch.grps

    sorted_idx, cmts_dct, spc_dct_ord = srt_mch.return_mech_df()
    rxn_param_dct_sort = reordered_mech(rxn_param_dct, sorted_idx)
    
    return rxn_param_dct_sort, spc_dct_ord, cmts_dct, pes_groups, rxns_filter


def _sort_objs(spc_str, mech_str, sort_lst, isolate_spc, stereo_optns=False):
    """ Build the sort-mech object
    """

    # Build mech information
    # spc_dct = sparser.build_spc_dct(spc_str, SPC_TYPE)
    spc_dct = new_sparser.parse_mech_spc_dct(
        spc_str, chk_ste=stereo_optns, chk_match=stereo_optns, verbose=stereo_optns, canon_ent=stereo_optns)

    rxn_param_dct = mparser.parse_mechanism(
        mech_str, MECH_TYPE)

    # Build the sorted mechanism and species objects
    srt_mch = sorting(rxn_param_dct, spc_dct, sort_lst, isolate_spc)
    # spc_dct_ord = sparser.reorder_by_atomcount(spc_dct)

    return srt_mch, rxn_param_dct


# Functions that perform the individual sorting process
def sorting(rxn_param_dct, spc_dct, sort_lst, isolate_species):
    """ Uses the SortMech class to sort mechanism info and
        returns the sorted indices and the corresponding comments.

    :param rxn_param_dct: reaction parameter dictionary
    :param spc_dct: species dictionary
    :param sort_lst: list with sorting criteria
    :param isolate_species: species you want to isolate in the final mechanism
    :type isolate_species: list()

    calls sorting functions in mechanalyzer/pes
    returns the rxn indices associated with the comments about sorting
    """

    srt_mch = sort_fct.SortMech(rxn_param_dct, spc_dct)
    srt_mch.sort(sort_lst, isolate_species)

    return srt_mch


def reordered_mech(rxn_param_dct, sorted_idx):
    """ Sort the reaction parameter dcitionary using the indices from
        sort functions.

        :param rxn_param_dct: non-sorted reaction parameter dictionary
        :type rxn_param_dct: dct
        :param sorted_idx: indices of the rxn_param_dct in the desired order
        :type sorted_idx: list
        :return rxn_param_dct_sorted: sorted reaction parameter dictionary
        :rtype: dct
    """

    sorted_val = list(map(rxn_param_dct.get, sorted_idx))
    rxn_param_dct_sorted = dict(zip(sorted_idx, sorted_val))

    return rxn_param_dct_sorted
