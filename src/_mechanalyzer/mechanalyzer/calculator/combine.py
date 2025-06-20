""" Combines objects pertaining to mechanisms
"""

import copy
from mechanalyzer.calculator import compare


def comb_mechs(rxn_param_dct1, rxn_param_dct2, spc_nasa7_dct1, spc_nasa7_dct2,
               mech_spc_dct1, mech_spc_dct2, ste_mech1_only=False,
               strip_ste=False):
    """ Receives two sets of objects describing two mechanisms (rates, thermo,
        and species) and outputs a single combined mechanism. In the case of
        duplicate values, the first mechanism takes precedence
    """

    # If there is stereo in mech1 ONLY, overwrite the strip_ste input because
    # must strip stereo when comparing
    if ste_mech1_only:
        strip_ste = True

    # Get the instructions for renaming the species
    print('Inside comb_mechs: getting rename instructions...')
    rename_instr = compare.get_rename_instr(
        mech_spc_dct1, mech_spc_dct2, strip_ste=strip_ste)

    print('first set of rename_instr:\n', rename_instr)

    # If indicated, do some stereo checks, etc.
    if ste_mech1_only:
        # Note the flipped order
        ste_instr = compare.get_rename_instr(mech_spc_dct2, mech_spc_dct1,
                                                strip_ste=True)
        print('second set of rename_instr (ste_instr):\n', ste_instr)
        _, ste_dct = compare.rename_species(rxn_param_dct1, ste_instr, 'rxn')
    else:
        ste_dct = None

    print('Inside comb_mechs: ste_dct:\n', ste_dct)

    # Combine each pair of dictionaries
    comb_rxn_param_dct = comb_dcts(rxn_param_dct1, rxn_param_dct2,
                                   rename_instr, target_type='rxn',
                                   ste_dct=ste_dct)
    comb_spc_nasa7_dct = comb_dcts(spc_nasa7_dct1, spc_nasa7_dct2,
                                   rename_instr, target_type='spc')
    comb_mech_spc_dct = comb_dcts(mech_spc_dct1, mech_spc_dct2, 
                                  rename_instr, target_type='spc')

    #print('Inside comb_mechs: comb_mech_spc_dct:\n', comb_mech_spc_dct)

    return comb_rxn_param_dct, comb_spc_nasa7_dct, comb_mech_spc_dct


def comb_mult_mechs(rxn_param_dcts, spc_nasa7_dcts, mech_spc_dcts):

    tot_rxn_param_dct = copy.deepcopy(rxn_param_dcts[0])
    tot_spc_nasa7_dct = copy.deepcopy(spc_nasa7_dcts[0])
    tot_mech_spc_dct = copy.deepcopy(mech_spc_dcts[0])
    ncombs = len(dcts) - 1  # n-1 combinations to do
    for idx in range(ncombs):
        tot_rxn_param_dct, tot_spc_nasa7_dct, tot_mech_spc_dct = comb_mechs(
            tot_rxn_param_dct, rxn_param_dcts[idx+1], 
            tot_spc_nasa7_dct, spc_nasa7_dcts[idx+1], 
            tot_mech_spc_dct, mech_spc_dcts[idx+1])

    return tot_rxn_param_dct, tot_spc_nasa7_dct, tot_mech_spc_dct
    

def comb_dcts(dct1, dct2, rename_instr, target_type='rxn', ste_dct=None):
    """ Combines two dictionaries; can be rxn_param_dcts, spc_nasa7_dcts, or
        mech_spc_dcts
    """

    ste_dct = ste_dct or {}
    # Rename the second dictionary to match the first
    renamed_dct2, _ = compare.rename_species(dct2, rename_instr, target_type)

    # Remove any instances in dct2 that are in dct1
    for key in dct1.keys():  # key is either spc or rxn
        if target_type == 'rxn':
            matching_rxn, _ = compare.assess_rxn_match(key, renamed_dct2)
            if key in ste_dct:
                renamed_dct2.pop(key)
            elif matching_rxn:
                renamed_dct2.pop(matching_rxn)
        else:  # 'spc'
            if key in renamed_dct2:
                renamed_dct2.pop(key)

    # Fix reactions that have one or more stereoisomers but are not in dct1
    #if target_type == 'rxn':
    #    for rxn, params in renamed_dct2.items():
    #         remap_dct = _remap_dct(rxn, rename_instr, ste_spc_dct)    

    # Add renamed and cleaned up dct2 to the combined dct
    comb_dct = copy.deepcopy(dct1)  # combined dct is initially just dct1
    for key, value in renamed_dct2.items():
        comb_dct[key] = value

    return comb_dct


def _remap_dct(rxn, rename_instr, ste_spc_dct):

    def check_spcs(spcs, remap_dct):
    
        remap_dct = copy.deepcopy(remap_dct)
        for spc in spcs:
            for non_ste_spc, ste_spc_list in ste_spc_dct.items():
                if spc in ste_spc_list:
                    remap_dct[spc] = ste_spc_list
                    break  # stop looking through the ste_spc_dct

        return remap_dct
        
    rcts, prds, _ = rxn
    remap_dct = {}
    remap_dct = check_spcs(rcts, remap_dct)   
    remap_dct = check_spcs(prds, remap_dct)   
    if remap_dct == {}:
        remap_dct = None

    return remap_dct
    




