""" functions for managing creck classes and dictionaries
"""
from collections import Counter
import copy
import numpy as np
import pandas as pd
from mechanalyzer.calculator import ktp_util
from mechanalyzer.calculator.rates import get_aligned_rxn_ratio_dct, merge_rxn_ktp_dcts, is_ktp_dct_withinboundaries
from mechanalyzer.calculator.formulas import sum_heavy_atoms_fmls
from mechanalyzer.calculator.thermo import extract_deltaX_therm, spc_therm_dct_df, spc_therm_df_dct
from mechanalyzer.parser._util import resort_ktp_labels
from chemkin_io.writer import format_rxn_name


R = ['R', 'H', 'OH', 'O2', 'O', 'CH3', 'HO2', 'HCO', 'C2H3']
RSR = ['C5H5', 'C7H7', 'C6H5CH2', 'C6H5O', 'A1O']


def build_creckclass_fromdct(rxn_creckclass_dct, classtype_dct = {}):
    """ build a dataframe with reactiontype, classtype, speciestype, bimoltype for each reaction

    Args:
        rxn_creckclass_dct (dict{rxn: {'speciestype':speciestype, 
        'reactiontype':reactiontype}}): assigned species and reaction type
        classtype_dct (dict{reactiontype(str): classtype(str)}, optional): _description_. Defaults to None.
    
    return creckclass_df: dataframe[['classtype','speciestype','reactiontype','bimoltype'][rxn]]
    """
    
    creckclass_df = pd.DataFrame.from_dict(rxn_creckclass_dct, orient='index', dtype=str,)
    creckclass_df = creckclass_df.replace([None, 'nan', 'NaN', np.nan, ''], 'UNSORTED')
    # add info on bimolecular type
    for speciestype, df_sptype in creckclass_df.groupby('speciestype'):
        for idx, reactiontype in df_sptype['reactiontype'].items():
            bimoltype = set_bimol_type(idx, speciestype, reactiontype)
            creckclass_df.at[idx, 'bimoltype'] = bimoltype
    for reactiontype, df_rxntype in creckclass_df.groupby('reactiontype'):
        rxns = df_rxntype.index
        if reactiontype in classtype_dct.keys():
            cltype = classtype_dct[reactiontype]
        else:
            cltype = 'UNSORTED'
        creckclass_df.loc[rxns, 'classtype'] = cltype
    
    return creckclass_df
        
def set_bimol_type(rxn, speciestype, reactiontype):
    """
    check if a bimolecular reaction is of type M+M, RSR+M, R+M, R+R, RSR+RSR
    from speciestype A-M/R/RSR and reactiontype
    NB for now, relies on rudimental classification. can be upgraded
    by analysing radical character in autochem.
    suggestion: check if molecule or radical, then for rsr check presence of multiple rad structures
    risk: for aromatics, resonance will always be detected
    """
    if speciestype == 'UNSORTED' or reactiontype == ['UNSORTED']:
        return 'UNSORTED'
    
    # return unimol for unimolecular reactions
    if len(rxn[0]) == 1:
        return 'UNIMOL'
    # CHECK SPECIES TYPE 1 - CLASSIFIED ACCORDING TO THE SPECIES TYPE FROM THE USER
    try:
        species_type_1 = speciestype.split('-')[-1]
    except AttributeError:
        return 'UNSORTED'

    # CHECK SPECIES TYPE 2
    try:
        species_type_2 = reactiontype.split('_')[1]
    except AttributeError:
        return 'UNSORTED'
    except IndexError:
        # print warning and classify as unimol if prod is unimol
        # might have been classified in from the backward direction
        if len(rxn[1]) == 1:
            print('rxn {} bimol classification according to backward direction. \
                check that it makes sense! class is {} '.format(rxn, speciestype + ' ' + reactiontype))
            return 'UNIMOL'
        else:
            print('rxn {} bimol classification did not work check class {}. \
                check that it makes sense! '.format(rxn, speciestype + ' ' + reactiontype))
            return 'UNSORTED'
        
    if species_type_2 in R:
        species_type_2 = 'R'
    elif species_type_2 in RSR:
        species_type_2 = 'RSR'
    elif '-' in species_type_2:
        species_type_2 = species_type_2.split('-')[1]
    else:
        species_type_2 = 'NA'

    type_list = [species_type_1, species_type_2]
    type_list.sort()

    return '+'.join(type_list)


def check_rxns_by_spctype(spc, spc_dct, creckclass_df, ftype = 'reactiontype'):
    """ check that a species has all reaction classes for 
    each of its species types

    Args:
        spc (str): species
        spc_dct: species dictionary with 'fct_grp' in keys
        how many functional groups of each type
        creckclass_df (dataframe): [['classtype','speciestype','reactiontype','bimoltype'][rxn]]
        ftype(str): missing type to check
        
    return missing_types: string containing missing ftypes and reaction example
    """
    if ftype not in list(creckclass_df.columns):
        raise ValueError('requested filter check ftype {} unavailable'.format(ftype))
    # extract species types from species dictionary
    spctypes = list(set(get_lumped_fct_grps(spc_dct)[spc].keys()))
    # remove unsorted rxns (reactiontype only)
    creckclass_df = creckclass_df[creckclass_df[ftype] != 'UNSORTED']
    #creckclass_df = creckclass_df[~(creckclass_df == 'UNSORTED').any(axis=1)]
    # find all rxns where the spc is involved in
    rxn_list = ktp_util.get_rxns_for_species(spc, list(creckclass_df.index))
    # get all reaction classes for the selected reactions
    spc_class_df = creckclass_df.loc[rxn_list]
    # find all class/rxn types of the selected species type
    missing_types = ''
    
    for spctype in spctypes:
        missing_types += '{}: \n'.format(spctype)
        # loop over the different species types and check for rxn classes
        spctype_df_spc = spc_class_df[spc_class_df['speciestype'] == spctype]
        spctype_df_all = creckclass_df[creckclass_df['speciestype'] == spctype]
        for ctype, ctype_df in spctype_df_all.groupby(ftype):
            if ctype not in spctype_df_spc[ftype].values:
                # get an example of the missing classtype
                missing_types += ctype + ' \t' + \
                    format_rxn_name(ctype_df.index[0]) + '\n'         
        missing_types += '\n'
    return missing_types
    
  
def check_unclassified(spc_dct, creckclass_df, ftype = 'reactiontype'):
    """ check unsorted reactions by reactants.
        if the reactant is classified according to a given species type,
        suggest possible classification 

    Args:
        spc_dct (dct): {spc: {'inchi':, 'fct_grp_dct':{group1: 1, group2: 1, ...},..}
        species dictionary identifying the functional groups
        creckclass_df (dataframe): [['classtype','speciestype','reactiontype','bimoltype'][rxn]]

    """
    
    unsorted_rxns = creckclass_df[creckclass_df[ftype] == 'UNSORTED'].index
    avail_sptypes = list(set(creckclass_df['speciestype'].values))
    
    # get spc_dct of functional groups available 
    avail_fct_grp_dct = {spc: [grp for grp in list(val['fct_grp'].keys()) 
                               if grp in avail_sptypes]
                         for spc, val in spc_dct.items()
                         } 

    # available types for subsets of species
    avail_ftypes_bysptypes = dict.fromkeys(avail_sptypes)
    for sptype, sptype_df in creckclass_df.groupby('speciestype'):
        avail_ftypes_bysptypes[sptype] = '\n'.join(list(set(sptype_df[ftype].values)))
        
    # check unsorted reactions
    unsorted_but_spavail = {spc: '' for spc in avail_sptypes}
    for rxn in unsorted_rxns:
        rxn_name = format_rxn_name(rxn)
        for rct in rxn[0]: #check reactants
            if rct in avail_fct_grp_dct:
                for spctype in avail_fct_grp_dct[rct]:
                    unsorted_but_spavail[spctype] += rxn_name + '\n'
    
    # sort output
    unsorted_but_spavail = {spc: val for spc, val in unsorted_but_spavail.items() 
                            if len(val) > 0}
    unsorted_but_avail = ''
    for spc, rxns in unsorted_but_spavail.items():
        unsorted_but_avail += '\nspeciestype: {}; potential reactions: \n'.format(spc)
        unsorted_but_avail += rxns + '\n'
        unsorted_but_avail += 'available {} for {} : \n'.format(ftype, spc)
        unsorted_but_avail += avail_ftypes_bysptypes[spc] + '\n'
    
    return unsorted_but_avail


def possible_additionalrxns_byrctypes(creckclass_df, reference_reaction, spc_dct, check_prod_overlap=True):
    """ check unsorted reactions by reactant types (corresponding to a specific class).
        if reactants are classified according to the given species types,
        suggests fit within class
        Useful to check if unclassified reactions in a mechanism might fit within a class

    Args:
        spc_dct (dct): {spc: {'inchi':, 'fct_grp_dct':{group1: 1, group2: 1, ...},..}
        species dictionary identifying the functional groups
        creckclass_df (dataframe): [['classtype','speciestype','reactiontype','bimoltype'][rxn]]
        
        check_prod_overlap (bool): check overlap between species types of products
    """

    unsorted_rxns = creckclass_df[creckclass_df['speciestype'] == 'UNSORTED'].index
    fct_grp_dct = get_lumped_fct_grps(spc_dct) # all functional groups including lumped
    # for each group: get ref reaction
    # reference_reaction(tuple): classified reaction to find compatibility with ((A, B), (C,), ((+M,)))
    ref_rcts, ref_prds = reference_reaction[0], reference_reaction[1]
    # check unsorted reactions
    possible_reactions = ''
    for rxn in unsorted_rxns:
        rxn_name = format_rxn_name(rxn)
        similarity_reactants = rcts_similarity_index(rxn[0], ref_rcts, fct_grp_dct, spc_dct, samenumber=True)
        if similarity_reactants > 0 and not check_prod_overlap:
            possible_reactions += rxn_name + '\n'
        elif similarity_reactants > 0 and check_prod_overlap:
            # also check products
            similarity_products = rcts_similarity_index(
                rxn[1], ref_prds, fct_grp_dct, spc_dct, samenumber=False) # ok for lumped prods
            if similarity_products > 0:
                possible_reactions += rxn_name + '\n'
                
    return possible_reactions

def rcts_similarity_index(rcts1, rcts2, fct_grp_dct, spc_dct, samenumber = True):
    """ are the two sets of rcts overlapping in terms of functional groups?

    Args:
        rcts1 (tuple): set 1 of reactants
        rcts2 (tuple): set 2 of reactants
        fct_grp_dct (dict): functional groups dictionary
        samenumber (bool): rcts1 and rcts2 must be have the same number of elements
    return similarity (int): index with number of similar groups; None if minimum is not reached
    
    """
    # TEST WITH PAH GROWTH REACTIONS
    if samenumber and len(rcts1) != len(rcts2):
        return 0
    
    for spc in rcts1 + rcts2:
        if spc not in spc_dct.keys():
            print('*Warning: species {} not available, \
                  but similarity checks will be computed with recognized species'.format(spc))
    # order based on total stoichiometry
    # allow also when species are not available in the formula list, but set a warning
    fml1, fml2 = [sum_heavy_atoms_fmls([spc_dct[rct]['fml'] for rct in rcts if rct in spc_dct.keys()]) for rcts in [rcts1, rcts2]]
    if fml2 >= fml1:
        rcts1, rcts2 = list(rcts1), list(rcts2)
    elif fml1 > fml2:
        rcts1, rcts2 = list(rcts2), list(rcts1)
    # minimum overlap
    required_overlap = min(len(rcts1), len(rcts2))
    # check that rcts1 is contained in rcts2
    similarity = 0
    # if exact same species overlap: add to index and remove from list
    intersection = Counter(rcts1) & Counter(rcts2)
    similarity += sum(intersection.values())
    
    # remove from lists
    for common_el in list(intersection.elements()):
        rcts1.remove(common_el)
        rcts2.remove(common_el)
    
    # if no more elements to check, it means the reaction is the same-
    # increase similarity to 10 to maximise it with respect to other potentially overlapping groups
    if len(rcts1) == 0 and len(rcts2) == 0:
        return 10
    
    # add functional groups
    rcts1_groups = {sp: list(fct_grp_dct[sp].keys())
                    for sp in rcts1 if sp in fct_grp_dct.keys()}
    rcts2_groups = {sp: list(fct_grp_dct[sp].keys())
                    for sp in rcts2 if sp in fct_grp_dct.keys()}

    # check if functional groups in 1 are contained in 2
    # first, reorder items to ensure largest sets are checked first
    # this ensures that for instance [A1-R, A1,CH3-R] is checked before [A1-R] to maximize overlap
    rcts1_groups = dict(
        sorted(rcts1_groups.items(), key=lambda item: len(item[1]), reverse=True))
    
    # loop
    for _, grps1 in rcts1_groups.items():
        #sp1: one species in set 1
        #grps1: groups corresponding to species in sp1
        matches = {}
        # find maximum intersection with groups in rcts2_groups2
        for sp2, grps2 in rcts2_groups.items():
            similarity_grp1grp2 = sum(((Counter(grps1) & Counter(grps2)).values()))
            if similarity_grp1grp2 > 0:
                matches[sp2] = similarity_grp1grp2
        # select maximum match and pop from dictionary 2
        if matches:
            matching_sp = max(matches, key=matches.get) # maximize matches; if equal N, takes first key
            del rcts2_groups[matching_sp] # delete corresponding instance so you don't check it again
            similarity += 1
    
    if similarity < required_overlap: return 0
    else: return similarity
   
    
def check_spctype_consistency(spc_dct, creckclass_df):
    """ check that the species type defining the reaction 
        is consistent with those extracted automatically for the reactants
    """
    inconsistent_spctype = pd.DataFrame(columns=['speciestype','rct/prd types'], dtype=object)
    # get spc_dct of functional groups 
    fct_grp_dct = get_lumped_fct_grps(spc_dct)
    # check sptypes
    creckclass_df = creckclass_df[creckclass_df['speciestype'] != 'UNSORTED']
    for rxn, sptype in creckclass_df['speciestype'].items():

        rxn_name = format_rxn_name(rxn)
        rcts_types = []
        spcs_to_check = list(rxn[0]) + list(rxn[1])*('=>' not in rxn_name)
        for spc in spcs_to_check:  # check reactants
            if spc in fct_grp_dct.keys():
                rcts_types.extend(list(fct_grp_dct[spc].keys()))

        if len(rcts_types) > 0:
            rcts_types = list(set(rcts_types))
            if sptype not in rcts_types and sptype != 'UNSORTED':
                inconsistent_spctype.loc[rxn_name, ['speciestype', 'rct/prd types']] = [sptype, ' '.join(rcts_types)]
                
    return inconsistent_spctype

def get_lumped_fct_grps(spc_dct):
    """ from a species dictionary, 
        get functional group dictionaries that include also lumped isomers

    Args:
        spc_dct (dct): full species dictionary
    """
    lumped_spcs_dct = {spc: [] for spc in spc_dct.keys() if 'LUMPED' not in spc}
    lumped_spcs = set(spc_dct.keys() - lumped_spcs_dct.keys())
    for lumped_spc in lumped_spcs:
        lumped_spcs_dct[lumped_spc.split('-LUMPED')[0]] += [lumped_spc]

    fct_grp_dct = dict.fromkeys(lumped_spcs_dct.keys())
    for spc in lumped_spcs_dct.keys():
        fct_grps_merged = Counter() #dictionary with number of occurrences of an element
        fct_grps_merged.update(spc_dct[spc]['fct_grp'])
        for lumped_spc in lumped_spcs_dct[spc]:
            fct_grps_merged.update(spc_dct[lumped_spc]['fct_grp'])
            
        fct_grp_dct[spc] = dict(fct_grps_merged)

    return fct_grp_dct

def rm_bw(reactions, sptype, spc_dct):
    """ remove backward reaction if present based on compatibility with selected species
    """
    fct_grp_dct = get_lumped_fct_grps(spc_dct)
    filtered_reactions = []
    for rxn in reactions:
        rxn_bw = (rxn[1], rxn[0], rxn[2])
        if rxn_bw not in reactions: # no bw in list, safe to add 
            filtered_reactions.append(rxn)
        elif rxn in filtered_reactions or rxn_bw in filtered_reactions:
            continue # already added
        elif rxn_bw in reactions: # to check
            rxn_iscompatible = any([sptype in
                                    fct_grp_dct[rct] for rct in rxn[0]])
            rxnbw_iscompatible = any([sptype in
                                        fct_grp_dct[rct] for rct in rxn[1]])
            if rxn_iscompatible and not rxnbw_iscompatible:
                filtered_reactions.append(rxn)
            elif rxnbw_iscompatible and not rxn_iscompatible:
                filtered_reactions.append(rxn_bw)
            else: # unsure case, put both
                print('*Warning: incompatibility between reaction and speciestype- add both fw and bw rxns to plot')
                filtered_reactions.extend([rxn, rxn_bw])

    return filtered_reactions


def compare_thermo_byclass(class_df, spc_dct, therm_dct, grouptype=['speciestype','reactiontype']):
    """ Check consistency of rate parameters within a class
        class_df (dataframe): dataframe of reactions (indexes) and type (reactiontype, speciestype,..)
        grouptype (list): group class_df by types
        therm_dct (dictionary): thermo dictionary to compare DH/DG of reactions
    """
    # remove unsorted rxns
    class_df = class_df[class_df[grouptype].ne('UNSORTED').all(axis=1)]
    # relabe ktp dictionary and class df keys

    # allocation for final dct
    compare_dct = {}
    # get therm df
    therm_df = spc_therm_dct_df(therm_dct)

    # group reactions by type
    for selected_type, df_stype in class_df.groupby(grouptype):
        therm_class_dct = {}
        therm_diff_class_dct = {}
        typename = ' '.join(selected_type)
        print('checking thermo for class {}'.format(selected_type))
        # rename reactions
        reactions = list(resort_ktp_labels(dict.fromkeys(df_stype.index)).keys())
        if grouptype[0] == 'speciestype':
            reactions = rm_bw(reactions, selected_type[0], spc_dct)
        # get ref rxn with non-merged prods
        ref_rxn = ktp_util.get_smallest_stoich_rxn(reactions, spc_dct)[0]
        if all(sp in therm_dct.keys() for sp in ref_rxn[0] + ref_rxn[1]):
            therm_rxn_ref = pd.concat([extract_deltaX_therm(therm_df, ref_rxn[0], ref_rxn[1], var) for var in ['H', 'Cp', 'S', 'G', 'lnQ']], axis=1)
            therm_rxn_ref.columns = ['H', 'Cp', 'S', 'G', 'lnQ']
        else:
            print('Thermo not checked for ref reaction {} of class {}'.format(format_rxn_name(ref_rxn), selected_type))
            continue
        # get thermo for the other rxns if available
        for rxn in reactions:
            if all(sp in therm_dct.keys() for sp in rxn[0] + rxn[1]):
                therm_rxn_dct = pd.concat([extract_deltaX_therm(therm_df, rxn[0], rxn[1], var) for var in ['H', 'Cp', 'S', 'G', 'lnQ']], axis=1)
                therm_rxn_dct.columns = ['H', 'Cp', 'S', 'G', 'lnQ']
                therm_class_dct[rxn] = therm_rxn_dct
                therm_diff_class_dct[rxn] = therm_rxn_dct - therm_rxn_ref
        # add to compare_dct
        compare_dct[typename] = dict(zip(['therm_rxn_dct', 'therm_rxn_diff_dct',], [spc_therm_df_dct(therm_class_dct), spc_therm_df_dct(therm_diff_class_dct)]))

    return compare_dct
    
def compare_byclass(class_df, spc_dct, rxn_ktp_dct, grouptype=['speciestype','reactiontype'],
                    diff_threshold=2, scaling_rules={}):
    """ Check consistency of rate parameters within a class
        class_df (dataframe): dataframe of reactions (indexes) and type (reactiontype, speciestype,..)
        rxn_ktp_dct (dictionary): reaction ktp dictionary at given T,P
        grouptype (list): group class_df by types
        diff_threshold (float): difference in kTP above which a warning is printed
        therm_dct (dictionary): thermo dictionary to compare DH/DG of reactions
        scaling_rules (dictionary): possibly in the future check for scaling parameters?
    """
    # remove unsorted rxns
    class_df = class_df[class_df[grouptype].ne('UNSORTED').all(axis=1)]
    # relabe ktp dictionary and class df keys
    rxn_ktp_dct = resort_ktp_labels(rxn_ktp_dct)
    # allocation for final dct
    compare_dct = {}

    # group reactions by type
    for selected_type, df_stype in class_df.groupby(grouptype):
        typename = ' '.join(selected_type)
        # rename reactions
        reactions = list(resort_ktp_labels(dict.fromkeys(df_stype.index)).keys())
        # remove backward reactions
        # ... how to recognize them? check correspondence of species type? not enough..
        if grouptype[0] == 'speciestype':
            reactions = rm_bw(reactions, selected_type[0], spc_dct)
        # get for each reaction : ktp dictionary
        sub_ktp_dct = {rxn: vals for rxn, vals in rxn_ktp_dct.items() if rxn in reactions}  
        # tot_stoich_dct = {reaction: sum_heavy_atoms_fmls([spc_dct[rct]['fml'] for rct in reaction[0]]) 
        #                   for reaction in reactions} 
        # ref_rxn = max(tot_stoich_dct, key=tot_stoich_dct.get) #if multiple keys?
        
        ###### merge multichannel, but work on different dictionaries
        # sum k(T,P) when reactants are the same (multichannel)
        merged_sub_ktp_dct = merge_rxn_ktp_dcts({}, sub_ktp_dct, sum_multichannel=True)
        reactions_merged = list(merged_sub_ktp_dct.keys())
    
        # if only one reaction left: continue, comparison does not make sense!
        if len(reactions_merged) < 2:
            continue
        ########## TO DO ############
        # normalization: in the future, rescale rate constants with appropriate scaling factors (scaling rules)
        #############################
        ref_rxn_merged = ktp_util.get_smallest_stoich_rxn(reactions_merged, spc_dct)[0]
        # comparison tables: ratio between reference rate constant at given T,P and all of the others
        ref_idx_merged = reactions_merged.index(ref_rxn_merged)
        ratio_list = get_aligned_rxn_ratio_dct(list(merged_sub_ktp_dct.values()),  ref_idx=ref_idx_merged)
        ratio_dcts = dict(zip(reactions_merged, ratio_list))

        warnings_dct = is_ktp_dct_withinboundaries(ratio_dcts, 1/diff_threshold, diff_threshold) # rxns where diff surpasses the threshold 

        warnings_dct_simplified = {rxn: warndct['n_instances'] for rxn, warndct in warnings_dct.items()}

        compare_dct[typename] = {'sub_ktp_dct': merged_sub_ktp_dct,
                                 'ratios_dct': ratio_dcts,
                                 'warnings': warnings_dct_simplified}

        ####### TO DO 
        # when a class has multiple channels like REC_C3.DD-RSR_BRANCH: split so that you sum in the class
        # or organize directly according to classtype ..
           
    return compare_dct


def align_byclass(rxn_param_dct, mech_spc_dct, scale_df, align_if_noscaleavail=False):
    """_summary_

    Args:
        mech_spc_dct (dict): species dictionary including also functional groups
        rxn_param_dct (dict): reaction parameter dictionary already filtered by reaction class
        scale_factors_df (dataframe): dataframe with scaling parameters in arrhenius format for given species
        align_if_noscaleavail (bool, optional): Align to reference reaction if no scaling factor found in scale_factors_df
    reutrn rxn_param_dct_scaled: rxn_param_dct scaled according to indications
    """
    
    def ref_reaction(reaction, allreactions, mech_spc_dct, refspecies=''):
        """ get ref reaction based on smallest stoichiometry or filter by specific species
        """
        rcts1 = reaction[0]
        prds1 = reaction[1]
        fct_grp_dct = get_lumped_fct_grps(mech_spc_dct)
        if refspecies != '':
            rxns = [rxn for rxn in allreactions if refspecies in rxn[0]]
        else:
            rxns = ktp_util.get_smallest_stoich_rxn(allreactions, mech_spc_dct)
            # check, dovrebbe darti una lista di reazioni per canali multipli?

        similarity = {rxn: rcts_similarity_index(rcts1, rxn[0], fct_grp_dct,
                                                 mech_spc_dct, samenumber=True)
                      + rcts_similarity_index(prds1, rxn[1], fct_grp_dct,
                                                mech_spc_dct, samenumber=False)
                      for rxn in rxns}
        # ref rxn is the one with largest similarity

        max_sim = max(similarity.values())
        # if equal similarity, at least pick always the same one
        ref_rxn = sorted([k for k, v in similarity.items() if v == max_sim])[0]
        return ref_rxn # aggiungi anche rxn param dct qui! fallo per ogni reazione e prendi direttamente i parametri

    # cases:
    # multiple channels for reference reaction
    # multiple channels for new product 
    refspecies = ''
    scale_A, scale_n, scale_EA = [1, 0, 0]
    
    # options: do not scale if the rate constant is already very close to ref?
    reactions = list(rxn_param_dct.keys())
    rxn_param_dct_aligned = {} # rewrite only scaled reactions
    cmt_dct_aligned = {}
    for reaction in reactions:
        flitered_scale = [scale_df[scale_df['SPECIES'] == rct]
                          for rct in reaction[0] if not scale_df[scale_df['SPECIES'] == rct].empty]

        if len(flitered_scale) > 0:
            scale_A, scale_n, scale_EA = flitered_scale[0][['A', 'n', 'EA']].values[0]
            refspecies = flitered_scale[0]['REFSPECIES'].values[0]
        elif len(flitered_scale) == 0 and align_if_noscaleavail == False:
            continue # do not scale if no filter found
        # get reference reaction
        try:
            ref_rxn = ref_reaction(reaction, reactions, mech_spc_dct, refspecies=refspecies)
        except ValueError:
            print('SKIPPING REACTION- NOT FOUND FOR REF SPECIES: {}, REF {}'.format(reaction, refspecies))
            continue
        # continue if ref rxn is the same as the current reaction

        if ref_rxn == reaction:
            print('skipping ref reaction')
            continue
        # scale reference parameters
        ref_params = rxn_param_dct[ref_rxn]
        # scale with functions of param dct class found in autoreact._params.RxnParams
        rxn_param_dct_aligned[reaction] = ref_params.__scale__(scale_A, scale_n, scale_EA)
        #print(reaction, rxn_param_dct[reaction])
        cmt_dct_aligned[reaction] = {'inline': '!AUTOALIGN; ref: {}, scale: Aref*{:.1f}, EA+{:d}'.format(
            format_rxn_name(ref_rxn), float(scale_A), int(scale_EA))}

    return rxn_param_dct_aligned, cmt_dct_aligned
        
        
