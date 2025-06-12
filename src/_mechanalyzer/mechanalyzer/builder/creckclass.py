""" functions for managing creck classes and dictionaries
"""
from collections import Counter
import pandas as pd
from mechanalyzer.calculator.ktp_util import get_rxns_for_species
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
    
    creckclass_df = pd.DataFrame.from_dict(rxn_creckclass_dct, orient='index', dtype=str)
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
    if speciestype == '' or reactiontype in ['', 'UNSORTED']:
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

def check_rxns_by_spctype(spc, spctypes, creckclass_df, ftype = 'reactiontype'):
    """ check that a species has all reaction classes for 
    each of its species types

    Args:
        spc (str): species
        spctypes (list(str)): grp1, grp2,..
        how many functional groups of each type
        creckclass_df (dataframe): [['classtype','speciestype','reactiontype','bimoltype'][rxn]]
        ftype(str): missing type to check
        
    return missing_types: string containing missing ftypes and reaction example
    """
    if ftype not in list(creckclass_df.columns):
        raise ValueError('requested filter check ftype {} unavailable'.format(ftype))
    # remove unsorted rxns (reactiontype only)
    creckclass_df = creckclass_df[creckclass_df[ftype] != 'UNSORTED']
    #creckclass_df = creckclass_df[~(creckclass_df == 'UNSORTED').any(axis=1)]
    # find all rxns where the spc is involved in
    rxn_list = get_rxns_for_species(spc, list(creckclass_df.index))
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
        
    return missing_types
    
  
def check_unclassified(mech_spc_dct, creckclass_df, ftype = 'reactiontype'):
    """ check unsorted reactions by reactants.
        if the reactant is classified according to a given species type,
        suggest possible classification 

    Args:
        mech_spc_dct (dct): {spc: {'inchi':, 'fct_grp_dct':{group1: 1, group2: 1, ...},..}
        species dictionary identifying the functional groups
        creckclass_df (dataframe): [['classtype','speciestype','reactiontype','bimoltype'][rxn]]

    """
    
    unsorted_rxns = creckclass_df[creckclass_df[ftype] == 'UNSORTED'].index
    avail_sptypes = list(set(creckclass_df['speciestype'].values))
    
    # get spc_dct of functional groups available 
    avail_fct_grp_dct = {spc: [grp for grp in list(val['fct_grp'].keys()) 
                               if grp in avail_sptypes]
                         for spc, val in mech_spc_dct.items()
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
        unsorted_but_avail += 'speciestype: {}; potential reactions: \n'.format(spc)
        unsorted_but_avail += rxns + '\n'
        unsorted_but_avail += 'available {} for {} : \n'.format(ftype, spc)
        unsorted_but_avail += avail_ftypes_bysptypes[spc]
    
    return unsorted_but_avail

def check_spctype_consistency(mech_spc_dct, creckclass_df):
    """ check that the species type defining the reaction 
        is consistent with those extracted automatically for the reactants
    """
    inconsistent_spctype = pd.DataFrame(columns=['speciestype','rct/prd types'], dtype=object)
    # get spc_dct of functional groups 
    fct_grp_dct = get_lumped_fct_grps(mech_spc_dct)
    print(fct_grp_dct['C5H6'])
    # check sptypes
    creckclass_df = creckclass_df[creckclass_df['speciestype'] != 'UNSORTED']
    for rxn, sptype in creckclass_df['speciestype'].items():

        rxn_name = format_rxn_name(rxn)
        rcts_types = []
        spcs_to_check = list(rxn[0]) + list(rxn[1])*('=>' not in rxn_name)

        for spc in spcs_to_check:  # check reactants
            if spc in fct_grp_dct.keys():
                rcts_types += list(fct_grp_dct[spc].keys())
                
        if len(rcts_types) > 0:
            if sptype not in rcts_types:
                print(rxn_name, sptype, rcts_types)
                # inconsistent_spctype.loc[rxn_name, ['speciestype', 'rct/prd types']] = [sptype, ' '.join(rcts_types)]
                # print(rxn_name, inconsistent_spctype.loc[rxn_name].values())
                
    return inconsistent_spctype

def get_lumped_fct_grps(mech_spc_dct):
    """ from a species dictionary, get functional group dictionaries that include also lumped isomers

    Args:
        mech_spc_dct (dct): full species dictionary
    """
    lumped_spcs_dct = {spc: [] for spc in mech_spc_dct.keys() if 'LUMPED' not in spc}
    lumped_spcs = set(mech_spc_dct.keys() - lumped_spcs_dct.keys())
    for lumped_spc in lumped_spcs:
        lumped_spcs_dct[lumped_spc.split('-LUMPED')[0]] += [lumped_spc]

    fct_grp_dct = dict.fromkeys(lumped_spcs_dct.keys())
    for spc in lumped_spcs_dct.keys():
        fct_grps_merged = Counter()
        fct_grps_merged.update(mech_spc_dct[spc]['fct_grp'])
        for lumped_spc in lumped_spcs_dct[spc]:
            fct_grps_merged.update(mech_spc_dct[lumped_spc]['fct_grp'])
            
        fct_grp_dct[spc] = dict(fct_grps_merged)

    return fct_grp_dct
        
