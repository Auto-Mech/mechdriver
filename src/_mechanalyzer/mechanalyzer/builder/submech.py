"""
Functions for
- submechanism identification
- broad rxn classification
"""

import copy
import pandas as pd
from mechanalyzer.calculator import formulas

# atoms order: C, H, O, N, S, CL
STOICH_DEFAULT = [0, 2, 2, 0, 0, 0] # RETAIL ALL BELOW THIS STOICHIOMETRY + THAT OF THE FUEL
STOICH_DEFAULT_DEL = [1, 0, 2, 0, 0, 0] # DELETE ALL ABOVE THESE STOICH + THAT OF THE FUEL - useful to delete LT above a certain limit
#STOICH_DEFAULT_DEL = [1, 0, 0, 0, 0, 0] # DELETE ALL ABOVE THESE STOICH simultaneous condition on C and O only
STOICH_DCT_ADD = {
    'FUEL': [0, 0, 0, 0, 0, 0],
    'FUEL_RAD': [0, -1, 0, 0, 0, 0],
    'FUEL_ADD_H': [0, 1, 0, 0, 0, 0],
    'FUEL_ADD_CH3': [1, 3, 0, 0, 0, 0],
    'FUEL_ADD_O': [0, 0, 1, 0, 0, 0],
    'FUEL_ADD_OH': [0, 1, 1, 0, 0, 0],
    'FUEL_ADD_O2': [0, 0, 2, 0, 0, 0],
    'R_CH3': [1, 2, 0, 0, 0, 0],
    'R_O': [0, -1, 1, 0, 0, 0],
    'R_O2': [0, -1, 2, 0, 0, 0],
    'R_O4': [0, -1, 4, 0, 0, 0],
    'R_O3-H': [0, -2, 3, 0, 0, 0]}


def prescreen_species_subset(spc_dct, fuel = None, stoich = None, stoich_def = None):
    """ get
        stoichiometry to filter the mech
        starting list of species and corresponding dataframe
        species dictionary and formulas of excluded species
    """
    
    if stoich is not None:
        # extract stoichiometry from string
        stoich = formulas.extract_fml_list_fromstr(stoich)
        
    if fuel is None and stoich is None:
        raise ValueError('submech called, but neither fuel nor stoich filter specified')
    
    elif fuel is not None:
        # extract fuel subset first
        species_list, species_subset_df, stoich_fuel = species_subset_fuel(fuel, spc_dct)
        fml_df_excluded, spc_dct_excluded = species_subset_excluded(species_list, spc_dct)
        # evaluate stoich if not in input
        if stoich is None:
            print('No stoich. filter specified: set as default as stoich {} of fuel + C H O N S Cl: {}'.format(stoich_fuel, stoich_def))
            stoich = stoich_fuel + stoich_def # add also other stoichiometries to the list
            
        return stoich, species_list, species_subset_df, fml_df_excluded, spc_dct_excluded
    
    elif fuel is None:
        # no filter applied, so extract regular fml_df
        fml_df_excluded = formulas.extract_fml_df(spc_dct)
        
        return stoich, [], pd.Series(), fml_df_excluded, spc_dct


def species_subset_excluded(species_list, spc_dct):
    """ fml_df and spc_dct with species not in species_list
    """
    # filter extracted species from spc_dct first
    spc_dct_excluded = copy.deepcopy(spc_dct)
    [spc_dct_excluded.pop(sp) for sp in species_list] #all remaining species
    fml_df_excluded = formulas.extract_fml_df(spc_dct_excluded) #all remaining formulas
    
    return fml_df_excluded, spc_dct_excluded

def species_subset_keep(spc_dct, fuel = None, stoich = None):
    """ call species_subset: extract all species with required stoichiometry
        below stoich: stoich_fuel + stoich_limit (if stoich_limit is negative, you can e.g., delete all core species with respect to the fuel)
        in the reaction, all of the reactants OR all of the products will have to be in the species list
    """
    stoich, species_list, species_subset_df, fml_df_excluded, _ = prescreen_species_subset(spc_dct, fuel = fuel, stoich = stoich, stoich_def = STOICH_DEFAULT)
        
    # extract sub-species : species with stoich. <= the given one (must apply to all species)
    species = formulas.extract_species_sub(stoich, fml_df_excluded) # extract all species with stoichiometry <= stoich
    print('Core (subfuel) species: \n')
    [print(sp) for sp in species]
    series = pd.Series('SUBFUEL', index=species) # these reactions will be called "subfuel"
    species_subset_df = species_subset_df._append(series)
    species_list += species

    return species_list, species_subset_df

def species_subset_del(spc_dct, fuel = None, stoich = None):
    """ extract fuel submech
        then delete from species everything else above a certain stoichiometry
        when sorting: at least all reactants or all products must be in the species list
    """
    stoich, species_list, species_subset_df, fml_df_excluded, spc_dct_excluded = prescreen_species_subset(spc_dct, fuel = fuel, stoich = stoich, stoich_def = STOICH_DEFAULT_DEL)

    # add to the species list all species above a certain stoichiometry
    # useful if you want to keep e.g. growth but not oxidation
    # extract all species with stoichiometry above the selected one
    species_del = formulas.extract_species_above(stoich, fml_df_excluded) 
    [spc_dct_excluded.pop(sp) for sp in species_del] # delete all species above that stoichiometry
    species_list += list(spc_dct_excluded.keys()) #consider all the other species
    
    # from excluded species, extract core rxns vs rxns of bigger species
    fml_df_excluded = formulas.extract_fml_df(spc_dct_excluded)
    if fuel:
        core_species = formulas.extract_species_core(stoich, fml_df_excluded)
        print('Core species: \n')
        [print(sp) for sp in core_species]
        series = pd.Series('CORE', index=core_species)
        species_subset_df = species_subset_df._append(series)
        # list all non-fuel species and non-core species as "supfuel"
        supfuel_list = list(spc_dct_excluded.keys())
        supfuel_list = [sp for sp in supfuel_list if sp not in core_species]
        print('Super-fuel species: \n')
        [print(sp) for sp in supfuel_list]
        series = pd.Series('SUPFUEL', index=supfuel_list)
        species_subset_df = species_subset_df._append(series)
    else:
        # integrate all species as 'core'
        print('Core species: \n')
        [print(sp) for sp in species_list]
        series = pd.Series('CORE', index=species_list)
        species_subset_df = species_subset_df._append(series)
    
    return species_list, species_subset_df


def species_subset_fuel(fuel, spc_dct):
    """ Given the name of a fuel in a spc_dct, it finds all species
        involved in its combustion mech by stoichiometry

        - fuel isomers
        - fuel radicals
        - main radical additions to fuel: fuel+H/CH3/OH/O2/O/HO2, R+O
        - low T mech: RO2, RO4, RO3-H

        Returns all species with corresponding stoichiometries with
        subset identifiers: 'fuel', 'fuel_rad', 'fuel_add_H', ...

        :param fuel: name of the fuel (one of the possible isomers)
        :type fuel: str
        :param spc_dct: species dict
        :type spc_dct: dict[]
        :rtype: pandas.dataframe
    """

    species_list = []
    species_subset_df = pd.Series()

    # extract formulas
    fml_df = formulas.extract_fml_df(spc_dct)

    # Generate list of stoichiometries to extract and the corresponding labels
    stoich_fuel = fml_df.loc[fuel][['nC', 'nH', 'nO', 'nN', 'nS', 'nCl']].values

    # extract desired species and assign labels

    for stoich_type, stoich in STOICH_DCT_ADD.items():
        species = formulas.extract_species(stoich + stoich_fuel, fml_df)
        series = pd.Series(stoich_type, index=species)
        species_subset_df = species_subset_df._append(series)
        species_list += species

    print('Fuel mech species')
    [print(sp) for sp in species_list]

    return species_list, species_subset_df, stoich_fuel

