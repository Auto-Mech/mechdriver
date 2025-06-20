""" extract formulas from spc dct and put in dataframes
"""

import pandas as pd
import automol

# da string a dct
#>>> automol.form.from_string('C6H4O2')
# {'C': 6, 'H': 4, 'O': 2}

def extract_fml_list_fromstr(fmlstr):
    """ 'C6H4O2' => [6, 4, 2, 0, 0, 0]
        from string to list; order is C H O N S Cl
    """
    fml_dct = automol.form.from_string(fmlstr)
    list_el = ['C', 'H', 'O', 'N', 'S', 'Cl']
    fml_list = [automol.form.element_count(fml_dct, X) for X in list_el]
    return fml_list

def extract_fml_df(spc_dct):
    """ Given species dictionary, builds a formula Pandas dataframe:
            index = species names; columns 'fml' (stoichiometry)
            'nC','nH','nO' (number of C/H/O atmoms in the formula)

        :param spc_dct:
        :type spc_dct: dict[]
    """

    fml_df = pd.DataFrame(index=list(spc_dct.keys()),
                          columns=['fml', 'nC', 'nH', 'nO', 'nN', 'nS', 'nCl'])
    for key in spc_dct.keys():
        if 'ts' not in key and 'global' not in key:
            ich = spc_dct[key]['inchi']
            fml_dct = automol.chi.formula(ich)
            fml_df.loc[key, 'fml'] = automol.form.string2(fml_dct)
            for X in ['C', 'H', 'O', 'N', 'S', 'Cl']:
                fml_df.at[key, 'n{}'.format(X)] = automol.form.element_count(fml_dct, X)

    return fml_df

def extract_species(n_at, fml_df):
    """ Extracts species corresponding to a given stoichiometry n_CHONSCl
        from the formulas dataframe.

        :param n_cho: numpy array [x,y,z,...] where
            x = n of C atoms, y = n of H atoms, ...
        :type: n_cho: numpy.ndarray
        :param fml_df: dataframe with
            index = species names; columns 'fml' (stoichiometry),
            'nC','nH','nO' (number of C/H/O atmoms in the formula)
        :returns: list of species with the corresponding n_CHONSCl
    """

    species_set = list(fml_df[(fml_df['nC'] == n_at[0]) & (
        fml_df['nH'] == n_at[1]) & (fml_df['nO'] == n_at[2]) & (
        fml_df['nN'] == n_at[3]) & (fml_df['nS'] == n_at[4]) & (
        fml_df['nCl'] == n_at[5])].index)

    return species_set

def extract_species_sub(n_at, fml_df):
    """ as above but leq
    """

    species_set = list(fml_df[(fml_df['nC'] <= n_at[0]) & (
        fml_df['nH'] <= n_at[1]) & (fml_df['nO'] <= n_at[2]) & (
        fml_df['nN'] <= n_at[3]) & (fml_df['nS'] <= n_at[4]) & (
        fml_df['nCl'] <= n_at[5])].index)

    return species_set

def extract_species_core(n_at, fml_df):
    """ extract any fml with less C/N/S/Cl atoms than indicated. exclude N containing species
        do not filter on oxygen assuming that small species can also oxidize
    """
    species_set = list(fml_df[(fml_df['nC'] <= n_at[0]) & (
        fml_df['nN'] <= n_at[3]) & (fml_df['nS'] <= n_at[4]) & (
        fml_df['nCl'] <= n_at[5])].index)
    
    if 'N2' not in species_set:
        species_set.append('N2')
    return species_set

def extract_species_above(n_at, fml_df):
    """ extract any fml above the indicated limits for C/O/N/S/Cl
    """
    species_set = list(fml_df[(fml_df['nC'] >= n_at[0] ) & (
        fml_df['nN'] >= n_at[3]) & (fml_df['nS'] >= n_at[4]) & (
        fml_df['nCl'] >= n_at[5]) & (fml_df['nO'] >= n_at[2])].index)
    
    if 'N2' in species_set:
        species_set.remove('N2')
    
    return species_set