"""
Useful functions for mechanalyzer.parser
- mechanism reading
- extract useful info from the mechanism (reactions, formulas..)
"""

import copy
from autoparse import find
from automol.form import element_count as n_el
from automol.form import atom_count
from automol.form import join
from automol.form import string2
from automol.chi import formula


def order_rct_bystoich(rct_names_lst, spc_dct=None):
    """ Reorder reactants and products based on the higher number of atoms
        If no species dictionary is given:
            reorder just according to name length,
        If length or stoichiometry is the same:
            reorder by alphabetical order
    """
    rct_names_lst_ordered = copy.deepcopy(rct_names_lst)
    ich_dct = {}
    if spc_dct:
        for key in spc_dct.keys():
            if 'ts' not in key and 'global' not in key:
                ich_dct[key] = spc_dct[key]['inchi']

        for key, val in enumerate(rct_names_lst_ordered):
            rct_names = val
            rct_ichs = list(map(ich_dct.__getitem__, rct_names))
            fml_rct = list(map(formula, rct_ichs))
            atoms_rct = list(map(atom_count, fml_rct))
            if len(rct_names) == 2:
                if atoms_rct[1] > atoms_rct[0]:
                    # swap places of reactants 1 and 2
                    rct_names_lst_ordered[key] = (rct_names[1], rct_names[0])
                elif atoms_rct[1] == atoms_rct[0]:
                    rct_names_lst_ordered[key] = order_names(rct_names)
#                    rct_names = list(rct_names)
#                    rct_names.sort()
#                    rct_names_lst_ordered[key] = tuple(rct_names)

    else:
        for key, val in enumerate(rct_names_lst_ordered):
            rct_names = val
            rct_names_lst_ordered[key] = order_names(rct_names)
#            if len(rct_names) == 2:
#                if len(rct_names[1]) > len(rct_names[0]):
#                    # swap places of reactants 1 and 2
#                    rct_names_lst_ordered[key] = (rct_names[1], rct_names[0])
#                elif len(rct_names[1]) == len(rct_names[0]):
#                    rct_names = list(rct_names)
#                    rct_names.sort()
#                    rct_names_lst_ordered[key] = tuple(rct_names)

    return rct_names_lst_ordered


def order_names(rct_names):
    """ order names according to higher N of atoms, or length;
        if same length, use alphabetical order
        returns a reordered tuple
    """
    if len(rct_names) == 2:
        
        atoms_rct = extract_numfromname(rct_names)
        if atoms_rct[1] > atoms_rct[0]:
            # order according to larger number of atoms
            rct_names = (rct_names[1], rct_names[0])
        elif atoms_rct[1] == atoms_rct[0]:
            # change criterion: names length
            if len(rct_names[1]) > len(rct_names[0]):
                # swap places of reactants 1 and 2
                rct_names = (rct_names[1], rct_names[0])
            elif len(rct_names[1]) == len(rct_names[0]):
                # alphabetical
                rct_names_lst = list(rct_names)
                rct_names_lst.sort()
                rct_names = tuple(rct_names_lst)
        # else: do nothing - original name is fine
                        
    return rct_names


def extract_numfromname(rct_names):
    """ extract numbers from a name string and return their sum
    """
    numbers_tuples = [find.all_captures(r'\d+', r) for r in rct_names]
    # more advanced for the future: from the reactant string, do like string.index('numberfound')
    # and check the character that preceeds it to understand if it is a number of C, N, H, .. atoms
    numbers_tuples = [nums if isinstance(nums, tuple) else (0,) for nums in numbers_tuples]
    sums = [sum(list(map(int, nums))) for nums in numbers_tuples]
    
    return sums


def remove_fw_rxns(rxn_dct, rxn_lst):
    """ remove all rxns of rxn_list from rxn_dct
    """
    for rxn in rxn_lst:
        if rxn in rxn_dct.keys():
            rxn_dct.pop(rxn)

    return rxn_dct


def remove_rev_rxns(rxn_dct, rxn_lst):
    """
    take a rxn dictionary
    consider rxn_list: for each rxn,
    also search for its backward rxn and remove if present
    """
    for rxn in rxn_lst:
        rev_rxn = (rxn[1], rxn[0], (None,))
        if rev_rxn in rxn_dct.keys():
            rxn_dct.pop(rev_rxn)

    return rxn_dct


def resort_ktp_labels(ktp_dct):
    """
    takes ktp dct
    resort labels of reactants and products according to N of atoms, name length or alphabetical order
    - only way to guarantee consistency in dct keys when applying prompt
    - older version: only sort product labels; now sort both for consistency
    """

    new_ktp_dct = {}
    for key, val in ktp_dct.items():
        # prods = list(key[1])
        # prods.sort()
        # new_key = (key[0], tuple(prods), key[2])
        reacs = order_names(key[0])
        prods = order_names(key[1])
        new_key = (reacs, prods, key[2]) 
        new_ktp_dct[new_key] = copy.deepcopy(val)

    return new_ktp_dct


def count_atoms(fml_list):
    """ Count Cl, S, N, C, O, H atoms in formula list

        :param fml_list: stoich chemical formula
        :type fml_list: list of dictionaries [dict[str:int], ]
        :rtype: list [int, ], int = nClnSnNnCnOnH
    """
    fml_num_list = []
    for fml in fml_list:

        fml_num = (1e7*n_el(fml, 'Cl')+1e6*n_el(fml, 'S')+1e5*n_el(fml, 'N')
                   + 1e3*n_el(fml, 'C')+1e2*n_el(fml, 'O')+n_el(fml, 'H'))
        fml_num_list.append(fml_num)

    return fml_num_list


def extract_spc(spc):
    """
    extract species 1 from tuple
    """
    _spc1, _spc2 = [], []
    for _spc in spc:
        _spc1.append(_spc[0])
        if len(_spc) > 1:
            # bimol species
            _spc2.append(_spc[1])
        else:
            # unimol species
            _spc2.append('')

    return _spc1, _spc2


def get_mult(spc_tuple, spc_dct):
    """
    extracts the total multiplicity of the set of species
    spc_tuple = (A,B,C,..)
    spc_dct: species dictionary
    returns integer of the multiplicity
    """
    mult = 1
    if isinstance(spc_tuple, str):
        spc_tuple = tuple([spc_tuple])

    for spc in spc_tuple:
        mult *= spc_dct[spc]['mult']

    return mult


def get_fml(rxn_ichs):
    '''
    rxn_icn: inchis of the species of one side of a reaction (ich1, ich2, ..)
    returns: formula dictionary, formula string
    '''
    formula_dct = ''
    for rct_ich in rxn_ichs:
        formula_i_dct = formula(rct_ich)
        formula_dct = join(formula_dct, formula_i_dct)
    formula_str = string2(formula_dct)

    return formula_dct, formula_str
