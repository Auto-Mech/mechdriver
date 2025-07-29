"""
functions called in the sorter
to be sorted . . .
"""
import numpy
import copy
import pandas as pd
import numpy as np
from mechanalyzer.calculator import formulas
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

def get_smallest_stoich_rxn(reactions, spc_dct):
    """ get reaction with smallest stoichiometry based on formulas in species dictionary

    Args:
        reactions (list(tuples)): [((rct1, rct2,),(prd1,),((+M),)), ...]
        spc_dct (dict): species dictionary {spc: {'fml':, 'inchi':,..}}
    """
    # choose reference reaction as that for the smallest species
    tot_stoich_dct = {reaction: formulas.sum_heavy_atoms_fmls([spc_dct[rct]['fml'] for rct in reaction[0]])
                             for reaction in reactions}
    min_val = min(tot_stoich_dct.values())
    return [k for k, v in tot_stoich_dct.items() if v == min_val]

# functions for ktp dictionaries moved to rates.py

