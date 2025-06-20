"""
Functions for mechanism reading and sorting
"""

import sys
import autoparse.pattern as app
from ioformat import ptt
from ioformat import remove_comment_lines
from mechanalyzer.parser import ckin_ as ckin


def parse_mechanism(mech_str, mech_type):
    """ Get the reactions and species from the mechanism input
    """

    # Parse the info from the chemkin file
    if mech_type == 'chemkin':
        rxn_param_dct = ckin.parse_rxn_param_dct(mech_str)
    else:
        raise NotImplementedError

    return rxn_param_dct

# Parse the auxiliary file used to sort a mechanism
def parse_sort(sort_str):
    """ Parse the string from the sort.dat input file that contains various
        parameters used to sort a mechanism.

        Returns the list of species to isolate from the mechanism as well
        as the criteria used to sort the mechanism.

        :param sort_str: string for the sort.dat file
        :type sort_str: str
        :rtype: (list(str), list(str), int)
    """
    # remove comments
    sort_str = remove_comment_lines(
                sort_str, delim_pattern=app.escape('#'))
    sort_str = remove_comment_lines(
            sort_str, delim_pattern=app.escape('!'))

    # Read and format information from the isolate_submech block
    spc_block = ptt.end_block(sort_str, 'isolate_submech')

    if not spc_block: # this should be checked because spc_block might be None - section not mandatory
        spc_lst = []
    elif not str.isspace(spc_block):
        spc_lst = list(ptt.values_from_block(
            spc_block, val_ptt=app.one_or_more(app.CKINSAFE_CHAR)))
    else:
        spc_lst = []

    # Read and format information from the sort_mech block
    isol_block = ptt.end_block(sort_str, 'sort_mech')
    # main criteria
    crit_block = ptt.paren_blocks(
        isol_block, key='criteria')
    
    if crit_block:
        crit_tup = ptt.values_from_block(
            crit_block[0][1], val_ptt=app.one_or_more(app.URLSAFE_CHAR))
    else:
        crit_tup = ()
        
    head_block = ptt.keyword_value_blocks(
        isol_block, key='n_criteria_headers')
    nhead = int(head_block[0][1]) if head_block is not None else 0
    
    keepbelow = ptt.keyword_value_blocks(
        isol_block, key='stoich_keepbelow')
    if keepbelow is not None:
        spc_lst += ['keepbelow ' + keepbelow[0][1].strip(),]
    deleteabove = ptt.keyword_value_blocks(
        isol_block, key='stoich_deleteabove')
    if deleteabove is not None:
        spc_lst += ['deleteabove ' + deleteabove[0][1].strip(),]
    singlespecies = ptt.keyword_value_blocks(
        isol_block, key='singlespecies')
    if singlespecies is not None:
        if singlespecies[0][1].strip() == 'True':
            spc_lst += ['singlespecies']
    
    if keepbelow is not None and deleteabove is not None:
        raise ValueError('Cannot have both keepbelow and deleteabove criteria - incompatible!')
    
    sort_tup = crit_tup + (nhead,)
    sort_lst = list(sort_tup)

    # prompt criteria
    prompt_block = ptt.end_block(sort_str, 'prompt_filter')
    prompt_filter_dct = {}
    if prompt_block is not None:
        prompt_block = ptt.keyword_value_blocks(
            prompt_block)        
        dct_0 = dict(prompt_block)
        for key in dct_0.keys():
            prompt_filter_dct[key] = float(dct_0[key])

        
    # Print an error message for an isol block
    if isol_block is None:
        print('*ERROR: sort_mech section is not defined')

    return spc_lst, sort_lst, prompt_filter_dct

def parse_classtype(classtype_str):
    """ parse the class type string. useful to process larger class groups within a (CRECK) mech

    Args:
        classtype_str (str): string containing reaction types classified in larger groups
                             inline comments should be removed in advance.
        reaction_typeclass_dct (dct{reactiontype(str): classtype(str)}):
        dictionary with correspondence between reaction type and class type
    """

    grp_blocks = ptt.named_end_blocks(classtype_str, 'classtype', footer='classtype')
    reaction_typeclass_dct = {}
    for classtype, block in grp_blocks.items():
        reactiontypes = block.split()
        if reactiontypes:
            for reactiontype in reactiontypes:
                reaction_typeclass_dct[reactiontype] = classtype
            
    return reaction_typeclass_dct