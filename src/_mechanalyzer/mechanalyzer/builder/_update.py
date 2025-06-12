""" Functions which handle updating the master and temporary
    objects describing the reactions and species of a mechanism

    Need to add code to remove unstable species as reactants
    (going to bimol prods?)
"""

import itertools
import automol
from autoreact.params import RxnParams
import thermfit
from mechanalyzer.builder._names import rxn_ich_to_name
from mechanalyzer.builder._names import ich_name_dct
from mechanalyzer.builder._names import functional_group_name
from mechanalyzer.builder._names import stereo_name_suffix


# Handles Species Object Updates
def update_spc_dct_from_reactions(rxns, spc_dct, rename=False,
                                  enant_label=True, spc_orig_name_dct=None):
    """ Update a species with species from a set of reactions

        :param enant_label: Include the enantiomer label?
        :type enant_label: bool
    """

    spc_lst = _spc_from_reactions(rxns)
    spc_dct = update_spc_dct(spc_lst, spc_dct, rename=rename,
                             enant_label=enant_label,
                             spc_orig_name_dct=spc_orig_name_dct)

    return spc_dct


def update_spc_dct(spc_infos, spc_dct, rename=False, enant_label=True,
                   spc_orig_name_dct=None):
    """ Update the species dictionary with a list of species

        :param enant_label: Include the enantiomer label?
        :type enant_label: bool
    """
    spc_orig_name_dct = {} if spc_orig_name_dct is None else spc_orig_name_dct

    print('\nAdding new unique species to mechanism by',
          'adding to mechanism spc_dct...\n')

    # determine if charges and multiplicities should be considered
    # when making name dictionary
    has_inf = False
    if spc_infos:
        if not isinstance(spc_infos[0], str):
            has_inf = True

    _info_name_dct = ich_name_dct(spc_dct, incl_mult=has_inf, incl_chg=has_inf)

    # Add species dict to mech dct if it is not already in mechanism
    # Build a lst of species that have been added to the mechanism
    i = 0
    for info in spc_infos:
        if info not in _info_name_dct:
            # Generate a functional group name
            ich = info[0] if has_inf else info
            if not rename and spc_orig_name_dct:
                orig_name = spc_orig_name_dct[info]
                ste_lbl = stereo_name_suffix(ich, enant_label=enant_label)
                name = f'{orig_name}-{ste_lbl}' if ste_lbl else orig_name
                print('original name')
            else:
                print('new name')
                name = functional_group_name(ich, name='',
                                             enant_label=enant_label)
            print(f"InChI {ich} is giving name {name}")

            # Generate the data dct
            rgt_dct = thermfit.create_spec(ich)
            if has_inf:
                rgt_dct['charge'] = info[1]
                rgt_dct['mult'] = info[2]
            # Add to the overall mechanism spc_dct and new species lst
            smi = automol.chi.smiles(ich)
            print(f'Adding species {name} = {smi} = {ich}')

            if name in spc_dct:
                print("WARNING: GENERATED NAME ALREADY IN DCT!!!")
                print(f" - generated: {name} {ich}")
                print(f" - in dct   : {name} {spc_dct[name]['inchi']}")
                if has_inf:
                    if spc_dct[name]['mult'] != rgt_dct['mult']:
                        print('But it has a diff mult, renaming...')
                        mult_str = {1: 's', 2: 'd', 3: 't', 4: 'q'}
                        name = mult_str[rgt_dct['mult']] + name
                    if spc_dct[name]['charge'] != rgt_dct['charge']:
                        print('But it has a diff mult, renaming...')
                        chg_str = {-1: 'ani', 1: 'cat'}
                        name = chg_str[rgt_dct['charge']] + name

            spc_dct.update({name: rgt_dct})

    return spc_dct


# Handles Reaction Object Updates
def update_rxn_dct(rxn_lst, rxn_dct, spc_dct):
    """ Update the reaction dictionary with a list of reactions
    """

    print('\nAdding new unique reactions to mechanism...\n')

    rxn_dct = rxn_dct if rxn_dct is not None else {}
    for rxn in rxn_lst:
        rxn_wname = rxn_ich_to_name(rxn, spc_dct)
        if _unique_reaction(rxn_wname, rxn_dct):

            # Convert to names and print message
            print(f'Adding reaction {rxn_wname} to param dct')

            rxn_dct[rxn_wname] = RxnParams(
                arr_dct={'arr_tuples': ((1.0, 0.0, 0.0),)})

    return rxn_dct


# Removal functions
def remove_spc_not_in_reactions(rxn_param_dct, mech_spc_dct):
    """ Remove species from the spc dct not currently
        in the list of reactions in the rxn_dct
    """

    # spc_in_rxns = _spc_from_reactions(rxns)
    spc_in_rxns = ()
    for rxn in rxn_param_dct:
        spc_in_rxns += rxn[0]
        spc_in_rxns += rxn[1]
    spc_in_rxns = set(spc_in_rxns)

    new_mech_spc_dct = {}
    for name, dct in mech_spc_dct.items():
        if name in spc_in_rxns:
            new_mech_spc_dct[name] = dct
        elif any(x in name for x in ('cbh0_', 'cbh1_', 'cbh2_', 'cbh_3')):
            new_mech_spc_dct[name] = dct
        else:
            print(f'Remove species: {name}')

    return new_mech_spc_dct


def remove_improper_reactions(rxn_param_dct, mech_spc_dct,
                              stereo=True, reverse=True):
    """ Remove reactions from the mechanism that do not correspond
        to proper, physical elementary step reactions.
    """
    print('call improper')
    ste_rxn_param_dct = {}
    for rxn, params in rxn_param_dct.items():
        print(f'Checking {rxn[0]}->{rxn[1]}')
        rcts_ich = tuple(mech_spc_dct[rct]['inchi'] for rct in rxn[0])
        prds_ich = tuple(mech_spc_dct[prd]['inchi'] for prd in rxn[1])

        rxn_objs = automol.reac.from_chis(
            rcts_ich, prds_ich, stereo=stereo)
        if rxn_objs:
            rxn_class = automol.reac.class_(rxn_objs[0])
            print(f' - Keep: IDd {rxn_class} for reaction {rxn[0]}->{rxn[1]}')
            ste_rxn_param_dct[rxn] = params
        else:
            # Check if the reverse reaction cannot be ID'd
            if reverse:
                rxn_objs = automol.reac.from_chis(
                    prds_ich, rcts_ich, stereo=stereo)
                if rxn_objs:
                    rev_rxn = (rxn[1], rxn[0], rxn[2])
                    ste_rxn_param_dct[rev_rxn] = params
                    print(
                        f' - Keep: (Reverse) IDd {rxn_class} for '
                        f'reaction {rxn[0]}->{rxn[1]}')
                else:
                    # print(f'Removing reaction {rcts_ich}->{prds_ich}')
                    print(
                        ' - Remove: No ID in either direction for '
                        f'reaction {rxn[0]}->{rxn[1]}')
            else:
                # print(f'Removing reaction {rcts_ich}->{prds_ich}')
                print(f' - Remove: No ID for reaction {rxn[0]}->{rxn[1]}')

    return ste_rxn_param_dct


def remove_unstable_reactions(rxn_param_dct, mech_spc_dct):
    """ Remove reaction A -> B + C where A is an unstable reactant.
        This is because these reactions are readded with other codes.
    """
    _rxn_param_dct = {}
    for rxn, params in rxn_param_dct.items():
        reacs, prods, _ = rxn
        if len(reacs) == 1 and len(prods) == 2:
            rct_geo = automol.chi.geometry(mech_spc_dct[reacs[0]]['inchi'])
            rct_zma = automol.geom.zmatrix(rct_geo)
            instab_zmas = automol.reac.instability_product_zmas(rct_zma)
            if instab_zmas:
                print('Removing reaction', reacs, prods)
                _rxn_param_dct[rxn] = params
            else:
                _rxn_param_dct[rxn] = params
        else:
            _rxn_param_dct[rxn] = params

    return _rxn_param_dct


def _unique_reaction(rxn, rxn_dct):
    """ Determine if a reaction is in the parameter_dictionary

        does not deal with the third body, so function does not work
    """
    return not any(rxn in rxn_dct for rxn in _make_reaction_permutations(rxn))


# Other helper functions
def _spc_from_reactions(rxns):
    """ Build a species dictionary from a list of reactions
        which define a reaction using the inchi strings
    """
    rgts = ()
    for rxn in rxns:
        rgts += tuple(itertools.chain(*rxn))
    uni_rgts = tuple(
        rgt for rgt in
        automol.util.remove_duplicates_with_order(rgts)
        if rgt not in (None, '+M', '(+M)'))
    return uni_rgts


def _make_reaction_permutations(rxn):
    """ reactions
    """
    # Get all reactions R=P including all permutations of R and P
    rct_perms = tuple(itertools.permutations(rxn[0]))
    prd_perms = tuple(itertools.permutations(rxn[1]))

    all_rxns = (
        tuple(itertools.product(rct_perms, prd_perms)) +
        tuple(itertools.product(prd_perms, rct_perms))
    )

    # Remove duplicates
    all_rxns = tuple(set(all_rxns))

    # Re-add the third body
    third_body = rxn[2]
    all_rxns = tuple([(*rxn, third_body) for rxn in all_rxns])
    return all_rxns
