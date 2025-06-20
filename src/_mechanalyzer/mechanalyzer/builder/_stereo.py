"""
  Code to expand the mechanism using stereochemistry

    maybe remove all inchis that are not incomplete?
"""
import os
import itertools as it
from copy import deepcopy

import automol
from chemkin_io.writer._util import format_rxn_name
from autorun import execute_function_in_parallel
import mechanalyzer.parser
from mechanalyzer.builder._update import update_spc_dct_from_reactions
from mechanalyzer.builder._update import update_rxn_dct
from mechanalyzer.builder._names import rxn_name_str


# MAIN CALLABLE
def expand_mech_stereo(inp_mech_rxn_dct, inp_mech_spc_dct, nprocs='auto',
                       enant=True):
    """ Build list of stereochemistry to reactions

        :param enant: Include all enantiomers? Otherwise, includes only
            canonical enantiomer species and reactions.
        :type enant: bool

        Currently, we assume that the species in them mech_spc_dct have
        stereochemistry already added to them.
    """

    def _expand(name_ich_dct, rxns, output_queue):
        """ Expand reactons
        """
        # First loop over all reactions to add stereo
        all_ste_rxns = []
        for rxn in rxns:

            log1 = f'\nExpanding Stereo for Reaction: {format_rxn_name(rxn)}\n'
            print(log1)

            # Reformat reaction to use InChI instead of mechanism name
            # Split thrdbdy off, not needed for stereo code, add back later
            rxn_ich = _rxn_name_to_ich(rxn, name_ich_dct)
            _rxn_ich = (rxn_ich[0], rxn_ich[1])
            thrdbdy = rxn_ich[2]

            # Build list of all stereochemically allowed versions of reaction
            ste_rxns_lst, log2 = _ste_rxn_lsts(_rxn_ich, enant=enant)
            print(log2)
            # Appropriately format the reactions with third body
            ste_rxns_lst = _add_third(ste_rxns_lst, thrdbdy)
            all_ste_rxns.append(ste_rxns_lst)
        output_queue.put(tuple(all_ste_rxns))
        print(f'Processor {os.getpid()} finished')

    # Dictionaries to map name to inchi
    name_ich_dct = mechanalyzer.parser.spc.name_inchi_dct(inp_mech_spc_dct)

    # Generate all stereo reactons from the initial set
    rxns = tuple(inp_mech_rxn_dct.keys())
    args = (name_ich_dct,)
    pes_noste_rxns_dct = _rxns_noste_pes_dct(rxns, name_ich_dct)
    all_ste_rxns = ()

    # Loop over the PES (stoich similar)
    forms = list(pes_noste_rxns_dct.keys())
    # ste_rxn_cnt = 0
    for formula in forms:
        noste_rxns_dct = pes_noste_rxns_dct[formula]
        print('PES: {} has {:g} reactions'.format(
            formula, len(noste_rxns_dct.keys())))
        pes_gra, ccs_dct = _pes_gra(noste_rxns_dct)
        # Loop over ccs (connected channels)
        for _, rxns in ccs_dct.items():
            args = (name_ich_dct,)
            ste_rxns = [
                key for key, val in noste_rxns_dct.items() if val in rxns]
            # expand each reaction to all valid stoich versions
            ste_rxns = execute_function_in_parallel(
                _expand, ste_rxns, args, nprocs=nprocs)
            for rxn_lst in ste_rxns:
                all_ste_rxns += rxn_lst
    return all_ste_rxns


def expand_mech_stereo_debug(inp_mech_rxn_dct, inp_mech_spc_dct, enant=True):
    """ Build list of stereochemistry to reactions

        :param enant: Include all enantiomers? Otherwise, includes only
            canonical enantiomer species and reactions.
        :type enant: bool

        Currently, we assume that the species in them mech_spc_dct have
        stereochemistry already added to them.
    """

    # Dictionaries to map name to inchi
    name_ich_dct = mechanalyzer.parser.spc.name_inchi_dct(inp_mech_spc_dct)

    # Generate all stereo reactons from the initial set
    rxns = tuple(inp_mech_rxn_dct.keys())
    pes_noste_rxns_dct = _rxns_noste_pes_dct(rxns, name_ich_dct)
    all_ste_rxns = []
    failed_rxns = []

    spc_orig_name_dct = {}

    # Loop over the PES (stoich similar)
    forms = list(pes_noste_rxns_dct.keys())
    for formula in forms:
        noste_rxns_dct = pes_noste_rxns_dct[formula]
        print('PES: {} has {:g} reactions'.format(
            formula, len(noste_rxns_dct.keys())))
        pes_gra, ccs_dct = _pes_gra(noste_rxns_dct)
        # Loop over ccs (connected channels)
        for _, rxns in ccs_dct.items():
            noste_rxns = [
                key for key, val in noste_rxns_dct.items() if val in rxns]
            for rxn in noste_rxns:
                print(f'\nExpanding Stereo for Reaction: '
                      f'{format_rxn_name(rxn)}\n')

                # Reformat reaction to use InChI instead of mechanism name
                # Split thrdbdy off, not needed for stereo code, add back later
                rgts = rxn[0] + rxn[1]
                rxn_ich = _rxn_name_to_ich(rxn, name_ich_dct)
                _rxn_ich = (rxn_ich[0], rxn_ich[1])
                thrdbdy = rxn_ich[2]

                # Build list of all stereochemically allowed versions of
                # reaction
                try:
                    ste_rxns_lst, log2 = _ste_rxn_lsts(_rxn_ich, enant=enant)
                    rgt_ichs_lst = [r+p for r, p in ste_rxns_lst]
                    for rgt_ichs in rgt_ichs_lst:
                        spc_orig_name_dct.update(dict(zip(rgt_ichs, rgts)))
                    print(log2)
                    # Appropriately format the reactions with third body
                    ste_rxns_lst = _add_third(ste_rxns_lst, thrdbdy)
                    all_ste_rxns.extend(ste_rxns_lst)
                except Exception as e:
                    failed_rxns.append(rxn)
                    print(f'Reaction {format_rxn_name(rxn)} failed\n'
                          f'Error: {e}\n')

    return all_ste_rxns, spc_orig_name_dct, failed_rxns


def valid_enantiomerically(ste_mech_spc_dct):
    """ are there NOT enantiomeric species in this reaction list
    """
    rule_out = False
    for spc_name_i, spc_name_j in list(
            it.combinations(ste_mech_spc_dct.keys(), 2)):
        no_ste_ich_i = automol.chi.standard_form(
            ste_mech_spc_dct[spc_name_i]['inchi'], stereo=False)
        no_ste_ich_j = automol.chi.standard_form(
            ste_mech_spc_dct[spc_name_j]['inchi'], stereo=False)
        if no_ste_ich_i == no_ste_ich_j:
            if automol.chi.are_enantiomers(
                    ste_mech_spc_dct[spc_name_i]['inchi'],
                    ste_mech_spc_dct[spc_name_j]['inchi']):
                print(
                    'ruling out because',
                    ste_mech_spc_dct[spc_name_i]['inchi'],
                    ste_mech_spc_dct[spc_name_j]['inchi'])
                rule_out = True
                break
    return rule_out is False


def _make_ste_rxn_dct(ste_rxns):
    """  A dictionary with a stereoless reaction ichs as the
         key and stereo expansions ich reaction list as value
    """
    ste_rxn_dct = {}
    for rxn_lst in ste_rxns:
        ichs1, ichs2, _ = rxn_lst[0]
        noste_rxn, _, _ = _noste_rxn((ichs1, ichs2))
        ste_rxn_dct[noste_rxn] = rxn_lst
    return ste_rxn_dct


def _make_ccs_rxn_gra(ste_rxn_dct, pes_gra):
    """ merges the pes_gra and the reaction list
    """
    ccs_rxn_gra = {}
    for noste_rxn, ste_exp_lst in ste_rxn_dct.items():
        ccs_rxn_gra[noste_rxn] = (
            pes_gra[noste_rxn],
            ste_exp_lst)
    return ccs_rxn_gra


def _is_ste_conn(rxn, rxn_lst, prev_rxns_lst):
    """ Are to rxns connected stereochemically?
    """
    all_pes_ichs = []
    ichsi, ichsj, _ = rxn
    i_is_conn = False
    j_is_conn = False
    for rxn_i in rxn_lst:
        if rxn_i in prev_rxns_lst or prev_rxns_lst == []:
            ichs1, ichs2, _ = rxn_i
            all_pes_ichs.append(ichs1)
            all_pes_ichs.append(ichs2)
    if ichsi in all_pes_ichs:
        i_is_conn = True
    if ichsj in all_pes_ichs:
        j_is_conn = True
    return i_is_conn or j_is_conn or prev_rxns_lst == ['fake']


def _missed_rxn(rxn, rxn_lst, ignore_rxns):
    """ Are to rxns connected stereochemically?
    """
    all_pes_ichs = []
    ichsi, ichsj, _ = rxn
    enant_rxns = []
    for rxn_i in ignore_rxns:
        enant_rxns.append(_enant_rxn(rxn_i))
    enant_spcs = [ich for ich in rxn_i[:2] for rxn_i in enant_rxns]
    i_is_conn = False
    j_is_conn = False
    for rxn_i in rxn_lst:
        if rxn_i not in ignore_rxns:
            ichs1, ichs2, _ = rxn_i
            if not any(ich in enant_spcs for ich in ichs1 + ichs2):
                all_pes_ichs.append(
                    [ichs for ichs in ichs1 + ichs2])# if automol.inchi.is_enantiomer(ichs1)])
    if ichsi in all_pes_ichs:
        i_is_conn = True
    if ichsj in all_pes_ichs:
        j_is_conn = True
    return i_is_conn or j_is_conn


def _enant_rxn(rxn_i):
    """ Returns the enantiomer of a reaction,
        uses None if there is no enantiomer
    """
    check_ent = [
        tuple(map(automol.chi.reflect, ichs)) for ichs in rxn_i[:-1]]
    if check_ent[0] == rxn_i[0]:
        check_ent[0] = None
    if check_ent[1] == rxn_i[1]:
        check_ent[1] = None
    check_ent.append(rxn_i[-1])
    return check_ent


def _split_ste_ccs(ccs_rxn_gra):
    """ seperates a ccs into a dictionary of sccs
        where a numeric index is a key and a
        list of stereo-reactions (in ichs) is val
    """
    def _get_best_combo(exp_rxns_lst, maxcombo=None):
        nonste_rxns_lst = ()
        ste_rxns_lst = ()
        for rxn in exp_rxns_lst:
            enantiomer = False
            for ich in rxn[0] + rxn[1]:
                if automol.inchi.is_enantiomer(ich):
                    enantiomer = True
            if enantiomer:
                ste_rxns_lst += (rxn,)
            else:
                nonste_rxns_lst += (rxn,)
        rxn_combos = (0)
        if maxcombo is None:
            maxcombo = len(ste_rxns_lst)
        if maxcombo > 10 and maxcombo < 15:
            maxcombo = 12
        elif maxcombo > 14:
            maxcombo = 5
        print('setting r max', maxcombo)
        #for r in range(1, maxcombo+1, 1):# if len(ste_rxns_lst) % 6 == 0 else 2):
        for r in range(1, maxcombo, 1):# if  len(ste_rxns_lst) % 6 == 0 or len(ste_rxns_lst)==16 else 2):
            if rxn_combos == ():
                break
            print('combining r=', r)
            last_combo = deepcopy(rxn_combos)
            rxn_combos = ()
            for rxn_combo in it.combinations(ste_rxns_lst, r):
                good_combo = True
                for rxn_i in rxn_combo:
                    check_ent = _enant_rxn(rxn_i)
                    if _is_ste_conn(check_ent, rxn_combo, []):
                        good_combo = False
                        break
                    if _missed_rxn(rxn_i, ste_rxns_lst, rxn_combo):
                        good_combo = False
                        break
                if good_combo:
                    rxn_combos += (rxn_combo,)
        if rxn_combos != () and r < 4:
            last_combo = deepcopy(rxn_combos)
        return (combo + nonste_rxns_lst for combo in last_combo)

    def _recursive_step(
            noste_rxn, ccs_rxn_gra, old_sccs_rxn_gra,
            considered_rxns, prev_rxn_lst):
        adj_noste_rxn_lst, exp_rxns_lst = ccs_rxn_gra[noste_rxn]
        considered_rxns.append(noste_rxn)

        # Loop through reaction's expansions and put them
        # in appropriate sccs graph
        exp_rxns_lst = list(set(exp_rxns_lst))
        print(
            'next reaction in ccs has {:g} expansions'.format(len(exp_rxns_lst)))
        sccs_rxn_gra = {}

        maxcombo=None
        best_combos = _get_best_combo(exp_rxns_lst, maxcombo=maxcombo)
        print("WHAT")
        for rxn_set in  best_combos:
            print('rxn set', rxn_set)
            is_new_stereo = True

            # Loop over propogating sccs graphs
            for idx in old_sccs_rxn_gra:
                if not all(
                        _is_ste_conn(
                            rxn_i, old_sccs_rxn_gra[idx], prev_rxn_lst)
                        for rxn_i in rxn_set):
                    continue
                combo_allowed = True
                for rxn_i in rxn_set:
                    # check if this stereo-rxn is connected to the sccs
                    # if _is_ste_conn(rxn_i, old_sccs_rxn_gra[idx], prev_rxn_lst):
                    if True:
                        check_ent = _enant_rxn(rxn_i)
                        # check if this rxns enantiomers are also in the sccs
                        if _is_ste_conn(check_ent, old_sccs_rxn_gra[idx], []):
                            combo_allowed = False
                            break
                    else:
                        combo_allowed = False
                        break
                if combo_allowed:
                    sccs_rxn_gra[len(sccs_rxn_gra)] = deepcopy(old_sccs_rxn_gra[idx]) + rxn_set
                    is_new_stereo = False
            if is_new_stereo:
                sccs_rxn_gra[len(sccs_rxn_gra)] = rxn_set
                # sccs_rxn_gra[len(sccs_rxn_gra)] = (rxn_i,)

        print()
        print('resulting in the following set of sccss')
        for idx in sccs_rxn_gra.keys():
            print(idx, len(sccs_rxn_gra[idx]))
            for p in sccs_rxn_gra[idx]:
                print(p)
        # walk through adjacent reactions and repeat the procedure recursively
        # until all reactions are considered
        for noste_rxn_i in sorted(
                adj_noste_rxn_lst, key=lambda d: len(ccs_rxn_gra[d][0]) + len(ccs_rxn_gra[d][1]), reverse=False):
            # for noste_rxn_i in adj_noste_rxn_lst:
            if noste_rxn_i not in considered_rxns:
                sccs_rxn_gra, considered_rxns = _recursive_step(
                    noste_rxn_i, ccs_rxn_gra, sccs_rxn_gra,
                    considered_rxns, exp_rxns_lst)
        return sccs_rxn_gra, considered_rxns

    def _find_allowed_sccs_combos(sccs_rxn_gra, spc_key_dct, idxs_lst, output_queue):
        combo_lst = ()
        for (idx_i, idx_j) in idxs_lst:
            rxns_i = sccs_rxn_gra[idx_i]
            rxns_j = sccs_rxn_gra[idx_j]
            connects = False
            enantiomeric = False
            spcs = spc_key_dct[idx_i] + spc_key_dct[idx_j]
            spcs = list(set(spcs))
            for spc_a, spc_b in it.combinations(spcs, 2):
                if automol.inchi.are_enantiomers(spc_a, spc_b):
                    enantiomeric = True
                    break
            if enantiomeric:
                continue
            connects = True
            # for rxn_i in rxns_i:
            #     if _is_ste_conn(rxn_i, rxns_j, []):
            #         connects = True
            #         break
            if connects:
                print('combo attempt', idx_i, idx_j)
                combo_lst += ((idx_i, idx_j,),)
        output_queue.put(combo_lst)
    def _check_combo_grps(combo_set, combo_idxs_grp_lst, output_queue):
        full_combo_lst = ()
        for combo_idxs_lst in combo_idxs_grp_lst:
            #all_idxs = sorted(list(set(sum(combo_idxs_lst, ()))))
            all_idxs = list(combo_idxs_lst)
            print(all_idxs)
            if any(check_idxs not in combo_set for check_idxs in it.combinations(all_idxs, 2)):
                continue
            else:
                full_combo_lst += (tuple(all_idxs),)
        output_queue.put(full_combo_lst)
    def _add_idx_to_combo(combo_set, grp_idxs_lst, idxs, output_queue): 
        added_idxs_lst = ()
        print('idxs', idxs)
        for idx_i in idxs:
            print('doing idx_i', idx_i)
            for grp_idxs in grp_idxs_lst:
                if idx_i in grp_idxs:
                    continue
                all_idxs = sorted([idx_i, *grp_idxs])
                if any(check_idxs not in combo_set for check_idxs in it.combinations(all_idxs, 2)):
                    continue
                else:
                    added_idxs_lst += (tuple(all_idxs),)
        output_queue.put(added_idxs_lst)
    def _combine_sccss(sccs_rxn_gra, idxs_lst):
        new_sccs_rxn_gra = {}
        for i, idxs in enumerate(idxs_lst):
            rxns = ()
            for idx in idxs:
                rxns += sccs_rxn_gra[idx]
            new_sccs_rxn_gra[i] = tuple(set(rxns))
        return new_sccs_rxn_gra
    # initialize
    considered_rxns = []

    # pass through first reaction in a ccs
    start_key = list(ccs_rxn_gra.keys())[0]
    adj_noste_rxn_lst, exp_rxns_lst = ccs_rxn_gra[start_key]
    exp_rxns_lst = list(set(exp_rxns_lst))
    for key_i in ccs_rxn_gra:
        adj_noste_rxn_lst_i, exp_rxns_lst_i = ccs_rxn_gra[key_i]
        exp_rxns_lst_i = list(set(exp_rxns_lst_i))
        #if len(exp_rxns_lst_i) == 16:
        if len(exp_rxns_lst_i) < len(exp_rxns_lst):
            start_key = key_i
            adj_noste_rxn_lst = adj_noste_rxn_lst_i
            exp_rxns_lst = exp_rxns_lst_i
    #    elif len(exp_rxns_lst_i) == len(exp_rxns_lst) and len(adj_noste_rxn_lst_i) > len(adj_noste_rxn_lst):
    #        start_key = key_i
    #        adj_noste_rxn_lst = adj_noste_rxn_lst_i
    #        exp_rxns_lst = exp_rxns_lst_i
    sccs_rxn_gra = {}
    considered_rxns = []


    print('first reaction in ccs has {:g} expansions'.format(len(exp_rxns_lst)))
    for erxn in exp_rxns_lst:
        print(erxn)
    considered_rxns.append(start_key)
    sccs_rxn_gra = {}
    for idx, rxn_set in enumerate(_get_best_combo(exp_rxns_lst)):#, maxcombo=3)):
        sccs_rxn_gra[idx] = rxn_set

    print('resulting in the following set of sccss')
    for idx in sccs_rxn_gra.keys():
        print(idx, len(sccs_rxn_gra[idx]))
        for p in sccs_rxn_gra[idx]:
            print(p)
    # old_sccs_rxn_gra = {0: ()}
    # sccs_rxn_gra = {}

    # # Loop through reaction's expansions and put them
    # # in appropriate sccs graph
    # for rxn_i in exp_rxns_lst[1:]:
    #     is_new_stereo = True

    #     # Loop over propogating sccs graphs
    #     for idx in old_sccs_rxn_gra:

    #         # check if this stereo-rxn is connected to the sccs
    #         if _is_ste_conn(rxn_i, old_sccs_rxn_gra[idx], ['hmm']):

    #             # check if this rxns enantiomers are also in the sccs
    #             check_ent = _enant_rxn(rxn_i)
    #             if _is_ste_conn(check_ent, old_sccs_rxn_gra[idx], []):
    #                 continue

    #             sccs_rxn_gra[len(sccs_rxn_gra)] = deepcopy(old_sccs_rxn_gra[idx]) + (rxn_i,)
    #             is_new_stereo = False
    #     if is_new_stereo:
    #         # if it was just totally independent, create a new sccs
    #         sccs_rxn_gra[len(sccs_rxn_gra)] = (rxn_i,)
    # # walk through adjacent reactions and repeat the procedure recursively
    # # until all reactions are considered
    # print('first sccs rxn gra', sccs_rxn_gra)
    for noste_rxn_i in sorted(adj_noste_rxn_lst, key=lambda d: len(ccs_rxn_gra[d][0]) + len(ccs_rxn_gra[d][1]), reverse=False):
        # adj_noste_rxn_lst:
        if noste_rxn_i not in considered_rxns:
            sccs_rxn_gra, considered_rxns = _recursive_step(
                noste_rxn_i, ccs_rxn_gra, sccs_rxn_gra, considered_rxns,
                exp_rxns_lst)

    # get the species of each sccs
    spc_key_dct = {}
    for idx, rxns in sccs_rxn_gra.items():
        print(idx, len(sccs_rxn_gra[idx]))
        spcs = []
        for rcts, prds, _ in rxns:
            spcs.extend([rct for rct in rcts if rct not in spcs])
            spcs.extend([prd for prd in prds if prd not in spcs])
        spc_key_dct[idx] = spcs

    # make all possible combinations of sub_sccs
    # print('len of r=10', sum(1 for i in it.combinations(sccs_idxs_i, 10)))
    # print('len of r=9', sum(1 for i in it.combinations(sccs_idxs_i, 9)))
    # print('len of r=8', sum(1 for i in it.combinations(sccs_idxs_i, 8)))
    # print('len of r=7', sum(1 for i in it.combinations(sccs_idxs_i, 7)))
    #print('len of r=4', sum(1 for i in it.combinations(sccs_idxs_i, 4)))
    #print('len of r=3', sum(1 for i in it.combinations(sccs_idxs_i, 3)))
    sccs_idxs_i = list(sccs_rxn_gra.keys())
    idxs_lst = list(it.combinations(sccs_idxs_i, 2))
    args = (sccs_rxn_gra, spc_key_dct)
    allowed_combo_lst = execute_function_in_parallel(
        _find_allowed_sccs_combos, idxs_lst, args, nprocs=30)
    #two_idxs_lst = ()
    #for (two_idxs, rxn_lst, spcs) in new_sccs_rxn_lst:
    #    two_idxs = tuple(sorted(two_idxs))
    #    two_idxs_lst += (two_idxs,)
    #    sccs_rxn_gra[two_idxs] = rxn_lst
    #    spc_key_dct[two_idxs] = spcs
    #enant_lst = [idxs for idxs in idxs_lst if idxs not in two_idxs_lst]
    # enant_lst = set(idxs_lst) - set(two_idxs_lst)
    # three_idxs_lst = ()
    # for idx_i in sccs_idxs_i:
    #     for (idx_j, idx_k) in two_idxs_lst:
    #         if idx_i == idx_j and idx_i == idx_k:
    #             continue
    #         if tuple(sorted((idx_i, idx_j,))) in enant_lst:
    #             continue
    #         if tuple(sorted((idx_i, idx_k,))) in enant_lst:
    #             continue
    #     three_idxs_lst += ((idx_i, (idx_j, idx_k,),),)
    # print('three idxs', len(three_idxs_lst), three_idxs_lst)
    # args = (sccs_rxn_gra, spc_key_dct)
    # new_sccs_rxn_lst = execute_function_in_parallel(
    #     _combine_sccss, three_idxs_lst, args, nprocs=30)
    print('all combo list len', len(allowed_combo_lst))
    full_combo_lst = ()
    combo_set = frozenset(allowed_combo_lst)
    allow_combo_dct = {}
    for idx in sccs_idxs_i:
       allow_combo_dct[idx] = []
       for idx_i, idx_j in combo_set:
            if idx_i == idx:
                allow_combo_dct[idx] += [idx_j]
            elif idx_j == idx:
                allow_combo_dct[idx] += [idx_i]
    print('hi', allow_combo_dct)
   
    all_best_combos = ()
    for idx in allow_combo_dct:
        allow_idxs = allow_combo_dct[idx]
        allow_idxs = [allow_idx for allow_idx in allow_idxs if allow_idx > idx]
        print(idx, allow_idxs)
        better_r =  True
        if idx < 10:
            new_best_combos = ()
            for r in range(2, 8):#len(allow_idxs)):
                if not better_r:
                    break
                better_r = False
                best_combos = deepcopy(new_best_combos)
                new_best_combos = ()
                print('trying r of', r)
                combo_idxs_lst = it.combinations(allow_idxs, r)
                for combo_idxs in combo_idxs_lst:
                    all_idxs = sorted(combo_idxs)
                    if any(check_idxs not in combo_set for check_idxs in it.combinations(all_idxs, 2)):
                        continue
                    else:
                        better_r = True
                        new_best_combos += ((idx,) + tuple(all_idxs),)
        else:
            best_combos = ()
            combo_idxs_lst = it.combinations(allow_idxs, r)
            for combo_idxs in combo_idxs_lst:
                all_idxs = sorted(combo_idxs)
                if any(check_idxs not in combo_set for check_idxs in it.combinations(all_idxs, 2)):
                    continue
                else:
                    better_r = True
                    best_combos += ((idx,) + tuple(all_idxs),)
        all_best_combos += best_combos    
        print('best idxs', idx, best_combos) 



    #combo_idxs_grp_lst = list(it.combinations(allowed_combo_lst, 3))
    # combo_idxs_grp_lst = list(it.combinations(sccs_idxs_i, 3))
    # print('checking combo groups for r=3', len(combo_idxs_grp_lst))
    # args = (combo_set,)
    # full_combo_lst = execute_function_in_parallel(
    #     _check_combo_grps, combo_idxs_grp_lst, args, nprocs=30)
    # three_idxs_lst = ()
    # args = (combo_set, combo_set,)
    # print(sccs_idxs_i)
    # three_idxs_lst = execute_function_in_parallel(
    #     _add_idx_to_combo, sccs_idxs_i, args, nprocs=30)
    # print('check length', len(three_idxs_lst))
    # args = (combo_set, three_idxs_lst,)
    # more_idxs_lst = execute_function_in_parallel(
    #     _add_idx_to_combo, sccs_idxs_i, args, nprocs=30)
    # print('check length', len(more_idxs_lst))
    # args = (combo_set, more_idxs_lst,)
    # more_idxs_lst = execute_function_in_parallel(
    #     _add_idx_to_combo, sccs_idxs_i, args, nprocs=30)
    # print('check length', len(more_idxs_lst))
    # args = (combo_set, more_idxs_lst,)
    # more_idxs_lst = execute_function_in_parallel(
    #     _add_idx_to_combo, sccs_idxs_i, args, nprocs=30)
    # print('check length', len(more_idxs_lst))
    # #for idx_k in sccs_idxs_i:
    # #    print('doing 3s for ', idx_k)
    # #    for idx_i, idx_j in combo_set:
    # #        if idx_i == idx_j and idx_i == idx_k:
    # #            continue
    # #        if not tuple(sorted((idx_i, idx_k,))) in combo_set:
    # #            continue
    # #        if not tuple(sorted((idx_j, idx_k,))) in combo_set:
    # #            continue
    # #        three_idxs_lst += ((idx_i, idx_j, idx_k,),)
    # #four_idxs_lst = ()
    # #for idx_l in sccs_idxs_i:
    # #    print('doing 4s for ', idx_l)
    # #    for idx_i, idx_j, idx_k in three_idxs_lst:
    # #        if idx_l in [idx_i, idx_j, idx_k]:
    # #            continue
    # #        all_idxs = sorted([idx_i, idx_j, idx_k, idx_l])
    # #        if any(check_idxs not in combo_set for check_idxs in it.combinations(all_idxs, 2)):
    # #            continue
    # #        else:
    # #            four_idxs_lst += (tuple(all_idxs),)
    # #print('check length', len(four_idxs_lst))
    #    
    # #print(len(full_combo_lst))
    more_idxs_lst = all_best_combos
    sccs_rxn_gra = _combine_sccss(
        sccs_rxn_gra, more_idxs_lst)
    sccs_idxs_i = list(sccs_rxn_gra.keys())
    idxs_lst = list(it.combinations(sccs_idxs_i, 2))
    spc_key_dct = {}
    for idx, rxns in sccs_rxn_gra.items():
        print(idx, len(sccs_rxn_gra[idx]))
        spcs = []
        for rcts, prds, _ in rxns:
            spcs.extend([rct for rct in rcts if rct not in spcs])
            spcs.extend([prd for prd in prds if prd not in spcs])
        spc_key_dct[idx] = spcs
    # make all possible combinations of sub_sccs
    args = (sccs_rxn_gra, spc_key_dct)
    allowed_combo_lst = execute_function_in_parallel(
        _find_allowed_sccs_combos, idxs_lst, args, nprocs=30)
    print('any further combos?', allowed_combo_lst)
    #two_idxs_lst = ()
    #sccs_rxn_gra, full_combo_lst)
    #print('parallel list')
    #print(len(new_sccs_rxn_lst))
    # for idx_i in sccs_idxs_i:
    #     sccs_idxs_j = list(sccs_rxn_gra.keys())
    #     for idx_j in sccs_idxs_j:
    #         if idx_i >= idx_j:
    #             continue
    #         print('combo attempt', idx_i, idx_j)
    #         connects = False
    #         rxns_i = sccs_rxn_gra[idx_i]
    #         rxns_j = sccs_rxn_gra[idx_j]
    #         enantiomeric = False
    #         spcs = spc_key_dct[idx_i] + spc_key_dct[idx_j]
    #         for spc_a, spc_b in it.combinations(spcs, 2):
    #             if automol.inchi.are_enantiomers(spc_a, spc_b):
    #                 enantiomeric = True
    #                 break
    #         if enantiomeric:
    #             continue
    #         for rxn_i in rxns_i:
    #             if _is_ste_conn(rxn_i, rxns_j, []):
    #                 connects = True
    #                 break
    #         if connects:
    #             print('new length', len(rxns_i + rxns_j))
    #             sccs_rxn_gra[len(sccs_rxn_gra)] = tuple(set(sorted(rxns_i + rxns_j)))
    #             spc_key_dct[len(sccs_rxn_gra)-1] = spc_key_dct[idx_i] + spc_key_dct[idx_j]
    print('BEST NUMBER OF SCCS')
    max_len = max([len(lst) for lst in sccs_rxn_gra.values()])
    keep_sccs_rxn_gra = {}
    new_idx = 0
    for idx in sccs_rxn_gra.keys():
        if len(sccs_rxn_gra[idx]) == max_len:
            if tuple(sorted(sccs_rxn_gra[idx])) not in keep_sccs_rxn_gra.values():
                keep_sccs_rxn_gra[new_idx] = tuple(sorted(sccs_rxn_gra[idx]))
                new_idx += 1
    for idx, rxns in keep_sccs_rxn_gra.items():
        print(idx, len(rxns))
    return keep_sccs_rxn_gra


def _pes_gra(noste_rxn_dct):
    """ seperates a list of reactions into graphs
        pes_gra: graph of a reaction (key)
            and a list (value) of the reactions with the same stoichiometry
        ccs_gra: graph of reaction (key
            and list (value) of reactions that are connected through wells
    """
    def _recursive_add(rxna, pes_gra, sub_pes_gra_i):
        for rxnb in pes_gra[rxna]:
            if rxnb in sub_pes_gra_i:
                continue
            sub_pes_gra_i += (rxnb,)
            sub_pes_gra_i = _recursive_add(rxnb, pes_gra, sub_pes_gra_i)
        return sub_pes_gra_i

    # initialize
    pes_gra = {}
    ccs_gra = {}
    idx = 0

    # loop over reactions
    for rxna in noste_rxn_dct:
        ichsi, ichsj = noste_rxn_dct[rxna]
        for rxnb in noste_rxn_dct:
            if rxna == rxnb:
                continue
            noste_rxnb = noste_rxn_dct[rxnb]
            if (
                    (len(ichsi) == 1 and ichsi in noste_rxnb) or
                    (len(ichsj) == 1 and ichsj in noste_rxnb)):
                if (ichsi, ichsj) not in pes_gra:
                    pes_gra[(ichsi, ichsj)] = (noste_rxnb,)
                else:
                    pes_gra[(ichsi, ichsj)] += (noste_rxnb,)
        if (ichsi, ichsj) not in pes_gra:
            print(rxna, ' IS DISCONNECTED')
            pes_gra[(ichsi, ichsj)] = ()

    # loops over PESes
    for rxna in pes_gra:
        if idx not in ccs_gra:
            ccs_gra[idx] = (rxna,)
        elif not any(
                rxna in ccs_gra[idx_j] for idx_j in ccs_gra):
            idx += 1
            ccs_gra[idx] = (rxna,)
        else:
            continue
        ccs_gra[idx] = _recursive_add(rxna, pes_gra, ccs_gra[idx])

    for idx in ccs_gra:
        print('CCS {:g}: {:g} reactions'.format(idx, len(ccs_gra[idx])))

    return pes_gra, ccs_gra


def _noste_rxn(rxn_ichs):
    ichs1, ichs2 = rxn_ichs
    noste_ichs1 = ()
    noste_ichs2 = ()
    form1 = ()
    form2 = ()
    for ich in ichs1:
        noste_ichs1 += (automol.chi.standard_form(ich, stereo=False),)
        form1 += (automol.chi.formula(ich),)
    for ich in ichs2:
        noste_ichs2 += (automol.chi.standard_form(ich, stereo=False),)
        form2 += (automol.chi.formula(ich),)
    form1 = (automol.form.join_sequence(form1) if len(form1) > 1
             else form1[0])
    form2 = (automol.form.join_sequence(form2) if len(form2) > 1
             else form2[0])
    noste_rxn = (automol.chi.sorted_(noste_ichs1),
                 automol.chi.sorted_(noste_ichs2))
    return noste_rxn, form1, form2


def _rxns_noste_pes_dct(rxns, name_ich_dct):
    """Dictionary of dictionaries for FORMULA: {origanalRXN: RXnwithoutstereo}
    """
    noste_dct = {}
    for rxn in rxns:
        ichs1, ichs2, _ = _rxn_name_to_ich(rxn, name_ich_dct)
        noste_rxn, form1, form2 = _noste_rxn((ichs1, ichs2,))
        if not automol.form.string(form1) in noste_dct:
            noste_dct[automol.form.string(form1)] = {rxn: noste_rxn}
        else:
            noste_dct[automol.form.string(form1)][rxn] = noste_rxn
    return noste_dct


def _sort_expansion(all_ste_rxns):
    rxn_ich_count = {}
    for rxn_lst in all_ste_rxns:
        for rxn in rxn_lst:
            ichs1, ichs2, _ = rxn
            for ich in ichs1 + ichs2:
                ich_no_ste = automol.chi.standard_form(ich, stereo=False)
                if ich_no_ste not in rxn_ich_count:
                    rxn_ich_count[ich_no_ste] = 1
                elif ich_no_ste != ich:
                    rxn_ich_count[ich_no_ste] += 1
    sort_val_lst = []
    for rxn_lst in all_ste_rxns:
        sort_val = 0
        for rxn in rxn_lst:
            ichs1, ichs2, _ = rxn
            for ich in ichs1 + ichs2:
                ich_no_ste = automol.chi.standard_form(ich, stereo=False)
                sort_val += rxn_ich_count[ich_no_ste]
        sort_val_lst.append(sort_val)

    return [x for _, x in sorted(zip(sort_val_lst, all_ste_rxns), reverse=True)]


def remove_stereochemistry(inp_mech_rxn_dct, inp_mech_spc_dct):
    """ Generate a mechanism with all stereochemistry removed
    """

    print('Removing stereochemistry from the species and reactions')

    # Loop over the reactions and generate the variants without stereo
    name_ich_dct = mechanalyzer.parser.spc.name_inchi_dct(inp_mech_spc_dct)

    noste_rxns = ()
    for rxn in inp_mech_rxn_dct:

        # Write rxn in terms of inchi, then remove the inchi strings
        rxn_ich = _rxn_name_to_ich(rxn, name_ich_dct)
        rxn_ich_noste = _remove_rxn_stereo(rxn_ich)

        if rxn_ich_noste not in noste_rxns:
            noste_rxns += (rxn_ich_noste,)

    # Update the mechanism objects with unique spc and rxns
    noste_spc_dct, noste_rxn_dct = {}, {}
    noste_spc_dct = update_spc_dct_from_reactions(noste_rxns, noste_spc_dct)
    noste_rxn_dct = update_rxn_dct(noste_rxns, noste_rxn_dct, noste_spc_dct)

    return noste_rxn_dct, noste_spc_dct


# Build reaction lists
def _ste_rxn_lsts(rxn_ich, enant=True):
    """ Build reaction onjects

        :param enant: Include all enantiomers? Otherwise, includes only
            canonical enantiomer species and reactions.
        :type enant: bool

    """
    # Build reaction objects
    rxn_obj_sets = automol.reac.from_chis(
        rxn_ich[0], rxn_ich[1], stereo=False)
    try:
        rxn_obj = rxn_obj_sets[0]  # expand just with rxn object
    except (TypeError, IndexError):
        print('No ID', rxn_ich, rxn_obj_sets)

    # Build a list of stereo reactions
    ste_rxn_ichs = ()
    for ste_rxn in automol.reac.expand_stereo(rxn_obj, enant=enant):
        attempt = 1
        while attempt < 4:
            try:
                rct_ichs, prd_ichs = automol.reac.chis(ste_rxn, stereo=True)
                ste_rxn_ichs += ((rct_ichs, prd_ichs),)
                break
            except:
                attempt += 1

            if attempt == 3:
                print('Fail to get stereo in 3 attempts', rxn_ich)

    # Set log message
    log = f' - Reaction identified as {automol.reac.class_(rxn_obj)}.\n'

    return ste_rxn_ichs, log


# Functions to check and sort the reactions by stereochemistry
def _remove_enantiomer_reactions(ste_rxn_lst, reacs_stereo_inchi=None):
    """ Take all reactions that occur from stereochemically
        and determine which reactions should be kept and
        whihc are unneccessary

        There are two reduction methods to reduce the set.
        If the reactant inchi is given, we grab reactions that use that
        stereo. Otherwise, we use internal logic in autochem to enfore
        m0 stereochemistry in the InChI strings.
    """

    # Convert reactants stero inchi to set for comparisons
    reacs_stereo_inchi = set(reacs_stereo_inchi)

    # Remove redundant sets and rebuild proper list
    if reacs_stereo_inchi is not None:
        log = ' - Reducing reactions to those with reactant stereochemistry\n'
        # Checks the InChI of the reactants in each reaction to see if they
        # match the input stereo inchi
        f_ste_rxn_lst = tuple(rxn for rxn in ste_rxn_lst
                              if set(rxn[0]) == reacs_stereo_inchi)
    else:
        log = ' - Reducing reactions to enforce InChI/m0 stereo throughout\n'
        f_ste_rxn_lst = automol.chi.filter_enantiomer_reactions(ste_rxn_lst)

    # Print the removed reactions
    removed_ste_rxn_lst = set(ste_rxn_lst) - set(f_ste_rxn_lst)

    return f_ste_rxn_lst, removed_ste_rxn_lst, log


# Diastereomer Abstraction Code
def diastereomer_abstractions(sccs_rxn_dct_lst, ccs_sccs_spc_dct,
                              chosen_idx_lst, all_chosen_ichs):
    """ Get additional (CCS, S-CCS) index pairs for diastereomer S-CCS
        that were missed in the initial selection process.

        Assuming this only occurs for abstractions of diastereomers.
    """

    dias_chosen_idx_lst = ()
    for ccs_idx, sccs_idx in chosen_idx_lst:
        if _ccs_is_abstraction(sccs_rxn_dct_lst[ccs_idx]):

            # Get indices for S-CCS containing diastereomers of chosen
            dias_sccs_idxs = _diastereomer_sccs_idxs(
                sccs_rxn_dct_lst, ccs_sccs_spc_dct,
                ccs_idx, sccs_idx,
                all_chosen_ichs)

            # Add indices to final list
            dias_chosen_idx_lst += dias_sccs_idxs
            # print('dias_sccs_idxs', ccs_idx, sccs_idx, dias_chosen_idx_lst)
            # print('is abs', ccs_idx, sccs_idx)
        else:
            print('not abs', ccs_idx, sccs_idx)

    return dias_chosen_idx_lst


def _ccs_is_abstraction(sccs_rxn_dct):
    """ ID if it is an abstraction by looking at reactions of each CCS.
    """

    # Check if each S-CCS is composed of a single reaction
    # Assume if first S-CCS is an abstraction, they all will be
    # Assumes the expansion is correct
    is_abstraction = False
    # print('a1', all(len(rxn_lst) == 1 for rxn_lst in sccs_rxn_dct.values()))
    # print('sccs rxn dct')
    # for x, y in sccs_rxn_dct.items():
    #     print(x, y)
    if all(len(rxn_lst) == 1 for rxn_lst in sccs_rxn_dct.values()):
        for sccs_rxn_lst in sccs_rxn_dct.values():
            # print('SCCS rxn lst TEST', sccs_rxn_lst)
            rxn_objs = automol.reac.from_chis(
                sccs_rxn_lst[0][0], sccs_rxn_lst[0][1])
            # print('rxn obj test', bool(rxn_obj))
            if rxn_objs:
                if automol.reac.class_(rxn_objs[0]) == 'hydrogen abstraction':
                    is_abstraction = True
                    break

    return is_abstraction


def _diastereomer_sccs_idxs(sccs_rxn_dct_lst, ccs_sccs_spc_dct,
                            chosen_ccs_idx, chosen_sccs_idx,
                            all_chosen_ichs):
    """ Find idxs for S-CCS that are diasteromer reactions to chosen S-CCS
    """

    # Get the S-CCSs for the CCS
    sccs_rxn_dct = sccs_rxn_dct_lst[chosen_ccs_idx]

    # Get reactants of reaction on the initially chosen S-CCS
    chosen_sccs_rxn = sccs_rxn_dct[chosen_sccs_idx][0]
    chosen_sccs_rxn = (chosen_sccs_rxn[0], chosen_sccs_rxn[1])

    # Find reactions on each S-CCS that diastereoisomeric with chosen S-CCS
    dias_sccs_idxs, dias_rxn_lst = (), ()
    for sccs_idx, sccs_rxn_lst in sccs_rxn_dct.items():
        if sccs_idx != chosen_sccs_idx:
            sccs_rxn = (sccs_rxn_lst[0][0], sccs_rxn_lst[0][1])
            if _possible_diastereoisomeric_reaction(chosen_sccs_rxn, sccs_rxn):
                dias_sccs_idxs += ((chosen_ccs_idx, sccs_idx),)
                dias_rxn_lst += (sccs_rxn[0] + sccs_rxn[1],)

    # Sort diastereoisomeric S-CCSs and corresponding species by overlap
    overlap = ()
    for (ccs_idx, sccs_idx) in dias_sccs_idxs:
        spc_dct = ccs_sccs_spc_dct[ccs_idx][sccs_idx]
        overlap += (len(set(all_chosen_ichs) & set(spc_dct)),)

    _sorter = sorted(zip(overlap, dias_sccs_idxs, dias_rxn_lst), reverse=True)
    dias_sccs_idxs = tuple(idx for _, idx, _ in _sorter)
    dias_rxn_lst = tuple(rxn for _, _, rxn in _sorter)

    # Loop through diastereomer S-CCSs (start w/ max overlap) & grab any S-CCS
    # that has a species that is not enantiomeric with chosen ichs
    # e.g., if R,R first, grab R,S dias then ignore S,S and S,R enantiomers
    final_dias_sccs_idxs, final_dias_ich_lst = (), ()
    for (ccs_idx, sccs_idx), dias_rxn in zip(dias_sccs_idxs, dias_rxn_lst):
        # Assess if any species in the reaction are enantiomers with the
        # the species we have chosen to maintain from all S-CCSs
        is_enant = False
        for dias_ich in dias_rxn:
            if any(automol.chi.are_enantiomers(dias_ich, ich)
                   for ich in final_dias_ich_lst):
                is_enant = True
            else:
                final_dias_ich_lst += (dias_ich,)

        if not is_enant:
            final_dias_sccs_idxs += ((ccs_idx, sccs_idx),)

    return final_dias_sccs_idxs


def _possible_diastereoisomeric_reaction(rxn_a, rxn_b):
    """ Assess if two reactions contain a diastereomer
    """

    has_dias = False
    for side_idx in (0, 1):
        for rgt_a, rgt_b in zip(rxn_a[side_idx], rxn_b[side_idx]):
            if automol.chi.are_diastereomers(rgt_a, rgt_b):
                has_dias = True
                break
        if has_dias:
            break

    return has_dias


# Formatters and printers
def _stereo_results(rxn, f_ste_rxns_lst, removed_ste_rxns_lst):
    """ Print the final filtered reactions and those removed
    """

    # Print final list of reactions
    log = f' - Stereochemical Versions of Reaction: {rxn_name_str(rxn)}\n'
    for ste_rxn in f_ste_rxns_lst:
        log += '    ' + rxn_name_str(ste_rxn, newline=True) + '\n'

    # Print removed reactions
    if removed_ste_rxns_lst:
        log += (' - Redundant, enantiomeric reactions '
                'precluded from final list\n')
        for ste_rxn in removed_ste_rxns_lst:
            log += '    ' + rxn_name_str(ste_rxn, newline=True) + '\n'

    return log


def _add_third(rxn_lst, thrdbdy):
    """ Format a rxn list to have the third-body added back
    """
    return tuple((rxn[0], rxn[1], thrdbdy) for rxn in rxn_lst)


def _rxn_name_to_ich(rxn, ich_dct):
    """ Convert a reacion written with spc names to spc inchis
        Third body list remains the same
    """

    # Convert reactant and product names to InChIs
    _rxn = (
        tuple(ich_dct.get(rgt) for rgt in rxn[0]),
        tuple(ich_dct.get(rgt) for rgt in rxn[1]),
        rxn[2]
    )

    # Set rxn_ich to None
    if (
        any(rgt is None for rgt in _rxn[0]) or
        any(rgt is None for rgt in _rxn[1])
    ):
        _rxn = None
    if _rxn is None:
        print('we got a none')
        print(rxn)
        print([ich_dct.get(rgt) for rgt in rxn[0]])
        print([ich_dct.get(rgt) for rgt in rxn[1]])
        print(ich_dct)
    return _rxn


def _remove_rxn_stereo(rxn):
    """ Generate rxn in inchi representation with no stereo
    """

    return (
        tuple(automol.chi.standard_form(ich, stereo=False) for ich in rxn[0]),
        tuple(automol.chi.standard_form(ich, stereo=False) for ich in rxn[1]),
        rxn[2]
    )


def _rxn_smiles(rxn):
    """ write a reaction into smles
    """
    return (
        tuple(automol.chi.smiles(rgt) for rgt in rxn[0]),
        tuple(automol.chi.smiles(rgt) for rgt in rxn[1]),
    )
