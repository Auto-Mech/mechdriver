import copy
from mechanalyzer.calculator import compare
from mechanalyzer.calculator import rates
from mechanalyzer.calculator.rates import check_p_t
from mechanalyzer.builder import _names as names
from ratefit.fit import _fit as fit
from automol import inchi
from automol import chi
from autoreact.params import RxnParams


def main(rxn_param_dct, spc_nasa7_dct, mech_spc_dct, temp_lst, pressures,
         dummy=False):

    # Get isomer sets, racemic_sets, and racemic rxn_param_dct
    iso_sets = find_iso_sets(mech_spc_dct)
    rac_sets, rac_names, rac_mech_spc_dct = get_rac_sets(
        iso_sets, mech_spc_dct)
    rac_rxn_param_dct = get_rac_rxn_param_dct(
        rac_sets, rac_names, rxn_param_dct)

    # Get the lumped rxn parameter dictionary
    lump_rxn_param_dct = lump(rac_rxn_param_dct, temp_lst, pressures,
                              dummy=dummy)

    # Get the thermo dct with only racemized species names
    rac_spc_nasa7_dct = get_rac_spc_nasa7_dct(rac_names, spc_nasa7_dct)

    return lump_rxn_param_dct, rac_spc_nasa7_dct, rac_mech_spc_dct


def lump(rac_rxn_param_dct, temps_lst, pressures, dummy=False):

    lump_rxn_param_dct = {}
    for rac_rxn, [rxns, params_lst] in rac_rxn_param_dct.items():
        if dummy:
            comb_params = RxnParams(arr_dct={'arr_tuples': ((1, 0, 0),)})
        else:
            lump_dct = {}
            for idx, rxn in enumerate(rxns):  # rxns in terms of spc names
                rcts, prds, _ = rxn
                params = params_lst[idx]
                sort_rcts = tuple(sorted(rcts))  # alphabetize
                if sort_rcts in lump_dct:
                    lump_dct[sort_rcts].append(params)
                else:
                    lump_dct[sort_rcts] = [params]
    
            # Get the list of ktp_dcts
            ktp_dct_lst = []
            for params_lst in lump_dct.values():
                for params in params_lst:
                    ktp_dct = rates.eval_params(params, temps_lst, pressures)
                    ktp_dct_lst.append(ktp_dct)
    
            # Add all the ktp_dcts together, then divide by averaging factor
            comb_ktp_dct = _sum_ktp_list(ktp_dct_lst)
            comb_ktp_dct = rates.mult_by_factor(comb_ktp_dct, 1 / len(lump_dct))
    
            # Refit the final ktp_dct
            try:
                comb_params, _ = fit.fit_ktp_dct(comb_ktp_dct, 'plog')
            except:
                comb_params = RxnParams(arr_dct={'arr_tuples': ((1e-5, 0, 0),)})
        # Store rates
        lump_rxn_param_dct[rac_rxn] = comb_params

    return lump_rxn_param_dct


def get_rac_spc_nasa7_dct(rac_names, spc_nasa7_dct):
    """ Takes a full, redundant NASA-7 thermo dict and simplifies to only the
        racemic species names
    """
    def swap_enant_name(spc):
        new_spc = spc
        if spc[-1] == '0':
            new_spc = spc[:-1] + '1'
        elif spc[-1] == '1':
            new_spc = spc[:-1] + '0'

        return new_spc

    print('Starting thermo selection')
    rac_spc_nasa7_dct = {}
    for rac_name in rac_names:
        if spc_nasa7_dct.get(rac_name) is not None:  # if spc in orig. thermo
            rac_spc_nasa7_dct[rac_name] = spc_nasa7_dct[rac_name]
        else:  # if spc not in original thermo, check for enantiomer
            other_rac = swap_enant_name(rac_name) 
            if spc_nasa7_dct.get(other_rac) is not None:
                # Note: store under the original racemic name
                rac_spc_nasa7_dct[rac_name] = spc_nasa7_dct[other_rac]
            else:
                print(f'neither {rac_name} nor its enant, {other_rac}, '
                      'in thermo')

    return rac_spc_nasa7_dct


def get_rac_rxn_param_dct(rac_sets, rac_names, rxn_param_dct):
    """ Note: this does not return a regular rxn_param_dct

        Has the form {rxn: [[rxn1, rxn2, ...], [params1, params2, ...]], ...}
        Keys are the racemized rxn names, while the values are the multiple
        reactions refering to the stereo-specific reactions
        couched under each racemized reaction name
    """

    def get_spc_names(rcts_or_prds, rac_sets, rac_names):
        spc_names = []
        error = False
        for rct_or_prd in rcts_or_prds:
            for rac_idx, rac_set in enumerate(rac_sets):
                if rct_or_prd in rac_set:
                    rac_name = rac_names[rac_idx]
                    spc_names.append(rac_name)
                    break
                # If on last rac_set and didn't break, the spc is missing
                if rac_idx == len(rac_sets) - 1:
                    print(f'Warning: species {rct_or_prd} not in rac_sets!')
                    error = True
        spc_names = sorted(spc_names)
        spc_names = tuple(spc_names)

        return spc_names, error

    rac_rxn_param_dct = {}
    for rxn, params in rxn_param_dct.items():
        rcts, prds, thrd_bod = rxn
        rct_names, rct_err = get_spc_names(rcts, rac_sets, rac_names)
        prd_names, prd_err = get_spc_names(prds, rac_sets, rac_names)
        # If any species missing from rac_sets, don't add current rxn
        if rct_err or prd_err:
            continue  # skip to next rxn
        rac_rxn_name = (rct_names, prd_names, thrd_bod)
        # Either append current vals or make new entry with them
        if rac_rxn_name in rac_rxn_param_dct:
            rac_rxn_param_dct[rac_rxn_name][0].append(rxn)
            rac_rxn_param_dct[rac_rxn_name][1].append(params)
        else:
            rac_rxn_param_dct[rac_rxn_name] = [[rxn], [params]]

    return rac_rxn_param_dct


def get_rac_sets(iso_sets, mech_spc_dct):

    rac_sets = []
    rac_names = []
    rac_mech_spc_dct = {}
    for iso_set in iso_sets:
        rac_dct = {}
        for iso in iso_set:
            spc_dct = mech_spc_dct[iso]
            orig_ich = spc_dct['canon_enant_ich']
            rac_ich = chi.racemic(orig_ich)
            if rac_ich in rac_dct:
                rac_dct[rac_ich].append(iso)
            else:  # if on the first time through
                # RENAMING OPTION 1: Keep things the same
                rac_name = iso

                # Leave this alone
                rac_dct[rac_ich] = [iso]
                new_spc_dct = copy.deepcopy(spc_dct)
                new_spc_dct['inchi'] = rac_ich
                rac_mech_spc_dct[rac_name] = new_spc_dct
                rac_names.append(rac_name)
        # Append to lists
        for rac_set in rac_dct.values():
           # RENAMING OPT. 2: Use _names (gives diff. names from Sarah's mech)
           # if len(rac_set) > 1:  # if an enantiomer
           #     # Use rac_ich (stored from last iso) to get rac_name
           #     rac_name = names.functional_group_name(rac_ich)
           # else:
           #     rac_name = iso
           # rac_names.append(rac_name)

            # Leave this alone
            rac_sets.append(rac_set)

    return rac_sets, rac_names, rac_mech_spc_dct


def find_iso_sets(mech_spc_dct, canon_ent=True):
    """ Finds all sets of stereoisomers in a mech_spc_dct

        Requires that the mech_spc_dct include canonical enantiomers!

        :return iso_sets: list of all sets of stereoisomers
        :rtype: [[iso1, iso2, ...], [iso1, iso2, ...], ...]
    """

    # Get two things: (i) species and (ii) inchis
    spcs = tuple(mech_spc_dct.keys())
    ichs = ()
    for spc, spc_dct in mech_spc_dct.items():
        assert spc_dct.get('canon_enant_ich') is not None, (
            f'The species {spc} is missing a canon_ent_ich!')
        ichs += (spc_dct['canon_enant_ich'],)

    # Get stereo sets, which are sets of species that are the same if one
    # ignores stereo (usually will be singles or doubles)
    iso_sets = []
    already_done = []
    for idx, ich in enumerate(ichs):
        # If this inchi has already been done, skip it
        if idx in already_done:
            continue
        # Otherwise, store ich and look through remaining ichs for match(es)
        iso_set = [spcs[idx]]
        for sub_idx in range(idx + 1, len(ichs)):
            curr_ich = ichs[sub_idx]
            if ich == curr_ich:
                iso_set.append(spcs[sub_idx])
                already_done.append(sub_idx)  # store that ich has been done
        iso_sets.append(iso_set)

    # Check that each iso_set has identical spc_dcts for all isomers
    for iso_set in iso_sets:
        for iso_idx, iso in enumerate(iso_set):
            if iso_idx == 0:  # if on first iso, store ref_spc_dct
                ref_spc_dct = mech_spc_dct[iso]
            else:  # if on later isos, check against reference
                spc_dct = mech_spc_dct[iso]
                same = _are_spc_dcts_same(ref_spc_dct, spc_dct, 
                                          canon_ent=canon_ent)
                assert same, (f'In the set of stereoisomes {iso_set}, the '
                               'species dcts are not the same')

    return iso_sets


def _are_spc_dcts_same(spc_dct1, spc_dct2, canon_ent=True):
    """ Checks if two spc dcts are the same
    """

    if canon_ent:
        ich1 = spc_dct1['canon_enant_ich']
    else:
        ich1 = spc_dct1['inchi']
    mlt1 = spc_dct1['mult']
    chg1 = spc_dct1['charge']
    exc1 = spc_dct1['exc_flag']
    fml1 = spc_dct1['fml']

    same = compare.are_spc_same(ich1, mlt1, chg1, exc1, fml1, spc_dct2, 
                                canon_ent=canon_ent)

    return same


def _sum_ktp_list(ktp_dct_lst):
    """ Sums all ktp_dcts in a list

        :param ktp_dct_list: list of ktp_dcts
        :type ktp_dct_list: list
        :return summed_ktp_dct: summed ktp_dct
        :rtype: dict
    """

    for idx, ktp_dct in enumerate(ktp_dct_lst):
        if idx == 0:
            summed_ktp_dct = copy.deepcopy(ktp_dct)
        else:
            summed_ktp_dct = rates.add_ktp_dcts(summed_ktp_dct, ktp_dct)

    return summed_ktp_dct


def check_bal_rxns(rxn_param_dct, mech_spc_dct):
    """ Checks that all rxns are chemically balanced
    """
    def _atom_counts(spcs, mech_spc_dct):
        """ Get atom counts for either rcts or prds 
        """
        atom_counts = {}
        for spc in spcs:
            spc_dct = mech_spc_dct[spc]
            fml = spc_dct['fml']
            for atom, count in fml.items():
                if atom in atom_counts:
                    atom_counts[atom] += count
                else:
                    atom_counts[atom] = count
        return atom_counts

    warn_msgs = []
    for rxn in rxn_param_dct:
        warn_msg = []
        rcts, prds, _ = rxn
        rct_atom_counts = _atom_counts(rcts, mech_spc_dct)
        prd_atom_counts = _atom_counts(prds, mech_spc_dct)
        if set(rct_atom_counts.keys()) != set(prd_atom_counts.keys()):
            warn_msg.append('Mismatching atom identities in rcts and prds '
                            f'for rxn {rxn}. Rcts: {rct_atom_counts.keys()} '
                            f'Prds: {prd_atom_counts.keys()}')
        else:
            for rct_atom, rct_count in rct_atom_counts.items():
                prd_count = prd_atom_counts[rct_atom]
                if rct_count != prd_count:
                    warn_msg.append(f'Different counts for {rct_atom}! '
                                    f'Rcts: {rct_count}, Prds: {prd_count}')
        if warn_msg != []:
            warn_msgs.append(warn_msg)

    assert warn_msgs == [], ('Unbalanced reactions exist! See below:\n'
                             f'{warn_msgs}')
    
