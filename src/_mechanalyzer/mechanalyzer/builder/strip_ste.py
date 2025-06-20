""" Removes stereo-specific species and rates from a mechanism by averaging or
    adding rate constants as needed
"""

import copy
from mechanalyzer.calculator import compare
from mechanalyzer.calculator import rates
from mechanalyzer.calculator.rates import check_p_t
from mechanalyzer.builder import _names as names
from ratefit.fit import _fit as fit
from automol import inchi
from automol import chi


def main(rxn_param_dct, mech_spc_dct, temps_lst, pressures):
    """ Main function; removes stereo from reactions and species and averages
        or adds rate constants as appropriate

        :param rxn_param_dct: rxn_param_dct with possible stereoisomers
        :type rxn_param_dct: {rxn: params, ...}
        :param mech_spc_dct: mech_spc_dct with possible stereoisomers
        :type mech_spc_dct: {spc: spc_dct, ...}
        :param temps_lst: temp arrays, one per pressure, for fitting (K)
        :type temps_lst: [np.array1, np.array2, ...]
        :param pressures: pressures for fitting (atm)
        :type pressures: list
        :return re_rxn_param_dct_comb: rxns with all stereo rates averaged
        :rtype: {rxn: params, ...}
        :return re_mech_spc_dct_comb: spcs with all stereo removed
        :rtype: {spc: spc_dct, ...}
    """

    print('rxn_param_dct:\n', rxn_param_dct)
    print('mech_spc_dct:\n', mech_spc_dct)

    # Check/reform temps list
    temps_lst = check_p_t(temps_lst, pressures)

    # Strip stereo layers from inchis and save non-stereo spcs for later
    mech_spc_dct_strpd, mech_spc_dct_no_ste = strip_mech_spc_dct(mech_spc_dct)
    print('mech_spc_dct_strpd:\n', mech_spc_dct_strpd)
    print('mech_spc_dct_no_ste:\n', mech_spc_dct_no_ste)

    # Get the iso_sets, i.e., the sets of species that are stereoisomers
    iso_sets = find_iso_sets(mech_spc_dct_strpd)
    print('iso_sets:\n', iso_sets)

    # Rename mech spc_dct to have the spc names be (stereo-stripped) inchis
    mech_spc_dct_strpd_ich = make_mech_spc_dct_ich(
        iso_sets, mech_spc_dct_strpd)
    print('mech_spc_dct_strpd_ich:\n', mech_spc_dct_strpd_ich)

    # Get the reactions and params for each iso in each iso set
    iso_sets_rxns = get_ste_rxns(rxn_param_dct, iso_sets)
    print('iso_sets_rxns:\n', iso_sets_rxns)

    # Rename the stereo species to stereo-stripped and get the aligned params
    algn_iso_sets_rxns = align_rxns(iso_sets_rxns, mech_spc_dct_strpd_ich,
                                    mech_spc_dct_strpd, rxn_param_dct)
    print('algn_iso_sets_rxns:\n', algn_iso_sets_rxns)

    # Combine (either add or average) rates for all stereo reactions
    iso_sets_par_comb = get_comb_params(algn_iso_sets_rxns, temps_lst, pressures)
    print('iso_sets_par_comb:\n', iso_sets_par_comb)

    # Combine all iso_sets back into a single rxn_param_dct
    rxn_param_dct_strpd = join_rxns(iso_sets_par_comb)
    print('rxn_param_dct_strpd:\n', rxn_param_dct_strpd)

    # Get rxn_param_dct with no stereoisomers
    rxn_param_dct_no_ste = get_no_ste_rxns(rxn_param_dct, iso_sets_rxns)
    print('rxn_param_dct_no_ste:\n', rxn_param_dct_no_ste)

    # Create stereo-free mechanism names using the stereo-stripped inchis and
    # rename all species in the mechanism with these names
    re_mech_spc_dct, re_rxn_param_dct = regenerate_names(
        mech_spc_dct_strpd_ich, rxn_param_dct_strpd)

    print('re_rxn_param_dct:\n', re_rxn_param_dct)

    # Reunite the non-stereo stuff with the newly stereo-stripped stuff
    re_rxn_param_dct_comb, re_mech_spc_dct_comb = comb_strpd_and_no_ste(
        re_mech_spc_dct, re_rxn_param_dct, mech_spc_dct_no_ste,
        rxn_param_dct_no_ste)

    return re_rxn_param_dct_comb, re_mech_spc_dct_comb


def comb_strpd_and_no_ste(re_mech_spc_dct, re_rxn_param_dct,
                          mech_spc_dct_no_ste, rxn_param_dct_no_ste):
    """ Reunites the long-lost non-stereo stuff with the now-stereo-free stuff

        :param re_mech_spc_dct: stereo-stripped spcs
        :type re_mech_spc_dct: dict
        :param re_rxn_param_dct: stereo-stripped reaction params
        :type re_rxn_param_dct: dict
        :param mech_spc_dct_no_ste: all non-stereo spcs, unchanged
        :type mech_spc_dct_no_ste: dict
        :param rxn_param_dct_no_ste: only reactions without stereoisomers
        :type rxn_param_dct_no_ste: dict
        :return mech_spc_dct_final: mech_spc_dct with all stereo stripped and
            all non-stereo species added back in
        :rtype: dict
        :return mech_spc_dct_final: mech_spc_dct with all stereo stripped and
            all non-stereo species added back in
        :rtype: dict
    """

    # Deepcopy to prevent changes
    mech_spc_dct_final = copy.deepcopy(re_mech_spc_dct)
    rxn_param_dct_final = copy.deepcopy(re_rxn_param_dct)

    # Combine dictionaries
    mech_spc_dct_final.update(mech_spc_dct_no_ste)
    rxn_param_dct_final.update(rxn_param_dct_no_ste)

    return rxn_param_dct_final, mech_spc_dct_final


def regenerate_names(mech_spc_dct_strpd_ich, rxn_param_dct_strpd_ich):
    """ Regenerates mechanism names according to inchis, then renames inchis

        :param mech_spc_dct_strpd_ich: dct with only stereo specific spcs,
            but with stereo stripped from the inchis and inchis as spcs names
        :type: dict
        :param rxn_param_dct_strpd_ich: a normal rxn_param_dct, where the
            stereo specificity has been removed (names in inchis)
        :type: rxn_param_dct
        :return re_mech_spc_dct: names converted back to mechanism version
        :rtype: dict
        :return re_rxn_param_dct: names converted back to mechanism version
        :rtype: dict
    """

    map_dct = names.functional_group_name_dct(mech_spc_dct_strpd_ich,
                                              force_rename=True)
    re_mech_spc_dct, re_rxn_param_dct = names.remap_mechanism_names(
        mech_spc_dct_strpd_ich, rxn_param_dct_strpd_ich, map_dct)

    return re_mech_spc_dct, re_rxn_param_dct


def join_rxns(iso_sets_par_comb):
    """ Takes iso_sets where each rxn has a single set of params and combines
        into a single rxn_param_dct

        :param iso_sets_par_comb: list of rxn_param_dcts, one per iso_set,
            where each rxn has a single, combined RxnParams object
        :type: [{rxn1: params, ...}, {rxn1: params, ...}, ...]
        :return rxn_param_dct_strpd_ich: a normal rxn_param_dct, where the
            stereo specificity has been removed (names in inchis)
        :rtype: rxn_param_dct
    """

    rxn_param_dct_strpd_ich = {}
    for iso_set in iso_sets_par_comb:
        for rxn, params in iso_set.items():
            assert rxn not in rxn_param_dct_strpd_ich, (
                f'The rxn {rxn} already exists!')
            rxn_param_dct_strpd_ich[rxn] = params

    return rxn_param_dct_strpd_ich


def get_comb_params(algn_iso_sets_rxns, temps_lst, pressures):
    """ For each aligned rxn_param_dct (one for each iso_set), combine the
        RxnParams by averaging or adding the rate constants

        :param algn_iso_sets_rxns: list of aligned rxn_param_dcts, one for
            each iso_set
        :type: [{rxn1: [params1, params2, ...], ...},
                {rxn1: [params1, params2, ...], ...}, ...]
        :param temps_lst: temperatures at which to evaluate rates
        :type temps_lst: [numpy.array, numpy.array, ...]
        :param pressures: pressures at which to evaluate rates
        :type pressures: list
        :return iso_sets_par_comb: list of rxn_param_dcts, one per iso_set,
            where each rxn has a single, combined RxnParams object
        :rtype: [{rxn1: params, ...}, {rxn1: params, ...}, ...]
    """

    iso_sets_par_comb = []
    for iso_set in algn_iso_sets_rxns:
        new_iso_set = {}
        for rxn, params_lst in iso_set.items():
            new_iso_set[rxn] = _combine_single_rxn(rxn, params_lst, temps_lst,
                                                   pressures)
        iso_sets_par_comb.append(new_iso_set)

    return iso_sets_par_comb


def align_rxns(iso_sets_rxns, mech_spc_dct_strpd_ich,
               mech_spc_dct_strpd, rxn_param_dct):
    """ Renames all stereo species to inchis and aligns RxnParams under
        singular reaction names

        Note to self: this function has the most complicated logic; errors
            are likely to originate here

        :param iso_sets_rxns: rxn_params_dcts for each isomer in each iso_set
        :type: [[rxn_param_dict_iso1, rxn_param_dct_iso2, ...],
        :param mech_spc_dct_strpd_ich: dct with only stereo specific spcs,
            but with stereo stripped from the inchis and inchis as spcs names
        :type: dict
        :param mech_spc_dct_strpd: dct with only stereo specific spcs, but
            with stereo stripped from the inchis
        :type mech_spc_dct_strpd: dict
        :param rxn_param_dct: keys are rxns, values are RxnParams
        :type rxn_param_dct: dict
        :return algn_iso_sets_rxns: list of aligned rxn_param_dcts, one for
            each iso_set
        :rtype: [{rxn1: [params1, params2, ...], ...},
                 {rxn1: [params1, params2, ...], ...}, ...]
    """

    # Get the rename instructions
    rename_instr = compare.get_rename_instr(mech_spc_dct_strpd_ich,
                                            mech_spc_dct_strpd)

    # Loop over each iso_set and rename all the rxns
    algn_iso_sets_rxns = []
    all_rxns = []
    print('inside new align_rxns')
    for iso_set_rxns in iso_sets_rxns:
        print('iso_set_rxns:\n', iso_set_rxns)
        algn_iso_set_rxns = {}
        for iso_rxns in iso_set_rxns:
            print('iso_rxns:\n', iso_rxns)
            renamed_dct, ste_dct = compare.rename_species(iso_rxns, rename_instr)
            print('ste_dct:\n', ste_dct)
            # Loop over the ste_dct first; this catches reactions that have
            # duplicates *within the current iso_rxns*
            for new_rxn, old_rxns in ste_dct.items():
                print('new_rxn:\n', new_rxn)
                params_lst = [rxn_param_dct[old_rxn] for old_rxn in old_rxns]
                # If new rxn already in current aligned dct, extend params
                if new_rxn in algn_iso_set_rxns:
                    algn_iso_set_rxns[new_rxn].extend(params_lst)
                # If new rxn not in current aligned dct, create w/params
                else:
                    algn_iso_set_rxns[new_rxn] = params_lst
            # Loop over the renamed dct; this catches reactions that do not
            # have duplicates *within the current iso_rxns*
            for new_rxn, params in renamed_dct.items():
                if new_rxn not in ste_dct:  # only do if rxn not in ste_dct
                    if new_rxn in algn_iso_set_rxns:
                        algn_iso_set_rxns[new_rxn].append(params)
                    else:
                        algn_iso_set_rxns[new_rxn] = [params]
        # If any of the reactions in the current aligned_dct showed up in any
        # previous iso_sets, delete the reaction in the current aligned_dct
        for new_rxn in copy.deepcopy(algn_iso_set_rxns).keys():
            if new_rxn in all_rxns:
                algn_iso_set_rxns.pop(new_rxn)  # remove
            else:
                all_rxns.append(new_rxn)  # add to list of existing rxns
        # Store the current aligned_dct
        algn_iso_sets_rxns.append(algn_iso_set_rxns)

    return algn_iso_sets_rxns


def get_no_ste_rxns(rxn_param_dct, iso_sets_rxns):
    """ Gets a rxn_param_dct with all reactions that have no stereoisomers

        :param rxn_param_dct: keys are rxns, values are RxnParams
        :type rxn_param_dct: dict
        :param iso_sets_rxns: rxn_params_dcts for each isomer in each iso_set
        :type: [[rxn_param_dict_iso1, rxn_param_dct_iso2, ...],
            [rxn_param_dict_iso1, rxn_param_dct_iso2, ...], ...]
        :return rxn_param_dct_no_ste: only reactions without stereoisomers
        :rtype: dict
    """

    # Get list of all rxns involving stereoisomers
    all_iso_rxns = []
    for iso_set_rxns in iso_sets_rxns:
        for iso_rxns in iso_set_rxns:
            all_iso_rxns.extend(list(iso_rxns.keys()))

    rxn_param_dct_no_ste = copy.deepcopy(rxn_param_dct)
    for rxn in all_iso_rxns:
        if rxn in rxn_param_dct_no_ste:
            rxn_param_dct_no_ste.pop(rxn)

    return rxn_param_dct_no_ste


def get_ste_rxns(rxn_param_dct, iso_sets):
    """ For each isomer, gets a rxn param dct with all rxns containing that
        isomer

        :param rxn_param_dct: keys are rxns, values are RxnParams
        :type rxn_param_dct: dict
        :param iso_sets: list of all sets of stereoisomers
        :type iso_sets: [[iso1, iso2, ...], [iso1, iso2, ...], ...]
        :return iso_sets_rxns: rxn_params_dcts for each isomer in each iso_set
        :rtype: [[rxn_param_dict_iso1, rxn_param_dct_iso2, ...],
            [rxn_param_dict_iso1, rxn_param_dct_iso2, ...], ...]
    """
    def search_rxns(rxn_param_dct, spc):
        """ Searches a rxn_param_dct for all rxns that contain a species
        """
        filt_rxn_param_dct = {}
        for rxn, params in rxn_param_dct.items():
            rcts, prds, _ = rxn
            if spc in rcts or spc in prds:
                filt_rxn_param_dct[rxn] = params

        return filt_rxn_param_dct

    iso_sets_rxns = []
    # Loop over each set of isomers
    for iso_set in iso_sets:
        iso_set_rxns = []
        for iso in iso_set:
            # Get all rxns that contain this isomer and store
            iso_rxn_params = search_rxns(rxn_param_dct, iso)
            iso_set_rxns.append(iso_rxn_params)
        iso_sets_rxns.append(iso_set_rxns)

    return iso_sets_rxns


def make_mech_spc_dct_ich(iso_sets, mech_spc_dct_strpd):
    """ Renames a mech_spc_dct_strpd (i.e., all stereo stripped) to have the
        species names be inchis

        :param iso_sets: list of all sets of stereoisomers
        :type iso_sets: [[iso1, iso2, ...], [iso1, iso2, ...], ...]
        :param mech_spc_dct_strpd: dct with only stereo specific spcs, but
            with stereo stripped from the inchis
        :type mech_spc_dct_strpd: dict
        :return mech_spc_dct_strpd_ich: dct with only stereo specific spcs,
            but with stereo stripped from the inchis and inchis as spcs names
        :rtype: dict
    """

    mech_spc_dct_strpd_ich = {}
    for iso_set in iso_sets:
        iso = iso_set[0]  # just take first iso b/c all spc_dcts are same
        spc_dct = mech_spc_dct_strpd[iso]
        ich = spc_dct['inchi']
        mech_spc_dct_strpd_ich[ich] = spc_dct

    return mech_spc_dct_strpd_ich


def find_iso_sets(mech_spc_dct_strpd, canon_ent=False):
    """ Finds all sets of isomers in a stripped mech_spc_dct that are now the
        exact same species (since stereo has been stripped)

        :param mech_spc_dct_strpd: dct with only stereo specific spcs, but
            with stereo stripped from the inchis
        :type mech_spc_dct_strpd: dict
        :return iso_sets: list of all sets of stereoisomers
        :rtype: [[iso1, iso2, ...], [iso1, iso2, ...], ...]
    """

    # Get two things: (i) species and (ii) inchis
    spcs = tuple(mech_spc_dct_strpd.keys())
    ichs = ()
    for spc_dct in mech_spc_dct_strpd.values():
        if canon_ent:
            ichs += (spc_dct['canon_enant_ich'],)
        else:
            ichs += (spc_dct['inchi'],)

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
                ref_spc_dct = mech_spc_dct_strpd[iso]
            else:  # if on later isos, check against reference
                spc_dct = mech_spc_dct_strpd[iso]
                same = _are_spc_dcts_same(ref_spc_dct, spc_dct)
                assert same, (f'In the set of stereoisomes {iso_set}, the '
                               'species dcts are not the same (even after '
                               'stripping stereo)')

    return iso_sets


def strip_mech_spc_dct(mech_spc_dct, canon_ent=False):
    """ Removes stereochemistry from all species in a mech_spc_dct. Returns
        a new mech_spc_dct with all the stereo-specific species, but now
        stripped of the stereo. Also returns a separate mech_spc_dct of
        species that originally had no stereo in them (unchanged)

        :param mech_spc_dct: input mech_spc_dct
        :type mech_spc_dct: dict
        :return mech_spc_dct_strpd: dct with only stereo specific spcs, but
            with stereo stripped from the inchis
        :rtype: dict
        :return mech_spc_dct_no_ste: all non-stereo spcs, unchanged
        :rtype: dict
    """

    mech_spc_dct_strpd = {}
    mech_spc_dct_no_ste = {}
    for spc, spc_dct in mech_spc_dct.items():
        if canon_ent:
            orig_ich = copy.copy(spc_dct['canon_enant_ich'])
        else:
            orig_ich = copy.copy(spc_dct['inchi'])
        print('orig_ich: ', orig_ich)
        strpd_ich = chi.without_stereo(orig_ich)
        #try:
            #strpd_ich = inchi.without_stereo(orig_ich)
        #except:
            #print('spc: ', spc)
        # If the species is stereo-free, save in the no_ste dct
        if orig_ich == strpd_ich:
            mech_spc_dct_no_ste[spc] = spc_dct
        # If the species had stereo, save stereo-stripped info in strpd dct
        else:
            # Get the smiles and inchikey without stereo
            strpd_smi = chi.smiles(strpd_ich)
            strpd_ichkey = chi.inchi_key(strpd_ich)
            # Store the stereo-stripped information
            #spc_dct['smiles'] = strpd_smi
            #spc_dct['inchikey'] = strpd_ichkey
            spc_dct['inchi'] = strpd_ich
            mech_spc_dct_strpd[spc] = spc_dct

    return mech_spc_dct_strpd, mech_spc_dct_no_ste


def _are_spc_dcts_same(spc_dct1, spc_dct2):
    """ Checks if two spc dcts are the same
    """

    ich1 = spc_dct1['inchi']
    mlt1 = spc_dct1['mult']
    chg1 = spc_dct1['charge']
    exc1 = spc_dct1['exc_flag']
    fml1 = spc_dct1['fml']

    same = compare.are_spc_same(ich1, mlt1, chg1, exc1, fml1, spc_dct2)

    return same


def _combine_single_rxn(rxn, params_lst, temps_lst, pressures):
    """ Takes a list of rate parameters, calculates the k(T,P) values, either
        averages or adds them, and then refits them
    """

    if len(params_lst) == 0:  # if none are given
        return None

    # Calculate rates for each params in the list
    ktp_dct_lst = []
    for params in params_lst:
        if params is not None:
            # Get rates
            ktp_dct = rates.eval_params(params, temps_lst, pressures)
            ktp_dct_lst.append(ktp_dct)
            # Get fit form (used later)
            fit_method = params.get_existing_forms()[0]  # first one

    # Combine rate constants
    comb_ktp_dct = _sum_ktp_list(ktp_dct_lst)  # first step: add all the rates
    # If there are two stereo rxns, check if the rates should be averaged
    if len(params_lst) == 2:
        should_avg = _should_avg(rxn)  # depends on where stereo spc is
        if should_avg:
            comb_ktp_dct = rates.mult_by_factor(comb_ktp_dct, 0.5)
    # If there are four stereo rxns, divide rates by two; this is the same as
    # averaging two sets of two and then adding. The exact averaging
    elif len(params_lst) == 4:
        comb_ktp_dct = rates.mult_by_factor(comb_ktp_dct, 0.5)

    # Fit the averaged rate constants
    comb_params, _ = fit.fit_ktp_dct(comb_ktp_dct, fit_method)

    return comb_params


def _should_avg(rxn_ich):
    """ Determines if rate constants should be averaged; this should be done
        if the stereo species is in the products

        :param rxn_ich: description of a reaction, where any stereo species
            are described using their (stereo-stripped) inchis
        :type rxn_ich: ((rct1, rct2, ...), (prd1, prd2, ...), ...)
        :return should_avg: whether the rate constants should be averaged
        :rtype: Bool
    """

    should_avg = False
    _, prds, _ = rxn_ich
    # If the stereo spc is in the products, average the rates
    for prd in prds:
        if 'ChI' in prd:
            should_avg = True

    return should_avg


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
