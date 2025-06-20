""" This deals with merging mechanisms, one of which has stereo and the other
    of which does not
"""

import copy
import numpy as np
from automol import amchi
from automol import inchi
from automol import chi
from mechanalyzer.builder import strip_ste
from mechanalyzer.calculator import compare
from autoreact import params as params_module


def expand_all_rxns(mech_spc_dct_ste, mech_spc_dct_noste,
                    rxn_param_dct_noste):

    # Get the chirality of the species in the stereo mechanism
    _, sing_chiral, mult_chiral = classify_chiral(mech_spc_dct_ste)

    # Strip the stereo mechanism of its stereo and get the name maps
    mech_spc_dct_strpd = strip_mech_spc_dct(mech_spc_dct_ste)
    name_maps = get_name_maps(mech_spc_dct_strpd, mech_spc_dct_noste)

    # Expand each rxn
    new_rxn_param_dct = {}
    only_ste_rxn_param_dct = {}
    num_ste_rxns = 0
    num_added_rxns_chiral = 0
    num_added_rxns = 0
    for rxn, params in rxn_param_dct_noste.items():
        new_rxns, new_params, has_stereo, has_mult_chiral = expand_one_rxn(
            rxn, params, name_maps, sing_chiral, mult_chiral)
        # Populate with new rate parameters
        for new_rxn in new_rxns:
            new_rxn_param_dct[new_rxn] = new_params
        # If has_stereo, populate only stereo with new rate parameters
        if has_stereo:
            num_ste_rxns += 1
            for new_rxn in new_rxns:
                only_ste_rxn_param_dct[new_rxn] = new_params
        if has_mult_chiral:
            num_added_rxns_chiral += len(new_rxns)
        num_added_rxns += len(new_rxns) - 1

    print(f'# of ste-containing rxns in orig., no-ste mech: {num_ste_rxns}')
    print('num_added_rxns_chiral: ', num_added_rxns_chiral)
    print('num_added_rxns: ', num_added_rxns)
    print('len of no_ste mech: ', len(rxn_param_dct_noste))
    print('len of new mech: ', len(new_rxn_param_dct))

    return new_rxn_param_dct, only_ste_rxn_param_dct


def rename_spc(mech_spc_dct_ste, mech_spc_dct_noste, spc_nasa7_dct_ste,
               spc_nasa7_dct_noste):
    """ Take the no-stereo mechanism and replace any names with the stereo
        mech versions (both in the mech_spc_dct and the spc_nasa7_dct)
    """

    # Strip the stereo mechanism of its stereo and get the name maps
    mech_spc_dct_strpd = strip_mech_spc_dct(mech_spc_dct_ste)
    name_maps = get_name_maps(mech_spc_dct_strpd, mech_spc_dct_noste)

    # Loop over each species and rename as needed
    new_mech_spc_dct = {}
    new_spc_nasa7_dct = {}
    for spc, spc_dct in mech_spc_dct_noste.items():
        nasa7 = spc_nasa7_dct_noste[spc]
        if spc in name_maps:  # if in the renaming instructions
            for new_name in name_maps[spc]:
                spc_dct_ste = mech_spc_dct_ste[new_name]
                new_mech_spc_dct[new_name] = spc_dct_ste
                new_spc_nasa7_dct[new_name] = nasa7  # use non-ste thermo!
        else:  # just use the original spc name
            new_mech_spc_dct[spc] = spc_dct
            new_spc_nasa7_dct[spc] = nasa7

    return new_mech_spc_dct, new_spc_nasa7_dct


def expand_one_rxn(rxn, params, name_maps, sing_chiral, mult_chiral):
    """
    """

    def _new_rcts_or_prds(rcts_or_prds, name_maps, sing_chiral, mult_chiral):
        """ Creates new rcts and prds, where each species is now a list that
            has multiple entries if stereo is present
        """
        new_rcts_or_prds = []
        has_mult_chiral = False
        for rct_or_prd in rcts_or_prds:
            if rct_or_prd in name_maps and rct_or_prd not in sing_chiral:
                new_rcts_or_prds.append(name_maps[rct_or_prd])
            else:
                new_rcts_or_prds.append([rct_or_prd])
            # Check if spc in mult_chiral; just for debugging purposes
            for new_rct_or_prd in new_rcts_or_prds:
                for iso in new_rct_or_prd:
                    if iso in mult_chiral:
                        has_mult_chiral = True

        return new_rcts_or_prds, has_mult_chiral

    # Get the renamed rcts and prds
    rcts, prds, tbody = rxn
    new_rcts, has_mult_chiral_rcts = _new_rcts_or_prds(
        rcts, name_maps, sing_chiral, mult_chiral)
    new_prds, has_mult_chiral_prds = _new_rcts_or_prds(
        prds, name_maps, sing_chiral, mult_chiral)

    # Get the reaction parts array
    shape = []
    has_stereo = False
    for rct in new_rcts:
        shape.append(len(rct))
        if len(rct) > 1:
            has_stereo = True
    factor = 1  # factor by which rates are divided
    for prd in new_prds:
        shape.append(len(prd))
        factor /= len(prd)  # depends on # of stereoisomers in prods
        if len(prd) > 1:
            has_stereo = True
    rxn_parts = np.ndarray(tuple(shape), dtype='object')

    # Now, get all possible reactions by permuting
    new_rxn_names = []
    for idxs, _ in np.ndenumerate(rxn_parts):
        # Populate reactants
        curr_rcts = []
        for rct_idx, rct in enumerate(new_rcts):
            curr_rct = rct[idxs[rct_idx]]
            curr_rcts.append(curr_rct)
        # Populate products
        curr_prds = []
        for prd_idx, prd in enumerate(new_prds):
            curr_prd = prd[idxs[prd_idx + rct_idx + 1]]
            curr_prds.append(curr_prd)
        # Store
        new_rxn = (tuple(curr_rcts), tuple(curr_prds), tbody)
        new_rxn_names.append(new_rxn)

    # Scale rate constants by the factor
    new_params = params_module.multiply_factor(params, factor)

    # If either rcts or prds had mult_chiral, set to True
    has_mult_chiral = False
    if has_mult_chiral_rcts or has_mult_chiral_prds:
        has_mult_chiral = True

    return new_rxn_names, new_params, has_stereo, has_mult_chiral


def strip_mech_spc_dct(mech_spc_dct):
    """ Removes stereochemistry from all species in a mech_spc_dct. Returns
        a new mech_spc_dct with all the stereo-specific species, but with the
        inchis now stripped of the stereo

        :param mech_spc_dct: input mech_spc_dct
        :type mech_spc_dct: dict
        :return mech_spc_dct_strpd: dct with only stereo specific spcs, but
            with stereo stripped from the inchis
        :rtype: dict
    """

    mech_spc_dct_strpd = {}
    for spc, spc_dct in copy.deepcopy(mech_spc_dct).items():
        orig_ich = spc_dct['inchi']
        # Remove stereo from the inchi (if present)  (false keeps it as AmChI)
        strpd_ich = chi.without_stereo(orig_ich, reassess_amchi=False)
        # Get the smiles and inchikey without stereo
        strpd_smi = chi.smiles(strpd_ich)
        strpd_ichkey = chi.inchi_key(strpd_ich)
        # Store the stereo-stripped information
        spc_dct['smiles'] = strpd_smi
        spc_dct['inchikey'] = strpd_ichkey
        spc_dct['inchi'] = strpd_ich
        mech_spc_dct_strpd[spc] = spc_dct

    return mech_spc_dct_strpd


def get_name_maps(mech_spc_dct_strpd, mech_spc_dct_noste):
    """ Gets rename instructions for taking a stereo-free mech and adding
        stereo to it

        :param mech_spc_dct_strpd: dct with only stereo specific spcs, but
            with stereo stripped from the inchis
        :type mech_spc_dct_strpd: dict
        :return
    """

    rename_str = '-zz'

    name_maps = {}
    already_done = []
    for spc1, spc_dct1 in mech_spc_dct_noste.items():
        name_map = []
        ich1, mlt1, chg1, exc1, fml1 = compare._read_spc_dct(spc_dct1)
        for spc2, spc_dct2 in mech_spc_dct_strpd.items():
            spc_same = compare.are_spc_same(
                ich1, mlt1, chg1, exc1, fml1, spc_dct2)
            # If species are identical
            if spc_same:
                if spc1 != spc2:  # if spc names different, add to name_map
                    name_map.append(spc2)
                    already_done.append(spc2)
            # If species are different but have same name
            elif spc1 == spc2:
                name_map.append(spc2 + rename_str)
                # Note: don't add to already_done; this works for now
        if name_map != []:
            name_maps[spc1] = name_map

    return name_maps


def classify_chiral(mech_spc_dct):
    """ Sorts a mech_spc_dct into three categories: no chirality, single
        chirality, and multiple chiralities
    """
    def _chiral_count(ich):
        """ Counts the number of chiral sites on a species
        """

        bonds = len(amchi.bond_stereo_parities(ich))
        if bonds != 0:  # I think this should be true if it's E/Z
            count = 2
        else:  # no E/Z, so just count chiral sites
            atoms = len(amchi.atom_stereo_parities(ich))
            count = atoms

        return count

    no_chiral = []
    sing_chiral = []
    mult_chiral = []
    for spc, spc_dct in mech_spc_dct.items():
        count = _chiral_count(spc_dct['inchi'])
        if count == 0:
            no_chiral.append(spc)
        elif count == 1:
            sing_chiral.append(spc)
        else:
            mult_chiral.append(spc)

    return no_chiral, sing_chiral, mult_chiral


def remove_all_ts(mech_spc_dct):

    no_ts_spcs = tuple([spc for spc in mech_spc_dct.keys()
                        if 'rxndirn' not in mech_spc_dct[spc]])
    no_ts_mech_spc_dct = {}
    for spc in no_ts_spcs:
        no_ts_mech_spc_dct[spc] = mech_spc_dct[spc]

    return no_ts_mech_spc_dct

