""" build objects strictly used for running things
"""

import copy
import itertools
from mechanalyzer.builder import strip_ste
from mechanalyzer.builder import _names as names
import sys

# Overall run lst for both reactions and species
def combine(pes_rlst, spc_rlst):
    """ Combine
    """
    if pes_rlst is not None:
        run_rlst = copy.deepcopy(pes_rlst)
        if spc_rlst is not None:
            run_rlst.update(spc_rlst)
    else:
        run_rlst = copy.deepcopy(spc_rlst)

    return run_rlst


def run_lst(pes_dct, spc_dct, pes_idxs, spc_idxs):
    """ Build the iterable lists of PESs and species that we
        run over in the drivers.
    """
    # need a message to say pes take if both pes and spc used
    if pes_idxs is not None:
        _pes_run_lst = _lst_for_pes(pes_dct, pes_idxs)
    else:
        _pes_run_lst = None

    if spc_idxs is not None:
        _spc_run_lst = _lst_for_spc(spc_dct, spc_idxs)
    else:
        _spc_run_lst = None

    return _pes_run_lst, _spc_run_lst


def _lst_for_spc(spc_dct, spc_idxs):
    """ Get a dictionary of requested species matching the PES_DCT format
    """
    idx_lst = tuple(spc_idxs.values())[0]
    spc_lst = tuple(spc_dct)
    run_spc_lst = tuple(spc_lst[idx] for idx in idx_lst)
    run_dct = {('SPC', 0, idx_lst): run_spc_lst}

    # for idx, spc in enumerate(spc_dct):

    # _lst = tuple(spc for idx, spc in enumerate(spc_dct) if idx in
    #              tuple(spc_idxs.values())[0])
    # run_dct = {('SPC', 0, 0): _lst}

    return run_dct


def _lst_for_pes(pes_dct, run_pes_idxs):
    """ Get a dictionary of requested species matching the PES_DCT format
    """

    # Loop over input pes dictionary of full mechanism and grab PES and
    # channels associated with the requested run idxs
    red_pes_dct = {}
    for (form, pidx, sidx), chnls in pes_dct.items():
        # Grab PES if idx in run_pes_idx dct
        run_chnl_idxs = run_pes_idxs.get(pidx, None)
        if run_chnl_idxs is not None:
            # Grab the channels if they are in run_chnl_idxs
            red_chnls = ()
            for chnl in chnls:
                cidx, rgts = chnl
                if cidx in run_chnl_idxs:
                    # Only grabbing the reactants and products
                    # ignoring the bath gas for now
                    red_chnls += (
                        (cidx, (rgts[0], rgts[1])),
                    )

            # Only add to reduced dct if any chnls found,
            # main for loop over pes AND subpes, could add
            # empty list of channels for subpes without check
            if red_chnls:
                red_pes_dct[(form, pidx, sidx)] = red_chnls

    return red_pes_dct


def pes_groups(pes_dct, pes_grp_dct):
    """ Group the PES into ordered lists to handle multiPES effects.

        pes_dct = {(fml, pes_idx, sub_pes_idx): (chnls,)}
        pes_grp_idxs = ((pes_idx, subpes_idx), (pes_idx, subpes_idx),)

        Tries to sort by all single-PES groups, then multi-PES groups. There
        is a subsort on the PES-SUBPES indices after that.

        :return: pes_dct lst (({pes_dct}, {par_dct}), ...),
        where each dictionary contains all PESs that are member of PES group
    """

    # Build pes grp idx lists to loop over an build master list
    # run_pes_idxs = tuple(frozenset({x, y}) for _, x, y in pes_dct.keys())
    run_pes_idxs = tuple((x, y) for _, x, y in pes_dct.keys())

    # Get the groupings specified by the user
    if pes_grp_dct is not None:
        pes_grp_idxs = tuple(pes_grp_dct.keys())
        flat_pes_grp_idxs = tuple(itertools.chain(*pes_grp_idxs))
    else:
        pes_grp_idxs, flat_pes_grp_idxs = (), ()

    # Get all of the PESs not grouped, then add the ones grouped by user

    grp_lst = tuple(((x, y),) for (x, y) in run_pes_idxs
                    if (x, y) not in flat_pes_grp_idxs)
    grp_lst_sort = tuple(sorted(grp_lst, key=lambda x: x[0]))

    pes_grp_idxs_eff = []
    for grp in pes_grp_idxs:
        toappend = [subpes for subpes in grp
                        if subpes in run_pes_idxs]
        pes_grp_idxs_eff.append(tuple(toappend))
    pes_grp_idxs_eff= tuple(pes_grp_idxs_eff)
    # grp_lst_sort += pes_grp_idxs
    grp_lst_sort += pes_grp_idxs_eff # only add grps you asked for

    # Need to build a group for PESs
    pes_grps = ()
    for grp_idxs in grp_lst_sort:
        pes_grp = {}
        for idxs in grp_idxs:
            for (form, pidx, sidx), chnls in pes_dct.items():
                if idxs == (pidx, sidx):
                    pes_grp.update({(form, pidx, sidx): chnls})
        if pes_grp_dct is None:
            pes_grps += ((pes_grp, None),)
        else:
            if grp_idxs in pes_grp_dct.keys():
                pes_grps += ((pes_grp, pes_grp_dct.get(grp_idxs)),) 
            else:
                # check for each set of indices if grp_idxs are a subset
                new_pes_grp_dct = None
                for pes_lst, grp_dct in pes_grp_dct.items():
                    if all([idx in pes_lst for idx in grp_idxs]) and len(grp_idxs) >= 2:
                        # recreate grp
                        new_pes_grp_dct = {k: [] for k in ['peds', 'hot', 'en_limit']}
                        for idx in grp_idxs:
                            i = list(pes_lst).index(idx)
                            [new_pes_grp_dct[k].append(grp_dct[k][i]) for k in ['hot', 'en_limit']]
                        # now add modeltype and bf threshold
                        new_pes_grp_dct['modeltype'] = grp_dct['modeltype']
                        new_pes_grp_dct['bf_threshold'] = grp_dct['bf_threshold']
                        # add ped only if corresponding hot species is present
                        hots = list(itertools.chain(*new_pes_grp_dct['hot']))
                        for idx in grp_idxs:
                            i = list(pes_lst).index(idx)
                            ped_i = grp_dct['peds'][i]
                            ped_new = []
                            for ped_j in ped_i:
                                frags = ped_j.split('=')[1].split('+')
                                if any([fr in hots for fr in frags]):
                                    ped_new.append(ped_j)
                            new_pes_grp_dct['peds'].append(tuple(ped_new))

                        # check that all hot species are produced
                        ped_prds = list(itertools.chain(*new_pes_grp_dct['peds']))
                        ped_frags = list(itertools.chain(*[[frags.split('=')[1].split('+')] for frags in ped_prds]))
                        ped_frags = list(itertools.chain(*ped_frags))

                        if not all([hot in ped_frags for hot in hots]):
                            # back to single pess
                            for key, val in pes_grp.items():
                                pes_grps += (({key: val}, None,),)
                        else:
                            pes_grps += ((pes_grp, new_pes_grp_dct,),)
                        break

    return pes_grps

def species_groups(pes_rlst, spc_rlst, mech_spc_dct):
    """ Group the species that the user requested to run (given in the
        spc_rlst) that can be grouped by some means.

        Builds a list of species groups where each group
        is comprised of species that are connected by some relationship
        (i.e., they are stereoisomers of each other)

        Any species not connected by desired by a relationship are placed
        into a group of 1.

        List is ordered by having all the groups first, then the remaining
        species.

        Each group is assigned a name, right now it is generic string.

        spc_grps = (
            (name1, (grp_lst1),
            (name2, (grp_lst2),
        )
    """

    # Initialize final groups list to return
    spc_grps = {}

    # Get rlst into a set for later comparisons
    spc_lst = ()
    if pes_rlst is not None:
        for rxn_lst in pes_rlst.values():
            for rxn in rxn_lst:
                spc_lst += rxn[1][0]
                spc_lst += rxn[1][1]
        spc_lst = tuple(i for n, i in enumerate(spc_lst)
                        if i not in spc_lst[:n])
    if spc_rlst is not None:
        spc_lst += list(spc_rlst.values())[0]
    spc_lst_set = set(spc_lst)

    # Get the groups of species grouped by isomer
    # Keep all of the groups composed of species in the spc_rlst
    # iso groups returns all species that COULD have stereoisomers
    # if only one stereoisomer is present, we get a list of 1
    spc_flat_lst = [spc for group in spc_rlst.values() for spc in group]
    mech_spc_dct_no_ts = {spc: dct.copy() for spc, dct in mech_spc_dct.items()
                          if 'ts_' not in spc and spc in spc_flat_lst}
    mech_spc_dct_strpd, _ = strip_ste.strip_mech_spc_dct(mech_spc_dct_no_ts)
    iso_grps = strip_ste.find_iso_sets(mech_spc_dct_strpd)

    # print('iso grp test')
    # for grp in iso_grps:
    #     print(grp)
    # print('---')

    new_iso_grp_names = make_new_names(iso_grps, mech_spc_dct_strpd)

    spc_in_iso_grps = ()
    for idx, iso_grp in enumerate(iso_grps, start=1):
        if set(iso_grp) <= spc_lst_set:
            # Add to generic list to be used in next step
            spc_in_iso_grps += tuple(iso_grp)

            # Rename and add group to master group list
            new_name = new_iso_grp_names[idx - 1]  # -1 b/c starting with 1
            spc_grps.update({new_name: tuple(iso_grp)})
            #niso = len(iso_grp)
            #spc_grps.update({f's{niso}-{iso_grp[0]}': tuple(iso_grp)})
            # iso_name = remove_stereo_name_suffix(
            #     iso_grp[0])
            # spc_grps.update({f'combined-{iso_name}': tuple(iso_grp)})

    # Now get the rest of the spc_rlst not in the iso_grps
    for spc in spc_lst:
        if spc not in spc_in_iso_grps:
            spc_grps.update({'s1-' + spc: (spc,)})

    # Print message saying the groups if there are any
    if any(len(grp) > 1 for grp in spc_grps):
        print('Identified relationships among species used for groupings')
        print('Here are the groups that will be used to get combined info:')
        print('number. group-name group-members')
        for idx, (name, grp) in enumerate(spc_grps.items(), start=1):
            print(f'{idx}. {name} {grp}')

    return spc_grps


def spc_queue(runlst, fml):
    """ Build spc queue from the reaction lst for the drivers.

        Use the formula to discern if run_lst is SPC or PES
        :return spc_queue: all the species and corresponding models in rxn
        :rtype: list[(species, model),...]
    """

    if fml == 'SPC':
        _queue = runlst
    else:
        _ini_queue = []
        for (_, chnl) in runlst:
            _ini_queue += [rgt for rgts in chnl for rgt in rgts]

        # Remove duplicates from the queue, perserving order
        _queue = tuple(i for n, i in enumerate(_ini_queue)
                       if i not in _ini_queue[:n])

    return _queue


def make_new_names(iso_grps, mech_spc_dct_strpd):
    """ Generates new stereo-free species names for each iso_grp
    """

    map_dct = names.functional_group_name_dct(mech_spc_dct_strpd)
    new_iso_grp_names = []
    for iso_grp in iso_grps:
        iso = iso_grp[0]  # grab first one; inchis are same for all isos
        new_iso_grp_names.append(map_dct[iso])

    print('new_iso_grp_names:\n', new_iso_grp_names)
        
    return new_iso_grp_names    


