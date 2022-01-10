""" build objects strictly used for running things
"""

import copy
import itertools


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
                cidx, _ = chnl
                if cidx in run_chnl_idxs:
                    red_chnls += (chnl,)

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
    grp_lst_sort += pes_grp_idxs

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
            pes_grps += ((pes_grp, pes_grp_dct.get(grp_idxs)),)

    return pes_grps


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
