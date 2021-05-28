""" build objects strictly used for running things
"""

import copy


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
    _lst = tuple(spc for idx, spc in enumerate(spc_dct) if idx in
                 tuple(spc_idxs.values())[0])
    run_dct = {('SPC', 0, 0): _lst}

    return run_dct


def _lst_for_pes(pes_dct, run_pes_idxs):
    """ Get a dictionary of requested species matching the PES_DCT format
    """

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

            red_pes_dct[(form, pidx, sidx)] = red_chnls

    return red_pes_dct


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
