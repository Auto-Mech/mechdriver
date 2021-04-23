""" build objects strictly used for running things
"""

import itertools


# Overall run lst for both reactions and species
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

    _lst = tuple(spc for idx, spc in enumerate(spc_dct) if idx+1 in
                 spc_idxs.values())
    run_dct = {('SPC', 1, 1): _lst}
    # run_dct[('SPC', 1, 1)] = {'species': _lst, 'reacs': (), 'prods': ()}

    return run_dct


def _lst_for_pes(pes_dct, run_pes_idxs):
    """ Get a dictionary of requested species matching the PES_DCT format
    """

    print('pes_idxs\n', run_pes_idxs)
    print('pes_dct\n', pes_dct)

    prev_idx = -1
    red_pes_dct = {}
    for (form, pidx, sidx), chnls in pes_dct.items():

        # Count number of channels on the previous SUB-PES (if there is one)
        nchnls = len(chnls) if pidx == prev_idx else 0
        prev_idx = pidx

        # Grab all of the run channels from run_pes_idxs if
        # the pes index from pes_dct in the run_pes_idxs
        chnl_idxs = run_pes_idxs.get(pidx+1, None)
        if chnl_idxs is not None:
            # Increment chnl idxs by the num of channles on prev. sub-PES
            _chnl_idxs = tuple(idx+nchnls for idx in chnl_idxs)

            achnls = tuple(chnl for idx, chnl in enumerate(chnls)
                           if idx+1 in _chnl_idxs)

            red_pes_dct[(form, pidx, sidx)] = achnls
    
    # old info that is passed around
    #    run_lst.append(
    #        {'species': spc_queue,
    #         'reacs': list(rct_names_lst[idx]),
    #         'prods': list(prd_names_lst[idx]),
    #         # 'model': rxn_model_lst[idx],
    #         'chn_idx': rxn_chn_idxs[idx],
    #         'dummy': []}

    return red_pes_dct


def spc_queue(rxn_lst):
    """ Build spc queue from the reaction lst for the drivers
        :return spc_queue: all the species and corresponding models in rxn
        :rtype: list[(species, model),...]
    """

    # pron just taking the keys from dct so this should be simpler

    if 'all' in rxn_lst:
        spc_queue = rxn_lst.values()[0]
    else:
        spc_queue = []
        for inf, chnls in rxn_lst:
            spc_queue.append(chnls)
        spc_queue = itertools.chain(*spc_queue)
        spc_queue = tuple(set(spc_queue))

    return spc_queue


def split_queue(spc_queue):
    """ split the queue for multiple pfs
    """

    new_queue = []
    op_dct = {'*': 'multiply', '+': 'add', '/': 'divide', '-': 'substract'}
    for (spc_name, (pes_model, spc_model)) in spc_queue:
        coeffs = []
        operators = []
        models = []
        coeff = ''
        model = ''
        for char in spc_model:
            if char == '.' or char.isdigit():
                coeff += char
            elif char.isalpha():
                model += char
            elif char in op_dct:
                operators.append(op_dct[char])
                if coeff:
                    coeffs.append(float(coeff))
                else:
                    coeffs.append(1)
                models.append(model)
                coeff = ''
                model = ''
        if coeff:
            coeffs.append(float(coeff))
        else:
            coeffs.append(1)
        models.append(model)
        new_queue.append((spc_name, (pes_model, models, coeffs, operators)))
    return new_queue

    # OLD
    if 'all' in rxn_lst:
        spc_queue = rxn_lst['all']['species']
    else:
        spc_queue = []
        for rxn in rxn_lst:
            model = rxn['model']
            spc_queue.extend(((reac, model) for reac in rxn['reacs']))
            spc_queue.extend(((prod, model) for prod in rxn['prods']))

    return spc_queue


def split_queue(spc_queue):
    """ split the queue for multiple pfs
    """

    new_queue = []
    op_dct = {'*': 'multiply', '+': 'add', '/': 'divide', '-': 'substract'}
    for (spc_name, (pes_model, spc_model)) in spc_queue:
        coeffs = []
        operators = []
        models = []
        coeff = ''
        model = ''
        for char in spc_model:
            if char == '.' or char.isdigit():
                coeff += char
            elif char.isalpha():
                model += char
            elif char in op_dct:
                operators.append(op_dct[char])
                if coeff:
                    coeffs.append(float(coeff))
                else:
                    coeffs.append(1)
                models.append(model)
                coeff = ''
                model = ''
        if coeff:
            coeffs.append(float(coeff))
        else:
            coeffs.append(1)
        models.append(model)
        new_queue.append((spc_name, (pes_model, models, coeffs, operators)))
    return new_queue
