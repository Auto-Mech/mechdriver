""" build objects strictly used for running things
"""

# Overall run lst for both reactions and species
def run_lst(pes_dct, spc_dct, pes_idxs, spc_idxs):
    """ a
    """

    if pes_idxs is not None:

    elif spc_idxs is not None:
        _run_lst = _lst_for_spc(spc_dct, spc_idxs)
    else:
        _run_lst = None

    return _run_lst


def _lst_for_spc(spc_dct, spc_idxs):
    """ Get a dictionary of requested species matching the PES_DCT format
    """

    _lst = tuple(spc for idx, spc in enumerate(spc_dct) if idx+1 in 
                 spc_idxs.values())
    run_dct = {}
    run_dct[('SPC', 1, 1)] = {'species': _lst, 'reacs': (), 'prods': ()}

    return run_dct
    

def _lst_for_pes(pes_dct, pes_idxs):
    """ Get a dictionary of requested species matching the PES_DCT format
    """
    


def build_queue(rxn_lst):
    """ Build spc queue from the reaction lst for the drivers
        :return spc_queue: all the species and corresponding models in rxn
        :rtype: list[(species, model),...]
    """
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


def build_run_spc_dct(spc_dct, run_obj_dct):
    """ Get a dictionary of requested species matching the PES_DCT format
    """
    spc_nums = run_obj_dct['spc']
    run_spc_lst = []
    for idx, spc in enumerate(spc_dct):
        if spc != 'global':
            if idx+1 in spc_nums:
                run_spc_lst.append((spc, run_obj_dct['spc'][idx+1]))

    # Build the run dct
    run_dct = {}
    run_dct['all'] = {
        'species': run_spc_lst,
        'reacs': [],
        'prods': []
    }

    return run_dct
