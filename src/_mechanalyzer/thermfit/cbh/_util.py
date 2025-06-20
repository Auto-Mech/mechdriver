""" Generally useful functions
"""

import numpy
import automol.chi
import automol.graph


# misc
def stoich_gra(gra):
    """ convert gra to formula. required for ts graph
        since they kill the automol converter
    """
    atms = automol.graph.atoms(gra)
    stoich_dct = {}
    hcount = 0
    for atm in atms:
        hcount += numpy.floor(atms[atm][1])
        if atms[atm][0] in stoich_dct:
            stoich_dct[atms[atm][0]] += 1
        else:
            stoich_dct[atms[atm][0]] = 1
    if 'H' not in stoich_dct:
        stoich_dct['H'] = hcount
    else:
        stoich_dct['H'] += hcount
    return stoich_dct


# Assess points of interest
def branch_point(adj_atms_i, adj_atms_j=(), adj_atms_k=()):
    """ Branch point?
    """
    return max(1, len(adj_atms_i) + len(adj_atms_j) + len(adj_atms_k) - 1)


def branch_point2(nonhyd_adj_atms_dct, extended_site):
    """ Branch point?
    """
    branches = -1
    for atm in list(set(extended_site)):
        branches += len(nonhyd_adj_atms_dct[atm])
    return max(1, branches)


def terminal_moiety2(nonhyd_adj_atms_dct, extended_site):
    branches = 0
    terminal = 1
    for atm in list(set(extended_site)):
        branches += len(nonhyd_adj_atms_dct[atm])
    if branches < 2:
        terminal = 0
    return terminal


def terminal_moiety(adj_atms_i, adj_atms_j=(), adj_atms_k=(), endisterm=True):
    """ Get moiety of termial groups?
    """

    ret = 1
    if len(adj_atms_i) + len(adj_atms_j) < 2:
        ret = 0
    if ret == 0 and len(adj_atms_k) > 0 and not endisterm:
        ret = 1

    return ret


# Balance functions
def add2dic(dic, key, val=1):
    """ helper function to add a key to dct
    """
    if key in dic:
        dic[key] += val
    else:
        dic[key] = val


def balance(ich, frags):
    """ balance the equation?
    """
    stoichs = {}
    for frag in frags:
        _stoich = automol.chi.formula(frag)
        for atm in _stoich:
            if atm in stoichs:
                stoichs[atm] += _stoich[atm] * frags[frag]
            else:
                stoichs[atm] = _stoich[atm] * frags[frag]
    _stoich = automol.chi.formula(ich)
    for atom in _stoich:
        if atom in stoichs:
            stoichs[atom] += -_stoich[atom]
        else:
            stoichs[atom] = -_stoich[atom]
    balance_ = {x: round(-y, 3) for x, y in stoichs.items() if y != 0}
    return balance_


def balance_ts(gra, frags):
    """ balance the equation using graphs
    """
    stoichs = {}
    for frag in frags:
        if 'exp_gra' in frags[frag]:
            _stoich = automol.graph.formula(frags[frag]['exp_gra'])
        elif 'ts_gra' in frags[frag]:
            _stoich = stoich_gra(frags[frag]['ts_gra'])
        for atm in _stoich:
            if atm in stoichs:
                stoichs[atm] += _stoich[atm] * frags[frag]['coeff']
            else:
                stoichs[atm] = _stoich[atm] * frags[frag]['coeff']
    balance_ = {}
    _stoich = automol.graph.formula(gra)
    for atom in _stoich:
        if atom in stoichs:
            balance_[atom] = _stoich[atom] - stoichs[atom]
        else:
            balance_[atom] = _stoich[atom]
    balance_ = {x: y for x, y in balance_.items() if y != 0}
    return balance_


def balance_frags(ich, frags):
    """ balance the equation?
    """
    balance_ = balance(ich, frags)
    methane = automol.smiles.chi('C')
    water = automol.smiles.chi('O')
    ammonm = automol.smiles.chi('N')
    hydrgn = automol.smiles.chi('[H][H]')
    if 'C' in balance_:
        add2dic(frags, methane, balance_['C'])
    if 'N' in balance_:
        add2dic(frags, ammonm, balance_['N'])
    if 'O' in balance_:
        add2dic(frags, water, balance_['O'])
    balance_ = balance(ich, frags)
    if 'H' in balance_:
        add2dic(frags, hydrgn, balance_['H']/2)
    return frags


def balance_frags_ts(gra, frags):
    """ balance the equation?
    """
    balance_ = balance_ts(gra, frags)
    methane = automol.smiles.chi('C')
    water = automol.smiles.chi('O')
    ammonm = automol.smiles.chi('N')
    hydrgn = automol.smiles.chi('[H][H]')
    methane = automol.chi.graph(methane)
    water = automol.chi.graph(water)
    ammonm = automol.chi.graph(ammonm)
    hydrgn = automol.chi.graph(hydrgn)
    idx_dct = []
    for spc in [methane, water, ammonm, hydrgn]:
        spc = automol.graph.explicit(spc)
        found = False
        for frag in frags:
            if 'exp_gra' in frags[frag]:
                if automol.graph.isomorphism(frags[frag]['exp_gra'], spc):
                    idx = frag
                    found = True
                    break
        if not found:
            idx = len(frags.keys())
            frags[idx] = {}
            frags[idx]['exp_gra'] = spc
            frags[idx]['coeff'] = 0.0
        idx_dct.append(idx)
    if 'C' in balance_:
        add2dic(frags[idx_dct[0]], 'coeff', balance_['C'])
    if 'N' in balance_:
        add2dic(frags[idx_dct[1]], 'coeff', balance_['N'])
    if 'O' in balance_:
        add2dic(frags[idx_dct[2]], 'coeff', balance_['O'])
    balance_ = balance_ts(gra, frags)
    if 'H' in balance_:
        add2dic(frags[idx_dct[3]], 'coeff', balance_['H']/2)
    return frags


# I/O
def _lhs_rhs(frags):
    """ Determine the left-hand side and right-hand side of reaction
    """
    rhs = {}
    lhs = {}
    for frag in frags:
        if frags[frag] > 0:
            rhs[frag] = frags[frag]
        elif frags[frag] < 0:
            lhs[frag] = - frags[frag]
    return lhs, rhs


def print_lhs_rhs(ich, frags):
    """ print the fragments from each side of the reaction
    """
    lhs, rhs = _lhs_rhs(frags)
    lhsprint = automol.chi.smiles(ich)
    rhsprint = ''
    for frag, side in rhs.values():
        if rhsprint:
            rhsprint += f' +  {side:.1f} {automol.chi.smiles(frag)} '
        else:
            rhsprint = f' {rhs[frag]:.1f} {automol.chi.smiles(frag)} '
    for frag, side in lhs.values():
        lhsprint += f' +  {side:.1f} {automol.chi.smiles(frag)} '

    return f'{lhsprint} --> {rhsprint}'



def cleave_group_and_saturate(gra, bnd_ords, atmi, atmj):
    """
    Fragments the bnd between atmi-atmj and adds hydrogens
    to atmi for each bond order that was removed. Returns
    the graph that atmi is in
    INPUT:
    gra -- graph
    atmi -- atm we are saturating
    atmj -- atm we are cleaving
    OUTPUT
    gra -- new graph
    """
    #def _update_atom_stereo(gra, atm_idx):
    #    atm_stereo_pars = automol.graph.atom_stereo_parities(gra)
    #    if atm_stereo_pars[atm_idx] is not None: 
    #        new_atm_stereo_pars = automol.graph.calculate_priorities_and_assign_parities(gra)
    #        if new_atm_stereo_pars[atm_idx] is None:
    #            atm_stereo_pars[atm_idx] = None
    #            gra = automol.graph.set_atom_stereo_parities(gra, atm_stereo_pars)
    #    return gra 
    # Graphical info about molecule
    if frozenset({atmi, atmj}) in automol.graph.bonds(gra):
        order = bnd_ords[frozenset({atmi, atmj})]
        for _ in range(order):
            gra = automol.graph.add_bonded_atom(gra, 'H', atmi)
            gra = automol.graph.add_bonded_atom(gra, 'H', atmj)
        # gra = _update_atom_stereo(gra, atmi)
        # gra = _update_atom_stereo(gra, atmj)
        gra = automol.graph.remove_bonds(gra, (frozenset({atmi, atmj}),))
        gras = automol.graph.connected_components(gra)
        new_gra = None
        for gra_k in gras:
            atmsk = automol.graph.atom_keys(gra_k)
            if atmi in atmsk:
                new_gra = gra_k
        # # AVC: I don't think the following code does anything...
        # # (Commenting out for now, because it does nothing)
        # atm_stereo_pars = automol.graph.atom_stereo_parities(new_gra)
        # _, new_atm_stereo_gra, _, _ = automol.graph.calculate_stereo(new_gra)
        # new_atm_stereo_pars = automol.graph.atom_stereo_parities(new_atm_stereo_gra)
        # if atm_stereo_pars != new_atm_stereo_pars:
        #     new_gra = automol.graph.set_atom_stereo_parities(new_gra, new_atm_stereo_pars)
    else:
        new_gra = gra
    return new_gra
