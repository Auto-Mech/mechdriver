""" Utility functions for handling ts graphs for CBH schemes
"""

import numpy as np
import automol.graph
import automol.chi


# Check TS keys
def intersec(lst1, lst2):
    """ Determine index shared by broken and formed keys , verifying that
        they intersect
    """

    ret = None
    for atm in lst1:
        if atm in lst2:
            ret = atm
    assert ret is not None, (
        f'brk_key {lst1} and frm_key {lst2} do not intersect')
    return ret


def xor(lst1, lst2):
    """ Check

        :param lst1:
        :type lst1:
        :param lst2:
        :type lst2:
    """

    ret = None
    for atm in lst1:
        if atm not in lst2:
            ret = atm
    assert ret is not None, (
        f'problem with bond_key {lst1}')

    return ret


def ts_graph(gra, site1, site2=None):
    rad_atms = list(automol.graph.radical_atom_keys(gra, sing_res=True))
    unsat_atms_dct = automol.graph.atom_unpaired_electrons(automol.graph.kekule(gra))
    atm_vals = automol.graph.atomic_valences(gra)
    print('gra herea', gra)
    atms = automol.graph.atoms(gra)
    orig_bnds = automol.graph.bond_orders(gra)
    print('orig', orig_bnds)
    bnds = automol.graph.kekule_bond_orders(gra)
    print('kekbnds', bnds)
    for bnd, order in orig_bnds.items():
        if bnd in bnds:
            if orig_bnds[bnd] != bnds[bnd]:
                if abs(orig_bnds[bnd] - bnds[bnd]) < 0.99:
                    bnds[bnd] = orig_bnds[bnd]
            if any([atm in site1 for atm in bnd]):
                if bnds[bnd] > orig_bnds[bnd]:
                    bnds[bnd] = orig_bnds[bnd]
                    for atm in bnd:
                        if atm not in site1:
                            unsat_atms_dct[atm] += 1

    unsat_atms = [atm for (atm, sat) in unsat_atms_dct.items() if sat > 0]
    print('unsat atms', unsat_atms_dct)
    adj_atms = automol.graph.atoms_neighbor_atom_keys(gra)
    sites_lst = [site1]
    sites = site1

    # add forming key to reaction bonds that are forming double bonds
    if site2 is not None:
        if site1[1] == site2[0]:
            # second site is the breaking O=O bond
            sites.extend(site2)
            sites_lst.append(site2)
            brk_bnd = frozenset({site2[1], site2[2]})
            bnd_ord = bnds[brk_bnd]
            #bnds[brk_bnd] = (bnd_ord + 0.9, None)
            bnds[brk_bnd] = bnd_ord + 0.9
        else:
            # second site is a forming pi bond
            sites.extend(site2)
            sites_lst.append(site2)
            frm_bnd = frozenset({site2[0], site2[1]})
            bnd_ord = bnds[frm_bnd]
            #bnds[frm_bnd] = (bnd_ord + 0.1, None)
            bnds[frm_bnd] = bnd_ord + 0.1

    # switch resonances so dbl bnd isn't in rction site
    if len(unsat_atms) > 2:
        for unsat_a in unsat_atms:
            if unsat_a in sites:
                for unsat_b in adj_atms[unsat_a]:
                    if unsat_b in unsat_atms and unsat_b not in sites:
                        for unsat_c in adj_atms[unsat_b]:
                            if (unsat_c in unsat_atms and unsat_c != unsat_a
                                    and unsat_c not in sites):
                                unsat_ab = frozenset({unsat_a, unsat_b})
                                order_ab = bnds[unsat_ab]
                                if order_ab == 2:
                                    unsat_bc = frozenset({unsat_b, unsat_c})
                                    order_bc = bnds[unsat_bc]
                                    bnds[unsat_ab] = order_ab - 1
                                    bnds[unsat_bc] = order_bc + 1
                                else:
                                    unsat_bc = frozenset({unsat_b, unsat_c})
                                    order_bc = bnds[unsat_bc]
                                    bnds[unsat_ab] = order_ab
                                    bnds[unsat_bc] = order_bc + 1
    # fix the hydrogen valence of radical atms
    for atm in unsat_atms_dct:
        atm_vals[atm] -= unsat_atms_dct[atm]

    for site in sites_lst:
        abs_atm = site[0]
        trans_atm = site[1]
        don_atm = site[2]
        if trans_atm not in adj_atms[abs_atm]:
            adj_atms[abs_atm] = frozenset(
                {*list(adj_atms[abs_atm]), trans_atm})
        if trans_atm not in adj_atms[don_atm]:
            adj_atms[don_atm] = frozenset(
                {*list(adj_atms[don_atm]), trans_atm})
    #    bnd_ords[brk_bnd] = frozenset({list(bnd_ords[frm_bnd])[0] - 0.1})
    return atms, bnds, atm_vals, adj_atms, unsat_atms


def remove_zero_order_bnds(gra):
    """ Remove bonds of zero-order from a molecular graph.
    """

    atms, bnds = gra
    new_bnds = {}
    for bnd in bnds:
        if bnds[bnd][0] > 0:
            new_bnds[bnd] = bnds[bnd]

    return (atms, new_bnds)


def remove_hyd_from_adj_atms(atms, adj_atms, othersite=(), other_adj=()):
    """ Removes H atoms from all atms adjacent to a set of atoms requested
        by the user
    """

    new_adj_atms = ()
    for atm in adj_atms:
        if atms[atm][0] != 'H' and atm not in othersite:
            if atm not in other_adj:
                new_adj_atms += (atm,)

    return new_adj_atms


def remove_hyd_from_adj_atms2(atms, adj_atms_dct, extended_site):
    """ Removes H atoms from all atms adjacent to a set of atoms requested
        by the user
    """

    new_adj_atms_dct = {}
    for atm in adj_atms_dct:
        new_adj_atms = ()
        for adj_atm in adj_atms_dct[atm]:
            if atms[adj_atm][0] != 'H' and adj_atm not in extended_site:
                new_adj_atms += (adj_atm,)

        new_adj_atms_dct[atm] = new_adj_atms
    return new_adj_atms_dct

# GRAPH SPLITTING FUNCTIONS
def split_radradabs_gras(gras):
    """ Split a graph from radical-radical abstraction TS into the constituent
        reactant/products graphs.
    """

    rct_ichs = []
    prd_ichs = []
    atms, bnd_ords = gras
    atms = atms.copy()
    bnd_ords = bnd_ords.copy()
    for bnd_ord in bnd_ords:
        order, tmp = bnd_ords[bnd_ord]
        if abs(np.floor(order) - (order - 0.1)) < 0.01:
            bnd_ords[bnd_ord] = (round(order + 0.9, 1), tmp)
            atmai, atmbi = bnd_ord
            if abs(np.floor(atms[atmbi][1]) - (atms[atmbi][1]-0.9)) < 0.01:
                atmbi, atmai = atmai, atmbi
            atma = list(atms[atmai])
            atmb = list(atms[atmbi])
            atmb[1] = np.floor(atmb[1])
            atma[1] = np.floor(atma[1])
            atms[atmai] = tuple(atma)
            atms[atmbi] = tuple(atmb)
            for atmi in atms:
                if frozenset({atmi, atmbi}) in bnd_ords and atmi != atmai:
                    order, tmp = bnd_ords[frozenset({atmbi, atmi})]
                    if abs(np.floor(order) - (order - 0.9)) < 0.01:
                        atm = list(atms[atmi])
                        atm[1] = np.floor(atm[1])
                        atms[atmi] = tuple(atm)
                        bnd_ords[frozenset({atmbi, atmi})] = (round(
                            order - 0.9, 1), tmp)
                if frozenset({atmi, atmai}) in bnd_ords and atmi != atmbi:
                    order, tmp = bnd_ords[frozenset({atmai, atmi})]
                    if abs(np.floor(order) - (order - 0.9)) < 0.01:
                        atm = list(atms[atmi])
                        atm[1] = np.floor(atm[1])
                        atms[atmi] = tuple(atm)
                        bnd_ords[frozenset({atmai, atmi})] = (round(
                            order - 0.9, 1), tmp)
            rct_gra = remove_zero_order_bnds((atms, bnd_ords))
            atms, bnd_ords = rct_gra
    rct_gras = automol.graph.connected_components(rct_gra)
    for rgra in rct_gras:
        rct_ichs.append(automol.graph.chi(rgra))
    if len(rct_ichs) > 1:
        rct_ichs = automol.chi.sorted_(rct_ichs)
    atms, bnd_ords = gras
    for bnd_ord in bnd_ords:
        order, tmp = bnd_ords[bnd_ord]
        if abs(np.floor(order) - (order - 0.1)) < 0.01:
            bnd_ords[bnd_ord] = (round(order - 0.1, 1), tmp)
            atmai, atmbi = bnd_ord
            if abs(np.floor(atms[atmbi][1]) - (atms[atmbi][1]-0.9)) < 0.01:
                atmbi, atmai = atmai, atmbi
            atma = list(atms[atmai])
            atmb = list(atms[atmbi])
            atmb[1] = np.floor(atmb[1])
            atma[1] = np.floor(atma[1])
            atms[atmai] = tuple(atma)
            atms[atmbi] = tuple(atmb)
            for atmi in atms:
                if frozenset({atmi, atmbi}) in bnd_ords and atmi != atmai:
                    order, tmp = bnd_ords[frozenset({atmbi, atmi})]
                    if abs(np.floor(order) - (order - 0.9)) < 0.01:
                        atm = list(atms[atmi])
                        atm[1] = np.floor(atm[1])
                        atms[atmi] = tuple(atm)
                        atms[atmi] = tuple(atm)
                        order, tmp = bnd_ords[frozenset({atmbi, atmi})]
                        bnd_ords[frozenset({atmbi, atmi})] = (round(
                            order + 0.1, 1), tmp)
                if frozenset({atmi, atmai}) in bnd_ords and atmi != atmbi:
                    order, tmp = bnd_ords[frozenset({atmai, atmi})]
                    if abs(np.floor(order) - (order - 0.9)) < 0.01:
                        atm = list(atms[atmi])
                        atm[1] = np.floor(atm[1])
                        atms[atmi] = tuple(atm)
                        order, tmp = bnd_ords[frozenset({atmai, atmi})]
            prd_gra = remove_zero_order_bnds((atms, bnd_ords))
            atms, bnd_ords = prd_gra
    prd_gras = automol.graph.connected_components(prd_gra)
    for pgra in prd_gras:
        prd_ichs.append(automol.graph.chi(pgra))
    if len(prd_ichs) > 1:
        prd_ichs = automol.chi.sorted_(prd_ichs)
    return (rct_ichs, prd_ichs)


def _fix_sig_fig_issues(bnd_ords, atms):
    for bnd_ord in bnd_ords:
        order, tmp = bnd_ords[bnd_ord]
        if abs(np.floor(order) - (order - 0.1)) < 0.01:
            order = np.floor(order) + 0.1
        elif abs(np.floor(order) - (order - 0.9)) < 0.01:
            order = np.floor(order) + 0.9
        bnd_ords[bnd_ord] = (order, tmp,)
    for atm_idx in atms:
        atm, val, tmp = atms[atm_idx]
        if abs(np.floor(val) - (val - 0.1)) < 0.01:
            val = np.floor(val) + 0.1
        elif abs(np.floor(val) - (val - 0.9)) < 0.01:
            val = np.floor(val) + 0.9
        atms[atm_idx] = (atm, val, tmp)
        return bnd_ords, atms


def split_gras(gras):
    atms = automol.graph.atoms(gras)
    bnds = automol.graph.bonds(gras)
    rct_gras = (atms.copy(), bnds.copy())
    prd_gras = (atms.copy(), bnds.copy())
    for bnd in bnds:
        order, _ = bnds[bnd]
        if abs(np.floor(order) - (order - 0.1)) < 0.01:
            new_ord = np.floor(order)
            if new_ord < 1:
                rct_gras = automol.graph.remove_bonds(
                    rct_gras, (bnd,))
            else:
                rct_gras = automol.graph.set_bond_orders(
                    rct_gras, {bnd: new_ord})
            new_ord = np.floor(order) + 1
            prd_gras = automol.graph.set_bond_orders(
                prd_gras, {bnd: new_ord})
        elif abs(np.floor(order) - (order - 0.9)) < 0.01:
            new_ord = np.floor(order)
            if new_ord < 1:
                prd_gras = automol.graph.remove_bonds(
                    prd_gras, (bnd,))
            else:
                prd_gras = automol.graph.set_bond_orders(
                    prd_gras, {bnd: new_ord})
            new_ord = np.floor(order) + 1
            rct_gras = automol.graph.set_bond_orders(
                rct_gras, {bnd: new_ord})
    rct_gras = automol.graph.connected_components(rct_gras)
    prd_gras = automol.graph.connected_components(prd_gras)
    rct_ichs = []
    prd_ichs = []
    for rgra in rct_gras:
        rct_ichs.append(automol.graph.chi(rgra))
    for pgra in prd_gras:
        prd_ichs.append(automol.graph.chi(pgra))
    if len(rct_ichs) > 1:
        rct_ichs = automol.chi.sorted_(rct_ichs)
    if len(prd_ichs) > 1:
        prd_ichs = automol.chi.sorted_(prd_ichs)
    return rct_ichs, prd_ichs



def simplify_gra_frags(frags):
    """ ?
    """

    new_frags = {}
    for i, frag in enumerate(frags.keys()):
        if abs(frags[frag]['coeff']) > 0.01:
            new_frags[i] = frags[frag]

    return new_frags


def remove_frm_bnd(gra, brk_key, frm_key):
    """ Remove the formning bond and add the breaking bond
        of a molecular graph.
    """

    bond_keys = automol.graph.bond_keys(gra)
    if brk_key and brk_key not in bond_keys:
        gra = automol.graph.add_bonds(gra, [brk_key])
    if frm_key and frm_key in bond_keys:
        gra = automol.graph.remove_bonds(gra, [frm_key])
    return gra


def add_appropriate_pi_bonds(gra, frm_key):
    """ Add pi bonds to graphs
    """

    adj_atms = automol.graph.atoms_neighbor_atom_keys(gra)
    unsat_atms_dct = automol.graph.atom_unpaired_electrons(gra)
    atms, bnd_ords = gra
    brk_key = frozenset({})
    unsat_atms = []
    for atm in unsat_atms_dct:
        if unsat_atms_dct[atm] > 0:
            unsat_atms.append(atm)
    for atmi in unsat_atms:
        for atmj in unsat_atms:
            if atmi > atmj:
                if atmi in adj_atms[atmj]:
                    key = frozenset({atmi, atmj})
                    if not key == frm_key:
                        brk_key = key
                        bnd, tmp = bnd_ords[key]
                        bnd_ords[key] = (bnd + .9, tmp)

    return (atms, bnd_ords), brk_key


def elimination_second_forming_bond(gra, brk_key1, brk_key2):
    """ a
    """

    frm_bnd2 = frozenset({})
    adj_atms = automol.graph.atoms_neighbor_atom_keys(gra)
    for atm1 in brk_key1:
        for atm2 in brk_key2:
            if atm2 in adj_atms[atm1]:
                frm_bnd2 = [atm1, atm2]
                frm_bnd2.sort()
                frm_bnd2 = frozenset(frm_bnd2)

    return frm_bnd2


def ring_forming_forming_bond(gra, brk_key):
    """ Add in missing forming bond for ring forming scission reactions
    """
    frm_key = frozenset({})
    adj_atms = automol.graph.atoms_neighbor_atom_keys(gra)
    rad_atms = list(automol.graph.radical_atom_keys(gra, sing_res=True))
    form_atm1 = rad_atms[0]
    for break_atm in brk_key:
        if adj_atms[break_atm] > 1:
            form_atm2 = break_atm
            frm_key = frozenset({form_atm1, form_atm2})
    return frm_key


def elimination_find_brk_bnds(gra, frm_key):
    """ a
    """

    brk_key1 = frozenset({})
    brk_key2 = frozenset({})
    adj_atms = automol.graph.atoms_neighbor_atom_keys(gra)
    atms, _ = gra
    atm1, atm2 = frm_key
    atm3, atm4 = list(adj_atms[atm1])[0], list(adj_atms[atm2])[0]
    if atms[atm1][0] == 'H':
        brk_key1 = [atm1, atm3]
    elif atms[atm1][0] == 'O':
        for atm5 in adj_atms[atm3]:
            if atm5 != atm1:
                brk_key1 = [atm3, atm5]
    if atms[atm2][0] == 'H':
        brk_key2 = [atm2, atm4]
    elif atms[atm2][0] == 'O':
        for atm6 in adj_atms[atm4]:
            if atm6 != atm2:
                brk_key2 = [atm4, atm6]
    brk_key1.sort()
    brk_key2.sort()

    return frozenset(brk_key1), frozenset(brk_key2)


def split_bnd_keys(bnd_keys):
    """ Obtain the indiviual keys of a forming/breaking bond keys.
    """

    bnd_key1 = None
    bnd_key2 = None
    bnd_keys = list(bnd_keys)
    if len(bnd_keys) > 0:
        bnd_key1 = bnd_keys[0]
        if len(bnd_keys) > 1:
            bnd_key2 = bnd_keys[1]

    return bnd_key1, bnd_key2


def remove_dummies(zma, frm_key, brk_key, geo=None):
    """get zma and bond key idxs without dummy atoms
    """

    zgeo = automol.zmat.geometry(zma)
    brk_key2 = None
    if isinstance(brk_key, list):
        brk_key, brk_key2 = brk_key
    dummy_idxs = automol.geom.dummy_atom_indices(zgeo)
    for idx in dummy_idxs:
        if frm_key:
            frm1, frm2 = frm_key
            if idx < frm1:
                frm1 -= 1
            if idx < frm2:
                frm2 -= 1
            frm_key = frozenset({frm1, frm2})
        if brk_key:
            brk1, brk2 = brk_key
            if idx < brk1:
                brk1 -= 1
            if idx < brk2:
                brk2 -= 1
            brk_key = frozenset({brk1, brk2})
        if brk_key2:
            brk3, brk4 = brk_key2
            if idx < brk3:
                brk3 -= 1
            if idx < brk4:
                brk4 -= 1
            brk_key2 = frozenset({brk3, brk4})
    if not geo:
        geo = automol.geom.without_dummy_atoms(zgeo)
    gra = automol.geom.graph(geo)

    return gra, frm_key, brk_key, brk_key2
