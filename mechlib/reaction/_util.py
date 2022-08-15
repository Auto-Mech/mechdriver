""" Util for rxn objects
"""

import automol


def reverse_ts_zmatrix(zrxn):
    """ reverses a TS zma to be the other direction of reaction
    """
    back_zma = None
    back_zrxn = None
    nonreversible = [
        automol.par.ReactionClass.Typ.ELIMINATION,
        automol.par.ReactionClass.Typ.HOMOLYT_SCISSION,
        automol.par.ReactionClass.Typ.RING_FORM_SCISSION,
        automol.par.ReactionClass.Typ.SUBSTITUTION]

    if zrxn.class_ not in nonreversible:
        back_zrxn = automol.reac.reverse(zrxn)
    if back_zrxn is not None:
        rct_gras = automol.reac.reactant_graphs(back_zrxn)
        rct_geos = tuple([automol.graph.geometry(rgra) for rgra in rct_gras])
        rct_idxs, _ = back_zrxn.sort_order()
        back_zrxn = automol.reac.standard_keys(back_zrxn)
        rct_geos = tuple(map(rct_geos.__getitem__, rct_idxs))
        ts_geo = automol.reac.ts_geometry(back_zrxn, rct_geos, log=False)
        back_zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(
            back_zrxn, ts_geo)
        back_zrxn = automol.reac.relabel_for_zmatrix(
            back_zrxn, zma_keys, dummy_key_dct)
    return back_zrxn, back_zma


def zmatrix_conversion_keys(gra_forw, gra_back):
    """ returns the dictionary that transforms the forward graph 
        into the backward graph
    """
    iso_gra = convert_breaking_to_forming(gra_back)
    iso_dct = automol.graph.isomorphism(gra_forw, iso_gra)
    return iso_dct


def convert_breaking_to_forming(gra):
    bond_orders = automol.graph.bond_orders(gra)
    update_dct = {}
    for bond, order in bond_orders.items():
        if isinstance(order, float):
            if str(order)[-2:] == '.9':
                order -= .8
                update_dct[bond] = order
            if str(order)[-2:] == '.1':
                order += .8
                update_dct[bond] = order
    gra = automol.graph.set_bond_orders(gra, update_dct)
    return gra
