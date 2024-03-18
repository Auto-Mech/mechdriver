""" Util for rxn objects
"""

import automol
from automol import ReactionClass


def reverse_ts_zmatrix(zrxn):
    """ reverses a TS zma to be the other direction of reaction
    """
    back_zma = None
    back_zrxn = None

    if ReactionClass.is_reversible(automol.reac.class_(zrxn)):
        try:
            back_zrxn = automol.reac.reverse(zrxn)
            back_zrxn = automol.reac.with_structures(zrxn, "zmat")
            back_zma = automol.reac.ts_structure(zrxn)
        except:
            back_zma = None
            back_zrxn = None
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
