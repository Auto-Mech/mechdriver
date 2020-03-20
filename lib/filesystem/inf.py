"""
Build paths and file systesm given species and theory
'info' objects
"""

import autofile
from automol.mult.ts import _low as tslow
from automol.mult.ts import _high as tshigh
from lib.phydat import phycon
from lib.filesystem import read as fsread


def get_thy_info(method, thy_dct):
    """ convert theory level dictionary to theory information array
    """
    method_dct = thy_dct[method]
    err_msg = ''
    info = ['program', 'method', 'basis', 'orb_res']
    for i, inf in enumerate(info):
        if inf in method_dct:
            info[i] = method_dct[inf]
        else:
            err_msg = inf
    if err_msg:
        print('ERROR: No {} found'.format(err_msg))
    return info


def get_es_info(method, thy_dct):
    """
    Turn es dictionary in theory info array
    """
    if method == 'input':
        ret = ['input_geom', None, None, None]
    else:
        ret = get_thy_info(method, thy_dct)
    return ret


def get_spc_info(spc_dct):
    """ convert species dictionary to species_info array
    """
    err_msg = ''
    props = ['ich', 'chg', 'mul']
    for i, prop in enumerate(props):
        if prop in spc_dct:
            props[i] = spc_dct[prop]
        else:
            err_msg = prop
    if err_msg:
        print('ERROR: No {} found'.format(err_msg))
    return props


def rxn_info(reacs, prods, spc_dct, ts_mul='low'):
    """ prepare rxn info and reverse the reactants and products
        if reaction is endothermic
    """
    rxn_ichs = [[], []]
    rxn_chgs = [[], []]
    rxn_muls = [[], []]
    for spc in reacs:
        rxn_ichs[0].append(spc_dct[spc]['ich'])
        rxn_chgs[0].append(spc_dct[spc]['chg'])
        rxn_muls[0].append(spc_dct[spc]['mul'])
    for spc in prods:
        rxn_ichs[1].append(spc_dct[spc]['ich'])
        rxn_chgs[1].append(spc_dct[spc]['chg'])
        rxn_muls[1].append(spc_dct[spc]['mul'])
    rxn_ichs, rxn_chgs, rxn_muls = autofile.system.sort_together(
        rxn_ichs, rxn_chgs, rxn_muls)
    low_mul, high_mul, mul, chg = rxn_chg_mult(
        rxn_muls, rxn_chgs, ts_mul=ts_mul)

    return [rxn_ichs, rxn_chgs, rxn_muls, mul]


def assess_rxn_exo(reacs, prods, spc_dct, thy_info, save_prefix):
    """ Check the directionality of the reaction
    """
    [rxn_ichs, rxn_chgs, rxn_muls, _] = rxn_info(reacs, prods, spc_dct)
    rxn_exo = fsread.reaction_energy(
        save_prefix, rxn_ichs, rxn_chgs, rxn_muls, thy_info)
    print('reaction is {:.2f} endothermic'.format(rxn_exo*phycon.EH2KCAL))
    if rxn_exo > 0:
        rxn_ichs = rxn_ichs[::-1]
        rxn_chgs = rxn_chgs[::-1]
        rxn_muls = rxn_muls[::-1]
        reacs, prods = prods, reacs

    return reacs, prods


def rxn_chg_mult(rxn_muls, rxn_chgs, ts_mul='low'):
    """ prepare rxn info and reverse the reactants and products
        if reaction is endothermic
    """
    assert ts_mul in ('low', 'high')
    low_mul = min(tslow(rxn_muls[0]), tslow(rxn_muls[1]))
    high_mul = max(tshigh(rxn_muls[0]), tshigh(rxn_muls[1]))
    mul = low_mul if ts_mul == 'low' else high_mul

    print('aa', rxn_chgs)
    chg = 0
    for rchg in rxn_chgs[0]:
        chg += rchg

    return low_mul, high_mul, mul, chg
