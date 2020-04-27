"""
Handle construction and manipulation various info objects
used through out moldriver
"""

import autofile
from automol.mult.ts import _low as tslow
from automol.mult.ts import _high as tshigh
from lib.phydat import phycon
from lib.filesys.read import reaction_energy


# Handle info objects for theory
def get_es_info(method, thy_dct):
    """
    Turn es dictionary in theory info array
    """
    if method == 'input':
        ret = ['input_geom', None, None, None]
    else:
        ret = get_thy_info(method, thy_dct)
    return ret


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


def modify_orb_restrict(spc_info, thy_info):
    """ append to the theory level the orb restricted stuff
    """
    orb_restr = _set_orbital_restriction_label(spc_info, thy_info)
    thy_info = thy_info[0:3]
    thy_info.append(orb_restr)

    return thy_info


def _set_orbital_restriction_label(spc_info, thy_info):
    """ orbital restriction logical
    """
    mul = spc_info[2]
    if thy_info[3] == 'RR':
        orb_restr = 'R'
    elif thy_info[3] == 'UU':
        orb_restr = 'U'
    elif thy_info[3] == 'RU':
        if mul == 1:
            orb_restr = 'R'
        else:
            orb_restr = 'U'
    return orb_restr


# Handle info objects for species
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


# Handle info objects for reactions
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
    _, _, mul, _ = rxn_chg_mult(
        rxn_muls, rxn_chgs, ts_mul=ts_mul)

    return [rxn_ichs, rxn_chgs, rxn_muls, mul]


def assess_rxn_ene(reacs, prods, spc_dct, thy_info, ini_thy_info, save_prefix):
    """ Check the directionality of the reaction
    """
    [rxn_ichs, rxn_chgs, rxn_muls, _] = rxn_info(reacs, prods, spc_dct)
    try:
        rxn_ene = reaction_energy(
            save_prefix, rxn_ichs, rxn_chgs, rxn_muls, thy_info)
        method = thy_info
    except TypeError:
        rxn_ene = reaction_energy(
            save_prefix, rxn_ichs, rxn_chgs, rxn_muls, ini_thy_info)
        method = ini_thy_info
    except AssertionError:
        rxn_ene = reaction_energy(
            save_prefix, rxn_ichs, rxn_chgs, rxn_muls, ini_thy_info)
        method = ini_thy_info
    except IOError:
        rxn_ene = reaction_energy(
            save_prefix, rxn_ichs, rxn_chgs, rxn_muls, ini_thy_info)
        method = ini_thy_info
    print('    Reaction energy is {:.2f} at {} level'.format(
        rxn_ene*phycon.EH2KCAL, method[1]))
    if rxn_ene > 0:
        rxn_ichs = rxn_ichs[::-1]
        rxn_chgs = rxn_chgs[::-1]
        rxn_muls = rxn_muls[::-1]
        reacs, prods = prods, reacs
        print('    Since reaction is endothermic, flipping reaction to')
        print('      {} = {}'.format('+'.join(reacs), '+'.join(prods)))

    return reacs, prods


def rxn_chg_mult(rxn_muls, rxn_chgs, ts_mul='low'):
    """ prepare rxn info and reverse the reactants and products
        if reaction is endothermic
    """
    assert ts_mul in ('low', 'high')
    low_mul = min(tslow(rxn_muls[0]), tslow(rxn_muls[1]))
    high_mul = max(tshigh(rxn_muls[0]), tshigh(rxn_muls[1]))
    mul = low_mul if ts_mul == 'low' else high_mul

    chg = 0
    for rchg in rxn_chgs[0]:
        chg += rchg

    return low_mul, high_mul, mul, chg
