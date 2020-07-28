"""
Handle construction and manipulation of various info objects
used throughout moldriver
"""

import sys
import autofile
from lib.filesys.build import zma_fs_from_prefix


# Handle info objects for theory
def get_es_info(method, thy_dct):
    """
    Turn es dictionary into theory info array
    """
    if method == 'input':
        ret = ['input_geom', None, None, None]
    else:
        ret = get_thy_info(method, thy_dct)
    return ret


def get_thy_info(method, thy_dct):
    """ convert theory level dictionary to theory information array
    """
    method_dct = thy_dct.get(method, None)
    err_msg = ''
    info = ['program', 'method', 'basis', 'orb_res']
    if method_dct is not None:
        for i, inf in enumerate(info):
            if inf in method_dct:
                info[i] = method_dct.get(inf, None)
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
    props = ['inchi', 'charge', 'mult']
    for i, prop in enumerate(props):
        if prop in spc_dct:
            props[i] = spc_dct[prop]
        else:
            err_msg = prop
    if err_msg:
        print('ERROR: No {} found'.format(err_msg))
    return props


# Handle info objects for reactions
def rxn_info(reacs, prods, spc_dct, rxn_mul='low'):
    """ prepare rxn info and reverse the reactants and products
        if reaction is endothermic
    """
    rxn_ichs = [[], []]
    rxn_chgs = [[], []]
    rxn_muls = [[], []]
    for spc in reacs:
        rxn_ichs[0].append(spc_dct[spc]['inchi'])
        rxn_chgs[0].append(spc_dct[spc]['charge'])
        rxn_muls[0].append(spc_dct[spc]['mult'])
    for spc in prods:
        rxn_ichs[1].append(spc_dct[spc]['inchi'])
        rxn_chgs[1].append(spc_dct[spc]['charge'])
        rxn_muls[1].append(spc_dct[spc]['mult'])
    rxn_ichs, rxn_chgs, rxn_muls = autofile.schema.sort_together(
        rxn_ichs, rxn_chgs, rxn_muls)
    _, ts_mul_low, ts_mul_high = rxn_chg_mult(rxn_muls, rxn_chgs)
    if rxn_mul == 'low':
        mul = ts_mul_low
    else:
        mul = ts_mul_high

    return [rxn_ichs, rxn_chgs, rxn_muls, mul]


def rxn_info2(reacs, prods, spc_dct, rxn_mul='low'):
    """ prepare rxn info and reverse the reactants and products
        if reaction is endothermic
    """
    rxn_ichs = [[], []]
    rxn_chgs = [[], []]
    rxn_muls = [[], []]
    for spc in reacs:
        rxn_ichs[0].append(spc_dct[spc]['inchi'])
        rxn_chgs[0].append(spc_dct[spc]['charge'])
        rxn_muls[0].append(spc_dct[spc]['mult'])
    for spc in prods:
        rxn_ichs[1].append(spc_dct[spc]['inchi'])
        rxn_chgs[1].append(spc_dct[spc]['charge'])
        rxn_muls[1].append(spc_dct[spc]['mult'])

    return rxn_ichs, rxn_chgs, rxn_muls


def rxn_chg_mult(rxn_muls, rxn_chgs):
    """ evaluate the ts multiplicity from the multiplicities
        of the reactants and products
    """
    nrcts, nprds = len(rxn_muls[0]), len(rxn_muls[1])

    # Set the multiplicities
    rct_spin_sum, prd_spin_sum = 0, 0
    rct_muls, prd_muls = [], []
    if nrcts == 1 and nprds == 1:
        ts_mul_low = max(rxn_muls[0][0], rxn_muls[1][0])
        ts_mul_high = ts_mul_low
        # rad_rad = False
    else:
        for rct_mul in rxn_muls[0]:
            rct_spin_sum += (rct_mul - 1.)/2.
            rct_muls.append(rct_mul)
        for prd_mul in rxn_muls[1]:
            prd_spin_sum += (prd_mul - 1.)/2.
            prd_muls.append(prd_mul)
        # rct_chk = bool(min(rct_muls) == 1 or nrcts == 1)
        # prd_chk = bool(min(prd_muls) == 1 or nprds == 1)
        # if rct_chk and prd_chk:
        #     rad_rad = False
        ts_mul_low = min(rct_spin_sum, prd_spin_sum)
        ts_mul_low = int(round(2*ts_mul_low + 1))
        ts_mul_high = max(rct_spin_sum, prd_spin_sum)
        ts_mul_high = int(round(2*ts_mul_high + 1))

    # Set the charges
    ts_chg = 0
    for rct_chg in rxn_chgs[0]:
        ts_chg += rct_chg

    return ts_chg, ts_mul_low, ts_mul_high


def cnf_fs_zma_geo(filesys, locs):
    """ Get the geometry and zmatrix from a filesystem
    """

    # Read the zma
    zma_fs, _ = zma_fs_from_prefix(
        filesys[-1].path(locs), zma_idxs=[0])
    if zma_fs[-1].file.zmatrix.exists([0]):
        zma = zma_fs[-1].file.zmatrix.read([0])
    else:
        zma = None

    # Read the geom
    if filesys[-1].file.geometry.exists(locs):
        geo = filesys[-1].file.geometry.read(locs)
    else:
        geo = None

    # Check
    if zma is None and geo is None:
        print('*ERROR: Neither a Z-Matrix or a Cartesian Geom exists level')
        sys.exit()

    return zma, geo
