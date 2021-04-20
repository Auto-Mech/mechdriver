"""
 New reaction ID code
"""

import os
import autofile
import automol
from ioformat import ptt
import chemkin_io
from mechanalyzer.inf import rxn as rinfo
from mechanalyzer.inf import thy as tinfo
from phydat import phycon
from mechlib import filesys
from mechlib.filesys import build_fs


CLA_INP = 'inp/class.csv'


def build_reaction(rxn_info, ini_thy_info, zma_locs, save_prefix):
    """ reactin
    """

    # Try to build the Z-Matrix reaction object or identify from scratch
    zrxn, zma = _read_from_filesys(
        rxn_info, ini_thy_info, zma_locs, save_prefix)
    if zrxn is None:
        print('    Identifying class')
        zrxns, zmas = _id_reaction(rxn_info)
    else:
        zrxns = (zrxn,)
        zmas = (zma,)
        print('    Reading from fileysystem')

    print('    Reaction class identified as: {}'.format(zrxns[0].class_))
    # print('    Reaction class identified as: {}'.format(zrxn.class_))

    return zrxns, zmas


def _read_from_filesys(rxn_info, ini_thy_info, zma_locs, save_prefix):
    """ Check if reaction exists in the filesystem and has been identified
    """

    zrxn, zma = None, None

    sort_rxn_info = rinfo.sort(rxn_info, scheme='autofile')
    ts_info = rinfo.ts_info(rxn_info)
    mod_ini_thy_info = tinfo.modify_orb_label(ini_thy_info, ts_info)

    rxn_fs = autofile.fs.reaction(save_prefix)
    if rxn_fs[-1].exists(sort_rxn_info):
        _, cnf_save_fs = build_fs(
            save_prefix, save_prefix, 'CONFORMER',
            rxn_locs=sort_rxn_info,
            thy_locs=mod_ini_thy_info[1:],
            # this needs to be fixed for any case with more than one TS
            ts_locs=(0,))

        _, ini_min_cnf_path = filesys.mincnf.min_energy_conformer_locators(
            cnf_save_fs, mod_ini_thy_info)
        if ini_min_cnf_path:
            zma_fs = autofile.fs.zmatrix(ini_min_cnf_path)
            if zma_fs[-1].file.reaction.exists(zma_locs):
                zrxn = zma_fs[-1].file.reaction.read(zma_locs)
                zma = zma_fs[-1].file.zmatrix.read(zma_locs)

    return zrxn, zma


def _id_reaction(rxn_info):
    """ Identify the reaction and build the object
    """

    [rxn_ichs, _, _, _] = rxn_info   # replace with mechanalyzer grab
    rct_ichs, prd_ichs = rxn_ichs[0], rxn_ichs[1]

    zrxn_objs = automol.reac.rxn_objs_from_inchi(
        rct_ichs, prd_ichs, indexing='zma')
    zrxns = tuple(obj[0] for obj in zrxn_objs)
    zmas = tuple(obj[1] for obj in zrxn_objs)

    return zrxns, zmas


# from direction
def set_reaction_direction(reacs, prods, rxn_info,
                           thy_info, ini_thy_info, save_prefix,
                           direction='forw'):
    """ Set the reaction of a direction
    """

    if direction == 'forw':
        print('    User requested forward direction.')
    elif direction == 'back':
        print('    User requested reverse direction, flipping reaction.')
        reacs, prods = prods, reacs
    elif direction == 'exo':
        print('    User requested exothermic direction.',
              'Checking energies...')
        reacs, prods = assess_rxn_ene(
            reacs, prods, rxn_info, thy_info, ini_thy_info, save_prefix)

    print('    Running reaction as:')
    print('      {} = {}'.format('+'.join(reacs), '+'.join(prods)))

    return reacs, prods


# Functions for the exothermicity check
def assess_rxn_ene(reacs, prods, rxn_info,
                   thy_info, ini_thy_info, save_prefix):
    """ Check the directionality of the reaction
    """

    rxn_ene = reaction_energy(rxn_info, thy_info, ini_thy_info, save_prefix)
    method1, method2 = thy_info, ini_thy_info
    if rxn_ene is None:
        rxn_ene = reaction_energy(
            rxn_info, ini_thy_info, ini_thy_info, save_prefix)
        method1, method2 = ini_thy_info, ini_thy_info

    print('    Reaction energy is {:.2f} at {}//{} level'.format(
        rxn_ene*phycon.EH2KCAL, method1[1], method2[1]))

    if rxn_ene > 0:
        reacs, prods = prods, reacs
        rxn_info = rinfo.reverse(rxn_info)
        print('    Reaction is endothermic, flipping reaction.')

    return reacs, prods


def reaction_energy(rxn_info, sp_thy_info, geo_thy_info, save_prefix):
    """ reaction energy """

    rct_enes = reagent_energies(
        'reacs', rxn_info, sp_thy_info, geo_thy_info, save_prefix)
    prd_enes = reagent_energies(
        'prods', rxn_info, sp_thy_info, geo_thy_info, save_prefix)

    if rct_enes is not None and prd_enes is not None:
        rxn_ene = sum(prd_enes) - sum(rct_enes)
    else:
        rxn_ene = None

    return rxn_ene


def reagent_energies(rgt, rxn_info, sp_thy_info, geo_thy_info, save_prefix):
    """ reagent energies """

    assert rgt in ('reacs', 'prods')

    enes = []
    for rgt_info in rinfo.rgts_info(rxn_info, rgt):

        # Set filesys
        spc_save_fs = autofile.fs.species(save_prefix)
        spc_save_path = spc_save_fs[-1].path(rgt_info)

        mod_geo_thy_info = tinfo.modify_orb_label(geo_thy_info, rgt_info)
        mod_sp_thy_info = tinfo.modify_orb_label(sp_thy_info, rgt_info)
        thy_save_fs = autofile.fs.theory(spc_save_path)
        thy_save_path = thy_save_fs[-1].path(mod_geo_thy_info[1:4])
        cnf_save_fs = autofile.fs.conformer(thy_save_path)
        min_locs, min_path = filesys.min_energy_conformer_locators(
            cnf_save_fs, mod_geo_thy_info)
        # Read energy
        ene = None
        if min_locs:
            # Create run fs if that directory has been deleted to run the jobs
            sp_fs = autofile.fs.single_point(min_path)
            if sp_fs[-1].file.energy.exists(mod_sp_thy_info[1:4]):
                ene = sp_fs[-1].file.energy.read(mod_sp_thy_info[1:4])
        enes.append(ene)

    if any(ene is None for ene in enes):
        enes = None

    return enes
