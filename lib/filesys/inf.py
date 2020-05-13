"""
Handle construction and manipulation various info objects
used through out moldriver
"""

import sys
import autofile
import automol
from automol.mult.ts import _low as tslow
from automol.mult.ts import _high as tshigh
from automol.zmatrix import shifted_standard_zmas_graphs as shift_gra
from routines.es._routines import geom
from lib.phydat import phycon
from lib.filesys.mincnf import min_energy_conformer_locators


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


def reaction_energy(save_prefix, rxn_ich, rxn_chg, rxn_mul, thy_level):
    """ reaction energy """
    rct_ichs, prd_ichs = rxn_ich
    rct_chgs, prd_chgs = rxn_chg
    rct_muls, prd_muls = rxn_mul
    rct_enes = reagent_energies(
        save_prefix, rct_ichs, rct_chgs, rct_muls, thy_level)
    prd_enes = reagent_energies(
        save_prefix, prd_ichs, prd_chgs, prd_muls, thy_level)
    return sum(prd_enes) - sum(rct_enes)


def reagent_energies(save_prefix, rgt_ichs, rgt_chgs, rgt_muls, thy_level):
    """ reagent energies """
    enes = []
    for rgt_ich, rgt_chg, rgt_mul in zip(rgt_ichs, rgt_chgs, rgt_muls):
        spc_save_fs = autofile.fs.species(save_prefix)
        rgt_info = [rgt_ich, rgt_chg, rgt_mul]
        spc_save_path = spc_save_fs[-1].path(rgt_info)

        thy_lvl = modify_orb_restrict(rgt_info, thy_level)
        thy_save_fs = autofile.fs.theory(spc_save_path)
        thy_save_path = thy_save_fs[-1].path(thy_lvl[1:4])
        cnf_save_fs = autofile.fs.conformer(thy_save_path)
        min_cnf_locs = min_energy_conformer_locators(cnf_save_fs)
        ene = cnf_save_fs[-1].file.energy.read(min_cnf_locs)
        enes.append(ene)
    return enes


def get_zmas(
        reacs, prods, spc_dct, ini_thy_info, save_prefix, run_prefix,
        kickoff_size, kickoff_backward):
    """get the zmats for reactants and products using the initial level of theory
    """
    if len(reacs) > 2:
        ich = spc_dct[reacs[-1]]['ich']
        ichgeo = automol.inchi.geometry(ich)
        ichzma = automol.geom.zmatrix(ichgeo)
        reacs = reacs[:-1]
    elif len(prods) > 2:
        ich = spc_dct[prods[-1]]['ich']
        ichgeo = automol.inchi.geometry(ich)
        ichzma = automol.geom.zmatrix(ichgeo)
        prods = prods[:-1]
    rct_geos, rct_cnf_save_fs_lst = get_geos(
        reacs, spc_dct, ini_thy_info, save_prefix, run_prefix, kickoff_size,
        kickoff_backward)
    prd_geos, prd_cnf_save_fs_lst = get_geos(
        prods, spc_dct, ini_thy_info, save_prefix, run_prefix, kickoff_size,
        kickoff_backward)
    rct_zmas = list(map(automol.geom.zmatrix, rct_geos))
    prd_zmas = list(map(automol.geom.zmatrix, prd_geos))
    if len(rct_zmas) > 2:
        rct_zmas.append(ichzma)
    if len(prd_zmas) > 2:
        prd_zmas.append(ichzma)
    return rct_zmas, prd_zmas, rct_cnf_save_fs_lst, prd_cnf_save_fs_lst


def get_geos(
        spcs, spc_dct, ini_thy_info, save_prefix, run_prefix, kickoff_size,
        kickoff_backward):
    """get geos for reactants and products using the initial level of theory
    """
    spc_geos = []
    cnf_save_fs_lst = []
    for spc in spcs:
        spc_info = [spc_dct[spc]['ich'],
                    spc_dct[spc]['chg'],
                    spc_dct[spc]['mul']]
        ini_thy_lvl = modify_orb_restrict(spc_info, ini_thy_info)
        spc_save_fs = autofile.fs.species(save_prefix)
        spc_save_fs[-1].create(spc_info)
        spc_save_path = spc_save_fs[-1].path(spc_info)
        ini_thy_save_fs = autofile.fs.theory(spc_save_path)
        ini_thy_save_path = ini_thy_save_fs[-1].path(ini_thy_lvl[1:4])
        cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)
        cnf_save_fs_lst.append(cnf_save_fs)
        min_cnf_locs = min_energy_conformer_locators(cnf_save_fs)
        if min_cnf_locs:
            geo = cnf_save_fs[-1].file.geometry.read(min_cnf_locs)
        else:
            spc_run_fs = autofile.fs.species(run_prefix)
            spc_run_fs[-1].create(spc_info)
            spc_run_path = spc_run_fs[-1].path(spc_info)
            ini_thy_run_fs = autofile.fs.theory(spc_run_path)
            ini_thy_run_path = ini_thy_run_fs[-1].path(ini_thy_lvl[1:4])
            cnf_run_fs = autofile.fs.conformer(ini_thy_run_path)
            run_fs = autofile.fs.run(ini_thy_run_path)
            run_fs[0].create()
            geo = geom.reference_geometry(
                spc_dct[spc], ini_thy_lvl, ini_thy_lvl,
                ini_thy_run_fs, ini_thy_save_fs,
                ini_thy_save_fs,
                cnf_run_fs, cnf_save_fs,
                run_fs,
                kickoff_size=kickoff_size,
                kickoff_backward=kickoff_backward,
                overwrite=False)
        spc_geos.append(geo)
    return spc_geos, cnf_save_fs_lst


def get_zma_geo(filesys, locs):
    """ Get the geometry and zmatrix from a filesystem
    """
    if filesys[-1].file.zmatrix.exists(locs):
        zma = filesys[-1].file.zmatrix.read(locs)
    else:
        zma = None

    if filesys[-1].file.geometry.exists(locs):
        geo = filesys[-1].file.geometry.read(locs)
    else:
        geo = None

    # Check
    if zma is None and geo is None:
        print('*ERROR: Neither a Z-Matrix or a Cartesian Geom exists level')
        sys.exit()

    return zma, geo


def min_dist_conformer_zma(dist_name, cnf_save_fs):
    """ locators for minimum energy conformer """
    cnf_locs_lst = cnf_save_fs[-1].existing()
    cnf_zmas = [cnf_save_fs[-1].file.zmatrix.read(locs)
                for locs in cnf_locs_lst]
    min_dist = 100.
    min_zma = []
    for zma in cnf_zmas:
        dist = automol.zmatrix.values(zma)[dist_name]
        if dist < min_dist:
            min_dist = dist
            min_zma = zma
    min_zma = [min_zma]
    return min_zma


def min_dist_conformer_zma_geo(dist_coords, cnf_save_fs):
    """ locators for minimum energy conformer """
    cnf_locs_lst = cnf_save_fs[-1].existing()
    cnf_zmas = [cnf_save_fs[-1].file.zmatrix.read(locs)
                for locs in cnf_locs_lst]
    min_dist = 100.
    min_zma = []
    for zma in cnf_zmas:
        zmas, _ = shift_gra([zma])
        zma = zmas[0]
        geo = automol.zmatrix.geometry(zma)
        dist = automol.geom.distance(geo, *list(dist_coords))
        if dist < min_dist:
            min_dist = dist
            min_zma = zma
    min_zma = [min_zma]
    return min_zma
