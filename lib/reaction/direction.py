"""
 Functions to handle direction of a reaction
"""

import os
import automol
import autofile
import chemkin_io
from ioformat import remove_whitespace
# from routines.es._routines import geom
from lib.phydat import phycon
from lib.filesys.mincnf import min_energy_conformer_locators
from lib.filesys.inf import modify_orb_restrict
from lib.amech_io.parser import ptt


CLA_INP = 'inp/class.csv'


# Main direction function
def set_reaction_direction(reacs, prods, spc_dct, cla_dct,
                           thy_info, ini_thy_info, save_prefix,
                           direction='forw'):
    """ Set the reaction of a direction
    """

    # Check if reaction is present in the class direction
    if cla_dct:
        given_class, flip_rxn = set_class_with_dct(
            cla_dct, reacs, prods)
        if flip_rxn:
            reacs, prods = prods, reacs
    else:
        given_class = None

    # If no class, given set direction to requested direction
    if given_class is not None:
        print('    Reaction present in class dct, Setting direction to that.')
    else:
        if direction == 'forw':
            print('    User requested forward direction.')
        elif direction == 'back':
            print('    User requested reverse direction, flipping reaction.')
            reacs, prods = prods, reacs
        elif direction == 'exo':
            print('    User requested exothermic direction.',
                  'Checking energies...')
            reacs, prods = assess_rxn_ene(
                reacs, prods, spc_dct, thy_info, ini_thy_info, save_prefix)

    print('    Running reaction as:')
    print('      {} = {}'.format('+'.join(reacs), '+'.join(prods)))

    return reacs, prods, given_class


# Handle setting reaction directions with the class dictionary
def set_class_with_dct(cla_dct, reacs, prods):
    """ set the class using the class dictionary
    """
    rxn = (reacs, prods)
    rxn_rev = (prods, reacs)
    if rxn in cla_dct:
        given_class = cla_dct[rxn]
        flip_rxn = False
    elif rxn_rev in cla_dct:
        given_class = cla_dct[rxn_rev]
        flip_rxn = True
    else:
        given_class = None
        flip_rxn = False

    return given_class, flip_rxn


def parse_rxn_class_file(job_path):
    """ Read the class dictionary
    """

    if os.path.exists(os.path.join(job_path, CLA_INP)):
        print('  class.dat found. Reading contents...')
        cla_str = ptt.read_inp_str(job_path, CLA_INP, remove_comments='#')
        cla_dct = _build_cla_dct(cla_str)
    else:
        print('  No class.dat found.')
        cla_dct = {}

    return cla_dct


def _build_cla_dct(cla_str):
    """ read file
    """
    cla_dct = {}
    cla_str = remove_whitespace(cla_str)
    for line in cla_str.splitlines():
        # try:
        [rxn_line, rclass] = line.split('||')
        reacs = chemkin_io.parser.reaction.reactant_names(rxn_line)
        prods = chemkin_io.parser.reaction.product_names(rxn_line)
        cla_dct[(reacs, prods)] = rclass
        # except:
        #     print('*ERROR: Error in formatting line')
        #     print(line)
        #     sys.exit()

    return cla_dct


# Functions for the exothermicity check
def assess_rxn_ene(reacs, prods, spc_dct, thy_info, ini_thy_info, save_prefix):
    """ Check the directionality of the reaction
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
        reacs, prods = prods, reacs
        print('    Reaction is endothermic, flipping reaction.')

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
        min_cnf_locs = min_energy_conformer_locators(cnf_save_fs, thy_lvl)
        ene = cnf_save_fs[-1].file.energy.read(min_cnf_locs)
        enes.append(ene)
    return enes


def get_zmas(
        reacs, prods, spc_dct, ini_thy_info, save_prefix, run_prefix,
        kickoff_size, kickoff_backward):
    """get the zmats for reactants and products using the initial level of theory
    """
    if len(reacs) > 2:
        ich = spc_dct[reacs[-1]]['inchi']
        ichgeo = automol.inchi.geometry(ich)
        ichzma = automol.geom.zmatrix(ichgeo)
        reacs = reacs[:-1]
    elif len(prods) > 2:
        ich = spc_dct[prods[-1]]['inchi']
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
        spc_info = [spc_dct[spc]['inchi'],
                    spc_dct[spc]['charge'],
                    spc_dct[spc]['mult']]
        ini_thy_lvl = modify_orb_restrict(spc_info, ini_thy_info)
        spc_save_fs = autofile.fs.species(save_prefix)
        spc_save_fs[-1].create(spc_info)
        spc_save_path = spc_save_fs[-1].path(spc_info)
        ini_thy_save_fs = autofile.fs.theory(spc_save_path)
        ini_thy_save_path = ini_thy_save_fs[-1].path(ini_thy_lvl[1:4])
        cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)
        cnf_save_fs_lst.append(cnf_save_fs)
        min_cnf_locs = min_energy_conformer_locators(cnf_save_fs, ini_thy_lvl)
        if min_cnf_locs:
            geo = cnf_save_fs[-1].file.geometry.read(min_cnf_locs)
        # else:
        #     spc_run_fs = autofile.fs.species(run_prefix)
        #     spc_run_fs[-1].create(spc_info)
        #     spc_run_path = spc_run_fs[-1].path(spc_info)
        #     ini_thy_run_fs = autofile.fs.theory(spc_run_path)
        #     ini_thy_run_path = ini_thy_run_fs[-1].path(ini_thy_lvl[1:4])
        #     cnf_run_fs = autofile.fs.conformer(ini_thy_run_path)
        #     run_fs = autofile.fs.run(ini_thy_run_path)
        #     run_fs[0].create()
        #     geo = geom.reference_geometry(
        #         spc_dct[spc], spc_info,
        #         ini_thy_lvl, ini_thy_lvl,
        #         ini_thy_run_fs, ini_thy_save_fs,
        #         ini_thy_save_fs,
        #         cnf_run_fs, cnf_save_fs,
        #         run_fs,
        #         opt_script_str, overwrite,
        #         kickoff_size=kickoff_size,
        #         kickoff_backward=kickoff_backward)
        spc_geos.append(geo)
    return spc_geos, cnf_save_fs_lst
