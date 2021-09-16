"""
 New reaction ID code
"""

import autofile
import automol
from mechanalyzer.inf import rxn as rinfo
from mechanalyzer.inf import thy as tinfo
from phydat import phycon
from mechlib import filesys


CLA_INP = 'inp/class.csv'

def build_reaction(rxn_info, ini_thy_info, zma_locs, save_prefix,
                   id_missing=True, re_id=False):
    """ For a given reaction, attempt to identify its reaction class obtain
        the corresponding Z-Matrix reaction object.

        Function will attempt to read the filesytem for the appropriate
        reaction object in the appropriate latyer
            save/RXN/THY/TS/CONFS/Z/

        If file is not located it will build Z-Matrices for the reactants
        and run the class identifier.

        :param rxn_info: Mechanalyzer reaction info object
        :type rxn_info: tuple(tuple(str/int))
        :param ini_thy_info: Mechanalyzer theory info object (input lvl)
        :type ini_thy_info: tuple(str)
        :param zma_locs: locs for Z filesys to search for reaction object
        :type zma_locs: tuple(int)
        :param save_prefix: root path to the save filesys
        :type save_prefix: str
    """

    zrxns, zmas, rclasses = None, (), ()

    # Try and read the reaction from filesys if requested
    if not re_id:
        zrxns, zmas = filesys.read.reactions(
            rxn_info, ini_thy_info, zma_locs, save_prefix)
        if zrxns is not None:
            # zrxns = (zrxn,)
            # zmas = (zma,)
            print('    Reading from fileysystem...')
    else:
        # unsafe without checking if zrxn id matches what is in save...
        print('    Requested Reidentification regardless of what is in SAVE')

    if zrxns is None:
        if id_missing:
            print('    Identifying class...')
            zrxns, zmas = _id_reaction(rxn_info)
            if zrxns is None:
                zrxns = 'MISSING-SKIP'
        else:
            zrxns = 'MISSING-ADD'

    if 'MISSING' not in zrxns:
        for zrxn in zrxns:
            rclasses += (_mod_class(zrxn.class_, rxn_info),)

        print('    Reaction class identified as: {}'.format(
            automol.par.string(rclasses[0])))

        print('    There are {} configuration(s) of transition state'.format(
            len(zrxns)))

    return zrxns, zmas, rclasses


def _id_reaction(rxn_info):
    """ Identify the reaction and build the object

        :param rxn_info: reaction info object
        :type rxn_info: mechanalyzer.inf.rxn object
        :rtype: (tuple(automol.Reaction object), tuple(automol.zmat object))
    """

    rxn_ichs = rinfo.value(rxn_info, 'inchi')
    rct_ichs, prd_ichs = rxn_ichs[0], rxn_ichs[1]

    zrxn_objs = automol.reac.rxn_objs_from_inchi(
        rct_ichs, prd_ichs, indexing='zma')

    zrxns, zmas = (), ()
    for objs in zrxn_objs:
        zrxn, zma, _, _ = objs
        zrxns += (zrxn,)
        zmas += (zma,)

    if not zrxns:
        zrxns, zmas = None, None

    return zrxns, zmas


def _mod_class(class_typ, rxn_info):
    """ Create the object containing the full description of the
        reaction class, including its type, spin state, and whether
        it is a radical-radical combination.

        :param class_typ: reaction class type
        :type class_typ: str
        :param rxn_info: ???
        :tpye rxn_info: ???
    """

    # Set the spin of the reaction to high/low
    _fake_class = (class_typ, '', '', False)
    if automol.par.need_spin_designation(_fake_class):
        ts_mul = rinfo.value(rxn_info, 'tsmult')
        high_mul = rinfo.ts_mult(rxn_info, rxn_mul='high')
        _spin = 'high-spin' if ts_mul == high_mul else 'low-spin'
    else:
        _spin = ''

    # Determine if it iss intersystem crossing
    # rxn_muls = rinfo.value(rxn_info, 'mult')
    # is_isc = automol.reac.intersystem_crossing(rxn_muls)

    return automol.par.reaction_class_from_data(
        class_typ, _spin, rinfo.radrad(rxn_info), False)


# from direction (loop over the reactions around split)
def set_reaction_direction(reacs, prods, rxn_info,
                           thy_info, ini_thy_info, save_prefix,
                           direction='forw'):
    """ Arrange the reactants and products in the order corresponding
        to the desired direction of the reaction. If the direction
        is the exothermic direction, than the species energies are read
        from the filesystem at the level of theory.

        :param reacs: list of names of the reactants
        :type reacs: tuple(str)
        :param prods: list of names of the products
        :type prods: tuple(str)
        :param rxn_info:
        :type rxn_info: tuple(tuple(str), tuple(int), tuple(int), int)
        :param thy_info: ???
        :type thy_info: ??
        :param ini_thy_info: ??
        :param direction: direction to set reaction to
        :type direction: str
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
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
    else:
        raise NotImplementedError

    print('    Running reaction as: {} = {}'.format(
        '+'.join(reacs), '+'.join(prods)))

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
            sp_fs = autofile.fs.single_point(min_path)
            if sp_fs[-1].file.energy.exists(mod_sp_thy_info[1:4]):
                ene = sp_fs[-1].file.energy.read(mod_sp_thy_info[1:4])
        enes.append(ene)

    if any(ene is None for ene in enes):
        enes = None

    return enes
