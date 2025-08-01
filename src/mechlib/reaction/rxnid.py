"""
 New reaction ID code
"""

import autofile
import automol
from mechanalyzer.inf import rxn as rinfo
from mechanalyzer.inf import thy as tinfo
from phydat import phycon
from mechlib import filesys


def build_reaction(rxn_info, ini_thy_info, zma_locs, save_prefix,
                   id_missing=True, re_id=False, canonical=True):
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
    status = 'MISSING'

    # Try and read the reaction from filesys if requested
    if not re_id:
        print('    Trying to find class info from fileysystem...')
        zrxns, zmas = filesys.read.reactions(
            rxn_info, ini_thy_info, zma_locs, save_prefix)
        status = 'FOUND' if zrxns is not None else 'MISSING'
        if status == 'FOUND':
            print('    Info found. Reading from fileysystem...')
    else:
        # unsafe without checking if zrxn id matches what is in save...
        print('    Requested Reidentification regardless of what is in SAVE')

    # Try and identify reaction if not rxn obj found
    if status == 'MISSING':
        if id_missing:
            print('    No info found. Identifying class...')
            zrxns, zmas = _id_reaction(
                rxn_info, ini_thy_info, save_prefix, canonical=canonical)
            status = 'FOUND' if zrxns is not None else 'MISSING-SKIP'
        else:
            status = 'MISSING-ADD'

    # Build a tuple with the full description of the reaction class, if ID'd
    if status not in ('MISSING-SKIP', 'MISSING-ADD'):
        for zrxn in zrxns:
            rclasses += (_mod_class(automol.reac.class_(zrxn), rxn_info),)

        print('    Reaction class identified as: '
              f'{automol.ReactionInfo.string(rclasses[0])}')
        print(f'    There are {len(zrxns)} '
              'configuration(s) of transition state')

    return zrxns, zmas, rclasses, status


def _id_reaction(rxn_info, thy_info, save_prefix, canonical=True):
    """ Identify the reaction and build the object

        :param rxn_info: reaction info object
        :type rxn_info: mechanalyzer.inf.rxn object
        :rtype: (tuple(automol.Reaction object), tuple(automol.zmat object))
    """

    # Check the save filesystem for the reactant and product geometries
    rct_geos, prd_geos, rct_paths, prd_paths = reagent_geometries(
        rxn_info, thy_info, save_prefix, canonical=canonical)
    # Identify reactants and products from geoms or InChIs, depending
    # We automatically assess and add stereo to the reaction object, as needed
    if any(rct_geos) and any(prd_geos):
        print('    Reaction ID from geometries from SAVE filesys')
        for i, path in enumerate(rct_paths):
            print(f'     - reactant {i+1}: {path}')
        for i, path in enumerate(prd_paths):
            print(f'     - product {i+1}: {path}')
        zrxn_objs = automol.reac.from_geometries(
            rct_geos, prd_geos, struc_typ="zmat", stereo=True, enant=False)

        # AVC comment: Why is this necessary? Is this a bug I need to fix?
        if not zrxn_objs and len(rxn_info[2][0]) == 2:
            if rxn_info[2][0][0] > 1 and rxn_info[2][0][1] > 1:
                zrxn_objs = automol.reac.from_geometries(
                    rct_geos, prd_geos, struc_typ="zmat", stereo=False)
                print('RXN UNIDENTIFIABLE FROM GEOS BUT ASSUMING IT IS R+O2: BUT GOING AHEAD')

    else:
        print('    Reaction ID from geometries from input InChIs')
        rxn_ichs = rinfo.value(rxn_info, 'inchi')
        rct_ichs, prd_ichs = rxn_ichs[0], rxn_ichs[1]

        zrxn_objs = automol.reac.from_chis(
            rct_ichs, prd_ichs, struc_typ='zmat', stereo=True, enant=False)

        # AVC comment: Why is this necessary? Is this a bug I need to fix?
        if not zrxn_objs and len(rxn_info[2][0]) == 2:
            if rxn_info[2][0][0] > 1 and rxn_info[2][0][1] > 1:
                zrxn_objs = automol.reac.from_chis(
                    rct_ichs, prd_ichs, struc_typ='zmat', stereo=False)
                print('RXN UNIDENTIFIABLE FROM ICHS BUT ASSUMING IT IS R+O2 AND GOING AHEAD')

    # Loop over the found reaction objects
    if zrxn_objs is not None:
        zrxns, zmas = (), ()
        for zrxn_obj in zrxn_objs:
            zrxn = zrxn_obj
            zma = automol.reac.ts_structure(zrxn_obj)
            zrxns += (zrxn,)
            zmas += (zma,)
        if len(zrxns) == 0 and len(zmas) == 0:
            zrxns, zmas = None, None
    else:
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
    _fake_class = automol.ReactionInfo(class_typ)
    if automol.ReactionInfo.requires_spin_designation(_fake_class):
        ts_mul = rinfo.value(rxn_info, 'tsmult')
        high_mul = rinfo.ts_mult(rxn_info, rxn_mul='high')
        _spin = automol.ReactionSpin.HIGH if ts_mul == high_mul else automol.ReactionSpin.LOW
    else:
        _spin = automol.ReactionSpin.NONE

    # Determine if it iss intersystem crossing
    # rxn_muls = rinfo.value(rxn_info, 'mult')
    # is_isc = automol.reac.intersystem_crossing(rxn_muls)

    return automol.ReactionInfo(
        class_=class_typ, spin=_spin, is_rad_rad=rinfo.radrad(rxn_info), is_isc=False
    )


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

    rct_str, prd_str = '+'.join(reacs), '+'.join(prods)
    print(f'    Running reaction as: {rct_str} = {prd_str}')

    return reacs, prods


# Functions for the exothermicity check
def reagent_geometries(rxn_info, thy_info, save_prefix, canonical=True):
    """ Identify the reaction and build the object

        :param rxn_info: reaction info object
        :type rxn_info: mechanalyzer.inf.rxn object
        :rtype: (tuple(automol.Reaction object), tuple(automol.zmat object))
    """
    def _canonical(rgt_info):
        canon_rgt_info = ()
        for spc_info in rgt_info:
            canon_ich = automol.chi.canonical_enantiomer(spc_info[0])
            # canon_ich = spc_info[0]
            canon_rgt_info += ((canon_ich, spc_info[1], spc_info[2]),)
        # canon_ichs = rgt_info[0]
        return canon_rgt_info

    # Check the save filesystem for the reactant and product geometries
    rct_info = rinfo.rgt_info(rxn_info, 'reacs')
    prd_info = rinfo.rgt_info(rxn_info, 'prods')
    if canonical:
        canon_enant_rct_info = _canonical(rct_info)
        canon_enant_prd_info = _canonical(prd_info)
        _rcts_cnf_fs = filesys.rcts_cnf_fs(
            canon_enant_rct_info, thy_info, None, save_prefix)
        _prds_cnf_fs = filesys.rcts_cnf_fs(
            canon_enant_prd_info, thy_info, None, save_prefix)
    else:
        _rcts_cnf_fs = filesys.rcts_cnf_fs(
            rct_info, thy_info, None, save_prefix)
        _prds_cnf_fs = filesys.rcts_cnf_fs(
            prd_info, thy_info, None, save_prefix)

    # If min cnfs found for all rcts and prds, read the geometries
    rct_geos, prd_geos = (), ()
    rct_paths, prd_paths = (), ()
    if (
        _rcts_cnf_fs.count(None) == 0 and _prds_cnf_fs.count(None) == 0
    ):
        for idx, (_, cnf_save_fs, min_locs, _) in enumerate(_rcts_cnf_fs):
            geo = cnf_save_fs[-1].file.geometry.read(min_locs)
            if canonical:
                if rct_info[idx] != canon_enant_rct_info[idx]:
                    geo = automol.geom.reflect_coordinates(geo)
            path = cnf_save_fs[-1].file.geometry.path(min_locs)
            rct_geos += (geo,)
            rct_paths += (path,)
        for idx, (_, cnf_save_fs, min_locs, _) in enumerate(_prds_cnf_fs):
            geo = cnf_save_fs[-1].file.geometry.read(min_locs)
            if canonical:
                if prd_info[idx] != canon_enant_prd_info[idx]:
                    geo = automol.geom.reflect_coordinates(geo)
            path = cnf_save_fs[-1].file.geometry.path(min_locs)
            prd_geos += (geo,)
            prd_paths += (path,)
    else:
        _rcts_cnf_fs = filesys.rcts_cnf_fs(
            rct_info, thy_info, None, save_prefix)
        _prds_cnf_fs = filesys.rcts_cnf_fs(
            prd_info, thy_info, None, save_prefix)
        if (
            _rcts_cnf_fs.count(None) == 0 and _prds_cnf_fs.count(None) == 0
        ):
            for idx, (_, cnf_save_fs, min_locs, _) in enumerate(_rcts_cnf_fs):
                geo = cnf_save_fs[-1].file.geometry.read(min_locs)
                path = cnf_save_fs[-1].file.geometry.path(min_locs)
                rct_geos += (geo,)
                rct_paths += (path,)
            for idx, (_, cnf_save_fs, min_locs, _) in enumerate(_prds_cnf_fs):
                geo = cnf_save_fs[-1].file.geometry.read(min_locs)
                path = cnf_save_fs[-1].file.geometry.path(min_locs)
                prd_geos += (geo,)
                prd_paths += (path,)

    return rct_geos, prd_geos, rct_paths, prd_paths


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

    print(f'    Reaction energy is {rxn_ene*phycon.EH2KCAL:.2f} '
          f'at {method1[1]}//{method2[1]} level')

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


def reagent_energies(
        rgt, rxn_info, sp_thy_info, geo_thy_info, save_prefix):
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
