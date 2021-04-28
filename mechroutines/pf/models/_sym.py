""" Handle symmetry factor stuff
"""

import automol
from autofile import fs
from mechlib.amech_io import printer as ioprinter


def symmetry_factor(pf_filesystems, spc_mod_dct_i, spc_dct_i, rotors,
                    grxn=None, zma=None):
    """ Calculate the symmetry factor for a species
        Note: ignoring for saddle pts the possibility that two configurations
        differ only in their torsional values.
        As a result, the symmetry factor is a lower bound of the true value
    """

    sym_factor = spc_dct_i.get('sym_factor')
    if sym_factor is not None:
        ioprinter.info_message(' - Reading symmetry number input by user:', sym_factor)
    else:

        zrxn = spc_dct_i.get('zrxn', None)
        if zrxn is not None:
            grxn = automol.reac.relabel_for_geometry(zrxn)
        else:
            grxn = None

        # if automol.geom.is_atom(geo):
        sym_model = spc_mod_dct_i['symm']['mod']

        # Obtain geometry, energy, and symmetry filesystem
        [cnf_fs, cnf_path, min_cnf_locs, _, _] = pf_filesystems['symm']
        geo = cnf_fs[-1].file.geometry.read(min_cnf_locs)

        # Obtain the external symssetry number
        ext_sym = automol.geom.external_symmetry_factor(geo)

        # Obtain the internal symmetry number using some routine
        if sym_model == 'sampling':

            # Set up the symmetry filesystem
            sym_fs = fs.symmetry(cnf_path)
            sym_geos = [geo]
            sym_geos += [sym_fs[-1].file.geometry.read(locs)
                         for locs in sym_fs[-1].existing()]

            # Obtain the internal
            if rotors:
                ioprinter.info_message(
                    ' - Determining internal sym number ',
                    'using sampling routine.')
                int_sym = int_sym_num_from_sampling(sym_geos, rotors, grxn=grxn, zma=zma)
            else:
                ioprinter.info_message(' - No torsions, internal sym is 1.0')
                int_sym = 1.0

        else:
            ioprinter.info_message(
                'No symmetry model requested, ',
                'setting internal sym factor to 1.0')
            int_sym = 1.0

        # Obtain overall number
        sym_factor = ext_sym * int_sym
        print('sym_factor test:', ext_sym, int_sym, sym_factor)

        # Reduce sym factor using rotor symmetries
        sym_factor = tors_reduced_sym_factor(sym_factor, rotors)

        # ioprinter.info_message('sym_factor test:', sym_factor)

    return sym_factor


def _modify_idxs(idxs_lst, removed_atms, dummy_atms):
    mod_idxs_lst = []
    no_dummy_idxs_lst = []
    for idxs in idxs_lst:
        mod_idxs = []
        for idx in idxs:
            mod_idx = idx
            for atm in dummy_atms:
                if atm < idx:
                    mod_idx -= 1
            mod_idxs.append(mod_idx)
        no_dummy_idxs_lst.append(mod_idxs)

    for idxs in no_dummy_idxs_lst:
        in_lst = True
        for idx in idxs:
            if idx in removed_atms:
                in_lst = False
        if in_lst:
            mod_idxs = []
            for idx in idxs:
                mod_idx = idx
                for atm in removed_atms:
                    if atm < idx:
                        mod_idx -= 1
                mod_idxs.append(mod_idx)
            mod_idxs_lst.append(mod_idxs)
    return mod_idxs_lst


def int_sym_num_from_sampling(sym_geos, rotors, grxn=None, zma=None):
    """ Determine the symmetry number for a given conformer geometry.
    (1) Explore the saved conformers to find the list of similar conformers -
        i.e. those with a coulomb matrix and energy that are equivalent
        to those for the reference geometry.
    (2) Expand each of those similar conformers by applying
        rotational permutations to each of the terminal groups.
    (3) Count how many distinct distance matrices there are in
        the fully expanded conformer list.
    """

    if grxn is not None:
        frm_bnd_keys = automol.reac.forming_bond_keys(grxn)
        brk_bnd_keys = automol.reac.breaking_bond_keys(grxn)
        tors_names = automol.rotor.names(rotors, flat=True)
        tors_idxs = [automol.zmat.coord_idxs(zma, name) for name in tors_names]
    else:
        frm_bnd_keys, brk_bnd_keys = frozenset({}), frozenset({})

    int_sym_num = 0
    # modify geometries to remove H's from rotatable XHn end group
    # this will be accounted for separately as multiplicative factor
    mod_sym_geos = []
    print('sym_geos test:', sym_geos)
    print('keys:', frm_bnd_keys, brk_bnd_keys)
    ts_bnds = ()
    if grxn is not None:
        ts_bnds = (frm_bnd_keys, brk_bnd_keys)
    for geo_sym_i in sym_geos:
        ret = automol.geom.end_group_symmetry_factor(
            geo_sym_i, frm_bnd_keys, brk_bnd_keys)
        mod_geo_sym_i, end_group_factor, removed_atms = ret
        if grxn is not None:
            mod_tors_idxs = _modify_idxs(tors_idxs, removed_atms, automol.zmat.dummy_keys(zma))
        # ioprinter.info_message('end_group_factor test:', end_group_factor)

        new_geom = True
        for mod_geo_sym_j in mod_sym_geos:
            if automol.geom.almost_equal_dist_matrix(
                    mod_geo_sym_i, mod_geo_sym_j, thresh=3e-1):
                if grxn is None:
                    tors_same = automol.geom.are_torsions_same(
                        mod_geo_sym_i, mod_geo_sym_j, ts_bnds=())
                else:
                    tors_same = automol.geom.are_torsions_same2(
                        mod_geo_sym_i, mod_geo_sym_j, mod_tors_idxs)
                if tors_same:
                    new_geom = False
                    break
        if new_geom:
            mod_sym_geos.append(mod_geo_sym_i)
            int_sym_num += 1
            print('sym_geo test:', mod_geo_sym_i, int_sym_num)

    int_sym_num *= end_group_factor
    print('end_group_factor:', end_group_factor)
    print('final int_sym_num:', int_sym_num)

    return int_sym_num


def tors_reduced_sym_factor(sym_factor, rotors):
    """ Decrease the overall molecular symmetry factor by the
        torsional mode symmetry numbers
    """
    tors_symms = automol.rotor.symmetries(rotors, flat=True)
    for symm in tors_symms:
        sym_factor /= symm

    return sym_factor
