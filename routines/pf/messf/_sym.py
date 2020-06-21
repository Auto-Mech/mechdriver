""" Handle symmetry factor stuff
"""

import automol
from autofile import fs
from lib import structure


def symmetry_factor(pf_filesystems, pf_models, spc_dct_i,
                    frm_bnd_keys=(), brk_bnd_keys=(), rotor_names=()):
    """ Calculate the symmetry factor for a species
        Note: ignoring for saddle pts the possibility that two configurations
        differ only in their torsional values.
        As a result, the symmetry factor is a lower bound of the true value
    """

    if 'sym_factor' in spc_dct_i:
        sym_factor = spc_dct_i['sym_factor']
        print(' - Reading symmetry number input by user:', sym_factor)
    else:

        # if automol.geom.is_atom(geo):
        sym_model = pf_models['sym']

        # Obtain geometry, energy, and symmetry filesystem
        [cnf_fs, cnf_path, min_cnf_locs, _, _] = pf_filesystems['sym']
        geo = cnf_fs[-1].file.geometry.read(min_cnf_locs)

        # Obtain the external symssetry number
        ext_sym = automol.geom.external_symmetry_factor(geo)

        # Obtain the internal symmetry number using some routine
        if sym_model == 'sampling':

            # Set up the symmetry filesystem
            sym_fs = fs.manager(cnf_path, 'SYMMETRY')
            sym_geos = [sym_fs[-1].file.geometry.read(locs)
                        for locs in sym_fs[-1].existing()]

            # Obtain the internal
            if rotor_names:
                print(' - Determining internal sym number ',
                      'using sampling routine.')
                int_sym = int_sym_num_from_sampling(
                    sym_geos,
                    frm_bnd_keys=frm_bnd_keys,
                    brk_bnd_keys=brk_bnd_keys)
                print('int', int_sym)
            else:
                print(' - No torsions, internal sym is 1.0')
                int_sym = 1.0

        else:
            print('No symmetry model requested, ',
                  'setting internal sym factor to 1.0')
            int_sym = 1.0

        # Obtain overall number
        sym_factor = ext_sym * int_sym

    return sym_factor


def int_sym_num_from_sampling(sym_geos, frm_bnd_keys=(), brk_bnd_keys=()):
    """ Determine the symmetry number for a given conformer geometry.
    (1) Explore the saved conformers to find the list of similar conformers -
        i.e. those with a coulomb matrix and energy that are equivalent
        to those for the reference geometry.
    (2) Expand each of those similar conformers by applying
        rotational permutations to each of the terminal groups.
    (3) Count how many distinct distance matrices there are in
        the fully expanded conformer list.
    """

    # Set saddle
    saddle = bool(frm_bnd_keys or brk_bnd_keys)

    int_sym_num = 0
    sym_geos2 = []
    for geo_sym_i in sym_geos:
        new_geos = automol.geom.rot_permutated_geoms(
            geo_sym_i, frm_bnd_keys, brk_bnd_keys)
        for new_geo in new_geos:
            new_geom = True
            for geo_sym_j in sym_geos2:
                if automol.geom.almost_equal_dist_mat(
                        new_geo, geo_sym_j, thresh=3e-1):
                    if saddle:
                        new_geom = False
                        break
                    if structure.geom.are_torsions_same(new_geo, geo_sym_j):
                        new_geom = False
                        break
            if new_geom:
                sym_geos2.append(new_geo)
            int_sym_num += 1

    return int_sym_num


def tors_reduced_sym_factor(sym_factor, tors_sym_nums):
    """ Decrease the overall molecular symmetry factor by the
        torsional mode symmetry numbers
    """
    for sym_num in tors_sym_nums:
        sym_factor /= sym_num

    return sym_factor
