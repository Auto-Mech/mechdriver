""" Handle symmetry factor stuff
"""

import automol
from routines.es._routines import conformer


def symmetry_factor(sym_model, spc_dct_i, spc_info, dist_names,
                    saddle, frm_bnd_key, brk_bnd_key, tors_names,
                    tors_cnf_save_fs, tors_min_cnf_locs,
                    sym_cnf_save_fs, sym_min_cnf_locs):
    """ Get the overall factor for a species
    """

    form_coords = []
    if 'sym_factor' in spc_dct_i:
        sym_factor = spc_dct_i['sym_factor']
        print('sym_factor from spc_dct_i:', sym_factor)
    else:

        # Obtain the external symmetry number
        ext_sym = automol.geom.external_symmetry_factor(geo)

        # Obtain the internal symmetry number using some routine
        if sym_model == 'sampling':

            if tors_names:
                int_sym = int_sym_num_from_sampling(
                geo, ene, cnf_save_fs, saddle,
                frm_bnd_key, brk_bnd_key, tors_names)
            else:
                int_sym = 1.0

            if not sym_min_cnf_locs:
                # Fix the return statement here
                print('ERROR: Reference geometry is missing for symmetry',
                      'for species {}'.format(spc_info[0]))
                return '', 0.
            sym_geo = sym_cnf_save_fs[-1].file.geometry.read(sym_min_cnf_locs)
            sym_ene = sym_cnf_save_fs[-1].file.energy.read(sym_min_cnf_locs)
            sym_factor = conformer.symmetry_factor(
                sym_geo, sym_ene, sym_cnf_save_fs, saddle,
                frm_bnd_key, brk_bnd_key, tors_names)
            print('sym_factor from conformer sampling:', sym_factor)

        elif sym_model == 'none':
            # print('Warning: no symmetry model requested,',
            #       'setting symmetry factor to 1.0')
            int_sym = 1.0

        # Obtain overall number
        sym_factor = ext_sym * int_sym

    return sym_factor


def symmetry_factor(
        geo, ene, cnf_save_fs, saddle=False, frm_bnd_key=(), brk_bnd_key=(),
        tors_names=()):
    """ obtain overall symmetry factor for a geometry as a product
        of the external symmetry factor and the internal symmetry number
        Note: ignoring for saddle pts the possibility that two configurations
        differ only in their torsional values.
        As a result, the symmetry factor is a lower bound of the true value
    """
    # Note: ignoring for saddle points the possibility that two configurations
    # differ only in their torsional values.
    # As a result, the symmetry factor is a lower bound of the true value
    ext_sym = automol.geom.external_symmetry_factor(geo)
    if not saddle:
        tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    if tors_names:
        int_sym = int_sym_num_from_sampling(
            geo, ene, cnf_save_fs, saddle,
            frm_bnd_key, brk_bnd_key, tors_names)
    else:
        int_sym = 1
    sym_fac = ext_sym * int_sym
    return sym_fac


def int_sym_num_from_sampling(
        geo, ene, cnf_save_fs, saddle=False, frm_bnd_key=(),
        brk_bnd_key=(), tors_names=()):
    """ Determine the symmetry number for a given conformer geometry.
    (1) Explore the saved conformers to find the list of similar conformers -
        i.e. those with a coulomb matrix and energy that are equivalent
        to those for the reference geometry.
    (2) Expand each of those similar conformers by applying
        rotational permutations to each of the terminal groups.
    (3) Count how many distinct distance matrices there are in
        the fully expanded conformer list.
    """

    # Note: ignoring for saddle points the possibility that two configurations
    # differ only in their torsional values.
    # As a result, the symmetry factor is a lower bound of the true value
    if automol.geom.is_atom(geo):
        int_sym_num = 1.
    else:
        if not saddle:
            tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
        if tors_names is None:
            int_sym_num = 1.
        else:
            ethrsh = 1.e-5
            locs_lst = cnf_save_fs[-1].existing()
            int_sym_num = 1.
            if locs_lst:
                enes = [cnf_save_fs[-1].file.energy.read(locs)
                        for locs in locs_lst]
                geos = [cnf_save_fs[-1].file.geometry.read(locs)
                        for locs in locs_lst]
                geo_sim = []
                geo_sim2 = []
                ene_sim = []
                for geoi, enei in zip(geos, enes):
                    if enei - enes[0] < ethrsh:
                        geo_lst = [geoi]
                        ene_lst = [enei]
                        unique = is_unique_coulomb_energy(
                            geo, ene, geo_lst, ene_lst)
                        if not unique:
                            geo_sim.append(geoi)
                            ene_sim.append(enei)

                int_sym_num = 0
                for geo_sim_i in geo_sim:
                    new_geos = automol.geom.rot_permutated_geoms(
                        geo_sim_i, saddle, frm_bnd_key, brk_bnd_key)
                    for new_geo in new_geos:
                        new_geom = True
                        for geo_sim_j in geo_sim2:
                            if automol.geom.almost_equal_dist_mat(
                                    new_geo, geo_sim_j, thresh=3e-1):
                                if saddle:
                                    new_geom = False
                                    break
                                if are_torsions_same(new_geo, geo_sim_j):
                                    new_geom = False
                                    break
                        if new_geom:
                            geo_sim2.append(new_geo)
                            int_sym_num += 1
    return int_sym_num


def tors_reduced_sym_factor(sym_factor, tors_sym_nums):
    """ Decrease the overall molecular symmetry factor by the
        torsional mode symmetry numbers
    """
    for sym_num in tors_sym_nums:
        sym_factor /= sym_num

    return sym_factor
