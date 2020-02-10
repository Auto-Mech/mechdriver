"""
 Handle stuff for fake species
"""

def set_fake_freqs(harm_min_cnf_locs_i, harm_min_cnf_locs_j,
                   harm_cnf_save_fs_i, harm_cnf_save_fs_j):
    """ Set fake frequencies
    """
    if harm_min_cnf_locs_i is not None:
        harm_geo_i = harm_cnf_save_fs_i.leaf.file.geometry.read(
            harm_min_cnf_locs_i)
        if harm_min_cnf_locs_j is not None:
            harm_geo_j = harm_cnf_save_fs_j.leaf.file.geometry.read(
                harm_min_cnf_locs_j)
    freqs = [30, 50, 70, 100, 200]
    ntrans = 5
    is_atom_i = automol.geom.is_atom(harm_geo_i)
    is_linear_i = automol.geom.is_linear(harm_geo_i)
    is_atom_j = automol.geom.is_atom(harm_geo_j)
    is_linear_j = automol.geom.is_linear(harm_geo_i)
    if is_atom_i:
        ntrans = ntrans - 3
    if is_atom_j:
        ntrans = ntrans - 3
    if is_linear_i:
        ntrans = ntrans - 2
    if is_linear_j:
        ntrans = ntrans - 2
    if is_atom_i and is_atom_j:
        ntrans = 0
    freqs = freqs[0:ntrans]

    return freqs


def combine_geos_in_fake_well(harm_min_cnf_locs_i, harm_min_cnf_locs_j,
                              harm_cnf_save_fs_i, harm_cnf_save_fs_j):
    """ put two geometries together in a fake well
    """
    if harm_min_cnf_locs_i is not None:
        harm_geo_i = harm_cnf_save_fs_i.leaf.file.geometry.read(
            harm_min_cnf_locs_i)
        if harm_min_cnf_locs_j is not None:
            harm_geo_j = harm_cnf_save_fs_j.leaf.file.geometry.read(
                harm_min_cnf_locs_j)
    max_z_i = max(atom[1][2] for atom in harm_geo_i)
    min_z_j = min(atom[1][2] for atom in harm_geo_j)
    harm_geo = harm_geo_i
    harm_geo_j = automol.geom.translated(
        harm_geo_j, [0., 0., max_z_i + 6. - min_z_j])
    harm_geo += harm_geo_j

    return harm_geo
