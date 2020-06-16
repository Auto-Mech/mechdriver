"""
 Handle stuff for fake species
"""

import automol


def set_fake_freqs(geom_i, geom_j):
    """ Set fake frequencies
    """

    freqs = (30, 50, 70, 100, 200)
    ntrans = 5
    is_atom_i = automol.geom.is_atom(geom_i)
    is_linear_i = automol.geom.is_linear(geom_i)
    is_atom_j = automol.geom.is_atom(geom_j)
    is_linear_j = automol.geom.is_linear(geom_i)
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


def combine_geos_in_fake_well(geom_i, geom_j):
    """ put two geometries together in a fake well
    """

    max_z_i = max(atom[1][2] for atom in geom_i)
    min_z_j = min(atom[1][2] for atom in geom_j)
    geom = geom_i
    geom_j = automol.geom.translated(
        geom_j, [0., 0., max_z_i + 6. - min_z_j])
    geom += geom_j

    return geom
