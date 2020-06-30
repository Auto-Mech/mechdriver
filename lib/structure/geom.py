"""
  Functions used for handling and comparing geometries
"""

import numpy
import automol


# Handle idx lists for zmas with dummys
def build_remdummy_shift_lst(zma):
    """ Assess the zma for dummys and build a list to shift values
        derived from zma
    """

    atom_symbols = automol.zmatrix.symbols(zma)
    dummy_idx = []
    for atm_idx, atm in enumerate(atom_symbols):
        if atm == 'X':
            dummy_idx.append(atm_idx)
    remdummy = numpy.zeros(len(zma[0]))
    for dummy in dummy_idx:
        for idx, _ in enumerate(remdummy):
            if dummy < idx:
                remdummy[idx] += 1

    return remdummy


def shift_vals_from_dummy(vals, zma):
    """ Shift a set of values using remdummy
        Shift requires indices be 1-indexed
    """

    dummy_idxs = automol.zmatrix.atom_indices(zma, 'X', match=True)

    shift_vals = []
    for val in vals:
        shift = 0
        for dummy in dummy_idxs:
            if val >= dummy:
                shift += 1
        shift_vals.append(val+shift)

    return shift_vals


# Various checks and assessments for geometries
def is_atom_closest_to_bond_atom(zma, idx_rad, bond_dist):
    """ Check to see whether the radical atom is still closest to the bond
        formation site.
    """
    geo = automol.zmatrix.geometry(zma)
    atom_closest = True
    for idx, _ in enumerate(geo):
        if idx < idx_rad:
            distance = automol.geom.distance(geo, idx, idx_rad)
            if distance < bond_dist-0.01:
                atom_closest = False
                # print('idx test:', idx, distance, bond_dist)
    return atom_closest


def calc_rxn_angle(ts_zma, frm_bnd_keys, brk_bnd_keys, rxn_class):
    """ Calculate the angle over a forming-and-breaking bond
    """

    angle = None
    if 'abstraction' in rxn_class or 'addition' in rxn_class:
        if frm_bnd_keys and brk_bnd_keys:

            ang_atms = [0, 0, 0]
            cent_atm = list(set(brk_bnd_keys) & set(frm_bnd_keys))
            if cent_atm:
                ang_atms[1] = cent_atm[0]
                for idx in brk_bnd_keys:
                    if idx != ang_atms[1]:
                        ang_atms[0] = idx
                for idx in frm_bnd_keys:
                    if idx != ang_atms[1]:
                        ang_atms[2] = idx

                geom = automol.zmatrix.geometry(ts_zma)
                angle = automol.geom.central_angle(
                    geom, *ang_atms)

    return angle


def is_unique_coulomb_energy(geo, ene, geo_list, ene_list):
    """ compare given geo with list of geos all to see if any have the same
    coulomb spectrum and energy
    """
    unique = True
    for idx, geoi in enumerate(geo_list):
        enei = ene_list[idx]
        etol = 2.e-5
        if abs(ene-enei) < etol:
            if automol.geom.almost_equal_coulomb_spectrum(
                    geo, geoi, rtol=1e-2):
                unique = False
    return unique


def is_unique_dist_mat_energy(geo, ene, geo_list, ene_list):
    """ compare given geo with list of geos all to see if any have the same
    distance matrix and energy
    """
    unique = True
    for idx, geoi in enumerate(geo_list):
        enei = ene_list[idx]
        etol = 2.e-5
        if abs(ene-enei) < etol:
            if automol.geom.almost_equal_distance_matrix(
                    geo, geoi, thresh=1e-1):
                unique = False
    return unique


def is_unique_stereo_dist_mat_energy(geo, ene, geo_list, ene_list):
    """ compare given geo with list of geos all to see if any have the same
    distance matrix and energy and stereo specific inchi
    """
    unique = True
    ich = automol.convert.geom.inchi(geo)
    for idx, geoi in enumerate(geo_list):
        enei = ene_list[idx]
        etol = 2.e-5
        ichi = automol.convert.geom.inchi(geoi)
        # check energy
        if abs(ene-enei) < etol:
            # check distance matrix
            if automol.geom.almost_equal_distance_matrix(
                    geo, geoi, thresh=1e-1):
                # check stereo by generates stero label
                ichi = automol.convert.geom.inchi(geoi)
                if ich == ichi:
                    unique = False
    return unique


def are_torsions_same(geo, geoi):
    """ compare all torsional angle values
    """
    dtol = 0.09
    same_dihed = True
    zma = automol.geom.zmatrix(geo)
    tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    zmai = automol.geom.zmatrix(geoi)
    tors_namesi = automol.geom.zmatrix_torsion_coordinate_names(geoi)
    for idx, tors_name in enumerate(tors_names):
        val = automol.zmatrix.values(zma)[tors_name]
        vali = automol.zmatrix.values(zmai)[tors_namesi[idx]]
        valip = vali+2.*numpy.pi
        valim = vali-2.*numpy.pi
        vchk1 = abs(val - vali)
        vchk2 = abs(val - valip)
        vchk3 = abs(val - valim)
        if vchk1 > dtol and vchk2 > dtol and vchk3 > dtol:
            same_dihed = False
    return same_dihed


def is_unique_tors_dist_mat_energy(geo, ene, geo_list, ene_list, saddle):
    """ compare given geo with list of geos all to see if any have the same
    coulomb spectrum and energy and stereo specific inchi
    """
    unique = True
    etol = 2.e-5
    for idx, geoi in enumerate(geo_list):
        enei = ene_list[idx]
        # check energy
        if abs(ene-enei) < etol:
            # check distance matrix
            if automol.geom.almost_equal_distance_matrix(
                    geo, geoi, thresh=3e-1):
                # check dihedrals
                # for now only do this for minima
                # but this can create problems for TSs as well
                # - e.g., CH2OH = CH2O + H
                if saddle:
                    unique = False
                elif are_torsions_same(geo, geoi):
                    unique = False
                # if are_torsions_same(geo, geoi):
                    # unique = False
    return unique
