""" Generate the information necessary to product the vrctst input files
"""

import os
import subprocess
import automol
import varecof_io


def input_prep(ts_zma, rct_zmas, dist_name):
    """ prepare all the input files for a vrc-tst calculation
    """

    # Set info on the form indices
    bnd_frm_idxs = automol.zmatrix.bond_idxs(ts_zma, dist_name)
    min_idx, max_idx = min(bnd_frm_idxs), max(bnd_frm_idxs)

    # Build geometries needed for the varecof run
    total_geom, frag_geoms, frag_geoms_wdummy = fragment_geometries(
        ts_zma, rct_zmas, min_idx, max_idx)

    # Set information for the pivot points needed in divsur.inp
    frames, npivots = build_pivot_frames(
        min_idx, max_idx, total_geom, frag_geoms, frag_geoms_wdummy)
    pivot_angles = calc_pivot_angles(frag_geoms, frag_geoms_wdummy, frames)
    pivot_xyzs = calc_pivot_xyzs(min_idx, max_idx, total_geom, frag_geoms)

    # Write the long- and short-range divsur input files
    r1dists_lr = [8., 6., 5., 4.5, 4.]
    lr_divsur_inp_str = varecof_io.writer.input_file.divsur(
        r1dists_lr, 1, 1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])

    # Write the short-range divsur files
    r1dists_sr = [4., 3.8, 3.6, 3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2]
    r2dists_sr = [4., 3.8, 3.6, 3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2]
    d1dists = [0.01, 0.5, 1.]
    d2dists = [0.01, 0.5, 1.]
    t1angs = [pivot_angles[0]] if pivot_angles[0] is not None else []
    t2angs = [pivot_angles[1]] if pivot_angles[0] is not None else []
    if automol.geom.is_atom(frag_geoms[0]):
        d1dists = []
        t1angs = []
    if automol.geom.is_atom(frag_geoms[1]):
        d2dists = []
        t2angs = []
    if automol.geom.is_linear(frag_geoms[0]):
        d1dists = [0.]
        t1angs = []
    if automol.geom.is_linear(frag_geoms[1]):
        d2dists = [0.]
        t2angs = []
    conditions = {'delta_r': 0}
    sr_divsur_inp_str = varecof_io.writer.input_file.divsur(
        r1dists_sr, npivots[0], npivots[1], pivot_xyzs[0], pivot_xyzs[1],
        frame1=frames[0],
        frame2=frames[1],
        d1dists=d1dists,
        d2dists=d2dists,
        t1angs=t1angs,
        t2angs=t2angs,
        r2dists=r2dists_sr,
        **conditions)

    # Write the structure input files
    struct_inp_str = varecof_io.writer.input_file.structure(
        frag_geoms_wdummy[0], frag_geoms_wdummy[1])

    # Write the tst.inp file
    nsamp_max = 2000
    nsamp_min = 50
    flux_err = 10
    pes_size = 2
    faces, faces_symm = assess_face_symmetries()
    tst_inp_str = varecof_io.writer.input_file.tst(
        nsamp_max, nsamp_min, flux_err, pes_size, faces, faces_symm)

    # Write the potential energy surface input string
    exe_path = '/blues/gpfs/home/sjklipp/bin/molpro'
    base_name = 'mol'
    els_inp_str = varecof_io.writer.input_file.elec_struct(
        exe_path, base_name)

    # Write the mc_flux.inp input string
    mc_flux_inp_str = varecof_io.writer.input_file.mc_flux()

    # Write the convert.inp input string
    conv_inp_str = varecof_io.writer.input_file.convert()

    input_strs = [
        struct_inp_str, lr_divsur_inp_str, sr_divsur_inp_str, tst_inp_str,
        els_inp_str, mc_flux_inp_str, conv_inp_str]

    return input_strs


def build_correction_potential(mep_distances,
                               potentials,
                               bnd_frm_idxs,
                               fortran_compiler,
                               fortran_make_path,
                               bnd_frm_syms=(),
                               species_name='',
                               pot_labels=(),
                               pot_file_names=()):
    """  use the MEP potentials to compile the correction potential .so file
    """

    # Write strings corresponding to each of the correction potential files
    species_corr_str = varecof_io.writer.corr_potentials.species(
        mep_distances,
        potentials,
        bnd_frm_idxs,
        bnd_frm_syms=bnd_frm_syms,
        species_name=species_name,
        pot_labels=pot_labels)
    dummy_corr_str = varecof_io.writer.corr_potentials.dummy()
    pot_aux_str = varecof_io.writer.corr_potentials.auxiliary()
    makefile_str = varecof_io.writer.corr_potentials.makefile(
        fortran_compiler,
        pot_file_names=pot_file_names)

    # Write all of the files needed to build the correction potential
    with open('mol_corr.f', 'w') as mol_corr_file:
        mol_corr_file.write(species_corr_str)
    with open('dummy_corr.f', 'w') as dummy_corr_file:
        dummy_corr_file.write(dummy_corr_str)
    with open('pot_aux.f', 'w') as pot_aux_file:
        pot_aux_file.write(pot_aux_str)
    with open('makefile', 'w') as makefile_file:
        makefile_file.write(makefile_str)

    # Compile the correction potential
    varecof_io.writer.corr_potentials.compile_corr_pot(
        fortran_make_path)


def fragment_geometries(ts_zma, rct_zmas, min_idx, max_idx):
    """ Generate the fragment geometries from the ts Z-matrix and the
        indices involved in the forming bond
    """

    # Get the MEP geometries
    mep_total_geo = automol.zmatrix.geometry(ts_zma)
    mep_fgeos = [mep_total_geo[:max_idx], mep_total_geo[max_idx:]]

    # Get the geometries of the isolated fragments (the reactants)
    iso_fgeos = [automol.zmatrix.geometry(zma) for zma in rct_zmas]

    # Get the geometries for the structure.inp file
    iso_fgeos_wdummy = []
    mol_data = zip(mep_fgeos, iso_fgeos, (max_idx, min_idx))
    for i, (mep_fgeo, iso_fgeo, idx) in enumerate(mol_data):

        if not automol.geom.is_atom(mep_fgeo):

            # Build MEPFragGeom+X coordinates using MEP geometry
            x_coord = mep_total_geo[idx][1]
            dummy_row = ('X', x_coord)
            if i == 0:
                mep_geo_wdummy = mep_fgeo + (dummy_row,)
                x_idx = len(mep_geo_wdummy) - 1
                a1_idx = 0
            else:
                mep_geo_wdummy = (dummy_row,) + mep_fgeo
                x_idx = 0
                a1_idx = 1

            # Calculate coords to define X position in IsoFragGeom structure
            xyz1 = iso_fgeo[a1_idx][1]
            xyz2 = iso_fgeo[a1_idx+1][1]
            xdistance = automol.geom.distance(
                mep_geo_wdummy, x_idx, a1_idx)
            xangle = automol.geom.central_angle(
                mep_geo_wdummy, x_idx, a1_idx, a1_idx+1)
            if len(mep_fgeo) > 2:
                xyz3 = iso_fgeo[a1_idx+2][1]
                xdihedral = automol.geom.dihedral_angle(
                    mep_geo_wdummy, x_idx, a1_idx, a1_idx+1, a1_idx+2)
            else:
                xyz3 = 0.0
                xdihedral = 0.0

            # Calculate the X Position for the IsoFrag structure
            xyzp = automol.geom.find_xyzp_using_internals(
                xyz1, xyz2, xyz3, xdistance, xangle, xdihedral)

            # Generate the IsoFragGeom+X coordinates for the structure.inp file
            if i == 0:
                iso_geo_wdummy = iso_fgeo + (('X', xyzp),)
            else:
                iso_geo_wdummy = (('X', xyzp),) + iso_fgeo

            # Append to final geoms
            iso_fgeos_wdummy.append(iso_geo_wdummy)

        else:
            # If atom, set IsoFragGeom+X coords equal to mep_geo
            iso_fgeos_wdummy.append(mep_fgeo)

    return mep_total_geo, iso_fgeos, iso_fgeos_wdummy


def assess_face_symmetries():
    """ check the symmetry of the faces for each fragment
    """

    devnull = open(os.devnull, 'w')
    subprocess.check_call(['./convert_struct', 'divsur.inp'],
                          stdout=devnull, stderr=devnull)

    # Read fragment geoms from divsur.out with coordinates in the divsur frame
    with open('divsur.out', 'r') as divsur_file:
        divsur_string = divsur_file.read()
    fgeo1, fgeo2 = varecof_io.reader.divsur.frag_geoms_divsur_frame(
        divsur_string)
    fgeo1 = automol.geom.from_string(fgeo1)
    fgeo2 = automol.geom.from_string(fgeo2)

    # Reflect the dummy atom (pivot position) about the xy plane
    f1_dummy_idx = len(fgeo1) - 1
    f2_dummy_idx = 0
    fgeo1_reflect = automol.geom.reflect_coordinates(
        fgeo1, [f1_dummy_idx], ['x', 'y'])
    fgeo2_reflect = automol.geom.reflect_coordinates(
        fgeo2, [f2_dummy_idx], ['x', 'y'])

    # Compate coloumb spectrum for each geom to its reflected version
    fgeo1_symm = automol.geom.almost_equal_coulomb_spectrum(
        fgeo1, fgeo1_reflect, rtol=5e-2)
    fgeo2_symm = automol.geom.almost_equal_coulomb_spectrum(
        fgeo2, fgeo2_reflect, rtol=5e-2)

    # Set the face and face_sym keywords based on the above tests
    if fgeo1_symm and fgeo2_symm:
        faces = [0, 1]
        face_symm = 4
    elif fgeo1_symm and not fgeo2_symm:
        faces = [0, 1]
        face_symm = 2
    elif fgeo1_symm and not fgeo2_symm:
        faces = [0, 1]
        face_symm = 2
    elif not fgeo1_symm and not fgeo2_symm:
        faces = [0]
        face_symm = 1

    return faces, face_symm


def build_pivot_frames(min_idx, max_idx,
                       total_geom, frag_geoms, frag_geoms_wdummy):
    """ use geometries to get pivot info
        only set up for 1 or 2 pivot points
    """

    frames, npivots, = [], []
    geom_data = zip([min_idx, max_idx], frag_geoms, frag_geoms_wdummy)
    for i, (rxn_idx, geom, geom_wdummy) in enumerate(geom_data):

        # Single pivot point centered on atom
        print('geom test in build_pivot:', geom)
        if automol.geom.is_atom(geom):
            npivot = 1
            frame = [0, 0, 0, 0]
        # For linear species we place the pivot point on radical
        # with no displacment, so no need to find coordinates
        elif automol.geom.is_linear(geom):
            npivot = 2
            frame = [0, 0, 0, 0]
        else:
            # else we build an xy frame to easily place pivot point
            npivot = 2

            # Find the idx in each fragment bonded to the atom at the pivot pt
            for j, coords in enumerate(geom):
                if coords == total_geom[rxn_idx]:
                    bond_idx = j
                    break

            # For each fragment, get indices for a
            # chain (up to three atoms, that terminates at the dummy atom)
            gra = automol.geom.graph(geom)
            gra_neighbor_dct = automol.graph.atom_neighbor_keys(gra)
            bond_neighbors = gra_neighbor_dct[bond_idx]

            # Find idx in each fragment geom that corresponds to the bond index
            for j, idx in enumerate(bond_neighbors):
                if geom[idx][0] != 'H':
                    bond_neighbor_idx = idx
                    break
                elif geom[idx][0] == 'H' and j == (len(bond_neighbors) - 1):
                    bond_neighbor_idx = idx

            # Set up the frame indices for the divsur file
            if i == 0:
                pivot_idx = len(geom)
                frame = [bond_idx, bond_neighbor_idx, pivot_idx, bond_idx]
            else:
                pivot_idx = 0
                bond_idx += 1
                bond_neighbor_idx += 1
                frame = [bond_idx, bond_neighbor_idx, pivot_idx, bond_idx]
            frame = [val+1 for val in frame]

        # Append to lists
        frames.append(frame)
        npivots.append(npivot)

    return frames, npivots


def calc_pivot_angles(frag_geoms, frag_geoms_wdummy, frames):
    """ get the angle for the three atoms definining the frame
    """
    angles = []
    for geom, geom_wdummy, frame in zip(frag_geoms, frag_geoms_wdummy, frames):
        if automol.geom.is_atom(geom) or automol.geom.is_linear(geom):
            angle = None
        else:
            frame = [val-1 for val in frame]
            angle = automol.geom.central_angle(
                geom_wdummy, frame[2], frame[0], frame[1])

        angles.append(angle)

    return angles


def calc_pivot_xyzs(min_idx, max_idx, total_geom, frag_geoms):
    """ figure out where pivot point will be centered
        only linear speces need to have non-zero xyz, as they
        will not have a frame set up for them like atoms and
        polyatomics
    """
    xyzs = []
    for rxn_idx, geom in zip([min_idx, max_idx], frag_geoms):
        if automol.geom.is_linear(geom):
            xyz = total_geom[rxn_idx][1]
        else:
            xyz = [0.0, 0.0, 0.0]

        xyzs.append(xyz)

    return xyzs
