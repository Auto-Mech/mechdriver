""" Generate the information necessary to product the vrctst input files
"""

import varecof_io
import automol

def input_prep(ts_zma, dist_name):
    """ prepare all the input files for a vrc-tst calculation
    """
    bnd_frm_idxs = automol.zmatrix.bond_idxs(ts_zma, dist_name)
    min_idx = min(bnd_frm_idxs)
    max_idx = max(bnd_frm_idxs)
    print('ts_zma and Bnd_frm_idxs test:', ts_zma, bnd_frm_idxs)
    total_geom, frag_geoms, frag_geoms_wdummy = fragment_geometries(
        ts_zma, min_idx, max_idx)
    frames, npivots = build_pivot_frames(
        min_idx, max_idx, total_geom, frag_geoms, frag_geoms_wdummy)
    angles = calc_pivot_angles(frag_geoms, frag_geoms_wdummy, frames)
    xyzs = calc_pivot_xyzs(min_idx, max_idx, total_geom, frag_geoms)

    # Write the long- and short-range divsur input files
    r1dists_lr = [8., 6., 5., 4.5, 4.]
    #rdists_lr = [15., 12., 10., 9., 8.]
    lr_divsur_inp_str = varecof_io.writer.input_file.divsur(
        r1dists_lr, 1, 1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
    print('\nlong-range divsur input file:')
    print(lr_divsur_inp_str)

    # Write the short-range divsur files
    r1dists_sr = [4., 3.8, 3.6, 3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2]
    r2dists_sr = [4., 3.8, 3.6, 3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2]
    #rdists_sr = [8., 7.5, 7., 6.5, 6., 5.5, 5.25, 5., 4.75, 4.5, 4.25, 4., 3.75, 3.5]
    d1dists = [0.01, 0.5, 1.]
    d2dists = [0.01, 0.5, 1.]
    t1angs = [angles[0]] if angles[0] is not None else []
    t2angs = [angles[1]] if angles[0] is not None else []
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
    #print('sr_divsur test:', rdists_sr, npivots, xyzs, frames, d1dists, d2dists, t1angs, t2angs)
    conditions = {'delta_r': 0}
    sr_divsur_inp_str = varecof_io.writer.input_file.divsur(
        r1dists_sr, npivots[0], npivots[1], xyzs[0], xyzs[1],
        frame1=frames[0],
        frame2=frames[1],
        d1dists=d1dists,
        d2dists=d2dists,
        t1angs=t1angs,
        t2angs=t2angs,
        r2dists=r2dists_sr,
        **conditions)
    print('\nshort-range divsur input file:')
    print(sr_divsur_inp_str)
                  # Write the structure input files
    struct_inp_str = varecof_io.writer.input_file.structure(
        frag_geoms_wdummy[0], frag_geoms_wdummy[1])
    print('\nstructure.inp:')
    print(struct_inp_str)

    # Write the tst.inp file
    nsamp_max = 2000
    nsamp_min = 50
    flux_err = 10
    pes_size = 2
    tst_inp_str = varecof_io.writer.input_file.tst(
        nsamp_max, nsamp_min, flux_err, pes_size)
    print('\ntst.inp:')
    print(tst_inp_str)

    # Write the potential energy surface input string
    exe_path = '/blues/gpfs/home/sjklipp/bin/molpro'
    base_name = 'mol'
    els_inp_str = varecof_io.writer.input_file.elec_struct(
        exe_path, base_name)
    print('\n\nels.inp:')
    print(els_inp_str)

    # Write the mc_flux.inp input string
    mc_flux_inp_str = varecof_io.writer.input_file.mc_flux()
    print('\n\nmc_flux.inp:')
    print(mc_flux_inp_str)

    # Write the convert.inp input string
    conv_inp_str = varecof_io.writer.input_file.convert()
    print('\n\nconvert.inp:')
    print(conv_inp_str)

    input_str = [
        struct_inp_str, lr_divsur_inp_str, sr_divsur_inp_str, tst_inp_str,
        els_inp_str, mc_flux_inp_str, conv_inp_str]

    return input_str


def fragment_geometries(ts_zma, min_idx, max_idx):
    """ Generate the fragment geometries from the ts Z-matrix and the indices involved in the forming bond
    """

    # Get the geometry from a zmat on the MEP (here is ZMA)
    total_geom = automol.zmatrix.geometry(ts_zma)

    # Build list of the fragment geometries:
    # (1) without dummy atoms and,
    # (2) with dummy atom assuming position of 2nd bond atom
    print('total_geom test:', total_geom)
    #bnd_frm_idxs = sorted(bnd_frm_idxs)
    #frag_geoms = [total_geom[:bnd_frm_idxs[1]], total_geom[bnd_frm_idxs[1]:]]
    frag_geoms = [total_geom[:max_idx], total_geom[max_idx:]]
    if not automol.geom.is_atom(frag_geoms[0]):
        dummy_row = ('X', total_geom[max_idx][1])
        geo1_wdummy = frag_geoms[0] + (dummy_row,)
    else:
        geo1_wdummy = frag_geoms[0]
    if not automol.geom.is_atom(frag_geoms[1]):
        dummy_row = ('X', total_geom[min_idx][1])
        geo2_wdummy = (dummy_row,) + frag_geoms[1]
    else:
        geo2_wdummy = frag_geoms[1]
    frag_geoms_wdummy = [geo1_wdummy, geo2_wdummy]

    return total_geom, frag_geoms, frag_geoms_wdummy


def build_pivot_frames(min_idx, max_idx, total_geom, frag_geoms, frag_geoms_wdummy):
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
