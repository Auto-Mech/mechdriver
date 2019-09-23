"""
 routine to do vrc tst
"""

import automol
import varecof_io


GEO = (('C', (-2.5191523897, -0.5633331563, -0.0362123265)),
       ('C', (0.3603207881, -0.5331062394, -0.0031081202)),
       ('H', (-3.2419651558, -0.1921939889, -1.9771748689)),
       ('H', (-3.2382988468, 0.9670747472, 1.2069612071)),
       ('H', (1.0576235240, 1.3283804620, -0.6737608211)),
       ('H', (1.1103827994, -2.0334799496, -1.2645938535)),
       ('H', (1.0582003205, -0.8619737454, 1.9486756979)),
       ('O', (-5.7269478282, -6.6386704961, 0.8614633075)),
       ('H', (-5.3970505334, -7.9045286243, -0.4840514744)))
BND_IDX = [0, 7]


def test__routine():
    """ run the routine for vrctst
    """

    # Get the fragment geoms, including a dummy
    # make switch symbol and switch coord fxn in automol.geom
    geo1 = GEO[:BND_IDX[1]]
    dummy_row = ('X', GEO[BND_IDX[1]][1])
    geo1_wdummy = geo1 + (dummy_row,)
    geo2 = GEO[BND_IDX[1]:]
    dummy_row = ('X', GEO[BND_IDX[0]][1])
    geo2_wdummy = (dummy_row,) + geo2

    print('GEO')
    for x in GEO:
        print(x)
    print('geo1')
    for x in geo1:
        print(x)
    print('geo1_wdummy')
    for x in geo1_wdummy:
        print(x)
    print('geo2')
    for x in geo2:
        print(x)
    print('geo2_wdummy')
    for x in geo2_wdummy:
        print(x)

    # For each fragment, get indices for a
    # chain (up to three atoms, that terminates at the dummy atom)
    gra1 = automol.geom.graph(geo1)
    gra2 = automol.geom.graph(geo2)
    gra1_neighbor_dct = automol.graph.atom_neighbor_keys(gra1)
    gra2_neighbor_dct = automol.graph.atom_neighbor_keys(gra2)

    # Find the index in each fragment geom that corresponds to the bond index
    for i, coords in enumerate(geo1):
        if coords == GEO[BND_IDX[0]]:
            f1_bnd_idx = i
    for i, coords in enumerate(geo2):
        if coords == GEO[BND_IDX[1]]:
            f2_bnd_idx = i
    f1_bnd_neighbors = gra1_neighbor_dct[f1_bnd_idx]
    f2_bnd_neighbors = gra2_neighbor_dct[f2_bnd_idx]

    print('f1_bnd_idx')
    print(f1_bnd_idx)
    print('f2_bnd_idx')
    print(f2_bnd_idx)

    print('f1 bond neighbors')
    print(f1_bnd_neighbors)
    print('f2 bond neighbors')
    print(f2_bnd_neighbors)

    # Get the bond atom keys of the neighbors
    for i, idx in enumerate(f1_bnd_neighbors):
        if geo1[idx][0] != 'H':
            f1_bnd_neighbor_idx = idx
            break
        elif geo1[idx][0] == 'H' and i == (len(f1_bnd_neighbors) - 1):
            f1_bnd_neighbor_idx = idx
    for i, idx in enumerate(f2_bnd_neighbors):
        if geo2[idx][0] != 'H':
            f2_bnd_neighbor_idx = idx
            break
        elif geo2[idx][0] == 'H' and i == (len(f2_bnd_neighbors) - 1):
            f2_bnd_neighbor_idx = idx

    print('f1 bond_neighbor idx')
    print(f1_bnd_neighbor_idx)
    print('f2 bond_neighbor idx')
    print(f2_bnd_neighbor_idx)

    # Build the list of indexes for fragments starting with pivot
    f1_idxs = [len(geo1), f1_bnd_idx, f1_bnd_neighbor_idx]
    f2_idxs = [0, f2_bnd_idx+1, f2_bnd_neighbor_idx+1]
    print('f1_idxs')
    print(f1_idxs)
    print('f2_idxs')
    print(f2_idxs)

    # Compute the needed coordinates for defining positions
    f1_dist = automol.geom.distance(
        geo1_wdummy, f1_idxs[0], f1_idxs[1])
    f1_angle = automol.geom.central_angle(
        geo1_wdummy, f1_idxs[0], f1_idxs[1], f1_idxs[2])
    f2_dist = automol.geom.distance(
        geo2_wdummy, f2_idxs[0], f2_idxs[1])
    f2_angle = automol.geom.central_angle(
        geo2_wdummy, f2_idxs[0], f2_idxs[1], f2_idxs[2])

    print('f1 coords')
    print(f1_dist)
    print(f1_angle)
    print('f2 coords')
    print(f2_dist)
    print(f2_angle)

    # Set up the frame indices
    f1_frame = [f1_idxs[1], f1_idxs[0], f1_idxs[2], f1_idxs[1]]
    f2_frame = [f2_idxs[1], f2_idxs[0], f2_idxs[2], f2_idxs[1]]
    print('f1_frame')
    print(f1_frame)
    print('f2_frame')
    print(f2_frame)

    # Write the long-range and short-range divsur files
    rdists_lr = [10.5, 9.0, 8.0, 7.5, 7.0, 6.5, 6.0,
                 5.5, 5.25, 5.0, 4.5, 4.25, 4.0]
    lr_divsur_inp_str = varecof_io.writer.input_file.divsur(
        rdists_lr, 1, 1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])

    # Write the short-range divsur files
    if automol.geom.is_atom(geo1):
        npivot1 = 1
    else:
        npivot1 = 2
    if automol.geom.is_atom(geo2):
        npivot2 = 1
    else:
        npivot2 = 2

    # xyz might just be 0.0,0.0,0.0
    # pivot_xyz1 = geo1[f1_bnd_idx][1]
    # pivot_xyz2 = geo2[f2_bnd_idx][1]
    rdists_sr = [4.25, 4.0, 3.75, 3.5, 3.25, 3.0, 2.75]
    d1dists = [0.05, 0.15, 0.25]
    d2dists = [0.05, 0.15, 0.25]
    sr_divsur_inp_str = varecof_io.writer.input_file.divsur(
        rdists_sr, npivot1, npivot2, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0],
        frame1=[val + 1 for val in f1_frame],
        frame2=[val + 1 for val in f2_frame],
        d1dists=d1dists,
        d2dists=d2dists,
        t1angs=[f1_angle],
        t2angs=[f2_angle])

    print('\n\n\ndivsur input files')
    print('\nlong-range divsur input file:')
    print(lr_divsur_inp_str)
    print('\n\nshort-range divsur input file:')
    print(sr_divsur_inp_str)


if __name__ == '__main__':
    test__routine()
