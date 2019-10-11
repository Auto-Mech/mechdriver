"""
 routine to do vrc tst
"""

import automol
import varecof_io


# GEOM INPUT FOR C2H5 + OH RXN
GEO = (('C', (-2.5191523897, -0.5633331563, -0.0362123265)),
       ('C', (0.3603207881, -0.5331062394, -0.0031081202)),
       ('H', (-3.2419651558, -0.1921939889, -1.9771748689)),
       ('H', (-3.2382988468, 0.9670747472, 1.2069612071)),
       ('H', (1.0576235240, 1.3283804620, -0.6737608211)),
       ('H', (1.1103827994, -2.0334799496, -1.2645938535)),
       ('H', (1.0582003205, -0.8619737454, 1.9486756979)),
       ('O', (-5.7269478282, -6.6386704961, 0.8614633075)),
       ('H', (-5.3970505334, -7.9045286243, -0.4840514744)))
RXN_IDXS = [0, 7]

# GEOM INPUT FOR CH3 + H RXN
ZMA = ((('C', (None, None, None), (None, None, None)),
        ('H', (0, None, None), ('R1', None, None)),
        ('H', (0, 1, None), ('R2', 'A2', None)),
        ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
        ('H', (0, 1, 2), ('R4', 'A4', 'D4'))),
       {'R1': 2.09646,
        'R2': 2.09646, 'A2': 1.9106293854507126,
        'R3': 2.09646, 'A3': 1.9106293854507126, 'D3': 4.1887902047863905,
        'R4': 12.09646, 'A4': 1.9106293854507126, 'D4': 2.0943951023931953})
RXN_IDXS2 = [0, 4]

# Grid for the distances for points in the divsur file (in angstrom)
RDISTS_LR = [10.5, 9.0, 8.0, 7.5, 7.0, 6.5, 6.0,
             5.5, 5.25, 5.0, 4.5, 4.25, 4.0]
RDISTS_SR = [4.25, 4.0, 3.75, 3.5, 3.25, 3.0, 2.75],
D1DISTS = [0.05, 0.15, 0.25],
D2DISTS = [0.05, 0.15, 0.25],


def test__build_pivot_frames():
    """ run the routine for vrctst
        ret params:
        ret npivots:
        ret angles
    """

    # Build fragment geometries
    geoms = [GEO[:RXN_IDXS[1]], GEO[RXN_IDXS[1]:]]

    # Build fragment geomeries w/ dummy atom assuming position of 2nd bond atom
    dummy_row = ('X', GEO[RXN_IDXS[1]][1])
    geo1_wdummy = geoms[0] + (dummy_row,)
    dummy_row = ('X', GEO[RXN_IDXS[0]][1])
    geo2_wdummy = (dummy_row,) + geoms[1]
    geoms_wdummy = [geo1_wdummy, geo2_wdummy]

    # Loop over fragment geometries to get the info needed
    frames, npivots, angles = [], [], []
    geom_data = zip(RXN_IDXS, geoms, geoms_wdummy)
    for i, (rxn_idx, geom, geom_wdummy) in enumerate(geom_data):

        if automol.geom.is_atom(geom):
            npivot = 1
            frame = [0, 0, 0, 0]
            angle = None
        else:
            npivot = 2

            # Set idx for the pivot dummy atom for each fragment (wdummy)
            if i == 0:
                pivot_idx = len(geom)
            else:
                pivot_idx = 0

            # Find the idx in each fragment bonded to the atom at the pivot pt
            for j, coords in enumerate(geom):
                if coords == GEO[rxn_idx]:
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
                # Compute needed coordinates for defining positions in divsur
                angle = automol.geom.central_angle(
                    geom_wdummy, pivot_idx, bond_idx, bond_neighbor_idx)
            else:
                pivot_idx = 0
                frame = [bond_idx+1, bond_neighbor_idx+1, pivot_idx, bond_idx+1]
                # Compute needed coordinates for defining positions in divsur
                angle = automol.geom.central_angle(
                    geom_wdummy, pivot_idx, bond_idx+1, bond_neighbor_idx+1)
            frame = [val+1 for val in frame]

        # Append to lists
        frames.append(frame)
        npivots.append(npivot)
        angles.append(angle)

    return frames, npivots, angles, geoms_wdummy


def test__write_lr_divsur():
    """ write a long-range divsur file
    """
    # Write the long-range divsur files
    rdists_lr = [10.5, 9.0, 8.0, 7.5, 7.0, 6.5, 6.0,
                 5.5, 5.25, 5.0, 4.5, 4.25, 4.0]
    lr_divsur_inp_str = varecof_io.writer.input_file.divsur(
        rdists_lr, 1, 1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])

    print('\nlong-range divsur input file:')
    print(lr_divsur_inp_str)


def test__write_sr_divsur(frames, npivots, angles,
                          rdists=[4.25, 4.0, 3.75, 3.5, 3.25, 3.0, 2.75],
                          d2dists=[0.05, 0.15, 0.25],
                          d1dists=[0.05, 0.15, 0.25]):
    """ write a short-range divsur file
    """
    sr_divsur_inp_str = varecof_io.writer.input_file.divsur(
        rdists, npivots[0], npivots[1], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0],
        frame1=frames[0],
        frame2=frames[1],
        d1dists=d1dists,
        d2dists=d2dists,
        t1angs=[angles[0]],
        t2angs=[angles[1]])

    print('\nshort-range divsur input file:')
    print(sr_divsur_inp_str)


if __name__ == '__main__':
    print('mep geom')
    print(automol.geom.string(GEO))
    FRAMES, NPIVOTS, ANGLES, GEOMS_WDUMMY = test__build_pivot_frames()
    print('\ngeo1')
    print(automol.geom.string(GEOMS_WDUMMY[0]))
    print('\ngeo2')
    print(automol.geom.string(GEOMS_WDUMMY[1]))
    test__write_lr_divsur()
    test__write_sr_divsur(FRAMES, NPIVOTS, ANGLES)
