"""
 routine to do vrc tst
"""

import automol
import varecof_io
from pivot import build_pivot_frames
from pivot import calc_pivot_angles
from pivot import calc_pivot_xyzs


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


def vrctst():
    """ run the routine for vrctst
    """

    # run mep scan (some code in scripts/es.py)

    # calculate the infinite seperation energy

    # Get the geometry from a zmat on the MEP (here is ZMA)
    total_geom = automol.zmatrix.geometry(ZMA)

    # Build list of the fragment geometries:
    # (1) without dummy atoms and,
    # (2) with dummy atom assuming position of 2nd bond atom
    frag_geoms = [GEO[:RXN_IDXS[1]], GEO[RXN_IDXS[1]:]]
    dummy_row = ('X', GEO[RXN_IDXS[1]][1])
    geo1_wdummy = frag_geoms[0] + (dummy_row,)
    dummy_row = ('X', GEO[RXN_IDXS[0]][1])
    geo2_wdummy = (dummy_row,) + frag_geoms[1]
    frag_geoms_wdummy = [geo1_wdummy, geo2_wdummy]

    # Get pivot point info
    frames, npivots = build_pivot_frames(
        RXN_IDXS2, total_geom, frag_geoms, frag_geoms_wdummy)
    angles = calc_pivot_angles(frag_geoms, frag_geoms_wdummy, frames)
    xyzs = calc_pivot_xyzs(RXN_IDXS2, total_geom, frag_geoms)

    # Write the long- and short-range divsur input files
    lr_divsur_inp_str = varecof_io.writer.input_file.divsur(
        RDISTS_LR, 1, 1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
    print('\nlong-range divsur input file:')
    print(lr_divsur_inp_str)

    # Write the short-range divsur files
    t1angs = [angles[0]] if angles[0] is not None else None
    t2angs = [angles[1]] if angles[0] is not None else None
    sr_divsur_inp_str = varecof_io.writer.input_file.divsur(
        RDISTS_SR, npivots[0], npivots[1], xyzs[0], xyzs[1],
        frame1=frames[0],
        frame2=frames[1],
        d1dists=D1DISTS,
        d2dists=D2DISTS,
        t1angs=t1angs,
        t2angs=t2angs)
    print('\nshort-range divsur input file:')
    print(sr_divsur_inp_str)

    # Write the structure input files
    struct_inp_str = varecof_io.writer.input_file.structure(
        frag_geoms_wdummy[0], frag_geoms_wdummy[1])
    print('\nstructure.inp:')
    print(struct_inp_str)

    # Write the tst.inp file
    nsamp_max = 2000
    nsamp_min = 500
    flux_err = 5
    pes_size = 1
    tst_inp_str = varecof_io.writer.input_file.tst(
        nsamp_max, nsamp_min, flux_err, pes_size)
    print('\ntst.inp:')
    print(tst_inp_str)


