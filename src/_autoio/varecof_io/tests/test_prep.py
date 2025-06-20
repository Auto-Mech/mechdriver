""" test automol.varecof_io.writer.prep
"""

import numpy
import automol
import varecof_io.writer


# C2H6 + H,OH, CH3
CH3CH2_ZMA = (
    ('C', (None, None, None), (None, None, None),
     (None, None, None)),
    ('C', (0, None, None), ('R1', None, None),
     (2.8034554750087315, None, None)),
    ('H', (0, 1, None), ('R2', 'A2', None),
     (2.0477351881268895, 2.0981132181844426, None)),
    ('H', (0, 1, 2), ('R3', 'A3', 'D3'),
     (2.047734870248051, 2.0981200780710814, 3.2899037073022956)),
    ('H', (1, 0, 2), ('R4', 'A4', 'D4'),
     (2.0658681830076175, 1.924883755832615, 1.4966025648868544)),
    ('H', (1, 0, 4), ('R5', 'A5', 'D5'),
     (2.066515824593951, 1.9255815274161476, 2.0975452731245743)),
    ('H', (1, 0, 4), ('R6', 'A6', 'D6'),
     (2.0665179597230052, 1.9255796786166188, 4.185639485867064)))
H_ZMA = (
    ('H', (None, None, None), (None, None, None), (None, None, None)),)
OH_ZMA = (
    ('O', (None, None, None), (None, None, None),
     (None, None, None)),
    ('H', (0, None, None), ('R1', None, None),
     (1.8477888590491485, None, None)))
CH3_ZMA = (
    ('C', (None, None, None), (None, None, None),
     (None, None, None)),
    ('H', (0, None, None), ('R1', None, None),
     (2.04572282078478, None, None)),
    ('H', (0, 1, None), ('R2', 'A2', None),
     (2.045723120545338, 2.0943969516247156, None)),
    ('H', (0, 1, 2), ('R3', 'A3', 'D3'),
     (2.0457228854610228, 2.0943921145330617, 3.1415961306090057)))
CH3CH2_H_ZMA = (
    ('C', (None, None, None), (None, None, None),
     (None, None, None)),
    ('C', (0, None, None), ('R1', None, None),
     (2.857375545788057, None, None)),
    ('H', (0, 1, None), ('R2', 'A2', None),
     (2.067518646663273, 1.929774680570445, None)),
    ('H', (0, 1, 2), ('R3', 'A3', 'D3'),
     (2.067519421597435, 1.92977503072225, 2.0943955639047283)),
    ('H', (0, 1, 2), ('R4', 'A4', 'D4'),
     (2.0675198156854693, 1.9297748255978684, 4.188790396308539)),
    ('H', (1, 0, 2), ('R5', 'A5', 'D5'),
     (2.0675195325608655, 1.9297751281189566, 3.141594868756542)),
    ('H', (1, 0, 5), ('R6', 'A6', 'D6'),
     (2.0675190726592323, 1.9297747755130488, 4.18878925757192)),
    ('H', (1, 0, 5), ('R7', 'A7', 'D7'),
     (10.000, 1.9297750256469979, 2.0943942776405065)))
FRM_KEYS = (1, 7)
CH3CH2_OH_ZMA = ()
CH3CH2_CH3_ZMA = ()


def test__constraint():
    """ test varecof_io.writer.intramolecular_constraint_dct
    """

    ref_const_dct = {
        'R1': 2.857375545788057,
        'R2': 2.067518646663273,
        'A2': 1.929774680570445,
        'R3': 2.067519421597435,
        'A3': 1.92977503072225,
        'D3': 2.0943955639047283,
        'R4': 2.0675198156854693,
        'A4': 1.9297748255978684,
        'D4': 4.188790396308539,
        'R5': 2.0675195325608655,
        'A5': 1.9297751281189566,
        'D5': 3.141594868756542,
        'R6': 2.0675190726592323,
        'A6': 1.9297747755130488,
        'D6': 4.18878925757192}

    const_dct = varecof_io.writer.intramolecular_constraint_dct(
        CH3CH2_H_ZMA, (CH3CH2_ZMA, H_ZMA))
    assert ref_const_dct == const_dct  # add numpy check


def test__pivot():
    """ test varecof_io.writer.fragment_geometries
        test varecof_io.writer.build_pivot_frames
        test varecof_io.writer.calc_pivot_angles
        test varecof_io.writer.calc_pivot_xyzs
    """

    ref_tot_geo = (
        ('C', (0.0, 0.0, 0.0)),
        ('C', (0.0, 0.0, 2.857375545788057)),
        ('H', (0.0, 1.935727223564657, -0.726356297040919)),
        ('H', (1.6763889117962993, -0.9678646210561269, -0.726357247087748)),
        ('H', (-1.6763899924114845, -0.9678637852978632, -0.726356988473018)),
        ('H', (4.287959697162433e-06, -1.9357277279070393, 3.583733020392789)),
        ('H', (-1.6763904632916056, 0.9678616511395963, 3.583732176271997)),
        ('H', (8.108209279316, 4.6812918417618885, 6.370557717988566)))

    ref_isol_geos = (
        (('C', (0.0, 0.0, 2.8034554750087315)),
         ('C', (0.0, 0.0, 0.0)),
         ('H', (0.8473387628772707, -1.7427409842207437, 3.5213402402323344)),
         ('H', (1.0956829577601233, 1.5983181400102573, 3.521337399315168)),
         ('H', (-1.9323777167961476, 0.1436340237618356, 3.5197633532363044)),
         ('H', (0.0, 1.7695715856390393, -1.0304541736496604)),
         ('H', (-0.26148486078388367, -1.750138072983014, -1.030466152722124)),
         ('X', (9.33894282675913, -0.6646035686457461, -3.513182172200509))),
        (('H', (8.108209279316, 4.6812918417618889, 6.370557717988566)),)
    )

    ref_a1_idxs = (1, 1)

    tot_geo, isol_fgeos, a1_idxs = varecof_io.writer.fragment_geometries(
        CH3CH2_H_ZMA, (CH3CH2_ZMA, H_ZMA), FRM_KEYS)

    assert automol.geom.almost_equal_dist_matrix(ref_tot_geo, tot_geo)
    assert automol.geom.almost_equal_dist_matrix(
        ref_isol_geos[0], isol_fgeos[0])
    assert automol.geom.almost_equal_dist_matrix(
        ref_isol_geos[1], isol_fgeos[1])
    assert ref_a1_idxs == a1_idxs

    ref_frames = ((2, 1, 8, 2), (0, 0, 0, 0))
    ref_npivots = (2, 1)
    ref_angles = (1.9297750256469979, None)
    ref_xyzs = ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))

    frames, npivots = varecof_io.writer.build_pivot_frames(
        isol_fgeos, a1_idxs)
    angles = varecof_io.writer.calc_pivot_angles(
        isol_fgeos, frames)
    xyzs = varecof_io.writer.calc_pivot_xyzs(
        tot_geo, isol_fgeos, FRM_KEYS)

    assert ref_frames == frames
    assert ref_npivots == npivots
    assert numpy.isclose(ref_angles[0], angles[0])
    assert angles[1] is None
    for ref_xyz, xyz in zip(ref_xyzs, xyzs):
        assert numpy.allclose(ref_xyz, xyz)


def test__face_symm():
    """ test varecof_io.writer.assess_face_symmetries(fgeo1, fgeo2)
    """

    fgeo1 = automol.zmat.geometry(CH3CH2_ZMA)
    fgeo2 = automol.zmat.geometry(H_ZMA)
    faces, face_symm = varecof_io.writer.assess_face_symmetries(fgeo1, fgeo2)

    print(faces)
    print(face_symm)

    # assert faces == ()
    # assert face_symm = 0


if __name__ == '__main__':
    test__pivot()
