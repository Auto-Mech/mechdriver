""" test autofile.system
"""
import os
import numbers
import tempfile
import numpy
import pytest
import automol
import autofile.info
import autofile.system

PREFIX = tempfile.mkdtemp()
print(PREFIX)

# create a dummy root DataSeries for testing
ROOT_SPEC_DFILE = autofile.system.file_.locator(
    file_prefix='dir',
    map_dct_={
        'loc1': lambda locs: locs[0],
        'loc2': lambda locs: locs[1],
        'other': lambda locs: 'something else',
    },
    loc_keys=['loc1', 'loc2'],
)


def root_data_series_directory(prefix):
    """ root DataSeries
    """
    return autofile.system.model.DataSeries(
        prefix,
        map_=lambda x: os.path.join(*map(str, x)),
        nlocs=2,
        depth=2,
        loc_dfile=ROOT_SPEC_DFILE,)


def test__file__input_file():
    """ test autofile.system.file_.input_file
    """
    ref_inp_str = '<input file contents>'

    inp_dfile = autofile.system.file_.input_file('test')

    assert not inp_dfile.exists(PREFIX)
    inp_dfile.write(ref_inp_str, PREFIX)
    assert inp_dfile.exists(PREFIX)

    inp_str = inp_dfile.read(PREFIX)
    assert inp_str == ref_inp_str
    print(inp_str)


def test__file__output_file():
    """ test autofile.system.file_.output_file
    """
    ref_out_str = '<output file contents>'

    out_dfile = autofile.system.file_.output_file('test')

    assert not out_dfile.exists(PREFIX)
    out_dfile.write(ref_out_str, PREFIX)
    assert out_dfile.exists(PREFIX)

    out_str = out_dfile.read(PREFIX)
    assert out_str == ref_out_str
    print(out_str)


def test__file__information():
    """ test autofile.system.file_.information
    """
    def information(nsamp, tors_ranges):
        """ base information object
        """
        tors_ranges = autofile.info.Info(**dict(tors_ranges))
        assert isinstance(nsamp, numbers.Integral)
        inf_obj = autofile.info.Info(nsamp=nsamp, tors_ranges=tors_ranges)
        assert autofile.info.matches_function_signature(inf_obj, information)
        return inf_obj

    ref_inf_obj = information(
        nsamp=4, tors_ranges={'d1': (0., 1.), 'd2': (0., 3.)})

    inf_dfile = autofile.system.file_.information('test', function=information)

    assert not inf_dfile.exists(PREFIX)
    inf_dfile.write(ref_inf_obj, PREFIX)
    assert inf_dfile.exists(PREFIX)

    inf_obj = inf_dfile.read(PREFIX)
    assert inf_obj == ref_inf_obj
    print(inf_obj)


def test__file__energy():
    """ test autofile.system.file_.energy
    """
    ref_ene = -187.38518070487598

    ene_dfile = autofile.system.file_.energy('test')

    assert not ene_dfile.exists(PREFIX)
    ene_dfile.write(ref_ene, PREFIX)
    assert ene_dfile.exists(PREFIX)

    ene = ene_dfile.read(PREFIX)
    assert numpy.isclose(ene, ref_ene)
    print(ene)


def test__file__geometry():
    """ test autofile.system.file_.geometry
    """
    ref_geo = (('C', (0.066541036329, -0.86543409422, -0.56994517889)),
               ('O', (0.066541036329, -0.86543409422, 2.13152981129)),
               ('O', (0.066541036329, 1.6165813318, -1.63686376233)),
               ('H', (-1.52331011945, -1.99731957213, -1.31521725797)),
               ('H', (1.84099386813, -1.76479255185, -1.16213243427)),
               ('H', (-1.61114836922, -0.17751142359, 2.6046492029)),
               ('H', (-1.61092727126, 2.32295906780, -1.19178601663)))

    geo_dfile = autofile.system.file_.geometry('test')

    assert not geo_dfile.exists(PREFIX)
    geo_dfile.write(ref_geo, PREFIX)
    assert geo_dfile.exists(PREFIX)

    geo = geo_dfile.read(PREFIX)
    assert automol.geom.almost_equal(geo, ref_geo)
    print(geo)


def test__file__gradient():
    """ test autofile.system.file_.gradient
    """
    ref_grad = ((0.00004103632, 0.00003409422, 0.00004517889),
                (0.00004103632, 0.00003409422, 0.00002981129),
                (0.00004103632, 0.00008133180, 0.00006376233),
                (0.00001011945, 0.00001957213, 0.00001725797),
                (0.00009386813, 0.00009255185, 0.00003243427),
                (0.00004836922, 0.00001142359, 0.00004920290),
                (0.00002727126, 0.00005906780, 0.00008601663))

    grad_dfile = autofile.system.file_.gradient('test')

    assert not grad_dfile.exists(PREFIX)
    grad_dfile.write(ref_grad, PREFIX)
    assert grad_dfile.exists(PREFIX)

    grad = grad_dfile.read(PREFIX)
    assert numpy.allclose(grad, ref_grad)
    print(grad)


def test__file__hessian():
    """ test autofile.system.file_.hessian
    """
    ref_hess = (
        (-0.21406, 0., 0., -0.06169, 0., 0., 0.27574, 0., 0.),
        (0., 2.05336, 0.12105, 0., -0.09598, 0.08316, 0., -1.95737, -0.20421),
        (0., 0.12105, 0.19177, 0., -0.05579, -0.38831, 0., -0.06525, 0.19654),
        (-0.06169, 0., 0., 0.0316, 0., 0., 0.03009, 0., 0.),
        (0., -0.09598, -0.05579, 0., 0.12501, -0.06487, 0., -0.02902,
         0.12066),
        (0., 0.08316, -0.38831, 0., -0.06487, 0.44623, 0., -0.01829,
         -0.05792),
        (0.27574, 0., 0., 0.03009, 0., 0., -0.30583, 0., 0.),
        (0., -1.95737, -0.06525, 0., -0.02902, -0.01829, 0., 1.9864,
         0.08354),
        (0., -0.20421, 0.19654, 0., 0.12066, -0.05792, 0., 0.08354,
         -0.13862))

    hess_dfile = autofile.system.file_.hessian('test')

    assert not hess_dfile.exists(PREFIX)
    hess_dfile.write(ref_hess, PREFIX)
    assert hess_dfile.exists(PREFIX)

    hess = hess_dfile.read(PREFIX)
    assert numpy.allclose(hess, ref_hess)
    print(hess)


def test__file__zmatrix():
    """ test autofile.system.file_.zmatrix
    """
    ref_zma = (
        (('C', (None, None, None), (None, None, None)),
         ('O', (0, None, None), ('r1', None, None)),
         ('O', (0, 1, None), ('r2', 'a1', None)),
         ('H', (0, 1, 2), ('r3', 'a2', 'd1')),
         ('H', (0, 1, 2), ('r4', 'a3', 'd2')),
         ('H', (1, 0, 2), ('r5', 'a4', 'd3')),
         ('H', (2, 0, 1), ('r6', 'a5', 'd4'))),
        {'r1': 2.65933,
         'r2': 2.65933, 'a1': 1.90743,
         'r3': 2.06844, 'a2': 1.93366, 'd1': 4.1477,
         'r4': 2.06548, 'a3': 1.89469, 'd2': 2.06369,
         'r5': 1.83126, 'a4': 1.86751, 'd3': 1.44253,
         'r6': 1.83126, 'a5': 1.86751, 'd4': 4.84065})

    zma_dfile = autofile.system.file_.zmatrix('test')

    assert not zma_dfile.exists(PREFIX)
    zma_dfile.write(ref_zma, PREFIX)
    assert zma_dfile.exists(PREFIX)

    zma = zma_dfile.read(PREFIX)
    assert automol.zmatrix.almost_equal(zma, ref_zma)
    print(zma)


def test__file__vmatrix():
    """ test autofile.system.file_.vmatrix
    """
    ref_vma = (('C', (None, None, None), (None, None, None)),
               ('O', (0, None, None), ('r1', None, None)),
               ('O', (0, 1, None), ('r2', 'a1', None)),
               ('H', (0, 1, 2), ('r3', 'a2', 'd1')),
               ('H', (0, 1, 2), ('r4', 'a3', 'd2')),
               ('H', (1, 0, 2), ('r5', 'a4', 'd3')),
               ('H', (2, 0, 1), ('r6', 'a5', 'd4')))

    vma_dfile = autofile.system.file_.vmatrix('test')

    assert not vma_dfile.exists(PREFIX)
    vma_dfile.write(ref_vma, PREFIX)
    assert vma_dfile.exists(PREFIX)

    vma = vma_dfile.read(PREFIX)
    assert vma == ref_vma
    print(vma)


def test__file__trajectory():
    """ test autofile.system.file_.trajectory
    """
    ref_geos = [
        (('C', (0.0, 0.0, 0.0)),
         ('O', (0.0, 0.0, 2.699694868173)),
         ('O', (0.0, 2.503038629201, -1.011409768236)),
         ('H', (-1.683942509299, -1.076047850358, -0.583313101501)),
         ('H', (1.684063451772, -0.943916309940, -0.779079279468)),
         ('H', (1.56980872050, 0.913848877557, 3.152002706027)),
         ('H', (-1.57051358834, 3.264399836517, -0.334901043405))),
        (('C', (0.0, 0.0, 0.0)),
         ('O', (0.0, 0.0, 2.70915105770)),
         ('O', (0.0, 2.55808068205, -0.83913477573)),
         ('H', (-1.660164085463, -1.04177010816, -0.73213470306)),
         ('H', (1.711679909369, -0.895873802652, -0.779058492481)),
         ('H', (0.0238181080852, -1.813377410537, 3.16912929390)),
         ('H', (-1.36240560905, 3.348313125118, 0.1732746576216)))]
    ref_comments = ['energy: -187.3894105487809',
                    'energy: -187.3850624381528']

    ref_traj = list(zip(ref_comments, ref_geos))

    traj_dfile = autofile.system.file_.trajectory('test')

    assert not traj_dfile.exists(PREFIX)
    traj_dfile.write(ref_traj, PREFIX)
    assert traj_dfile.exists(PREFIX)

    # I'm not going to bother implementing a reader, since the trajectory files
    # are for human use only -- we aren't going to use this for data storage


def test__file__lennard_jones_epsilon():
    """ test autofile.system.file_.lennard_jones_epsilon
    """
    ref_eps = 247.880866746988

    eps_dfile = autofile.system.file_.lennard_jones_epsilon('test')

    assert not eps_dfile.exists(PREFIX)
    eps_dfile.write(ref_eps, PREFIX)
    assert eps_dfile.exists(PREFIX)

    eps = eps_dfile.read(PREFIX)
    assert numpy.isclose(eps, ref_eps)
    print(eps)


def test__file__lennard_jones_sigma():
    """ test autofile.system.file_.lennard_jones_sigma
    """
    ref_sig = 3.55018590361446

    sig_dfile = autofile.system.file_.lennard_jones_sigma('test')

    assert not sig_dfile.exists(PREFIX)
    sig_dfile.write(ref_sig, PREFIX)
    assert sig_dfile.exists(PREFIX)

    sig = sig_dfile.read(PREFIX)
    assert numpy.isclose(sig, ref_sig)
    print(sig)


def test__dir__run_trunk():
    """ test dir_.run_trunk
    """
    prefix = os.path.join(PREFIX, 'run_trunk')
    os.mkdir(prefix)

    root_ds = root_data_series_directory(prefix)
    ds_ = autofile.system.dir_.run_trunk(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]

    for root_locs in root_locs_lst:
        locs = root_locs

        assert not ds_.exists(locs)
        ds_.create(locs)
        assert ds_.exists(locs)

    assert sorted(root_ds.existing()) == sorted(root_locs_lst)


def test__dir__run_leaf():
    """ test dir_.run_leaf
    """
    prefix = os.path.join(PREFIX, 'run_leaf')
    os.mkdir(prefix)

    root_ds = root_data_series_directory(prefix)
    ds_ = autofile.system.dir_.run_leaf(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]
    leaf_locs_lst = [
        ['energy'],
        ['gradient'],
        ['hessian'],
        ['optimization'],
    ]

    for root_locs in root_locs_lst:
        for leaf_locs in leaf_locs_lst:
            locs = root_locs + leaf_locs

            assert not ds_.exists(locs)
            ds_.create(locs)
            assert ds_.exists(locs)

    assert sorted(root_ds.existing()) == sorted(root_locs_lst)

    print(ds_.existing(root_locs_lst[-1]))
    for root_locs in root_locs_lst:
        assert (sorted(ds_.existing(root_locs, relative=True)) ==
                sorted(leaf_locs_lst))


def test__dir__subrun_leaf():
    """ test dir_.subrun_leaf
    """
    prefix = os.path.join(PREFIX, 'subrun_leaf')
    os.mkdir(prefix)

    root_ds = root_data_series_directory(prefix)
    ds_ = autofile.system.dir_.subrun_leaf(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]
    leaf_locs_lst = [
        [0, 0],
        [0, 1],
        [0, 2],
        [1, 0],
        [1, 1],
        [2, 0],
    ]

    for root_locs in root_locs_lst:
        for leaf_locs in leaf_locs_lst:
            locs = root_locs + leaf_locs

            assert not ds_.exists(locs)
            ds_.create(locs)
            assert ds_.exists(locs)

    assert sorted(root_ds.existing()) == sorted(root_locs_lst)

    print(ds_.existing(root_locs_lst[-1]))
    for root_locs in root_locs_lst:
        assert (sorted(ds_.existing(root_locs, relative=True)) ==
                sorted(leaf_locs_lst))


def test__dir__species_trunk():
    """ test dir_.species_trunk
    """
    prefix = os.path.join(PREFIX, 'species_trunk')
    os.mkdir(prefix)

    # without a root directory
    ds_ = autofile.system.dir_.species_trunk(prefix)

    assert not ds_.exists()
    ds_.create()
    assert ds_.exists()

    # with a root directory
    root_ds = root_data_series_directory(prefix)
    ds_ = autofile.system.dir_.species_trunk(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]

    for root_locs in root_locs_lst:
        locs = root_locs

        assert not ds_.exists(locs)
        ds_.create(locs)
        assert ds_.exists(locs)

    assert sorted(root_ds.existing()) == sorted(root_locs_lst)


def test__dir__species_leaf():
    """ test dir_.species_leaf
    """
    prefix = os.path.join(PREFIX, 'species_leaf')
    os.mkdir(prefix)

    root_ds = root_data_series_directory(prefix)
    ds_ = autofile.system.dir_.species_leaf(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]
    branch_locs_lst = [
        ['InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+', 0, 1],
        ['InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1-', 0, 1],
        ['InChI=1S/C5H5O/c1-2-3-4-5-6/h1-5H/b4-3-', 0, 2],
        ['InChI=1S/O', 0, 1],
        ['InChI=1S/O', 0, 3],
    ]

    for root_locs in root_locs_lst:
        for branch_locs in branch_locs_lst:
            locs = root_locs + branch_locs

            assert not ds_.exists(locs)
            ds_.create(locs)
            assert ds_.exists(locs)

    assert sorted(root_ds.existing()) == sorted(root_locs_lst)

    print(ds_.existing(root_locs_lst[-1]))
    for root_locs in root_locs_lst:
        assert (sorted(ds_.existing(root_locs, relative=True)) ==
                sorted(branch_locs_lst))


def test__dir__reaction_trunk():
    """ test dir_.reaction_trunk
    """
    prefix = os.path.join(PREFIX, 'reaction_trunk')
    os.mkdir(prefix)

    # without a root directory
    ds_ = autofile.system.dir_.reaction_trunk(prefix)

    assert not ds_.exists()
    ds_.create()
    assert ds_.exists()

    # with a root directory
    root_ds = root_data_series_directory(prefix)
    ds_ = autofile.system.dir_.reaction_trunk(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]

    for root_locs in root_locs_lst:
        locs = root_locs

        assert not ds_.exists(locs)
        ds_.create(locs)
        assert ds_.exists(locs)

    assert sorted(root_ds.existing()) == sorted(root_locs_lst)


def test__dir__reaction_leaf():
    """ test dir_.reaction_leaf
    """
    prefix = os.path.join(PREFIX, 'reaction_leaf')
    os.mkdir(prefix)

    root_ds = root_data_series_directory(prefix)
    ds_ = autofile.system.dir_.reaction_leaf(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]

    branch_locs_lst = [
        [
            [['InChI=1S/C2H5O2/c1-2-4-3/h3H,1-2H2'],
             ['InChI=1S/C2H4/c1-2/h1-2H2', 'InChI=1S/HO2/c1-2/h1H']],
            [[0], [0, 0]],
            [[2], [1, 2]],
            2,
        ],
        [
            [['InChI=1S/CH2/h1H2', 'InChI=1S/H2O/h1H2'],
             ['InChI=1S/CH3/h1H3', 'InChI=1S/HO/h1H']],
            [[0, 0], [0, 0]],
            [[1, 1], [2, 2]],
            1,
        ],
        [
            [['InChI=1S/CH3O3/c2-1-4-3/h2H,1H2'],
             ['InChI=1S/CH3O3/c2-1-4-3/h3H,1H2']],
            [[0], [0]],
            [[2], [2]],
            2,
        ],
    ]

    for root_locs in root_locs_lst:
        for branch_locs in branch_locs_lst:
            locs = root_locs + branch_locs

            assert not ds_.exists(locs)
            ds_.create(locs)
            assert ds_.exists(locs)

    assert sorted(root_ds.existing()) == sorted(root_locs_lst)

    print(ds_.existing(root_locs_lst[-1]))
    for root_locs in root_locs_lst:
        assert (sorted(ds_.existing(root_locs, relative=True)) ==
                sorted(branch_locs_lst))


def test__dir__theory_leaf():
    """ test dir_.theory_leaf
    """
    prefix = os.path.join(PREFIX, 'theory_leaf')
    os.mkdir(prefix)

    root_ds = root_data_series_directory(prefix)
    ds_ = autofile.system.dir_.theory_leaf(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]
    branch_locs_lst = [
        ['hf', 'sto-3g', True],
        ['hf', 'sto-3g', False],
        ['b3lyp', 'sto-3g', False],
        ['b3lyp', '6-31g*', False],
    ]

    for root_locs in root_locs_lst:
        for branch_locs in branch_locs_lst:
            locs = root_locs + branch_locs

            assert not ds_.exists(locs)
            ds_.create(locs)
            assert ds_.exists(locs)

    assert sorted(root_ds.existing()) == sorted(root_locs_lst)

    print(ds_.existing(root_locs_lst[-1]))
    for root_locs in root_locs_lst:
        assert (sorted(ds_.existing(root_locs, relative=True)) ==
                sorted(branch_locs_lst))


def test__dir__conformer_trunk():
    """ test dir_.conformer_trunk
    """
    prefix = os.path.join(PREFIX, 'conformer_trunk')
    os.mkdir(prefix)

    # without a root directory
    ds_ = autofile.system.dir_.conformer_trunk(prefix)

    assert not ds_.exists()
    ds_.create()
    assert ds_.exists()

    # with a root directory
    root_ds = root_data_series_directory(prefix)
    ds_ = autofile.system.dir_.conformer_trunk(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]

    for root_locs in root_locs_lst:
        locs = root_locs

        assert not ds_.exists(locs)
        ds_.create(locs)
        assert ds_.exists(locs)

    assert sorted(root_ds.existing()) == sorted(root_locs_lst)


def test__dir__conformer_leaf():
    """ test dir_.conformer_leaf
    """
    prefix = os.path.join(PREFIX, 'conformer_leaf')
    os.mkdir(prefix)

    root_ds = root_data_series_directory(prefix)
    ds_ = autofile.system.dir_.conformer_leaf(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]

    nconfs = 10
    branch_locs_lst = [
        [autofile.system.generate_new_conformer_id()] for _ in range(nconfs)]

    for root_locs in root_locs_lst:
        for branch_locs in branch_locs_lst:
            locs = root_locs + branch_locs

            assert not ds_.exists(locs)
            ds_.create(locs)
            assert ds_.exists(locs)

    assert sorted(root_ds.existing()) == sorted(root_locs_lst)

    print(ds_.existing(root_locs_lst[-1]))
    for root_locs in root_locs_lst:
        assert (sorted(ds_.existing(root_locs, relative=True)) ==
                sorted(branch_locs_lst))


def test__dir__single_point_trunk():
    """ test dir_.single_point_trunk
    """
    prefix = os.path.join(PREFIX, 'single_point_trunk')
    os.mkdir(prefix)

    root_ds = root_data_series_directory(prefix)
    ds_ = autofile.system.dir_.single_point_trunk(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]

    for root_locs in root_locs_lst:
        locs = root_locs

        assert not ds_.exists(locs)
        ds_.create(locs)
        assert ds_.exists(locs)

    assert sorted(root_ds.existing()) == sorted(root_locs_lst)


def test__dir__scan_trunk():
    """ test dir_.scan_trunk
    """
    prefix = os.path.join(PREFIX, 'scan_trunk')
    os.mkdir(prefix)

    # without a root directory
    ds_ = autofile.system.dir_.scan_trunk(prefix)

    assert not ds_.exists()
    ds_.create()
    assert ds_.exists()

    # with a root directory
    root_ds = root_data_series_directory(prefix)
    ds_ = autofile.system.dir_.scan_trunk(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]

    for root_locs in root_locs_lst:
        locs = root_locs

        assert not ds_.exists(locs)
        ds_.create(locs)
        assert ds_.exists(locs)

    assert sorted(root_ds.existing()) == sorted(root_locs_lst)


def test__dir__scan_branch():
    """ test dir_.scan_branch
    """
    prefix = os.path.join(PREFIX, 'scan_branch')
    os.mkdir(prefix)

    root_ds = root_data_series_directory(prefix)
    ds_ = autofile.system.dir_.scan_branch(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]
    branch_locs_lst = [
        [['d3']],
        [['d3', 'd4']],
    ]

    for root_locs in root_locs_lst:
        for branch_locs in branch_locs_lst:
            locs = root_locs + branch_locs

            assert not ds_.exists(locs)
            ds_.create(locs)
            assert ds_.exists(locs)

    assert sorted(root_ds.existing()) == sorted(root_locs_lst)

    print(ds_.existing(root_locs_lst[-1]))
    for root_locs in root_locs_lst:
        assert (sorted(ds_.existing(root_locs, relative=True)) ==
                sorted(branch_locs_lst))


def test__dir__scan_leaf():
    """ test dir_.scan_leaf
    """
    prefix = os.path.join(PREFIX, 'scan_leaf')
    os.mkdir(prefix)

    root_ds = root_data_series_directory(prefix)
    ds_ = autofile.system.dir_.scan_leaf(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]
    leaf_locs_lst = [
        [[0, 0]],
        [[1, 0]],
        [[2, 0]],
        [[0, 1]],
        [[1, 1]],
        [[2, 1]],
        [[0, 2]],
        [[1, 2]],
        [[2, 2]],
        [[0, 3]],
        [[1, 3]],
        [[2, 3]],
        [[0, 4]],
        [[1, 4]],
        [[2, 4]],
    ]

    for root_locs in root_locs_lst:
        for leaf_locs in leaf_locs_lst:
            locs = root_locs + leaf_locs

            assert not ds_.exists(locs)
            ds_.create(locs)
            assert ds_.exists(locs)

    assert sorted(root_ds.existing()) == sorted(root_locs_lst)

    print(ds_.existing(root_locs_lst[-1]))
    for root_locs in root_locs_lst:
        assert (sorted(ds_.existing(root_locs, relative=True)) ==
                sorted(leaf_locs_lst))


def test__dir__tau_trunk():
    """ test dir_.tau_trunk
    """
    prefix = os.path.join(PREFIX, 'tau_trunk')
    os.mkdir(prefix)

    # without a root directory
    ds_ = autofile.system.dir_.tau_trunk(prefix)

    assert not ds_.exists()
    ds_.create()
    assert ds_.exists()

    # with a root directory
    root_ds = root_data_series_directory(prefix)
    ds_ = autofile.system.dir_.tau_trunk(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]

    for root_locs in root_locs_lst:
        locs = root_locs

        assert not ds_.exists(locs)
        ds_.create(locs)
        assert ds_.exists(locs)

    assert sorted(root_ds.existing()) == sorted(root_locs_lst)


def test__dir__tau_leaf():
    """ test dir_.tau_leaf
    """
    prefix = os.path.join(PREFIX, 'tau_leaf')
    os.mkdir(prefix)

    root_ds = root_data_series_directory(prefix)
    ds_ = autofile.system.dir_.tau_leaf(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]

    nconfs = 10
    branch_locs_lst = [
        [autofile.system.generate_new_tau_id()] for _ in range(nconfs)]

    for root_locs in root_locs_lst:
        for branch_locs in branch_locs_lst:
            locs = root_locs + branch_locs

            assert not ds_.exists(locs)
            ds_.create(locs)
            assert ds_.exists(locs)

    assert sorted(root_ds.existing()) == sorted(root_locs_lst)

    print(ds_.existing(root_locs_lst[-1]))
    for root_locs in root_locs_lst:
        assert (sorted(ds_.existing(root_locs, relative=True)) ==
                sorted(branch_locs_lst))

    with pytest.raises(ValueError):
        ds_.remove(locs)
    assert ds_.exists(locs)


def test__dir__build_trunk():
    """ test dir_.build_trunk
    """
    prefix = os.path.join(PREFIX, 'build_trunk')
    os.mkdir(prefix)

    root_ds = root_data_series_directory(prefix)
    ds_ = autofile.system.dir_.build_trunk(prefix, root_ds=root_ds)

    root_alocs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]
    rlocs_lst = [
        ['MESS'],
    ]

    for root_alocs in root_alocs_lst:
        for rlocs in rlocs_lst:
            alocs = root_alocs + rlocs

            assert not ds_.exists(alocs)
            ds_.create(alocs)
            assert ds_.exists(alocs)

    assert sorted(root_ds.existing()) == sorted(root_alocs_lst)

    print(ds_.existing(root_alocs_lst[-1]))
    for root_alocs in root_alocs_lst:
        assert (sorted(ds_.existing(root_alocs, relative=True)) ==
                sorted(rlocs_lst))


if __name__ == '__main__':
    # test__file__input_file()
    # test__file__information()
    # test__file__energy()
    # test__file__geometry()
    # test__file__gradient()
    # test__file__hessian()
    # test__file__zmatrix()
    # test__file__vmatrix()
    # test__file__trajectory()
    # test__file__lennard_jones_epsilon()
    # test__file__lennard_jones_sigma()
    test__dir__run_trunk()
    test__dir__run_leaf()
    test__dir__subrun_leaf()
    test__dir__species_trunk()
    test__dir__species_leaf()
    test__dir__reaction_trunk()
    test__dir__reaction_leaf()
    test__dir__theory_leaf()
    test__dir__conformer_trunk()
    test__dir__conformer_leaf()
    test__dir__single_point_trunk()
    test__dir__scan_trunk()
    test__dir__scan_branch()
    test__dir__scan_leaf()
    test__dir__tau_trunk()
    test__dir__tau_leaf()
    test__dir__build_trunk()
