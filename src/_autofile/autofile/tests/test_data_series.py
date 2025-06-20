""" test autofile.schema.data_series
"""

import os
import tempfile
import pytest
import autofile.info
import autofile.schema


PREFIX = tempfile.mkdtemp()
print(PREFIX)

# create a dummy root DataSeries for testing
ROOT_SPEC_DFILE = autofile.schema.data_files.locator(
    file_prefix='dir',
    map_dct_={
        'loc1': lambda locs: locs[0],
        'loc2': lambda locs: locs[1],
        'other': lambda locs: 'something else',
    },
    loc_keys=['loc1', 'loc2'],
)


def root_data_series(prefix):
    """ root DataSeries
    """
    return autofile.model.DataSeries(
        prefix,
        map_=lambda x: os.path.join(*map(str, x)),
        nlocs=2,
        depth=2,
        loc_dfile=ROOT_SPEC_DFILE,)


def test__data_series__species_trunk():
    """ test data_series.species_trunk
    """
    prefix = os.path.join(PREFIX, 'species_trunk')
    os.mkdir(prefix)

    # without a root directory
    ds_ = autofile.schema.data_series.species_trunk(prefix)

    assert not ds_.exists()
    ds_.create()
    assert ds_.exists()

    # with a root directory
    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.species_trunk(prefix, root_ds=root_ds)

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


def test__data_series__species_leaf():
    """ test data_series.species_leaf
    """
    prefix = os.path.join(PREFIX, 'species_leaf')
    os.mkdir(prefix)

    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.species_leaf(prefix, root_ds=root_ds)

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


def test__data_series__reaction_trunk():
    """ test data_series.reaction_trunk
    """
    prefix = os.path.join(PREFIX, 'reaction_trunk')
    os.mkdir(prefix)

    # without a root directory
    ds_ = autofile.schema.data_series.reaction_trunk(prefix)

    assert not ds_.exists()
    ds_.create()
    assert ds_.exists()

    # with a root directory
    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.reaction_trunk(prefix, root_ds=root_ds)

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


def test__data_series__reaction_leaf():
    """ test data_series.reaction_leaf
    """
    prefix = os.path.join(PREFIX, 'reaction_leaf')
    os.mkdir(prefix)

    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.reaction_leaf(prefix, root_ds=root_ds)

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


def test__data_series__transition_state_trunk():
    """ test data_series.transition_state_trunk
    """
    prefix = os.path.join(PREFIX, 'transition_state_trunk')
    os.mkdir(prefix)

    # without a root directory
    ds_ = autofile.schema.data_series.transition_state_trunk(prefix)

    assert not ds_.exists()
    ds_.create()
    assert ds_.exists()

    # with a root directory
    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.transition_state_trunk(
        prefix, root_ds=root_ds)

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


def test__data_series__theory_leaf():
    """ test data_series.theory_leaf
    """
    prefix = os.path.join(PREFIX, 'theory_leaf')
    os.mkdir(prefix)

    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.theory_leaf(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]
    branch_locs_lst = [
        ['hf', 'sto-3g', 'R'],
        ['hf', 'sto-3g', 'U'],
        ['b3lyp', 'sto-3g', 'R'],
        ['b3lyp', '6-31g*', 'U'],
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


def test__data_series__conformer_trunk():
    """ test data_series.conformer_trunk
    """
    prefix = os.path.join(PREFIX, 'conformer_trunk')
    os.mkdir(prefix)

    # without a root directory
    ds_ = autofile.schema.data_series.conformer_trunk(prefix)

    assert not ds_.exists()
    ds_.create()
    assert ds_.exists()

    # with a root directory
    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.conformer_trunk(prefix, root_ds=root_ds)

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


def test__data_series__conformer_leaf():
    """ test data_series.conformer_leaf
    """
    prefix = os.path.join(PREFIX, 'conformer_leaf')
    os.mkdir(prefix)

    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.conformer_leaf(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]

    nconfs = 10
    branch_locs_lst = [
        [autofile.schema.generate_new_conformer_id()] for _ in range(nconfs)]

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


def test__data_series__symmetry_trunk():
    """ test data_series.symmetry_trunk
    """
    prefix = os.path.join(PREFIX, 'symmetry_trunk')
    os.mkdir(prefix)

    # without a root directory
    ds_ = autofile.schema.data_series.symmetry_trunk(prefix)

    assert not ds_.exists()
    ds_.create()
    assert ds_.exists()

    # with a root directory
    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.symmetry_trunk(prefix, root_ds=root_ds)

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


def test__data_series__zmatrix_trunk():
    """ test data_series.zmatrix_trunk
    """
    prefix = os.path.join(PREFIX, 'zmatrix_trunk')
    os.mkdir(prefix)

    # without a root directory
    ds_ = autofile.schema.data_series.zmatrix_trunk(prefix)

    assert not ds_.exists()
    ds_.create()
    assert ds_.exists()

    # with a root directory
    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.zmatrix_trunk(prefix, root_ds=root_ds)

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


def test__data_series__zmatrix_leaf():
    """ test data_series.zmatrix_leaf
    """
    prefix = os.path.join(PREFIX, 'zmatrix_leaf')
    os.mkdir(prefix)

    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.zmatrix_leaf(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]
    branch_locs_lst = [
        [0],
        [1],
        [2],
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


def test__data_series__single_point_trunk():
    """ test data_series.single_point_trunk
    """
    prefix = os.path.join(PREFIX, 'single_point_trunk')
    os.mkdir(prefix)

    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.single_point_trunk(prefix,
                                                         root_ds=root_ds)

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


def test__data_series__scan_trunk():
    """ test data_series.scan_trunk
    """
    prefix = os.path.join(PREFIX, 'scan_trunk')
    os.mkdir(prefix)

    # without a root directory
    ds_ = autofile.schema.data_series.scan_trunk(prefix)

    assert not ds_.exists()
    ds_.create()
    assert ds_.exists()

    # with a root directory
    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.scan_trunk(prefix, root_ds=root_ds)

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


def test__data_series__scan_branch():
    """ test data_series.scan_branch
    """
    prefix = os.path.join(PREFIX, 'scan_branch')
    os.mkdir(prefix)

    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.scan_branch(prefix, root_ds=root_ds)

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


def test__data_series__scan_leaf():
    """ test data_series.scan_leaf
    """
    prefix = os.path.join(PREFIX, 'scan_leaf')
    os.mkdir(prefix)

    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.scan_leaf(prefix, root_ds=root_ds)

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


def test__data_series__tau_trunk():
    """ test data_series.tau_trunk
    """
    prefix = os.path.join(PREFIX, 'tau_trunk')
    os.mkdir(prefix)

    # without a root directory
    ds_ = autofile.schema.data_series.tau_trunk(prefix)

    assert not ds_.exists()
    ds_.create()
    assert ds_.exists()

    # with a root directory
    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.tau_trunk(prefix, root_ds=root_ds)

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


def test__data_series__tau_leaf():
    """ test data_series.tau_leaf
    """
    prefix = os.path.join(PREFIX, 'tau_leaf')
    os.mkdir(prefix)

    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.tau_leaf(prefix, root_ds=root_ds)

    root_locs_lst = [
        [1, 'a'],
        [1, 'b'],
        [2, 'a'],
        [2, 'b'],
        [2, 'c'],
    ]

    nconfs = 10
    branch_locs_lst = [
        [autofile.schema.generate_new_tau_id()] for _ in range(nconfs)]

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


def test__data_series__run_trunk():
    """ test data_series.run_trunk
    """
    prefix = os.path.join(PREFIX, 'run_trunk')
    os.mkdir(prefix)

    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.run_trunk(prefix, root_ds=root_ds)

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


def test__data_series__run_leaf():
    """ test data_series.run_leaf
    """
    prefix = os.path.join(PREFIX, 'run_leaf')
    os.mkdir(prefix)

    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.run_leaf(prefix, root_ds=root_ds)

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


def test__data_series__subrun_leaf():
    """ test data_series.subrun_leaf
    """
    prefix = os.path.join(PREFIX, 'subrun_leaf')
    os.mkdir(prefix)

    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.subrun_leaf(prefix, root_ds=root_ds)

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


def test__data_series__build_trunk():
    """ test data_series.build_trunk
    """
    prefix = os.path.join(PREFIX, 'build_trunk')
    os.mkdir(prefix)

    root_ds = root_data_series(prefix)
    ds_ = autofile.schema.data_series.build_trunk(prefix, root_ds=root_ds)

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
