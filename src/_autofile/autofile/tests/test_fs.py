""" test autofile.fs
"""

import os
import tempfile
import numpy
import pytest

import automol
import autofile.fs

PREFIX = tempfile.mkdtemp()
print(PREFIX)


def test__species():
    """ test autofile.fs.species
    """
    prefix = os.path.join(PREFIX, 'species')
    os.mkdir(prefix)

    locs = ['InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+', 0, 1]
    spc_fs = autofile.fs.species(prefix)

    assert not spc_fs[-1].exists(locs)
    spc_fs[-1].create(locs)
    assert spc_fs[-1].exists(locs)


def test__reaction():
    """ test autofile.fs.reaction
    """
    prefix = os.path.join(PREFIX, 'reaction')
    os.mkdir(prefix)

    locs = [
        [['InChI=1S/C2H5O2/c1-2-4-3/h3H,1-2H2'],
         ['InChI=1S/C2H4/c1-2/h1-2H2', 'InChI=1S/HO2/c1-2/h1H']],
        [[0], [0, 0]],
        [[2], [1, 2]],
        2
    ]
    rxn_fs = autofile.fs.reaction(prefix)

    assert not rxn_fs[-1].exists(locs)
    rxn_fs[-1].create(locs)
    assert rxn_fs[-1].exists(locs)


def test__transition_state():
    """ test autofile.fs.transition_state
    """
    prefix = os.path.join(PREFIX, 'ts')
    os.mkdir(prefix)

    ts_fs = autofile.fs.transition_state(prefix)

    locs = [0]

    assert not ts_fs[-1].exists(locs)
    ts_fs[-1].create(locs)
    assert ts_fs[-1].exists(locs)

    # ref_vma = ()
    # ts_fs[-1].file.vmatrix.write(ref_vma, locs)
    # assert ts_fs[-1].file.vmatrix.read(locs) == ref_vma


def test__theory():
    """ test autofile.fs.theory
    """
    prefix = os.path.join(PREFIX, 'theory')
    os.mkdir(prefix)

    thy_fs = autofile.fs.theory(prefix)
    locs = ['hf', 'sto-3g', 'U']

    thy_fs[-1].create(locs)
    thy_path = thy_fs[-1].path(locs)
    print(thy_path)
    assert os.path.exists(thy_path)

    # active space test
    locs = ['casscf', '6-31g', [1104, 1103]]

    thy_fs[-1].create(locs)
    thy_path = thy_fs[-1].path(locs)
    print(thy_path)
    assert os.path.exists(thy_path)


def test__conformer():
    """ test autofile.fs.conformer
    """
    prefix = os.path.join(PREFIX, 'conformer')
    os.mkdir(prefix)

    rid = autofile.schema.generate_new_ring_id()
    cid = autofile.schema.generate_new_conformer_id()
    locs = [rid, cid]

    cnf_fs = autofile.fs.conformer(prefix)
    assert not cnf_fs[-1].exists(locs)
    cnf_fs[-1].create(locs)
    assert cnf_fs[-1].exists(locs)

    ref_inf_obj1 = autofile.schema.info_objects.conformer_trunk(1)
    ref_inf_obj2 = autofile.schema.info_objects.conformer_branch(1)
    cnf_fs[0].file.info.write(ref_inf_obj1)
    cnf_fs[1].file.info.write(ref_inf_obj2, [rid])

    inf_obj1 = cnf_fs[0].file.info.read()
    inf_obj2 = cnf_fs[1].file.info.read([rid])
    assert inf_obj1 == ref_inf_obj1
    assert inf_obj2 == ref_inf_obj2

    ref_geo = (
        ('O', (-1.1707387949811348, -0.8186819289555977, 0.1847602826946391)),
        ('O', (1.2517867512044956, 0.33310347836345394, 0.7805302053715455)),
        ('H', (-2.2059402840308078, 0.6739299449565304, 0.505211869486996)),
        ('H', (2.124892327807454, -0.18835149436438584, -0.7582977119813665)))
    cnf_fs[-1].file.geometry.write(ref_geo, locs)
    geo = cnf_fs[-1].file.geometry.read(locs)
    assert automol.geom.almost_equal_dist_matrix(ref_geo, geo)

    ref_vpt2_inf_obj = autofile.schema.info_objects.vpt2(
        fermi_treatment='gaussian_default')
    cnf_fs[-1].file.vpt2_info.write(ref_vpt2_inf_obj, locs)
    assert cnf_fs[-1].file.vpt2_info.read(locs) == ref_vpt2_inf_obj


def test__symmetry():
    """ test autofile.fs.symmetry
    """

    prefix = os.path.join(PREFIX, 'symmetry')
    os.mkdir(prefix)

    symm_locs = [autofile.schema.generate_new_conformer_id()]

    symm_fs = autofile.fs.symmetry(prefix)
    symm_fs[-1].create(symm_locs)

    ref_geo = (
        ('O', (-1.1707387949811348, -0.8186819289555977, 0.1847602826946391)),
        ('O', (1.2517867512044956, 0.33310347836345394, 0.7805302053715455)),
        ('H', (-2.2059402840308078, 0.6739299449565304, 0.505211869486996)),
        ('H', (2.124892327807454, -0.18835149436438584, -0.7582977119813665)))
    symm_fs[-1].file.geometry.write(ref_geo, symm_locs)
    geo = symm_fs[-1].file.geometry.read(symm_locs)
    assert automol.geom.almost_equal_dist_matrix(ref_geo, geo)


def test__single_point():
    """ test autofile.fs.single_point
    """
    prefix = os.path.join(PREFIX, 'single_point')
    os.mkdir(prefix)

    sp_fs = autofile.fs.single_point(prefix)
    locs = ['hf', 'sto-3g', 'U']

    ref_ene = 5.7
    sp_fs[-1].create(locs)
    sp_fs[-1].file.energy.write(ref_ene, locs)
    assert sp_fs[-1].file.energy.read(locs) == ref_ene


def test__high_spin():
    """ test autofile.fs.high_spin
    """
    prefix = os.path.join(PREFIX, 'high_spin')
    os.mkdir(prefix)

    sp_fs = autofile.fs.high_spin(prefix)
    locs = ['hf', 'sto-3g', 'U']

    ref_ene = 5.7
    sp_fs[-1].create(locs)
    sp_fs[-1].file.energy.write(ref_ene, locs)
    assert sp_fs[-1].file.energy.read(locs) == ref_ene


def test__zmatrix():
    """ test autofile.fs.zmatrix
    """
    prefix = os.path.join(PREFIX, 'zmatrix')
    os.mkdir(prefix)

    zma_fs = autofile.fs.zmatrix(prefix)
    locs = [0]

    ref_zma = (
        ('O', (None, None, None), (None, None, None), (None, None, None)),
        ('H', (0, None, None), ('R1', None, None), (1.84779, None, None)))

    zma_fs[-1].create(locs)
    zma_fs[-1].file.zmatrix.write(ref_zma, locs)
    assert automol.zmat.almost_equal(
        zma_fs[-1].file.zmatrix.read(locs), ref_zma)

    ref_rxn = automol.reac.from_forward_reverse(
        cla=automol.ReactionClass.HYDROGEN_ABSTRACTION,
        ftsg=(
            {0: ('C', 0, None), 1: ('H', 0, None), 2: ('H', 0, None),
             3: ('H', 0, None), 4: ('H', 0, None), 5: ('O', 0, None),
             6: ('H', 0, None)},
            {frozenset({5, 6}): (1, None), frozenset({0, 3}): (1, None),
             frozenset({4, 5}): (0.1, None),
             frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
             frozenset({0, 4}): (0.9, None)}),
        rtsg=(
            {0: ('O', 0, None), 1: ('H', 0, None), 2: ('H', 0, None),
             3: ('C', 0, None), 4: ('H', 0, None), 5: ('H', 0, None),
             6: ('H', 0, None)},
            {frozenset({0, 2}): (0.9, None),
             frozenset({2, 3}): (0.1, None),
             frozenset({3, 6}): (1, None), frozenset({0, 1}): (1, None),
             frozenset({3, 4}): (1, None), frozenset({3, 5}): (1, None)}),
        rcts_keys=((0, 1, 2, 3, 4), (5, 6)),
        prds_keys=((0, 1, 2), (3, 4, 5, 6)),
    )

    zma_fs[-1].file.reaction.write(ref_rxn, locs)
    rxn = zma_fs[-1].file.reaction.read(locs)

    assert rxn == ref_rxn


def test__scan():
    """ test autofile.fs.scan
    """
    prefix = os.path.join(PREFIX, 'scan')
    os.mkdir(prefix)

    scn_fs = autofile.fs.scan(prefix)
    locs = [['d3', 'd4'], [2, 3]]

    ref_inp_str = '<input string>'
    scn_fs[-1].create(locs)
    scn_fs[-1].file.geometry_input.write(ref_inp_str, locs)
    assert scn_fs[-1].file.geometry_input.read(locs) == ref_inp_str

    inf_obj = autofile.schema.info_objects.scan_branch(
        {'d4': numpy.array([0, 2])})
    scn_fs[1].create(['d4'])
    scn_fs[1].file.info.write(inf_obj, ['d4'])
    assert scn_fs[1].file.info.read(['d4']) == inf_obj


def test__cscan():
    """ test autofile.fs.cscan
    """
    prefix = os.path.join(PREFIX, 'cscan')
    os.mkdir(prefix)

    scn_fs = autofile.fs.cscan(prefix)
    locs = [{'R1': 1., 'A2': 2.3}, ['D3', 'D4'], [1.2, 2.9]]

    ref_inp_str = '<input string>'
    scn_fs[-1].create(locs)
    scn_fs[-1].file.geometry_input.write(ref_inp_str, locs)
    assert scn_fs[-1].file.geometry_input.read(locs) == ref_inp_str


def test__tau():
    """ test autofile.fs.tau
    """
    prefix = os.path.join(PREFIX, 'tau')
    os.mkdir(prefix)

    locs = [autofile.schema.generate_new_tau_id()]
    tau_fs = autofile.fs.tau(prefix)

    assert not tau_fs[-1].exists(locs)
    tau_fs[-1].create(locs)
    assert tau_fs[-1].exists(locs)

    ref_inf_obj = autofile.schema.info_objects.tau_trunk(
        2, {'D1': (3.14159, 6.28319)})
    tau_fs[0].file.info.write(ref_inf_obj)
    inf_obj = tau_fs[0].file.info.read()
    assert inf_obj.nsamp == ref_inf_obj.nsamp
    assert numpy.allclose(dict(ref_inf_obj.tors_ranges)['D1'],
                          dict(inf_obj.tors_ranges)['D1'])


def test__vrctst():
    """ test autofile.fs.vrctst
    """

    prefix = os.path.join(PREFIX, 'vrctst')
    os.mkdir(prefix)

    vrc_fs = autofile.fs.vrctst(prefix)
    locs = [0]

    ref_tst_str = '<TST STR>'
    ref_divsur_str = '<DIVSUR STR>'
    ref_molpro_str = '<MOLPRO STR>'
    ref_tml_str = '<TML STR>'
    ref_struct_str = '<STRUCT STR>'
    ref_pot_str = '<POT STR>'
    ref_flux_str = '<FLUX STR>'

    vrc_fs[-1].create(locs)
    vrc_fs[-1].file.vrctst_tst.write(ref_tst_str, locs)
    vrc_fs[-1].file.vrctst_divsur.write(ref_divsur_str, locs)
    vrc_fs[-1].file.vrctst_molpro.write(ref_molpro_str, locs)
    vrc_fs[-1].file.vrctst_tml.write(ref_tml_str, locs)
    vrc_fs[-1].file.vrctst_struct.write(ref_struct_str, locs)
    vrc_fs[-1].file.vrctst_pot.write(ref_pot_str, locs)
    vrc_fs[-1].file.vrctst_flux.write(ref_flux_str, locs)

    tst_str = vrc_fs[-1].file.vrctst_tst.read(locs)
    divsur_str = vrc_fs[-1].file.vrctst_divsur.read(locs)
    molpro_str = vrc_fs[-1].file.vrctst_molpro.read(locs)
    tml_str = vrc_fs[-1].file.vrctst_tml.read(locs)
    struct_str = vrc_fs[-1].file.vrctst_struct.read(locs)
    pot_str = vrc_fs[-1].file.vrctst_pot.read(locs)
    flux_str = vrc_fs[-1].file.vrctst_flux.read(locs)

    assert ref_tst_str == tst_str
    assert ref_divsur_str == divsur_str
    assert ref_molpro_str == molpro_str
    assert ref_tml_str == tml_str
    assert ref_struct_str == struct_str
    assert ref_pot_str == pot_str
    assert ref_flux_str == flux_str


def test__energy_transfer():
    """ test autofile.fs.energy_transfer
    """
    prefix = os.path.join(PREFIX, 'ETRANS')
    os.mkdir(prefix)

    etrans_fs = autofile.fs.energy_transfer(prefix)

    bath_locs = ['InChI=1S/N2/c1-2', 0, 1]
    thy_locs = ['hf', 'sto-3g', 'R']
    locs = bath_locs + thy_locs

    etrans_fs[-1].create(locs)

    ref_eps = 300.0
    ref_sig = 3.50
    ref_inp_str = '<INP STR>'
    ref_temp_str = '<TEMP STR>'
    ref_traj = (
        ((('H', (0.0, 0.0, 0.0)),), 'comment 1'),
        ((('H', (0.0, 0.0, 0.0)),), 'comment 2')
    )
    ref_inf_obj = autofile.schema.info_objects.lennard_jones(
        nsamp=10, program='OneDMin', version='1.0')

    etrans_fs[-1].file.lennard_jones_epsilon.write(ref_eps, locs)
    etrans_fs[-1].file.lennard_jones_sigma.write(ref_sig, locs)
    etrans_fs[-1].file.lennard_jones_input.write(ref_inp_str, locs)
    etrans_fs[-1].file.lennard_jones_elstruct.write(ref_temp_str, locs)
    etrans_fs[-1].file.trajectory.write(ref_traj, locs)
    etrans_fs[-1].file.lennard_jones_info.write(ref_inf_obj, locs)
    assert numpy.isclose(
        etrans_fs[-1].file.lennard_jones_epsilon.read(locs), ref_eps)
    assert numpy.isclose(
        etrans_fs[-1].file.lennard_jones_sigma.read(locs), ref_sig)
    assert etrans_fs[-1].file.lennard_jones_input.read(locs) == ref_inp_str
    assert etrans_fs[-1].file.lennard_jones_elstruct.read(locs) == ref_temp_str
    assert etrans_fs[-1].file.trajectory.read(locs) == ref_traj
    assert etrans_fs[-1].file.lennard_jones_info.read(locs) == ref_inf_obj


def test__run():
    """ test autofile.fs.run
    """
    prefix = os.path.join(PREFIX, 'run')
    os.mkdir(prefix)

    run_fs = autofile.fs.run(prefix)
    run_locs = ['gradient']

    run_fs[-1].create(run_locs)

    ref_inf_obj = autofile.schema.info_objects.run(
        job='gradient', prog='psi4', version='1.0',
        method='ccsd(t)', basis='cc-pvdz', status='succeeded')
    ref_inf_obj.utc_start_time = autofile.schema.utc_time()

    run_fs[-1].create(run_locs)
    run_fs[-1].file.info.write(ref_inf_obj, run_locs)
    assert run_fs[-1].file.info.read(run_locs) == ref_inf_obj


def test__subrun():
    """ test autofile.fs.subrun
    """
    prefix = os.path.join(PREFIX, 'subrun')
    os.mkdir(prefix)

    subrun_fs = autofile.fs.subrun(prefix)
    subrun_locs = [0, 1]

    ref_inp_str = '<input string>'
    subrun_fs[-1].create(subrun_locs)
    subrun_fs[-1].file.input.write(ref_inp_str, subrun_locs)
    assert subrun_fs[-1].file.input.read(subrun_locs) == ref_inp_str


def test__build():
    """ test autofile.fs.build
    """
    prefix = os.path.join(PREFIX, 'build')
    os.mkdir(prefix)

    build_fs = autofile.fs.build(prefix)

    ref_inp_str = '<input string>'
    build_fs[-1].create(['MESS', 'C2H5O', 0])
    build_fs[-1].file.input.write(ref_inp_str, ['MESS', 'C2H5O', 0])
    assert build_fs[-1].file.input.read(['MESS', 'C2H5O', 0]) == ref_inp_str


def test__json_io():
    """ test autofile.json_.read_json and write_json
    """
    prefix = os.path.join(PREFIX, 'tau')
    file_path = os.path.join(prefix, 'test.json')
    jsona = {'a': 1}
    jsonb = {'b': 2}
    autofile.json_.write_json(jsona, file_path)
    autofile.json_.write_json(jsonb, file_path)
    assert autofile.json_.read_json(file_path) == jsonb
    with pytest.raises(IOError):
        json_err = {'eek'}
        autofile.json_.write_json(json_err, file_path)
        # assert autofile.json_.read_json(file_path) == jsonb
    # with pytest.raises(IOError):
    # json_err = 'not_a_json'
    # autofile.io_.write_file(file_path, json_err)
    # assert autofile.json_.read_json(file_path) == jsonb


def test__json_tau_save():
    """ test <fs>.json.<property>.write
    """
    prefix = os.path.join(PREFIX, 'tau')
    tau_fs = autofile.fs.tau(prefix)
    save_path = tau_fs[-1].root.path()
    print(tau_fs[-1].json_path())
    tau_fs[-1].root.create()
    tau_fs[-1].json_create()
    print(tau_fs[-1].json_existing())
    print(tau_fs[-1].json_path())
    locs = [autofile.schema.generate_new_tau_id()]
    print(tau_fs[-1].exists(locs))
    print(tau_fs[-1].json_path(locs))
    sp_locs = ['hf', 'sto-3g', 'U']
    ref_ene = 5.7
    print(tau_fs[-1].json.energy.exists(locs))
    tau_fs[-1].json.energy.write(ref_ene, locs)
    assert tau_fs[-1].json.energy.exists(locs)
    print(tau_fs[-1].json.geometry.exists(locs))
    print(tau_fs[-1].json.gradient.exists(locs))
    print(tau_fs[-1].json.geometry_info.exists(locs))
    assert tau_fs[-1].json.energy.read(locs) == ref_ene

    tau_fs[-1].json_existing(locs=locs)
    tau_fs[-1].json.energy.write_all([ref_ene], [locs])
    tau_fs[-1].json.energy.read_all([locs])
    sp_save_fs = autofile.fs.single_point(
        save_path, json_layer=locs)
    sp_save_fs[-1].json_existing()
    sp_save_fs[-1].json_existing(locs=sp_locs)
    tau_fs[-1].json_existing(locs=sp_locs, json_layer=locs)
    sp_ene = 5.5
    print(sp_save_fs[-1].json.energy.exists(sp_locs))
    sp_save_fs[-1].json.energy.write(sp_ene, sp_locs)
    assert sp_save_fs[-1].json.energy.exists(sp_locs)
    assert sp_save_fs[-1].json.energy.read(sp_locs) == sp_ene
    sp_save_fs[-1].json_existing()
    sp_save_fs[-1].json_existing(locs=sp_locs)
    tau_fs[-1].json_existing(locs=sp_locs, json_layer=locs)


def test__manager():
    """ test
    """

    prefix = os.path.join(PREFIX, 'data1')
    _build_fs(prefix)

    # Set locs for this test
    spc_locs = ['InChI=1S/H2O/h1H2', 0, 1]
    thy_locs = ['hf', 'sto-3g', 'U']
    cnf_locs_1 = ['rQ5VxakIXDkDp', 'cdgZx6pwjFtcX']
    cnf_locs_2 = ['rTW2L_4wU-PVV', 'caI7Oi5445ya2']

    # Set cnf fs using the manager
    cnf_fs = autofile.fs.manager(
        prefix, [['SPECIES', spc_locs], ['THEORY', thy_locs]], 'CONFORMER')

    # Test cnf fs with removing files
    # Removing files/layers withut setting removable
    rem1, rem2 = False, False
    try:
        cnf_fs[-1].remove(cnf_locs_1)
    except ValueError:
        rem1 = True
    try:
        cnf_fs[-1].file.geometry_input.remove(cnf_locs_1)
    except ValueError:
        rem2 = True
    assert rem1 and rem2

    cnf_fs[-1].removable = True
    if cnf_fs[-1].file.geometry_input.exists(cnf_locs_2):
        cnf_fs[-1].file.geometry_input.removable = True
        cnf_fs[-1].file.geometry_input.remove(cnf_locs_2)
    assert not cnf_fs[-1].file.geometry_input.exists(cnf_locs_2)
    cnf_fs[-1].remove(cnf_locs_2)


def test__iterate_managers():
    """ test autofile.fs.iterate_managers
    """

    def _get_set(locs):
        """ Figure out set name from cnf locs
        """
        for key, val in FAKE_LOCS_SETS_DCT.items():
            if locs == val[2]:
                set_name = key
                break
        return set_name

    # Build fs
    prefix = os.path.join(PREFIX, 'data2')
    _build_fs(prefix)

    # Generate managers for grapping parts of fake fileystem
    cnf_managers = tuple(autofile.fs.iterate_managers(
        prefix, ['SPECIES', 'THEORY'], 'CONFORMER'))

    # Loop over all cnf_fs and read file for correct data
    for cnf_fs in cnf_managers:
        for locs in cnf_fs[-1].existing():
            inp_str = cnf_fs[-1].file.geometry_input.read(locs)
            assert inp_str == FAKE_NAMES_DCT[_get_set(list(locs))]


def test__iterate_locators():
    """ test autofile.fs.iterate_locators
    """

    # Build fs
    prefix = os.path.join(PREFIX, 'data3')
    _build_fs(prefix)

    # Generate managers for grapping parts of fake fileystem
    spc_locators = tuple(autofile.fs.iterate_locators(
        prefix, ['SPECIES']))

    # Loop over species locs (generate geo from ich and match filesys)
    for spc_locs, in spc_locators:
        ich, _, _ = spc_locs
        ref_geo = automol.chi.geometry(ich)

        spc_fs = autofile.fs.species(prefix)
        thy_fs = autofile.fs.theory(spc_fs[-1].path(spc_locs))
        for thy_locs in thy_fs[-1].existing():
            cnf_fs = autofile.fs.conformer(thy_fs[-1].path(thy_locs))
            for cnf_locs in cnf_fs[-1].existing():
                geo = cnf_fs[-1].file.geometry.read(cnf_locs)
                assert automol.geom.almost_equal_dist_matrix(ref_geo, geo)


def test__path_prefix():
    """ test autofile.fs.path_prefix
    """
    prefix = os.path.join(PREFIX, 'data4')
    _build_fs(prefix)

    # Set locs for this test
    spc_locs = ['InChI=1S/H2O/h1H2', 0, 1]
    thy_locs = ['hf', 'sto-3g', 'U']
    cnf_locs_1 = ['rQ5VxakIXDkDp', 'cdgZx6pwjFtcX']

    # Set cnf fs using the manager
    cnf_fs = autofile.fs.manager(
        prefix, [['SPECIES', spc_locs], ['THEORY', thy_locs]], 'CONFORMER')
    spc_fs = autofile.fs.species(prefix)
    spc_prefix = spc_fs[-1].path(spc_locs)
    path_prefix = autofile.fs.path_prefix(
        cnf_fs[-1].path(cnf_locs_1), ['THEORY', 'CONFORMER'])
    assert spc_prefix == path_prefix


def _build_fs(prefix):
    """ Construct a filesystem to do stuff
    """

    os.mkdir(prefix)

    spc_fs = autofile.fs.species(prefix)
    for _set in ('set1', 'set2', 'set3', 'set4'):
        locs_set = FAKE_LOCS_SETS_DCT[_set]
        inp_name = FAKE_NAMES_DCT[_set]
        # spc layer
        spc_fs[-1].create(locs_set[0])
        spc_path = spc_fs[-1].path(locs_set[0])
        # thy layer
        thy_fs = autofile.fs.theory(spc_path)
        thy_fs[-1].create(locs_set[1])
        thy_path = thy_fs[-1].path(locs_set[1])
        # cnf layer
        cnf_fs = autofile.fs.conformer(thy_path)
        cnf_fs[-1].create(locs_set[2])
        cnf_fs[-1].file.geometry_input.write(inp_name, locs_set[2])
        geo = automol.chi.geometry(locs_set[0][0])
        cnf_fs[-1].file.geometry.write(geo, locs_set[2])


FAKE_LOCS_SETS_DCT = {
    'set1': [['InChI=1S/H2O/h1H2', 0, 1],
             ['hf', 'sto-3g', 'U'],
             ['rQ5VxakIXDkDp', 'cdgZx6pwjFtcX']],
    'set2': [['InChI=1S/H2O/h1H2', 0, 1],
             ['hf', 'sto-3g', 'U'],
             ['rTW2L_4wU-PVV', 'caI7Oi5445ya2']],
    'set3': [['InChI=1S/CH4/h1H4', 0, 1],
             ['hf', 'cc-pvdz', 'R'],
             ['r5-8BlffIKsw_', 'cWgRuV5qGsBrK']],
    'set4': [['InChI=1S/CH4/h1H4', 0, 1],
             ['hf', 'cc-pvdz', 'R'],
             ['ry0F_raMFxFnv', 'ckyd_uZJ-TmM7']]
}
FAKE_NAMES_DCT = {
    'set1': 'inp1',
    'set2': 'inp2',
    'set3': 'inp3',
    'set4': 'inp4'
}


if __name__ == '__main__':
    test__theory()
