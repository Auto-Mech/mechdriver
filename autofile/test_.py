""" test autofile.fs
"""
import os
import tempfile
import autofile.fs

PREFIX = tempfile.mkdtemp()
print(PREFIX)


def test__direction():
    """ test autofile.fs.direction
    """
    prefix = os.path.join(PREFIX, 'direction')
    os.mkdir(prefix)

    dir_fs = autofile.fs.direction(prefix)
    print(dir_fs.leaf.path([True]))

    ref_inp_str = '<input string>'
    print(dir_fs.leaf.file.geometry_input.path([True]))
    dir_fs.leaf.create([True])
    dir_fs.leaf.file.geometry_input.write(ref_inp_str, [True])
    assert dir_fs.leaf.file.geometry_input.read([True]) == ref_inp_str


def test__species():
    """ test autofile.fs.species
    """
    prefix = os.path.join(PREFIX, 'species')
    os.mkdir(prefix)

    locs = ['InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+', 0, 1]
    spc_fs = autofile.fs.species(prefix)
    print(spc_fs.leaf.path(locs))

    assert not spc_fs.leaf.exists(locs)
    spc_fs.leaf.create(locs)
    assert spc_fs.leaf.exists(locs)


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
    print(rxn_fs.leaf.path(locs))

    assert not rxn_fs.leaf.exists(locs)
    rxn_fs.leaf.create(locs)
    assert rxn_fs.leaf.exists(locs)


def test__ts():
    """ test autofile.fs.ts
    """
    prefix = os.path.join(PREFIX, 'ts')
    os.mkdir(prefix)

    ts_fs = autofile.fs.ts(prefix)
    print(ts_fs.trunk.path())

    assert not ts_fs.trunk.exists()
    ts_fs.trunk.create()
    assert ts_fs.trunk.exists()


def test__theory():
    """ test autofile.fs.theory
    """
    prefix = os.path.join(PREFIX, 'theory')
    os.mkdir(prefix)

    thy_fs = autofile.fs.theory(prefix)
    locs = ['hf', 'sto-3g', True]
    print(thy_fs.leaf.path(locs))

    ref_ene = 5.7
    print(thy_fs.leaf.file.energy.path(locs))
    thy_fs.leaf.create(locs)
    thy_fs.leaf.file.energy.write(ref_ene, locs)
    assert thy_fs.leaf.file.energy.read(locs) == ref_ene


def test__conformer():
    """ test autofile.fs.conformer
    """
    prefix = os.path.join(PREFIX, 'conformer')
    os.mkdir(prefix)

    locs = [autofile.system.generate_new_conformer_id()]
    cnf_fs = autofile.fs.conformer(prefix)
    print(cnf_fs.leaf.path(locs))

    assert not cnf_fs.leaf.exists(locs)
    cnf_fs.leaf.create(locs)
    assert cnf_fs.leaf.exists(locs)


def test__tau():
    """ test autofile.fs.tau
    """
    prefix = os.path.join(PREFIX, 'tau')
    os.mkdir(prefix)

    locs = [autofile.system.generate_new_tau_id()]
    tau_fs = autofile.fs.tau(prefix)
    print(tau_fs.leaf.path(locs))

    assert not tau_fs.leaf.exists(locs)
    tau_fs.leaf.create(locs)
    assert tau_fs.leaf.exists(locs)


def test__single_point():
    """ test autofile.fs.single_point
    """
    prefix = os.path.join(PREFIX, 'single_point')
    os.mkdir(prefix)

    sp_fs = autofile.fs.single_point(prefix)
    locs = ['hf', 'sto-3g', True]
    print(sp_fs.leaf.path(locs))

    ref_ene = 5.7
    print(sp_fs.leaf.file.energy.path(locs))
    sp_fs.leaf.create(locs)
    sp_fs.leaf.file.energy.write(ref_ene, locs)
    assert sp_fs.leaf.file.energy.read(locs) == ref_ene


# def test__energy_transfer():
#     """ test autofile.fs.energy_transfer
#     """
#     prefix = os.path.join(PREFIX, 'etrans')
#     os.mkdir(prefix)
#
#  #   spc_fs = autofile.fs.species(prefix)
#  #   spc_locs = ['InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+', 0, 1]
#  #   spc_fs.leaf.create(spc_locs)
#  #   spc_path = spc_fs.leaf.path(spc_locs)
#  #   print(spc_fs.leaf.path(spc_locs))
#
#     etrans_fs = autofile.fs.energy_transfer(prefix)
#     spc_locs = [['InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+', 0, 1],[],[]]
#     etrans_fs.leaf.create(spc_locs)
#     etrans_spc_path = etrans_fs.leaf.path(spc_locs)
#     print(etrans_spc_path)
#     thy_locs = [ 'hf', 'sto-3g', True]
#     etrans_thy_fs = autofile.fs.theory(etrans_spc_path)
#     etrans_thy_fs.leaf.create(thy_locs)
#     etrans_thy_path = etrans_thy_fs.leaf.path(thy_locs)
#     print(etrans_thy_path)
#     spc_path = print(etrans_spc_s.leaf.path(spc_locs))
#     print(etrans_spc_s.leaf.path(spc_locs))
#     locs = ['InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+', 0, 1]
#     etrans_fs.leaf.create(locs)
#     print(etrans_fs.leaf.path(locs))
#
#     ref_eps = 300.0
#     ref_sig = 3.50
#     print(etrans_fs.leaf.file.lennard_jones_epsilon.path(locs))
#     print(etrans_fs.leaf.file.lennard_jones_sigma.path(locs))
#     # etrans_fs.leaf.create(locs)
#     etrans_fs.leaf.file.lennard_jones_epsilon.write(ref_eps, locs)
#     etrans_fs.leaf.file.lennard_jones_sigma.write(ref_sig, locs)
#     assert etrans_fs.leaf.file.lennard_jones_epsilon.read(locs) == ref_eps
#     assert etrans_fs.leaf.file.lennard_jones_sigma.read(locs) == ref_sig


def test__scan():
    """ test autofile.fs.scan
    """
    prefix = os.path.join(PREFIX, 'scan')
    os.mkdir(prefix)

    scn_fs = autofile.fs.scan(prefix)
    locs = [['d3', 'd4'], [2, 3]]
    print(scn_fs.leaf.path(locs))

    ref_inp_str = '<input string>'
    print(scn_fs.leaf.file.geometry_input.path(locs))
    scn_fs.leaf.create(locs)
    scn_fs.leaf.file.geometry_input.write(ref_inp_str, locs)
    assert scn_fs.leaf.file.geometry_input.read(locs) == ref_inp_str


def test__cscan():
    """ test autofile.fs.cscan
    """
    prefix = os.path.join(PREFIX, 'cscan')
    os.mkdir(prefix)

    scn_fs = autofile.fs.cscan(prefix)
    # the dictionary at the end specifies the constraint values
    locs = [['d3', 'd4'], [1.2, 2.9], {'r1': 1., 'a2': 2.3}]
    print(scn_fs.leaf.path(locs))

    ref_inp_str = '<input string>'
    print(scn_fs.leaf.file.geometry_input.path(locs))
    scn_fs.leaf.create(locs)
    scn_fs.leaf.file.geometry_input.write(ref_inp_str, locs)
    assert scn_fs.leaf.file.geometry_input.read(locs) == ref_inp_str


def test__run():
    """ test autofile.fs.run
    """
    prefix = os.path.join(PREFIX, 'run')
    os.mkdir(prefix)

    run_fs = autofile.fs.run(prefix)
    print(run_fs.leaf.path(['gradient']))

    ref_inp_str = '<input string>'
    print(run_fs.leaf.file.input.path(['gradient']))
    run_fs.leaf.create(['gradient'])
    run_fs.leaf.file.input.write(ref_inp_str, ['gradient'])
    assert run_fs.leaf.file.input.read(['gradient']) == ref_inp_str


def test__build():
    """ test autofile.fs.build
    """
    prefix = os.path.join(PREFIX, 'build')
    os.mkdir(prefix)

    build_fs = autofile.fs.build(prefix)
    print(build_fs.leaf.path(['MESS', 0]))

    ref_inp_str = '<input string>'
    print(build_fs.leaf.file.input.path(['MESS', 0]))
    build_fs.leaf.create(['MESS', 0])
    build_fs.leaf.file.input.write(ref_inp_str, ['MESS', 0])
    assert build_fs.leaf.file.input.read(['MESS', 0]) == ref_inp_str


if __name__ == '__main__':
    # test__direction()
    # test__species()
    # test__reaction()
    # test__ts()
    # test__theory()
    # test__conformer()
    # test__tau()
    # test__single_point()
    # test__energy_transfer()
    # test__scan()
    # test__run()
    # test__build()
    test__cscan()
