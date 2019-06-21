""" test fs
"""
import os
import tempfile
import numpy
import automol
import autofile.system
from autofile import SFS
from autofile import RFS

PREFIX = tempfile.mkdtemp()
print(PREFIX)


def test__species():
    """ tets fsys.species
    """
    prefix = os.path.join(PREFIX, 'species')
    os.mkdir(prefix)

    specs_lst = [
        ('InChI=1S/O', 0, 3),
        ('InChI=1S/O', 0, 1),
        ('InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+', 0, 1),
        ('InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1-', 0, 1),
        ('InChI=1S/C5H5O/c1-2-3-4-5-6/h1-5H/b4-3-', 0, 2),
    ]

    for specs in specs_lst:
        assert not SFS.species.dir.exists(prefix, specs)
        SFS.species.dir.create(prefix, specs)
        assert SFS.species.dir.exists(prefix, specs)

    for specs in SFS.species.dir.existing(prefix):
        print(specs)


def test__reaction():
    """ tets fsys.reaction
    """
    prefix = os.path.join(PREFIX, 'reaction')
    os.mkdir(prefix)

    specs_lst = [
        (
            (('InChI=1S/C2H5O2/c1-2-4-3/h3H,1-2H2',),
             ('InChI=1S/C2H4/c1-2/h1-2H2', 'InChI=1S/HO2/c1-2/h1H')),
            ((0,), (0, 0)),
            ((2,), (1, 2)),
        ),
        (
            (('InChI=1S/CH3/h1H3', 'InChI=1S/HO/h1H'),
             ('InChI=1S/CH2/h1H2', 'InChI=1S/H2O/h1H2')),
            ((0, 0), (0, 0)),
            ((2, 2), (1, 1)),
        ),
        (
            (('InChI=1S/CH3O3/c2-1-4-3/h3H,1H2',),
             ('InChI=1S/CH3O3/c2-1-4-3/h2H,1H2',)),
            ((0,), (0,)),
            ((2,), (2,)),
        ),
    ]

    for specs in specs_lst:
        assert not RFS.reaction.dir.exists(prefix, specs)
        RFS.reaction.dir.create(prefix, specs)
        assert RFS.reaction.dir.exists(prefix, specs)

    for specs in RFS.reaction.dir.existing(prefix):
        print(specs)


def test__theory():
    """ tets fsys.theory
    """
    prefix = os.path.join(PREFIX, 'theory')
    os.mkdir(prefix)

    root_specs = ('InChI=1S/CH3/h1H3', 0, 2)
    specs_lst = [
        root_specs + ('hf', 'sto-3g', True),
        root_specs + ('hf', 'sto-3g', False),
        root_specs + ('b3lyp', 'sto-3g', False),
        root_specs + ('b3lyp', '6-31g*', False),
    ]

    for specs in specs_lst:
        assert not SFS.theory.dir.exists(prefix, specs)
        SFS.theory.dir.create(prefix, specs)
        assert SFS.theory.dir.exists(prefix, specs)

    for specs in SFS.theory.dir.existing(prefix, root_specs):
        print(specs)


def test__conformer():
    """ tets fsys.conformer
    """
    prefix = os.path.join(PREFIX, 'conformer')
    os.mkdir(prefix)

    nconfs = 10

    # generate IDs for each conformer
    cids = (autofile.system.generate_new_conformer_id()
            for _ in range(nconfs))

    root_specs = (
        'InChI=1S/C5H5O/c1-2-3-4-5-6/h1-5H/b4-3-', 0, 2,
        'hf', 'sto-3g', False)
    specs_lst = [root_specs + (cid,) for cid in cids]

    assert not SFS.conf_trunk.dir.exists(prefix, root_specs)
    SFS.conf_trunk.dir.create(prefix, root_specs)
    assert SFS.conf_trunk.dir.exists(prefix, root_specs)

    # create the trunk vmatrix file
    ref_trunk_vma = (('C', (None, None, None), (None, None, None)),
                     ('O', (0, None, None), ('r1', None, None)),
                     ('O', (0, 1, None), ('r2', 'a1', None)),
                     ('H', (0, 1, 2), ('r3', 'a2', 'd1')),
                     ('H', (0, 1, 2), ('r4', 'a3', 'd2')),
                     ('H', (1, 0, 2), ('r5', 'a4', 'd3')),
                     ('H', (2, 0, 1), ('r6', 'a5', 'd4')))
    SFS.conf_trunk.file.vmatrix.write(ref_trunk_vma, prefix, root_specs)
    trunk_vma = SFS.conf_trunk.file.vmatrix.read(prefix, root_specs)
    assert trunk_vma == ref_trunk_vma

    # create the trunk information file
    ref_trunk_inf_obj = autofile.system.info.conformer_trunk(
        nsamp=7, tors_ranges={'d3': (0, 6.283185307179586),
                              'd4': (0, 6.283185307179586)})
    SFS.conf_trunk.file.info.write(ref_trunk_inf_obj, prefix, root_specs)
    trunk_inf_obj = SFS.conf_trunk.file.info.read(prefix, root_specs)
    assert trunk_inf_obj == ref_trunk_inf_obj

    for specs in specs_lst:
        SFS.conf.dir.create(prefix, specs)

        ref_geom_inf_obj = autofile.system.info.run(
            job='optimization', prog='psi4', method='mp2', basis='sto-3g',
            status="succeeded")
        ref_grad_inf_obj = autofile.system.info.run(
            job='gradient', prog='psi4', method='mp2', basis='sto-3g',
            status="succeeded")
        ref_hess_inf_obj = autofile.system.info.run(
            job='hessian', prog='psi4', method='mp2', basis='sto-3g',
            status="succeeded")
        ref_geom_inp_str = '<geometry input file>'
        ref_grad_inp_str = '<gradient input file>'
        ref_hess_inp_str = '<hessian input file>'
        ref_ene = -187.38518070487598
        ref_geo = (('C', (0.066541036329, -0.86543409422, -0.56994517889)),
                   ('O', (0.066541036329, -0.86543409422, 2.13152981129)),
                   ('O', (0.066541036329, 1.6165813318, -1.63686376233)),
                   ('H', (-1.52331011945, -1.99731957213, -1.31521725797)),
                   ('H', (1.84099386813, -1.76479255185, -1.16213243427)),
                   ('H', (-1.61114836922, -0.17751142359, 2.6046492029)),
                   ('H', (-1.61092727126, 2.32295906780, -1.19178601663)))
        ref_grad = ((0.00004103632, 0.00003409422, 0.00004517889),
                    (0.00004103632, 0.00003409422, 0.00002981129),
                    (0.00004103632, 0.00008133180, 0.00006376233),
                    (0.00001011945, 0.00001957213, 0.00001725797),
                    (0.00009386813, 0.00009255185, 0.00003243427),
                    (0.00004836922, 0.00001142359, 0.00004920290),
                    (0.00002727126, 0.00005906780, 0.00008601663))
        # (I'm not bothering with the hessian for now)

        # writes
        SFS.conf.file.geometry_info.write(
            ref_geom_inf_obj, prefix, specs)
        SFS.conf.file.gradient_info.write(
            ref_grad_inf_obj, prefix, specs)
        SFS.conf.file.hessian_info.write(
            ref_hess_inf_obj, prefix, specs)
        SFS.conf.file.geometry_input.write(ref_geom_inp_str, prefix, specs)
        SFS.conf.file.gradient_input.write(ref_grad_inp_str, prefix, specs)
        SFS.conf.file.hessian_input.write(ref_hess_inp_str, prefix, specs)
        SFS.conf.file.energy.write(ref_ene, prefix, specs)
        SFS.conf.file.geometry.write(ref_geo, prefix, specs)
        SFS.conf.file.gradient.write(ref_grad, prefix, specs)

        # reads
        geom_inf_obj = SFS.conf.file.geometry_info.read(prefix, specs)
        grad_inf_obj = SFS.conf.file.gradient_info.read(prefix, specs)
        hess_inf_obj = SFS.conf.file.hessian_info.read(prefix, specs)
        geom_inp_str = SFS.conf.file.geometry_input.read(prefix, specs)
        grad_inp_str = SFS.conf.file.gradient_input.read(prefix, specs)
        hess_inp_str = SFS.conf.file.hessian_input.read(prefix, specs)
        ene = SFS.conf.file.energy.read(prefix, specs)
        geo = SFS.conf.file.geometry.read(prefix, specs)
        grad = SFS.conf.file.gradient.read(prefix, specs)

        # check read values
        assert geom_inf_obj == ref_geom_inf_obj
        assert grad_inf_obj == ref_grad_inf_obj
        assert hess_inf_obj == ref_hess_inf_obj
        assert geom_inp_str == ref_geom_inp_str
        assert grad_inp_str == ref_grad_inp_str
        assert hess_inp_str == ref_hess_inp_str
        assert numpy.isclose(ene, ref_ene)
        assert automol.geom.almost_equal(geo, ref_geo)
        assert numpy.allclose(grad, ref_grad)

    print(SFS.conf.dir.existing(prefix, root_specs))
    assert len(SFS.conf.dir.existing(prefix, root_specs)) == nconfs


def test__conformer_run():
    """ tets fsys.conformer_run
    """
    prefix = os.path.join(PREFIX, 'conformer_run')
    os.mkdir(prefix)

    cid = autofile.system.generate_new_conformer_id()
    root_specs = (
        'InChI=1S/C5H5O/c1-2-3-4-5-6/h1-5H/b4-3-', 0, 2,
        'hf', 'sto-3g', False, cid)
    specs_lst = (
        root_specs + ('energy',),
        root_specs + ('gradient',),
        root_specs + ('hessian',),
        root_specs + ('optimization',),
    )

    for specs in specs_lst:
        SFS.conf_run.dir.create(prefix, specs)
        job = specs[-1]
        ref_inf_obj = autofile.system.info.run(
            job=job, prog='psi4', method='hf', basis='sto-3g',
            status="succeeded")
        ref_inp_str = '<input file>'
        ref_out_str = '<output file>'

        SFS.conf_run.file.info.write(ref_inf_obj, prefix, specs)
        SFS.conf_run.file.input.write(ref_inp_str, prefix, specs)
        SFS.conf_run.file.output.write(ref_out_str, prefix, specs)


def test__single_point():
    """ tets fsys.single_point
    """
    prefix = os.path.join(PREFIX, 'single_point')
    os.mkdir(prefix)

    root_specs = (
        'InChI=1S/C5H5O/c1-2-3-4-5-6/h1-5H/b4-3-', 0, 2,
        'hf', 'sto-3g', False, 'HS2I1PDureBE')

    specs_lst = [
        root_specs + ('hf', 'sto-3g', True),
        root_specs + ('hf', 'sto-3g', False),
        root_specs + ('b3lyp', 'sto-3g', False),
        root_specs + ('b3lyp', '6-31g*', False),
    ]

    for specs in specs_lst:
        ref_ene = -187.38518070487598

        SFS.conf_sp.dir.create(prefix, specs)
        SFS.conf_sp.file.energy.write(ref_ene, prefix, specs)
        ene = SFS.conf_sp.file.energy.read(prefix, specs)
        assert numpy.isclose(ene, ref_ene)

    for specs in SFS.conf_sp.dir.existing(prefix, root_specs):
        print(specs)


def test__scan():
    """ tets fsys.scan
    """
    prefix = os.path.join(PREFIX, 'scan')
    os.mkdir(prefix)

    root_specs = (
        'InChI=1S/C5H5O/c1-2-3-4-5-6/h1-5H/b4-3-', 0, 2,
        'hf', 'sto-3g', False, 'HS2I1PDureBE')
    branch_specs = root_specs + (['d4', 'd8'],)

    SFS.scan_branch.dir.create(prefix, branch_specs)

    # create the trunk vmatrix file
    ref_trunk_vma = (('C', (None, None, None), (None, None, None)),
                     ('X', (0, None, None), ('r1', None, None)),
                     ('C', (0, 1, None), ('r2', 'a1', None)),
                     ('H', (0, 1, 2), ('r3', 'a2', 'd1')),
                     ('C', (2, 0, 1), ('r4', 'a3', 'd2')),
                     ('H', (2, 0, 4), ('r5', 'a4', 'd3')),
                     ('C', (4, 2, 0), ('r6', 'a5', 'd4')),
                     ('H', (4, 2, 6), ('r7', 'a6', 'd5')),
                     ('C', (6, 4, 2), ('r8', 'a7', 'd6')),
                     ('H', (6, 4, 8), ('r9', 'a8', 'd7')),
                     ('O', (8, 6, 4), ('r10', 'a9', 'd8')),
                     ('H', (8, 6, 10), ('r11', 'a10', 'd9')))
    SFS.scan_trunk.file.vmatrix.write(ref_trunk_vma, prefix, root_specs)
    trunk_vma = SFS.scan_trunk.file.vmatrix.read(prefix, root_specs)
    assert trunk_vma == ref_trunk_vma

    # create the branch information file
    ref_scan_inf_obj = autofile.system.info.scan_branch(
        tors_linspaces={'d3': (0, 6.283185307179586, 3),
                        'd4': (0, 6.283185307179586, 2)})
    SFS.scan_branch.file.info.write(ref_scan_inf_obj, prefix, branch_specs)
    scan_inf_obj = SFS.scan_branch.file.info.read(prefix, branch_specs)
    assert scan_inf_obj == ref_scan_inf_obj

    print(SFS.scan_branch.dir.existing(prefix, root_specs))

    main_specs_lst = (
        ([0, 0],),
        ([1, 0],),
        ([2, 0],),
        ([0, 1],),
        ([1, 1],),
        ([2, 1],),
    )
    specs_lst = tuple(branch_specs + main_specs
                      for main_specs in main_specs_lst)

    for specs in specs_lst:
        SFS.scan.dir.create(prefix, specs)

        ref_geom_inf_obj = autofile.system.info.run(
            job='optimization', prog='psi4', method='mp2', basis='sto-3g',
            status="succeeded")
        ref_geom_inp_str = '<geometry input file>'
        ref_ene = -187.38518070487598
        ref_geo = (('C', (0.066541036329, -0.86543409422, -0.56994517889)),
                   ('O', (0.066541036329, -0.86543409422, 2.13152981129)),
                   ('O', (0.066541036329, 1.6165813318, -1.63686376233)),
                   ('H', (-1.52331011945, -1.99731957213, -1.31521725797)),
                   ('H', (1.84099386813, -1.76479255185, -1.16213243427)),
                   ('H', (-1.61114836922, -0.17751142359, 2.6046492029)),
                   ('H', (-1.61092727126, 2.32295906780, -1.19178601663)))

        # writes
        SFS.scan.file.geometry_info.write(ref_geom_inf_obj, prefix, specs)
        SFS.scan.file.geometry_input.write(ref_geom_inp_str, prefix, specs)
        SFS.scan.file.energy.write(ref_ene, prefix, specs)
        SFS.scan.file.geometry.write(ref_geo, prefix, specs)

        # reads
        geom_inf_obj = SFS.scan.file.geometry_info.read(prefix, specs)
        geom_inp_str = SFS.scan.file.geometry_input.read(prefix, specs)
        ene = SFS.scan.file.energy.read(prefix, specs)
        geo = SFS.scan.file.geometry.read(prefix, specs)

        # check read values
        assert geom_inf_obj == ref_geom_inf_obj
        assert geom_inp_str == ref_geom_inp_str
        assert numpy.isclose(ene, ref_ene)
        assert automol.geom.almost_equal(geo, ref_geo)

    print(SFS.scan.dir.existing(prefix, branch_specs))


def test__scan_run():
    """ tets fsys.scan_run
    """
    prefix = os.path.join(PREFIX, 'scan_run')
    os.mkdir(prefix)

    cid = autofile.system.generate_new_conformer_id()
    root_specs = (
        'InChI=1S/C5H5O/c1-2-3-4-5-6/h1-5H/b4-3-', 0, 2,
        'hf', 'sto-3g', False, cid, ['d4', 'd8'], [0, 0])
    specs_lst = (
        root_specs + ('energy',),
        root_specs + ('gradient',),
        root_specs + ('hessian',),
        root_specs + ('optimization',),
    )

    for specs in specs_lst:
        SFS.scan_run.dir.create(prefix, specs)
        job = specs[-1]
        ref_inf_obj = autofile.system.info.run(
            job=job, prog='psi4', method='hf', basis='sto-3g',
            status="succeeded")
        ref_inp_str = '<input file>'
        ref_out_str = '<output file>'

        SFS.scan_run.file.info.write(ref_inf_obj, prefix, specs)
        SFS.scan_run.file.input.write(ref_inp_str, prefix, specs)
        SFS.scan_run.file.output.write(ref_out_str, prefix, specs)


def test__tau():
    """ test fsys.tau
    """
    prefix = os.path.join(PREFIX, 'tau')
    os.mkdir(prefix)

    ntaus = 10

    # generate IDs for each tau
    cids = (autofile.system.generate_new_conformer_id()
            for _ in range(ntaus))

    root_specs = (
        'InChI=1S/C5H5O/c1-2-3-4-5-6/h1-5H/b4-3-', 0, 2,
        'hf', 'sto-3g', False)
    specs_lst = [root_specs + (cid,) for cid in cids]

    assert not SFS.tau_trunk.dir.exists(prefix, root_specs)
    SFS.tau_trunk.dir.create(prefix, root_specs)
    assert SFS.tau_trunk.dir.exists(prefix, root_specs)

    # create the trunk vmatrix file
    ref_trunk_vma = (('C', (None, None, None), (None, None, None)),
                     ('O', (0, None, None), ('r1', None, None)),
                     ('O', (0, 1, None), ('r2', 'a1', None)),
                     ('H', (0, 1, 2), ('r3', 'a2', 'd1')),
                     ('H', (0, 1, 2), ('r4', 'a3', 'd2')),
                     ('H', (1, 0, 2), ('r5', 'a4', 'd3')),
                     ('H', (2, 0, 1), ('r6', 'a5', 'd4')))
    SFS.tau_trunk.file.vmatrix.write(ref_trunk_vma, prefix, root_specs)
    trunk_vma = SFS.tau_trunk.file.vmatrix.read(prefix, root_specs)
    assert trunk_vma == ref_trunk_vma

    # create the trunk information file
    ref_trunk_inf_obj = autofile.system.info.tau_trunk(
        nsamp=7, tors_ranges={'d3': (0, 6.283185307179586),
                              'd4': (0, 6.283185307179586)})
    SFS.tau_trunk.file.info.write(ref_trunk_inf_obj, prefix, root_specs)
    trunk_inf_obj = SFS.tau_trunk.file.info.read(prefix, root_specs)
    assert trunk_inf_obj == ref_trunk_inf_obj

    for specs in specs_lst:
        SFS.tau.dir.create(prefix, specs)

        ref_geom_inf_obj = autofile.system.info.run(
            job='optimization', prog='psi4', method='mp2', basis='sto-3g',
            status="succeeded")
        ref_grad_inf_obj = autofile.system.info.run(
            job='gradient', prog='psi4', method='mp2', basis='sto-3g',
            status="succeeded")
        ref_hess_inf_obj = autofile.system.info.run(
            job='hessian', prog='psi4', method='mp2', basis='sto-3g',
            status="succeeded")
        ref_geom_inp_str = '<geometry input file>'
        ref_grad_inp_str = '<gradient input file>'
        ref_hess_inp_str = '<hessian input file>'
        ref_ene = -187.38518070487598
        ref_geo = (('C', (0.066541036329, -0.86543409422, -0.56994517889)),
                   ('O', (0.066541036329, -0.86543409422, 2.13152981129)),
                   ('O', (0.066541036329, 1.6165813318, -1.63686376233)),
                   ('H', (-1.52331011945, -1.99731957213, -1.31521725797)),
                   ('H', (1.84099386813, -1.76479255185, -1.16213243427)),
                   ('H', (-1.61114836922, -0.17751142359, 2.6046492029)),
                   ('H', (-1.61092727126, 2.32295906780, -1.19178601663)))
        ref_grad = ((0.00004103632, 0.00003409422, 0.00004517889),
                    (0.00004103632, 0.00003409422, 0.00002981129),
                    (0.00004103632, 0.00008133180, 0.00006376233),
                    (0.00001011945, 0.00001957213, 0.00001725797),
                    (0.00009386813, 0.00009255185, 0.00003243427),
                    (0.00004836922, 0.00001142359, 0.00004920290),
                    (0.00002727126, 0.00005906780, 0.00008601663))
        # (I'm not bothering with the hessian for now)

        # writes
        SFS.tau.file.geometry_info.write(
            ref_geom_inf_obj, prefix, specs)
        SFS.tau.file.gradient_info.write(
            ref_grad_inf_obj, prefix, specs)
        SFS.tau.file.hessian_info.write(
            ref_hess_inf_obj, prefix, specs)
        SFS.tau.file.geometry_input.write(ref_geom_inp_str, prefix, specs)
        SFS.tau.file.gradient_input.write(ref_grad_inp_str, prefix, specs)
        SFS.tau.file.hessian_input.write(ref_hess_inp_str, prefix, specs)
        SFS.tau.file.energy.write(ref_ene, prefix, specs)
        SFS.tau.file.geometry.write(ref_geo, prefix, specs)
        SFS.tau.file.gradient.write(ref_grad, prefix, specs)

        # reads
        geom_inf_obj = SFS.tau.file.geometry_info.read(prefix, specs)
        grad_inf_obj = SFS.tau.file.gradient_info.read(prefix, specs)
        hess_inf_obj = SFS.tau.file.hessian_info.read(prefix, specs)
        geom_inp_str = SFS.tau.file.geometry_input.read(prefix, specs)
        grad_inp_str = SFS.tau.file.gradient_input.read(prefix, specs)
        hess_inp_str = SFS.tau.file.hessian_input.read(prefix, specs)
        ene = SFS.tau.file.energy.read(prefix, specs)
        geo = SFS.tau.file.geometry.read(prefix, specs)
        grad = SFS.tau.file.gradient.read(prefix, specs)

        # check read values
        assert geom_inf_obj == ref_geom_inf_obj
        assert grad_inf_obj == ref_grad_inf_obj
        assert hess_inf_obj == ref_hess_inf_obj
        assert geom_inp_str == ref_geom_inp_str
        assert grad_inp_str == ref_grad_inp_str
        assert hess_inp_str == ref_hess_inp_str
        assert numpy.isclose(ene, ref_ene)
        assert automol.geom.almost_equal(geo, ref_geo)
        assert numpy.allclose(grad, ref_grad)

    print(SFS.tau.dir.existing(prefix, root_specs))
    assert len(SFS.tau.dir.existing(prefix, root_specs)) == ntaus


def test__tau_run():
    """ tets fsys.tau_run
    """
    prefix = os.path.join(PREFIX, 'tau_run')
    os.mkdir(prefix)

    cid = autofile.system.generate_new_conformer_id()
    root_specs = (
        'InChI=1S/C5H5O/c1-2-3-4-5-6/h1-5H/b4-3-', 0, 2,
        'hf', 'sto-3g', False, cid)
    specs_lst = (
        root_specs + ('energy',),
        root_specs + ('gradient',),
        root_specs + ('hessian',),
        root_specs + ('optimization',),
    )

    for specs in specs_lst:
        SFS.tau_run.dir.create(prefix, specs)
        job = specs[-1]
        ref_inf_obj = autofile.system.info.run(
            job=job, prog='psi4', method='hf', basis='sto-3g',
            status="succeeded")
        ref_inp_str = '<input file>'
        ref_out_str = '<output file>'

        SFS.tau_run.file.info.write(ref_inf_obj, prefix, specs)
        SFS.tau_run.file.input.write(ref_inp_str, prefix, specs)
        SFS.tau_run.file.output.write(ref_out_str, prefix, specs)


if __name__ == '__main__':
    # test__species()
    # test__theory()
    # test__conformer()
    # test__scan()
    # test__conformer_run()
    # test__scan_run()
    # test__tau()
    # test__tau_run()
    # test__single_point()
    test__reaction()
