""" test fs
"""
import os
import tempfile
import numpy
import automol
import autofile.system
from autofile import fs

PREFIX = tempfile.mkdtemp()
print(PREFIX)


def test__species():
    """ tets fsys.species
    """
    prefix = os.path.join(PREFIX, 'species')
    os.mkdir(prefix)

    specs_lst = [
        ('InChI=1S/O', 3),
        ('InChI=1S/O', 1),
        ('InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+', 1),
        ('InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1-', 1),
        ('InChI=1S/C5H5O/c1-2-3-4-5-6/h1-5H/b4-3-', 2),
    ]

    for specs in specs_lst:
        assert not fs.species.dir.exists(prefix, specs)
        fs.species.dir.create(prefix, specs)
        assert fs.species.dir.exists(prefix, specs)

    for specs in fs.species.dir.existing(prefix):
        print(specs)


def test__theory():
    """ tets fsys.theory
    """
    prefix = os.path.join(PREFIX, 'theory')
    os.mkdir(prefix)

    root_specs = ('InChI=1S/CH3/h1H3', 2)
    specs_lst = [
        root_specs + ('hf', 'sto-3g', True),
        root_specs + ('hf', 'sto-3g', False),
        root_specs + ('b3lyp', 'sto-3g', False),
        root_specs + ('b3lyp', '6-31g*', False),
    ]

    for specs in specs_lst:
        assert not fs.theory.dir.exists(prefix, specs)
        fs.theory.dir.create(prefix, specs)
        assert fs.theory.dir.exists(prefix, specs)

    for specs in fs.theory.dir.existing(prefix, root_specs):
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
        'InChI=1S/C5H5O/c1-2-3-4-5-6/h1-5H/b4-3-', 2,
        'hf', 'sto-3g', False)
    specs_lst = [root_specs + (cid,) for cid in cids]

    assert not fs.conformer_trunk.dir.exists(prefix, root_specs)
    fs.conformer_trunk.dir.create(prefix, root_specs)
    assert fs.conformer_trunk.dir.exists(prefix, root_specs)

    # create the trunk vmatrix file
    ref_trunk_vma = (('C', (None, None, None), (None, None, None)),
                     ('O', (0, None, None), ('r1', None, None)),
                     ('O', (0, 1, None), ('r2', 'a1', None)),
                     ('H', (0, 1, 2), ('r3', 'a2', 'd1')),
                     ('H', (0, 1, 2), ('r4', 'a3', 'd2')),
                     ('H', (1, 0, 2), ('r5', 'a4', 'd3')),
                     ('H', (2, 0, 1), ('r6', 'a5', 'd4')))
    fs.conformer_trunk.file.vmatrix.write(ref_trunk_vma, prefix, root_specs)
    trunk_vma = fs.conformer_trunk.file.vmatrix.read(prefix, root_specs)
    assert trunk_vma == ref_trunk_vma

    # create the trunk information file
    ref_trunk_inf_obj = autofile.system.info.torsion_sampling(
        nsamp=7, tors_ranges={'d3': (0, 6.283185307179586),
                              'd4': (0, 6.283185307179586)})
    fs.conformer_trunk.file.info.write(ref_trunk_inf_obj, prefix, root_specs)
    trunk_inf_obj = fs.conformer_trunk.file.info.read(prefix, root_specs)
    assert trunk_inf_obj == ref_trunk_inf_obj

    for specs in specs_lst:
        fs.conformer.dir.create(prefix, specs)

        ref_geom_inf_obj = autofile.system.info.run(
            job='optimization', prog='psi4', method='mp2', basis='sto-3g')
        ref_grad_inf_obj = autofile.system.info.run(
            job='gradient', prog='psi4', method='mp2', basis='sto-3g')
        ref_hess_inf_obj = autofile.system.info.run(
            job='hessian', prog='psi4', method='mp2', basis='sto-3g')
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
        fs.conformer.file.geometry_information.write(
            ref_geom_inf_obj, prefix, specs)
        fs.conformer.file.gradient_information.write(
            ref_grad_inf_obj, prefix, specs)
        fs.conformer.file.hessian_information.write(
            ref_hess_inf_obj, prefix, specs)
        fs.conformer.file.geometry_input.write(ref_geom_inp_str, prefix, specs)
        fs.conformer.file.gradient_input.write(ref_grad_inp_str, prefix, specs)
        fs.conformer.file.hessian_input.write(ref_hess_inp_str, prefix, specs)
        fs.conformer.file.energy.write(ref_ene, prefix, specs)
        fs.conformer.file.geometry.write(ref_geo, prefix, specs)
        fs.conformer.file.gradient.write(ref_grad, prefix, specs)

        # reads
        geom_inf_obj = fs.conformer.file.geometry_information.read(prefix,
                                                                   specs)
        grad_inf_obj = fs.conformer.file.gradient_information.read(prefix,
                                                                   specs)
        hess_inf_obj = fs.conformer.file.hessian_information.read(prefix,
                                                                  specs)
        geom_inp_str = fs.conformer.file.geometry_input.read(prefix, specs)
        grad_inp_str = fs.conformer.file.gradient_input.read(prefix, specs)
        hess_inp_str = fs.conformer.file.hessian_input.read(prefix, specs)
        ene = fs.conformer.file.energy.read(prefix, specs)
        geo = fs.conformer.file.geometry.read(prefix, specs)
        grad = fs.conformer.file.gradient.read(prefix, specs)

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

    print(fs.conformer.dir.existing(prefix, root_specs))
    assert len(fs.conformer.dir.existing(prefix, root_specs)) == nconfs


def test__scan():
    """ tets fsys.scan
    """
    prefix = os.path.join(PREFIX, 'scan')
    os.mkdir(prefix)

    root_specs = (
        'InChI=1S/C5H5O/c1-2-3-4-5-6/h1-5H/b4-3-', 2,
        'hf', 'sto-3g', False, 'HS2I1PDureBE')
    branch_specs = root_specs + (['d4', 'd8'],)

    fs.scan_branch.dir.create(prefix, branch_specs)

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
    fs.scan_trunk.file.vmatrix.write(ref_trunk_vma, prefix, root_specs)
    trunk_vma = fs.scan_trunk.file.vmatrix.read(prefix, root_specs)
    assert trunk_vma == ref_trunk_vma

    # create the branch information file
    ref_scan_inf_obj = autofile.system.info.scan(
        tors_linspaces={'d3': (0, 6.283185307179586, 3),
                        'd4': (0, 6.283185307179586, 2)})
    fs.scan_branch.file.info.write(ref_scan_inf_obj, prefix, branch_specs)
    scan_inf_obj = fs.scan_branch.file.info.read(prefix, branch_specs)
    assert scan_inf_obj == ref_scan_inf_obj

    print(fs.scan_branch.dir.existing(prefix, root_specs))

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
        fs.scan.dir.create(prefix, specs)

        ref_geom_inf_obj = autofile.system.info.run(
            job='optimization', prog='psi4', method='mp2', basis='sto-3g')
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
        fs.scan.file.geometry_information.write(
            ref_geom_inf_obj, prefix, specs)
        fs.scan.file.geometry_input.write(ref_geom_inp_str, prefix, specs)
        fs.scan.file.energy.write(ref_ene, prefix, specs)
        fs.scan.file.geometry.write(ref_geo, prefix, specs)

        # reads
        geom_inf_obj = fs.scan.file.geometry_information.read(prefix, specs)
        geom_inp_str = fs.scan.file.geometry_input.read(prefix, specs)
        ene = fs.scan.file.energy.read(prefix, specs)
        geo = fs.scan.file.geometry.read(prefix, specs)

        # check read values
        assert geom_inf_obj == ref_geom_inf_obj
        assert geom_inp_str == ref_geom_inp_str
        assert numpy.isclose(ene, ref_ene)
        assert automol.geom.almost_equal(geo, ref_geo)

    print(fs.scan.dir.existing(prefix, branch_specs))


if __name__ == '__main__':
    # test__species()
    # test__theory()
    # test__conformer()
    test__scan()
