""" test the autofile.file module
"""
import os
import tempfile
import numpy
import automol
import autofile.info
import autofile.file

TMP_DIR = tempfile.mkdtemp()
print(TMP_DIR)


def test__information():
    """ test the information read/write functions
    """
    ref_inf_obj = autofile.info.Info(a=['b', 'c', 'd', 'e'],
                                     x=autofile.info.Info(y=1, z=2))

    inf_file_name = autofile.file.name.information('test')
    inf_file_path = os.path.join(TMP_DIR, inf_file_name)
    inf_str = autofile.file.write.information(ref_inf_obj)

    assert not os.path.isfile(inf_file_path)
    autofile.file.write_file(inf_file_path, inf_str)
    assert os.path.isfile(inf_file_path)

    inf_str = autofile.file.read_file(inf_file_path)
    inf_obj = autofile.file.read.information(inf_str)
    assert inf_obj == ref_inf_obj


def test__file():
    """ test the file read/write functions
    """
    inp_str = '<input file contents>'
    out_str = '<output file contents>'
    scr_str = '<shell script contents>'

    inp_file_name = autofile.file.name.input_file('test')
    out_file_name = autofile.file.name.output_file('test')
    scr_file_name = autofile.file.name.run_script('test')

    inp_file_path = os.path.join(TMP_DIR, inp_file_name)
    out_file_path = os.path.join(TMP_DIR, out_file_name)
    scr_file_path = os.path.join(TMP_DIR, scr_file_name)

    assert not os.path.isfile(inp_file_path)
    assert not os.path.isfile(out_file_path)
    assert not os.path.isfile(scr_file_path)
    autofile.file.write_file(inp_file_path, inp_str)
    autofile.file.write_file(out_file_path, out_str)
    autofile.file.write_file(scr_file_path, scr_str)


def test__energy():
    """ test the energy read/write functions
    """
    ref_ene = -75.00613628303537

    ene_file_name = autofile.file.name.energy('test')
    ene_file_path = os.path.join(TMP_DIR, ene_file_name)
    ene_str = autofile.file.write.energy(ref_ene)

    assert not os.path.isfile(ene_file_path)
    autofile.file.write_file(ene_file_path, ene_str)
    assert os.path.isfile(ene_file_path)

    ene_str = autofile.file.read_file(ene_file_path)
    ene = autofile.file.read.energy(ene_str)
    assert numpy.isclose(ref_ene, ene)


def test__geometry():
    """ test the geometry read/write functions
    """
    ref_geo = (('C', (-0.70116587131, 0.0146227007587, -0.016166607003)),
               ('O', (1.7323365056, -0.9538524899, -0.5617192010)),
               ('H', (-0.9827048283, 0.061897979239, 2.02901783816)),
               ('H', (-0.8787925682, 1.91673409124, -0.80019507919)),
               ('H', (-2.12093033745, -1.21447973767, -0.87411360631)),
               ('H', (2.9512589894, 0.17507745634, 0.22317665541)))

    geo_file_name = autofile.file.name.geometry('test')
    geo_file_path = os.path.join(TMP_DIR, geo_file_name)
    geo_str = autofile.file.write.geometry(ref_geo)

    assert not os.path.isfile(geo_file_path)
    autofile.file.write_file(geo_file_path, geo_str)
    assert os.path.isfile(geo_file_path)

    geo_str = autofile.file.read_file(geo_file_path)
    geo = autofile.file.read.geometry(geo_str)
    assert automol.geom.almost_equal(ref_geo, geo)


def test__trajectory():
    """ test the trajectory read/rwrite functions
    """
    ref_geo_lst = [
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
    ref_comments = [
        'energy: -187.38941054878092',
        'energy: -187.3850624381528']
    traj = tuple(zip(ref_comments, ref_geo_lst))

    traj_file_name = autofile.file.name.trajectory('test')
    traj_file_path = os.path.join(TMP_DIR, traj_file_name)
    traj_str = autofile.file.write.trajectory(traj)

    assert not os.path.isfile(traj_file_path)
    autofile.file.write_file(traj_file_path, traj_str)
    assert os.path.isfile(traj_file_path)


def test__zmatrix():
    """ test the zmatrix read/write functions
    """
    ref_zma = (
        (('C', (None, None, None), (None, None, None)),
         ('O', (0, None, None), ('r1', None, None)),
         ('H', (0, 1, None), ('r2', 'a1', None)),
         ('H', (0, 1, 2), ('r3', 'a2', 'd1')),
         ('H', (0, 1, 2), ('r4', 'a3', 'd2')),
         ('H', (1, 0, 2), ('r5', 'a4', 'd3'))),
        {'r1': 2.67535,
         'r2': 2.06501, 'a1': 1.9116242,
         'r3': 2.06501, 'a2': 1.9116242, 'd1': 2.108497362,
         'r4': 2.06458, 'a3': 1.9020947, 'd2': 4.195841334,
         'r5': 1.83748, 'a4': 1.8690905, 'd3': 5.228936625})

    zma_file_name = autofile.file.name.zmatrix('test')
    zma_file_path = os.path.join(TMP_DIR, zma_file_name)
    zma_str = autofile.file.write.zmatrix(ref_zma)

    assert not os.path.isfile(zma_file_path)
    autofile.file.write_file(zma_file_path, zma_str)
    assert os.path.isfile(zma_file_path)

    zma_str = autofile.file.read_file(zma_file_path)
    zma = autofile.file.read.zmatrix(zma_str)
    assert automol.zmatrix.almost_equal(ref_zma, zma)


def test__vmatrix():
    """ test the vmatrix read/write functions
    """
    ref_vma = (('C', (None, None, None), (None, None, None)),
               ('O', (0, None, None), ('r1', None, None)),
               ('H', (0, 1, None), ('r2', 'a1', None)),
               ('H', (0, 1, 2), ('r3', 'a2', 'd1')),
               ('H', (0, 1, 2), ('r4', 'a3', 'd2')),
               ('H', (1, 0, 2), ('r5', 'a4', 'd3')))

    vma_file_name = autofile.file.name.vmatrix('test')
    vma_file_path = os.path.join(TMP_DIR, vma_file_name)
    vma_str = autofile.file.write.vmatrix(ref_vma)

    assert not os.path.isfile(vma_file_path)
    autofile.file.write_file(vma_file_path, vma_str)
    assert os.path.isfile(vma_file_path)

    vma_str = autofile.file.read_file(vma_file_path)
    vma = autofile.file.read.vmatrix(vma_str)
    assert vma == ref_vma


def test__gradient():
    """ test the gradient read/write functions
    """
    ref_grad = (
        (-6.687065494e-05, -2.087360475e-05, -2.900518464e-05),
        (2.091815312e-05, -2.162539565e-05, 2.097302071e-05),
        (-9.17712798e-06, -2.97043797e-05, -1.378883952e-05),
        (2.244858355e-05, 6.0265176e-06, 1.00733398e-06),
        (4.676018254e-05, 5.473937899e-05, 1.744549e-05),
        (-1.407913705e-05, 1.143748381e-05, 3.3681804e-06))

    grad_file_name = autofile.file.name.gradient('test')
    grad_file_path = os.path.join(TMP_DIR, grad_file_name)
    grad_str = autofile.file.write.gradient(ref_grad)

    assert not os.path.isfile(grad_file_path)
    autofile.file.write_file(grad_file_path, grad_str)
    assert os.path.isfile(grad_file_path)

    grad_str = autofile.file.read_file(grad_file_path)
    grad = autofile.file.read.gradient(grad_str)
    assert numpy.allclose(ref_grad, grad)


def test__hessian():
    """ test the hessian read/write functions
    """
    ref_hess = (
        (0.70455411073336, -0.08287472697212, 0.0, -0.0798647434566,
         0.05385932364624, 0.0, -0.62468936727659, 0.02901540332588, 0.0),
        (-0.08287472697212, 0.73374800936533, 0.0, -0.06694564915575,
         -0.63928635896955, 0.0, 0.14982037612785, -0.09446165039567, 0.0),
        (0.0, 0.0, 0.00019470176975, 0.0, 0.0, -9.734161519e-05, 0.0, 0.0,
         -9.736015438e-05),
        (-0.0798647434566, -0.06694564915575, 0.0, 0.08499297871932,
         -0.02022391947876, 0.0, -0.00512823526272, 0.0871695686345, 0.0),
        (0.05385932364624, -0.63928635896955, 0.0, -0.02022391947876,
         0.65384367318451, 0.0, -0.03363540416748, -0.01455731421496, 0.0),
        (0.0, 0.0, -9.734161519e-05, 0.0, 0.0, 4.370597581e-05, 0.0, 0.0,
         5.363563938e-05),
        (-0.62468936727659, 0.14982037612785, 0.0, -0.00512823526272,
         -0.03363540416748, 0.0, 0.6298176025393, -0.11618497196038, 0.0),
        (0.02901540332588, -0.09446165039567, 0.0, 0.0871695686345,
         -0.01455731421496, 0.0, -0.11618497196038, 0.10901896461063, 0.0),
        (0.0, 0.0, -9.736015438e-05, 0.0, 0.0, 5.363563938e-05, 0.0, 0.0,
         4.3724515e-05)
    )

    hess_file_name = autofile.file.name.hessian('test')
    hess_file_path = os.path.join(TMP_DIR, hess_file_name)
    hess_str = autofile.file.write.hessian(ref_hess)

    assert not os.path.isfile(hess_file_path)
    autofile.file.write_file(hess_file_path, hess_str)
    assert os.path.isfile(hess_file_path)

    hess_str = autofile.file.read_file(hess_file_path)
    hess = autofile.file.read.hessian(hess_str)
    assert numpy.allclose(ref_hess, hess)


def test__lennard_jones_epsilon():
    """ test the epsilon read/write functions
    """
    ref_eps = 247.880866746988

    eps_file_name = autofile.file.name.lennard_jones_epsilon('test')
    eps_file_path = os.path.join(TMP_DIR, eps_file_name)
    eps_str = autofile.file.write.lennard_jones_epsilon(ref_eps)

    assert not os.path.isfile(eps_file_path)
    autofile.file.write_file(eps_file_path, eps_str)
    assert os.path.isfile(eps_file_path)

    eps_str = autofile.file.read_file(eps_file_path)
    eps = autofile.file.read.lennard_jones_epsilon(eps_str)
    assert numpy.isclose(ref_eps, eps)


def test__lennard_jones_sigma():
    """ test the sigma read/write functions
    """
    ref_sig = 3.55018590361446

    sig_file_name = autofile.file.name.lennard_jones_sigma('test')
    sig_file_path = os.path.join(TMP_DIR, sig_file_name)
    sig_str = autofile.file.write.lennard_jones_sigma(ref_sig)

    assert not os.path.isfile(sig_file_path)
    autofile.file.write_file(sig_file_path, sig_str)
    assert os.path.isfile(sig_file_path)

    sig_str = autofile.file.read_file(sig_file_path)
    sig = autofile.file.read.lennard_jones_sigma(sig_str)
    assert numpy.isclose(ref_sig, sig)


def test__anharmonicity_matrix():
    """ test the anharmonicity matrix
    """

    ref_xmat = ((1.0, 1.0), (1.0, 1.0))

    xmat_file_name = autofile.file.name.anharmonicity_matrix('test')
    xmat_file_path = os.path.join(TMP_DIR, xmat_file_name)
    xmat_str = autofile.file.write.anharmonicity_matrix(ref_xmat)

    print(xmat_str)

if __name__ == '__main__':
    # test__energy()
    # test__geometry()
    # test__zmatrix()
    # test__lennard_jones_epsilon()
    # test__lennard_jones_sigma()
    # test__file()
    # test__information()
    # test__gradient()
    # test__hessian()
    # test__trajectory()
    # test__vmatrix()
    test__anharmonicity_matrix()
