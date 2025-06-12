""" test autofile.schema.data_files
"""

import numbers
import tempfile
import numpy
import automol
import autofile.info
import autofile.schema


PREFIX = tempfile.mkdtemp()
print(PREFIX)


def test__data_files__input_file():
    """ test autofile.schema.data_files.input_file
    """
    ref_inp_str = '<input file contents>'

    inp_dfile = autofile.schema.data_files.input_file('test')

    assert not inp_dfile.exists(PREFIX)
    inp_dfile.write(ref_inp_str, PREFIX)
    assert inp_dfile.exists(PREFIX)

    inp_str = inp_dfile.read(PREFIX)
    assert inp_str == ref_inp_str


def test__data_files__output_file():
    """ test autofile.schema.data_files.output_file
    """
    ref_out_str = '<output file contents>'

    out_dfile = autofile.schema.data_files.output_file('test')

    assert not out_dfile.exists(PREFIX)
    out_dfile.write(ref_out_str, PREFIX)
    assert out_dfile.exists(PREFIX)

    out_str = out_dfile.read(PREFIX)
    assert out_str == ref_out_str


def test__data_files__information():
    """ test autofile.schema.data_files.information
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

    inf_dfile = autofile.schema.data_files.information(
        'test', function=information)

    assert not inf_dfile.exists(PREFIX)
    inf_dfile.write(ref_inf_obj, PREFIX)
    assert inf_dfile.exists(PREFIX)

    inf_obj = inf_dfile.read(PREFIX)
    assert inf_obj == ref_inf_obj


def test__data_files__energy():
    """ test autofile.schema.data_files.energy
    """
    ref_ene = -187.38518070487598

    ene_dfile = autofile.schema.data_files.energy('test')

    assert not ene_dfile.exists(PREFIX)
    ene_dfile.write(ref_ene, PREFIX)
    assert ene_dfile.exists(PREFIX)

    ene = ene_dfile.read(PREFIX)
    assert numpy.isclose(ene, ref_ene)


def test__data_files__geometry():
    """ test autofile.schema.data_files.geometry
    """
    ref_geo = (('C', (0.066541036329, -0.86543409422, -0.56994517889)),
               ('O', (0.066541036329, -0.86543409422, 2.13152981129)),
               ('O', (0.066541036329, 1.6165813318, -1.63686376233)),
               ('H', (-1.52331011945, -1.99731957213, -1.31521725797)),
               ('H', (1.84099386813, -1.76479255185, -1.16213243427)),
               ('H', (-1.61114836922, -0.17751142359, 2.6046492029)),
               ('H', (-1.61092727126, 2.32295906780, -1.19178601663)))

    geo_dfile = autofile.schema.data_files.geometry('test')

    assert not geo_dfile.exists(PREFIX)
    geo_dfile.write(ref_geo, PREFIX)
    assert geo_dfile.exists(PREFIX)

    geo = geo_dfile.read(PREFIX)
    assert automol.geom.almost_equal(geo, ref_geo)


def test__data_files__gradient():
    """ test autofile.schema.data_files.gradient
    """
    ref_grad = ((0.00004103632, 0.00003409422, 0.00004517889),
                (0.00004103632, 0.00003409422, 0.00002981129),
                (0.00004103632, 0.00008133180, 0.00006376233),
                (0.00001011945, 0.00001957213, 0.00001725797),
                (0.00009386813, 0.00009255185, 0.00003243427),
                (0.00004836922, 0.00001142359, 0.00004920290),
                (0.00002727126, 0.00005906780, 0.00008601663))

    grad_dfile = autofile.schema.data_files.gradient('test')

    assert not grad_dfile.exists(PREFIX)
    grad_dfile.write(ref_grad, PREFIX)
    assert grad_dfile.exists(PREFIX)

    grad = grad_dfile.read(PREFIX)
    assert numpy.allclose(grad, ref_grad)


def test__data_files__hessian():
    """ test autofile.schema.data_files.hessian
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

    hess_dfile = autofile.schema.data_files.hessian('test')

    assert not hess_dfile.exists(PREFIX)
    hess_dfile.write(ref_hess, PREFIX)
    assert hess_dfile.exists(PREFIX)

    hess = hess_dfile.read(PREFIX)
    assert numpy.allclose(hess, ref_hess)


def test__data_files__zmatrix():
    """ test autofile.schema.data_files.zmatrix
    """
    ref_zma = (
        ('C', (None, None, None), (None, None, None),
         (None, None, None)),
        ('O', (0, None, None), ('r1', None, None),
         (2.65933, None, None)),
        ('O', (0, 1, None), ('r2', 'a1', None),
         (2.65933, 1.90743, None)),
        ('H', (0, 1, 2), ('r3', 'a2', 'd1'),
         (2.06844, 1.93366, 4.1477)),
        ('H', (0, 1, 2), ('r4', 'a3', 'd2'),
         (2.06548, 1.89469, 2.06369)),
        ('H', (1, 0, 2), ('r5', 'a4', 'd3'),
         (1.83126, 1.86751, 1.44253)),
        ('H', (2, 0, 1), ('r6', 'a5', 'd4'),
         (1.83126, 1.86751, 4.84065)))

    zma_dfile = autofile.schema.data_files.zmatrix('test')

    assert not zma_dfile.exists(PREFIX)
    zma_dfile.write(ref_zma, PREFIX)
    assert zma_dfile.exists(PREFIX)

    zma = zma_dfile.read(PREFIX)
    print(zma)
    assert automol.zmat.almost_equal(zma, ref_zma)


def test__data_files__vmatrix():
    """ test autofile.schema.data_files.vmatrix
    """
    ref_vma = (('C', (None, None, None), (None, None, None)),
               ('O', (0, None, None), ('r1', None, None)),
               ('O', (0, 1, None), ('r2', 'a1', None)),
               ('H', (0, 1, 2), ('r3', 'a2', 'd1')),
               ('H', (0, 1, 2), ('r4', 'a3', 'd2')),
               ('H', (1, 0, 2), ('r5', 'a4', 'd3')),
               ('H', (2, 0, 1), ('r6', 'a5', 'd4')))

    vma_dfile = autofile.schema.data_files.vmatrix('test')

    assert not vma_dfile.exists(PREFIX)
    vma_dfile.write(ref_vma, PREFIX)
    assert vma_dfile.exists(PREFIX)

    vma = vma_dfile.read(PREFIX)
    assert vma == ref_vma


def test__data_files__trajectory():
    """ test autofile.schema.data_files.trajectory
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

    ref_traj = list(zip(ref_geos, ref_comments))

    traj_dfile = autofile.schema.data_files.trajectory('test')

    assert not traj_dfile.exists(PREFIX)
    traj_dfile.write(ref_traj, PREFIX)
    assert traj_dfile.exists(PREFIX)

    # I'm not going to bother implementing a reader, since the trajectory files
    # are for human use only -- we aren't going to use this for data storage


def test__data_files__lennard_jones_epsilon():
    """ test autofile.schema.data_files.lennard_jones_epsilon
    """
    ref_eps = 247.880866746988

    eps_dfile = autofile.schema.data_files.lennard_jones_epsilon('test')

    assert not eps_dfile.exists(PREFIX)
    eps_dfile.write(ref_eps, PREFIX)
    assert eps_dfile.exists(PREFIX)

    eps = eps_dfile.read(PREFIX)
    assert numpy.isclose(eps, ref_eps)


def test__data_files__lennard_jones_sigma():
    """ test autofile.schema.data_files.lennard_jones_sigma
    """
    ref_sig = 3.55018590361446

    sig_dfile = autofile.schema.data_files.lennard_jones_sigma('test')

    assert not sig_dfile.exists(PREFIX)
    sig_dfile.write(ref_sig, PREFIX)
    assert sig_dfile.exists(PREFIX)

    sig = sig_dfile.read(PREFIX)
    assert numpy.isclose(sig, ref_sig)


def test__data_files__external_symmetry_number():
    """ test autofile.schema.data_files.external_symmetry_number
    """
    ref_num = 1.55

    num_dfile = autofile.schema.data_files.external_symmetry_number('test')

    assert not num_dfile.exists(PREFIX)
    num_dfile.write(ref_num, PREFIX)
    assert num_dfile.exists(PREFIX)

    num = num_dfile.read(PREFIX)
    assert numpy.isclose(num, ref_num)


def test__data_files__internal_symmetry_number():
    """ test autofile.schema.data_files.internal_symmetry_number
    """
    ref_num = 2.00

    num_dfile = autofile.schema.data_files.internal_symmetry_number('test')

    assert not num_dfile.exists(PREFIX)
    num_dfile.write(ref_num, PREFIX)
    assert num_dfile.exists(PREFIX)

    num = num_dfile.read(PREFIX)
    assert numpy.isclose(num, ref_num)


def test__data_files__dipole_moment():
    """ test autofile.schema.data_files.dipole_moment
    """
    ref_vec = (1.00, 2.00, 3.00)

    vec_dfile = autofile.schema.data_files.dipole_moment('test')

    assert not vec_dfile.exists(PREFIX)
    vec_dfile.write(ref_vec, PREFIX)
    assert vec_dfile.exists(PREFIX)

    vec = vec_dfile.read(PREFIX)
    assert numpy.allclose(vec, ref_vec)


def test__data_files__polarizability():
    """ test autofile.schema.data_files.polarizability
    """

    ref_tensor = ((1.00, 2.00, 3.00),
                  (4.00, 5.00, 6.00),
                  (7.00, 8.00, 9.00))

    tensor_dfile = autofile.schema.data_files.polarizability('test')

    assert not tensor_dfile.exists(PREFIX)
    tensor_dfile.write(ref_tensor, PREFIX)
    assert tensor_dfile.exists(PREFIX)

    tensor = tensor_dfile.read(PREFIX)
    assert numpy.allclose(tensor, ref_tensor)
