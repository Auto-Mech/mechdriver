"""
 tests reading of projrot output
"""

import os
import numpy
from ioformat import pathtools
import projrot_io.reader


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')


def test__frequencies():
    """ test projrot_io.reader.rpht_output
    """

    ref_real1 = (1000.0, 2000.0, 3000.0, 4000.0, 5000.0)

    out_str1 = pathtools.read_file(DAT_PATH, 'min.out')
    real1, imag1 = projrot_io.reader.rpht_output(out_str1)
    assert numpy.allclose(real1, ref_real1)
    assert not imag1

    ref_real2 = (2000.0, 3000.0, 4000.0, 5000.0)
    ref_imag2 = (1111.11,)

    out_str2 = pathtools.read_file(DAT_PATH, 'one_imag.out')
    real2, imag2 = projrot_io.reader.rpht_output(out_str2)
    assert numpy.allclose(real2, ref_real2)
    assert numpy.allclose(imag2, ref_imag2)

    ref_real3 = (3000.0, 4000.0, 5000.0)
    ref_imag3 = (1111.11, 2222.22)

    out_str3 = pathtools.read_file(DAT_PATH, 'two_imag.out')
    real3, imag3 = projrot_io.reader.rpht_output(out_str3)
    assert numpy.allclose(real3, ref_real3)
    assert numpy.allclose(imag3, ref_imag3)


def test__normal_coordinates():
    """ test projrot_io.reader.normal_coordinates
    """

    out_str1 = pathtools.read_file(DAT_PATH, 'disp1.out')
    out_str2 = pathtools.read_file(DAT_PATH, 'disp2.out')
    out_str3 = pathtools.read_file(DAT_PATH, 'disp3.out')

    nc1 = projrot_io.reader.normal_coordinates(out_str1)
    nc2 = projrot_io.reader.normal_coordinates(out_str2)
    nc3 = projrot_io.reader.normal_coordinates(out_str3)

    ref_nc1 = (
        numpy.array([[0.32125344,  0.,  0.64250688],
                     [-0.3401507, -0., -0.18897261],
                     [0.5858151, -0.,  1.5117809],
                     [-0.05669178,  0., -0.56691784]]),
        numpy.array([[-0.,  0.17007535, -0.],
                     [0., -0.18897261, -0.],
                     [-0., -1.87082886, -0.],
                     [0.,  0.09448631,  0.]]),
        numpy.array([[0.18897261,  0., -0.54802058],
                     [-0.30235618, -0.,  0.94486306],
                     [0.62360962, -0.,  1.32280829],
                     [0.03779452,  0., -0.35904796]]),
        numpy.array([[-0.17007535,  0., -0.11338357],
                     [-0.,  0., -0.17007535],
                     [0.37794523,  0.,  1.81413708],
                     [0.15117809, -0.,  0.15117809]])
    )
    ref_nc2 = (
        numpy.array([[0.32125344,  0.,  0.64250688],
                     [-0.3401507, -0., -0.18897261],
                     [0.5858151, -0.,  1.5117809],
                     [-0.05669178,  0., -0.56691784]]),
        numpy.array([[-0.,  0.17007535, -0.],
                     [0., -0.18897261, -0.],
                     [-0., -1.87082886, -0.],
                     [0.,  0.09448631,  0.]]),
        numpy.array([[0.18897261,  0., -0.54802058],
                     [-0.30235618, -0.,  0.94486306],
                     [0.62360962, -0.,  1.32280829],
                     [0.03779452,  0., -0.35904796]]),
        numpy.array([[-0.17007535,  0., -0.11338357],
                     [-0.,  0., -0.17007535],
                     [0.37794523,  0.,  1.81413708],
                     [0.15117809, -0.,  0.15117809]]),
        numpy.array([[-0.13228083,  0.,  0.05669178],
                     [0.98265759, -0.,  0.15117809],
                     [0.2456644, -0.,  1.36060281],
                     [-0.75589045,  0., -0.26456166]])
    )
    ref_nc3 = (
        numpy.array([[0.32125344,  0.,  0.64250688],
                     [-0.3401507, -0., -0.18897261],
                     [0.5858151, -0.,  1.5117809],
                     [-0.05669178,  0., -0.56691784]]),
        numpy.array([[-0.,  0.17007535, -0.],
                     [0., -0.18897261, -0.],
                     [-0., -1.87082886, -0.],
                     [0.,  0.09448631,  0.]]),
        numpy.array([[0.18897261,  0., -0.54802058],
                     [-0.30235618, -0.,  0.94486306],
                     [0.62360962, -0.,  1.32280829],
                     [0.03779452,  0., -0.35904796]]),
        numpy.array([[-0.17007535,  0., -0.11338357],
                     [-0.,  0., -0.17007535],
                     [0.37794523,  0.,  1.81413708],
                     [0.15117809, -0.,  0.15117809]]),
        numpy.array([[-0.13228083,  0.,  0.05669178],
                     [0.98265759, -0.,  0.15117809],
                     [0.2456644, -0.,  1.36060281],
                     [-0.75589045,  0., -0.26456166]]),
        numpy.array([[0.11338357,  0., -0.01889726],
                     [0., -0., -0.],
                     [-1.81413708, -0.,  0.51022605],
                     [0., -0.,  0.]])
    )

    out_str1 = pathtools.read_file(DAT_PATH, 'disp1.out')
    out_str2 = pathtools.read_file(DAT_PATH, 'disp2.out')
    out_str3 = pathtools.read_file(DAT_PATH, 'disp3.out')

    nc1 = projrot_io.reader.normal_coordinates(out_str1)
    nc2 = projrot_io.reader.normal_coordinates(out_str2)
    nc3 = projrot_io.reader.normal_coordinates(out_str3)

    assert numpy.allclose(ref_nc1, nc1)
    assert numpy.allclose(ref_nc2, nc2)
    assert numpy.allclose(ref_nc3, nc3)
