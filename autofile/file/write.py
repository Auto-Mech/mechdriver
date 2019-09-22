""" string writers
"""
from io import StringIO as _StringIO
from numbers import Real as _Real
import numpy
import automol
import autofile.info


def information(inf_obj):
    """ write information (any dict/list combination) to a string
    """
    assert isinstance(inf_obj, autofile.info.Info)
    inf_str = autofile.info.string(inf_obj)
    return inf_str


def energy(ene):
    """ write an energy (hartree) to a string (hartree)
    """
    ene_str = _float(ene)
    return ene_str


def geometry(geo):
    """ write a geometry (bohr) to a string (angstrom)
    """
    assert automol.geom.is_valid(geo)
    xyz_str = automol.geom.xyz_string(geo)
    return xyz_str


def trajectory(traj):
    """ write a series of geometries (bohr) to a string (angstrom)

    (trajectory is given by a sequence of comment-line, geometry pairs)
    """
    comments, geo_lst = zip(*traj)
    assert all(isinstance(comment, str) and len(comment.splitlines()) == 1
               for comment in comments)
    assert all(map(automol.geom.is_valid, geo_lst))
    xyz_traj_str = automol.geom.xyz_trajectory_string(geo_lst,
                                                      comments=comments)
    return xyz_traj_str


def zmatrix(zma):
    """ write a zmatrix (bohr/radian) to a string (angstroms/degree)
    """
    assert automol.zmatrix.is_valid(zma)
    zma_str = automol.zmatrix.string(zma)
    return zma_str


def vmatrix(vma):
    """ write a variable zmatrix (bohr/radian) to a string (angstroms/degree)
    """
    assert automol.vmatrix.is_valid(vma)
    vma_str = automol.vmatrix.string(vma)
    return vma_str


def gradient(grad):
    """ write a gradient (hartree bohr^-1) to a string (hartree bohr^-1)
    """
    grad = numpy.array(grad)
    assert grad.ndim == 2 and grad.shape[1] == 3

    grad_str_io = _StringIO()
    numpy.savetxt(grad_str_io, grad)
    grad_str = grad_str_io.getvalue()
    grad_str_io.close()
    return grad_str


def hessian(hess):
    """ write a hessian (hartree bohr^-2) to a string (hartree bohr^-2)
    """
    hess = numpy.array(hess)
    assert hess.ndim == 2
    assert hess.shape[0] % 3 == 0 and hess.shape[0] == hess.shape[1]

    hess_str_io = _StringIO()
    numpy.savetxt(hess_str_io, hess)
    hess_str = hess_str_io.getvalue()
    hess_str_io.close()
    return hess_str


def harmonic_frequencies(freqs):
    """ write harmonic frequencies (cm^-1) to a string (cm^-1)
    """
    assert list(freqs) == sorted(freqs)
    return _frequencies(freqs)


def anharmonic_frequencies(freqs):
    """ write anharmonic frequencies (cm^-1) to a string (cm^-1)
    """
    assert list(freqs) == sorted(freqs)
    return _frequencies(freqs)


def projected_frequencies(freq):
    """ write projected frequencies (cm^-1) to a string (cm^-1)
    """
    assert list(freq) == sorted(freq)
    return _frequencies(freq)


def anharmonic_zpve(zpve):
    """ write the anharmonic ZPVE (hartree) to a string (hartree)
    """
    anh_zpve_str = _float(zpve)
    return anh_zpve_str


def anharmonicity_matrix(xmat):
    """ write anharmonicity matrix (cm^-1) to a string (cm^-1)
    """
    return _2d_square_matrix(xmat)


def vibro_rot_alpha_matrix(vibro_rot_mat):
    """ write vibro-rot alph matrix (cm^-1) to a string (cm^-1)
    """
    vibro_rot_mat = numpy.array(vibro_rot_mat)
    assert vibro_rot_mat.ndim == 2

    mat_str_io = _StringIO()
    numpy.savetxt(mat_str_io, vibro_rot_mat)
    mat_str = mat_str_io.getvalue()
    mat_str_io.close()
    return mat_str


def quartic_centrifugal_dist_consts(qcd_consts):
    """ write the quartic centrifugal distortion constant
        labels and values (cm^-1) to a string (cm^-1)
    """
    qcd_consts_str = ''
    for const in qcd_consts:
        qcd_consts_str += "{0:<6s}{1:>16.12f}\n".format(
            const[0], const[1])
    return qcd_consts_str


def lennard_jones_epsilon(eps):
    """ write a lennard-jones epsilon (waveunmbers) to a string (wavenumbers)
    """
    eps_str = _float(eps)
    return eps_str


def lennard_jones_sigma(sig):
    """ write a lennard-jones sigma (angstrom) to a string (angstrom)
    """
    sig_str = _float(sig)
    return sig_str


def external_symmetry_factor(esf):
    """ read an external symmetry factor from a string
    """
    esf_str = _float(esf)
    return esf_str


def _float(val):
    assert isinstance(val, _Real)
    val_str = str(val)
    return val_str


def _frequencies(freq):
    freq = numpy.array(freq)
    assert freq.ndim == 1
    freq_str = ""
    for val in freq:
        freq_str += "{:>8.1f}\n".format(val)
    return freq_str


def _2d_square_matrix(mat):
    mat = numpy.array(mat)
    assert mat.ndim == 2
    assert mat.shape[0] == mat.shape[1]

    mat_str_io = _StringIO()
    numpy.savetxt(mat_str_io, mat)
    mat_str = mat_str_io.getvalue()
    mat_str_io.close()
    return mat_str
