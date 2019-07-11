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


def harmonic_frequencies(freq):
    """ write harmonic frequencies (cm^-1) to a string (cm^-1)
    """
    assert list(freq) == sorted(freq)
    return _frequencies(freq)


def anharmonic_frequencies(freq):
    """ write anharmonic frequencies (cm^-1) to a string (cm^-1)
    """
    assert list(freq) == sorted(freq)
    return _frequencies(freq)


def anharmonicity_matrix(xmat):
    """ write anharmonicity matrix (cm^-1) to a string (cm^-1)
    """
    xmat = numpy.array(xmat)
    assert xmat.ndim == 2
    assert xmat.shape[0] == xmat.shape[1]

    xmat_str_io = _StringIO()
    numpy.savetxt(xmat_str_io, xmat)
    xmat_str = xmat_str_io.getvalue()
    xmat_str_io.close()
    return xmat_str


def projected_frequencies(freq):
    """ write projected frequencies (cm^-1) to a string (cm^-1)
    """
    assert list(freq) == sorted(freq)
    return _frequencies(freq)


def _frequencies(freq):
    freq = numpy.array(freq)
    assert freq.ndim == 1
    freq_str = ""
    for val in freq:
        freq_str += "{:>8.1f}\n".format(val)
    return freq_str


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
