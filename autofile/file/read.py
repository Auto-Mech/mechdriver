""" string readers
"""
from io import StringIO as _StringIO
import numpy
import automol
import autoparse.find as apf
import autofile.info


def information(inf_str):
    """ read information (any dict/list combination) from a string
    """
    inf_obj = autofile.info.from_string(inf_str)
    return inf_obj


def energy(ene_str):
    """ read an energy (hartree) from a string (hartree)
    """
    ene = _float(ene_str)
    return ene


def geometry(xyz_str):
    """ read a geometry (bohr) from a string (angstrom)
    """
    geo = automol.geom.from_xyz_string(xyz_str)
    return geo


def zmatrix(zma_str):
    """ read a zmatrix (bohr/radian) from a string (angstrom/degree)
    """
    zma = automol.zmatrix.from_string(zma_str)
    return zma


def vmatrix(vma_str):
    """ read a variable zmatrix (bohr/radian) from a string (angstrom/degree)
    """
    vma = automol.vmatrix.from_string(vma_str)
    return vma


def gradient(grad_str):
    """ read a gradient (hartree bohr^-1) from a string (hartree bohr^-1)
    """
    grad_str_io = _StringIO(grad_str)
    grad = numpy.loadtxt(grad_str_io)
    assert grad.ndim == 2 and grad.shape[1] == 3
    return tuple(map(tuple, grad))


def hessian(hess_str):
    """ read a hessian (hartree bohr^-2) from a string (hartree bohr^-2)
    """
    hess_str_io = _StringIO(hess_str)
    hess = numpy.loadtxt(hess_str_io)
    assert hess.ndim == 2
    assert hess.shape[0] % 3 == 0 and hess.shape[0] == hess.shape[1]
    return tuple(map(tuple, hess))


def harmonic_frequencies(freq_str):
    """ read harmonic frequencies (cm^-1) from a string (cm^-1)
    """
    return _frequencies(freq_str)


def anharmonic_frequencies(freq_str):
    """ read anharmonic frequencies (cm^-1) from a string (cm^-1)
    """
    return _frequencies(freq_str)


def projected_frequencies(freq_str):
    """ read projected frequencies (cm^-1) from a string (cm^-1)
    """
    return _frequencies(freq_str)


def anharmonic_zpve(anh_zpve_str):
    """ read the anharmonic zpve (hartree) from a string (hartree)
    """
    anh_zpve = _float(anh_zpve_str)
    return anh_zpve


def anharmonicity_matrix(xmat_str):
    """ read an anharmonicity matrix (cm^-1)
        from a string (cm^-1)
    """
    return _2d_square_matrix(xmat_str)


def vibro_rot_alpha_matrix(vibro_rot_str):
    """ read an vibro-rot alpha matrix (cm^-1)
        from a string (cm^-1)
    """
    return _2d_square_matrix(vibro_rot_str)


def quartic_centrifugal_dist_consts(qcd_consts_str):
    """ write the quartic centrifugal distortion constant
        labels and values (cm^-1) to a string (cm^-1)
    """
    qcd_consts_lines = qcd_consts_str.splitlines()
    qcd_consts = []
    for line in qcd_consts_lines:
        const = line.strip().split()
        qcd_consts.append([const[0], float(const[1])])
    qcd_consts = tuple(tuple(x) for x in qcd_consts)
    return qcd_consts


def lennard_jones_epsilon(eps_str):
    """ read a lennard-jones epsilon (waveunmbers) from a string (wavenumbers)
    """
    eps = _float(eps_str)
    return eps


def lennard_jones_sigma(sig_str):
    """ read a lennard-jones sigma (angstrom) from a string (angstrom)
    """
    sig = _float(sig_str)
    return sig


def external_symmetry_factor(esf_str):
    """ read an external symmetry factor from a string
    """
    esf = _float(esf_str)
    return esf


def _float(val_str):
    assert apf.is_number(val_str)
    val = float(val_str)
    return val


def _frequencies(freq_str):
    freq_str_io = _StringIO(freq_str)
    freqs = numpy.loadtxt(freq_str_io)
    assert freqs.ndim == 1
    return tuple(freqs)


def _2d_square_matrix(mat_str):
    mat_str_io = _StringIO(mat_str)
    mat = numpy.loadtxt(mat_str_io)
    assert mat.ndim == 2
    assert mat.shape[0] == mat.shape[1]
    return tuple(map(tuple, mat))
