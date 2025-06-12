""" string writers
    converts properties in internally used units and used data type
    to properties formatted in strings of externally preferred units dum
"""

import os
from io import StringIO as _StringIO
from numbers import Real as _Real
import numpy
import yaml
import automol
from phydat import phycon
import autofile.info


def information(inf_obj):
    """ write information (any dict/list combination) to a string

    :param inf: info yaml information
    :type inf: Info
    :return: info object as a string
    :rtype: str
    """
    assert isinstance(inf_obj, autofile.info.Info)
    inf_str = autofile.info.string(inf_obj)
    return inf_str


def instability(instab_rxn):
    """ write a reaction to a string
    :param instab_rxn: an automol Reaction object
    :type instab_rxn: automol.reac.Reaction
    :return: reaction string
    :rtype: str
    """
    rxn_str = automol.reac.string(instab_rxn)
    return rxn_str


def energy(ene):
    """ write an energy (hartree) to a string (hartree)

    :param ene: energy
    :type inf: float
    :return: energy
    :rtype: str
    """
    ene_str = _float(ene)
    return ene_str


def geometry(geo):
    """ write a geometry (bohr) to a string (angstrom)

    :param geo: geometry in autochem tuple format
    :type geo: tuple
    :return: goemetry as xyz format
    :rtype: str
    """
    assert automol.geom.is_valid(geo)
    xyz_str = automol.geom.xyz_string(geo)
    return xyz_str


def trajectory(traj):
    """ write a series of geometries (bohr) to a string (angstrom)


    (trajectory is given by a sequence of comment-line, geometry pairs)
    :param traj: traj object
    :type traj: tuple
    :return: trajectory
    :rtype: str
    """
    geo_lst, comments = zip(*traj)

    assert all(isinstance(comment, str) and len(comment.splitlines()) == 1
               for comment in comments)
    assert all(map(automol.geom.is_valid, geo_lst))
    xyz_traj_str = automol.geom.xyz_trajectory_string(geo_lst,
                                                      comments=comments)
    return xyz_traj_str


def zmatrix(zma):
    """ write a zmatrix (bohr/radian) to a string (angstroms/degree)

    :param zma: zmatrix in autochem tuple format
    :type zma: tuple
    :return: zmatrix as string
    :rtype: str
    """
    assert automol.zmat.is_valid(zma), (
        f'invalid zmatrix: {zma}'
    )
    zma_str = automol.zmat.string(zma)
    return zma_str


def vmatrix(vma):
    """ write a variable zmatrix (bohr/radian) to a string (angstroms/degree)

    :param vma: vmatrix in autochem tuple format
    :type vma: tuple
    :return: vmatrix string
    :rtype: str
    """
    assert automol.vmat.is_valid(vma)
    vma_str = automol.vmat.string(vma)
    return vma_str


def ring_torsions(ring_tors_dct):
    """ Write the torsions and their ranges (radian) to a string (degree).


        :param tors_lst: list of torsion objects
        :type tors_lst: tuple(automol torsion objects)
        :rtype: str
    """

    ring_tors_str = []
    for ring, tors_dct in ring_tors_dct.items():

        assert all(isinstance(key, str) and len(rng) == 2
                   and all(isinstance(x, _Real) for x in rng)
                   for key, rng in tors_dct.items())

        tors_names = list(tors_dct.keys())
        tors_names.sort(key=lambda x: int(x.split('D')[1]))
        sorted_dct = {}
        for name in tors_names:
            sorted_dct[name] = [
                tors_dct[name][0]*phycon.RAD2DEG,
                tors_dct[name][1]*phycon.RAD2DEG
            ]

        tors_str = yaml.dump(
            sorted_dct, default_flow_style=True, sort_keys=False)

        ring_tors_str.append(f'ring: {ring}\n{tors_str}')

    return os.linesep.join(ring_tors_str)


def torsions(tors_lst):
    """ Write the torsions and their ranges (radian) to a string (degree).


        :param tors_lst: list of torsion objects
        :type tors_lst: tuple(automol torsion objects)
        :rtype: str
    """
    tors_str = automol.data.tors.torsions_string(tors_lst)
    return tors_str


def gradient(grad):
    """ write a gradient (hartree bohr^-1) to a string (hartree bohr^-1)

    :param grad: gradient tuple
    :type grad: tuple
    :return: gradient string
    :rtype: str
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

    :param hess: hessian 3nx3n tuple of floats
    :type hess: tuple
    :return: hessian string
    :rtype: str
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

    :param freqs: freq tuple of floats
    :type freqs: tuple
    :return: frequencies as string
    :rtype: str
    """
    assert list(freqs) == sorted(freqs)
    return _frequencies(freqs)


def anharmonic_frequencies(freqs):
    """ write anharmonic frequencies (cm^-1) to a string (cm^-1)

    :param freqs: freq tuple of floats
    :type freqs: tuple
    :return: frequencies as string
    :rtype: str
    """
    assert list(freqs) == sorted(freqs)
    return _frequencies(freqs)


def cubic_force_constants(cfcs):
    """ Writes cubic force constants () to a string ().


        :param cfcs: cubic force constants
        :type cfcs: numpy.ndarray
        :rtype: str
    """
    return automol.util.tensor.string(cfcs, val_format='{0:>14.6f}')


def quartic_force_constants(qfcs):
    """ Writes quartic force constants () to a string ().


        :param cfcs: quartic force constants
        :type cfcs: numpy.ndarray
        :rtype: str
    """
    return automol.util.tensor.string(qfcs, val_format='{0:>14.6f}')


def harmonic_zpve(zpve):
    """ write the harmonic ZPVE (hartree) to a string (hartree)

    :param harm_zpve: zpve float
    :type harm_zpve: float
    :return: zpve as string
    :rtype: str
    """
    harm_zpve_str = _float(zpve)
    return harm_zpve_str


def anharmonic_zpve(zpve):
    """ write the anharmonic ZPVE (hartree) to a string (hartree)

    :param anh_zpve: zpve float
    :type anh_zpve: float
    :return: zpve as string
    :rtype: str
    """
    anh_zpve_str = _float(zpve)
    return anh_zpve_str


def anharmonicity_matrix(xmat):
    """ write anharmonicity matrix (cm^-1) to a string (cm^-1)

    :param xmat: anharmonicity xmatrix as nfreqxnfreq tuple
    :type xmat: tuple
    :return: xmat string
    :rtype: str
    """
    mat = numpy.array(xmat)
    assert mat.ndim in (0, 2)
    if mat.ndim == 2:
        assert mat.shape[0] == mat.shape[1]

    mat_str_io = _StringIO()
    numpy.savetxt(mat_str_io, mat)
    mat_str = mat_str_io.getvalue()
    mat_str_io.close()
    return mat_str


def vibro_rot_alpha_matrix(vibro_rot_mat):
    """ write vibro-rot alph matrix (cm^-1) to a string (cm^-1)

    :param vibro_rot: matrix as tuple
    :type vibro_rot: tuple
    :return: vibro-rot alpha matrix string
    :rtype: str
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
    :param qcd_consts: constants in a tuple
    :type qcd_consts: tuple
    :return: quartic centrifugal dist const string
    :rtype: str
    """
    qcd_consts_str = ''
    for const in qcd_consts:
        qcd_consts_str += f'{const[0]:<6s}{const[1]:>16.12f}\n'
    return qcd_consts_str


def lennard_jones_epsilon(eps):
    """ write a lennard-jones epsilon (waveunmbers) to a string (wavenumbers)

    :param eps: epsilon float
    :type eps_consts: float
    :return: epsilon string
    :rtype: str
    """
    eps_str = _float(eps)
    return eps_str


def lennard_jones_sigma(sig):
    """ write a lennard-jones sigma (angstrom) to a string (angstrom)

    :param sig: sigma float
    :type sig_consts: float
    :return: sigma string
    :rtype: str
    """
    sig_str = _float(sig)
    return sig_str


def external_symmetry_factor(esf):
    """ write an external symmetry factor to a string

    :param esf: external symmetry factor float
    :type esf_consts: float
    :return: external symmetry factor string
    :rtype: str
    """
    esf_str = _float(esf)
    return esf_str


def internal_symmetry_factor(isf):
    """ write an internal symmetry factor to a string

    :param isf: internal symmetry factor float
    :type isf_consts: float
    :return: internal symmetry factor string
    :rtype: str
    """
    isf_str = _float(isf)
    return isf_str


def dipole_moment(dip_mom):
    """ write a dipole moment vector to a string

    :param dip_mom: x,y,z dipole moment tuple
    :type dip_mom: tuple
    :return: x, y, z dipole moment vector string
    :rtype: str
    """
    dip_mom = numpy.array(dip_mom)
    assert dip_mom.ndim == 1
    assert dip_mom.shape[0] == 3

    dip_mom_str_io = _StringIO()
    numpy.savetxt(dip_mom_str_io, dip_mom)
    dip_mom_str = dip_mom_str_io.getvalue()
    dip_mom_str_io.close()
    return dip_mom_str


def polarizability(polar):
    """ write a polarizability tensor to a string

    :param polar: polarizability tensor
    :type polar: tuple
    :return: polarizability tensor
    :rtype: str
    """
    polar = numpy.array(polar)
    assert polar.ndim == 2
    assert polar.shape[0] == polar.shape[1] == 3

    polar_str_io = _StringIO()
    numpy.savetxt(polar_str_io, polar)
    polar_str = polar_str_io.getvalue()
    polar_str_io.close()
    return polar_str


def reaction(rxn):
    """ write a reaction to a string

    :param rxn: an automol Reaction object
    :type rxn: automol.reac.Reaction
    :return: reaction string
    :rtype: str
    """
    rxn_str = automol.reac.string(rxn)
    return rxn_str


def gradient_array(ndarray):
    """ transform numpy array into python list
    """
    out = ndarray
    if isinstance(out, numpy.ndarray):
        out = ndarray.tolist()
    return out


def _float(val):
    """ float to string
    """
    assert isinstance(val, _Real)
    val_str = str(val)
    return val_str


def _frequencies(freq):
    """ tuple of floats to string
    """
    freq = numpy.array(freq)
    assert freq.ndim == 1
    freq_str = ""
    for val in freq:
        freq_str += f'{val:>8.1f}\n'
    return freq_str
