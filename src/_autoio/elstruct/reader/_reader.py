""" Electronic structure program output reading module.

    Calls functions from the various program modules. Each module must provide
    a function that matches one in the module template --
    both the function name and signature are checked
    before calling the function.

    The resulting function signatures are exactly those in pm.Job.py
    with `prog` inserted as the first argument.
"""

import numpy
import automol
from elstruct.reader import program_modules as pm


def programs():
    """ Constructs a list of available electronic structure programs.
        At minimum, each program must have an energy reader to be enumerated.
    """
    return pm.program_modules_with_function(pm.Job.ENERGY)


# Molecular Info
def energy(prog, method, output_str):
    """ Reads the electronic energy from the output string.
        Returns the energy in Hartrees.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param method: electronic structure method
        :type method: str
        :param output_str: string of the program's output file
        :type output_str: str
    """

    ene = pm.call_module_function(
        prog, pm.Job.ENERGY,
        # *args
        method, output_str)
    if ene is not None:
        assert isinstance(ene, float)

    return ene


def gradient_programs():
    """ Constructs a list of program modules implementing
        gradient output readers.
    """
    return pm.program_modules_with_function(pm.Job.GRADIENT)


def gradient(prog, output_str):
    """ Reads the gradient from the output string.
        Returns the gradient in atomic units.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param output_str: string of the program's output file
        :type output_str: str
    """

    grad = pm.call_module_function(
        prog, pm.Job.GRADIENT,
        # *args
        output_str)

    # if grad is not None:
    #     assert all(isinstance(val, float) for val in grad)
    #     assert numpy.shape(grad)[1] == 3

    return grad


def hessian_programs():
    """ Constructs a list of program modules implementing
        Hessian output readers.
    """
    return pm.program_modules_with_function(pm.Job.HESSIAN)


def hessian(prog, output_str):
    """ Reads the Hessian from the output string.
        Returns the Hessian in atomic units.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param output_str: string of the program's output file
        :type output_str: str
    """

    hess = pm.call_module_function(
        prog, pm.Job.HESSIAN,
        # *args
        output_str)

    if hess is not None:
        assert numpy.allclose(hess, numpy.transpose(hess))

    return hess


def harmonic_frequencies_programs():
    """ Constructs a list of program modules implementing
        harmonic vibrarional frequency output readers.
    """
    return pm.program_modules_with_function(pm.Job.HARM_FREQS)


def harmonic_frequencies(prog, output_str):
    """ Reads the harmonic vibrational frequencies from the output string.
        Returns the frequencies in cm-1.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param output_str: string of the program's output file
        :type output_str: str
    """

    harm_freqs = pm.call_module_function(
        prog, pm.Job.HARM_FREQS,
        # *args
        output_str)

    if harm_freqs is not None:
        assert all(isinstance(val, float) for val in harm_freqs)
        # assert numpy.shape(harm_freqs)[1] == 3

    return harm_freqs


def normal_coordinates_programs():
    """ Constructs a list of program modules implementing
        normal coordinate output readers.
    """
    return pm.program_modules_with_function(pm.Job.NORM_COORDS)


def normal_coordinates(prog, output_str):
    """ Reads the displacement along the normal modes from the output string.
        Returns the coordinates in Bohr.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param output_str: string of the program's output file
        :type output_str: str
    """
    return pm.call_module_function(
        prog, pm.Job.NORM_COORDS,
        # *args
        output_str)


def irc_programs():
    """ Constructs a list of program modules implementing
        Intrinsic Reaction Coordinate output readers.
    """
    return pm.program_modules_with_function(pm.Job.IRC_PTS)


def irc_points(prog, output_str):
    """ Reads the geometries, gradients, and Hessians at each point along the
        Intrinsic Reaction Coordinate from the output string.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param output_str: string of the program's output file
        :type output_str: str
    """
    return pm.call_module_function(
        prog, pm.Job.IRC_PTS,
        # *args
        output_str)


def irc_path(prog, output_str):
    """ read irc_path from the output string

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param output_str: string of the program's output file
        :type output_str: str
    """
    return pm.call_module_function(
        prog, pm.Job.IRC_PATH,
        # *args
        output_str)


def opt_geometry_programs():
    """ Constructs a list of program modules implementing
        optimized geometry output readers.
    """
    return pm.program_modules_with_function(pm.Job.OPT_GEO)


def opt_geometry(prog, output_str):
    """ Reads the optimized geometry from the output string.
        Returns the geometry in Bohr.

        For robustness: if geometry read fails, try reading
        the z-matrix and converting.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param output_str: string of the program's output file
        :type output_str: str
    """
    try:
        geo = _opt_geometry(prog, output_str)
        geo = automol.geom.without_dummy_atoms(geo)
    except TypeError:
        zma = opt_zmatrix(prog, output_str)
        geo = automol.zmat.geometry(zma)
    return geo


def _opt_geometry(prog, output_str):
    """ Reads the optimized geometry from the output string.
        Returns the geometry in Bohr.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param output_str: string of the program's output file
        :type output_str: str
    """
    return pm.call_module_function(
        prog, pm.Job.OPT_GEO,
        # *args
        output_str)


def opt_zmatrix_programs():
    """ Contucts a list of program modules implementing
        optimized Z-Matrix output readers.
    """
    return pm.program_modules_with_function(pm.Job.OPT_ZMA)


def opt_zmatrices_programs():
    """ Contucts a list of program modules implementing
        optimized Z-Matrix output readers.
    """
    return pm.program_modules_with_function(pm.Job.OPT_ZMAS)


def opt_zmatrix(prog, output_str):
    """ Reads the optimized Z-Matrix from the output string.
        Returns the geometry in Bohr+Radians.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param output_str: string of the program's output file
        :type output_str: str
    """
    return pm.call_module_function(
        prog, pm.Job.OPT_ZMA,
        # *args
        output_str)


def opt_zmatrices(prog, output_str):
    """ Reads the optimized Z-Matrix from the output string.
        Returns the geometry in Bohr+Radians.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param output_str: string of the program's output file
        :type output_str: str
    """
    return pm.call_module_function(
        prog, pm.Job.OPT_ZMAS,
        # *args
        output_str)


def inp_zmatrix_programs():
    """ Contucts a list of program modules implementing
        optimized Z-Matrix output readers.
    """
    return pm.program_modules_with_function(pm.Job.INP_ZMA)


def inp_zmatrix(prog, input_str):
    """ Reads the optimized Z-Matrix from the output string.
        Returns the geometry in Bohr+Radians.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param output_str: string of the program's output file
        :type output_str: str
    """
    return pm.call_module_function(
        prog, pm.Job.INP_ZMA,
        # *args
        input_str)


def vpt2_programs():
    """ Constructs a list of program modules implementing
        2nd-order vibrational perturbation theory (VPT2) output readers.
    """
    return pm.program_modules_with_function(pm.Job.VPT2)


def vpt2(prog, output_str):
    """ Reads VPT2 information from the output string.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param output_str: string of the program's output file
        :type output_str: str
    """
    return pm.call_module_function(
        prog, pm.Job.VPT2,
        # *args
        output_str)


def dipole_moment_programs():
    """ Constructs a list of program modules implementing
        static dipole moment output readers.
    """
    return pm.program_modules_with_function(pm.Job.DIP_MOM)


def dipole_moment(prog, output_str):
    """ Reads the static dipole moment from the output string.
        Returns the dipole moment in Debye.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param output_str: string of the program's output file
        :type output_str: str
    """
    return pm.call_module_function(
        prog, pm.Job.DIP_MOM,
        # *args
        output_str)


def polarizability_programs():
    """ Constructs a list of program modules implementing
        polarizability tensor output readers.
    """
    return pm.program_modules_with_function(pm.Job.POLAR)


def polarizability(prog, output_str):
    """ Reads the polarizability tensor from the output string.
        Returns the dipole moment in _.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param output_str: string of the program's output file
        :type output_str: str
    """
    return pm.call_module_function(
        prog, pm.Job.POLAR,
        # *args
        output_str)


# Status
def has_normal_exit_message(prog, output_str):
    """ Assess whether the output file string contains the
        normal program exit message.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param output_str: string of the program's output file
        :type output_str: str
    """
    return pm.call_module_function(
        prog, pm.Job.EXIT_MSG,
        # *args
        output_str)


def error_list(prog):
    """ Constructs a list of errors that be identified from the output file.

        :param prog: the electronic structure program to use as a backend
        :type prog: str
    """
    return pm.call_module_function(
        prog, pm.Job.ERR_LST)


def success_list(prog):
    """ Constructs a list of successes that be identified from the output file.

        :param prog: the electronic structure program to use as a backend
        :type prog: str
    """
    return pm.call_module_function(
        prog, pm.Job.SUCCESS_LST)


def has_error_message(prog, error, output_str):
    """ Assess whether the output file string contains error messages
        for any of the procedures in the job.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param error: a key indicating the type of error message
        :type error: str
        :param output_str: string of the program's output file
        :type output_str: str
    """
    return pm.call_module_function(
        prog, pm.Job.ERR_MSG,
        # *args
        error, output_str)


def check_convergence_messages(prog, error, success, output_str):
    """ Assess whether the output file string contains messages
        denoting all of the requested procedures in the job have converged.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param error: a key indicating the type of error message
        :type error: str
        :param success: a key indicating the type of success message
        :type success: str
        :param output_str: string of the program's output file
        :type output_str: str
    """
    return pm.call_module_function(
        prog, pm.Job.CONV_MSG,
        # *args
        error, success, output_str)


# Versions
def program_name(prog, output_str):
    """ Reads the name of the electronic structure code from the output file.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param output_str: string of the program's output file
        :type output_str: str
    """
    return pm.call_module_function(
        prog, pm.Job.PROG_NAME,
        # *args
        output_str)


def program_version(prog, output_str):
    """ Reads the version of the electronic structure code from the output file.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param output_str: string of the program's output file
        :type output_str: str
    """
    return pm.call_module_function(
        prog, pm.Job.PROG_VERS,
        # *args
        output_str)
