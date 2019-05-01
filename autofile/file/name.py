""" file namers
"""


class Extension():
    """ file extensions """
    INFORMATION = '.yaml'
    INPUT_LOG = '.inp'
    OUTPUT_LOG = '.out'
    SHELL_SCRIPT = '.sh'
    ENERGY = '.ene'
    GEOMETRY = '.xyz'
    TRAJECTORY = '.t.xyz'
    ZMATRIX = '.zmat'
    VMATRIX = '.vmat'
    GRADIENT = '.grad'
    HESSIAN = '.hess'
    LJ_EPSILON = '.eps'
    LJ_SIGMA = '.sig'


def information(file_name):
    """ adds information extension, if missing
    """
    return _add_extension(file_name, Extension.INFORMATION)


def input_file(file_name):
    """ adds input file extension, if missing
    """
    return _add_extension(file_name, Extension.INPUT_LOG)


def output_file(file_name):
    """ adds output file extension, if missing
    """
    return _add_extension(file_name, Extension.OUTPUT_LOG)


def run_script(file_name):
    """ adds run script extension, if missing
    """
    return _add_extension(file_name, Extension.SHELL_SCRIPT)


def energy(file_name):
    """ adds energy extension, if missing
    """
    return _add_extension(file_name, Extension.ENERGY)


def geometry(file_name):
    """ adds geometry extension, if missing
    """
    return _add_extension(file_name, Extension.GEOMETRY)


def trajectory(file_name):
    """ adds trajectory extension, if missing
    """
    return _add_extension(file_name, Extension.TRAJECTORY)


def zmatrix(file_name):
    """ adds zmatrix extension, if missing
    """
    return _add_extension(file_name, Extension.ZMATRIX)


def vmatrix(file_name):
    """ adds variable zmatrix extension, if missing
    """
    return _add_extension(file_name, Extension.VMATRIX)


def gradient(file_name):
    """ adds gradient extension, if missing
    """
    return _add_extension(file_name, Extension.GRADIENT)


def hessian(file_name):
    """ adds hessian extension, if missing
    """
    return _add_extension(file_name, Extension.HESSIAN)


def lennard_jones_epsilon(file_name):
    """ adds lennard-jones epsilon extension, if missing
    """
    return _add_extension(file_name, Extension.LJ_EPSILON)


def lennard_jones_sigma(file_name):
    """ adds lennard-jones sigma extension, if missing
    """
    return _add_extension(file_name, Extension.LJ_SIGMA)


def _add_extension(file_name, ext):
    if not str(file_name).endswith(ext):
        file_name = '{}{}'.format(file_name, ext)
    return file_name
