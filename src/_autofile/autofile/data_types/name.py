""" file namers
"""


class Extension():
    """ file extensions """
    INFORMATION = '.yaml'
    INPUT_LOG = '.inp'
    OUTPUT_LOG = '.out'
    PROJROT_LOG = '.prot'
    TEMPLATE = '.temp'
    SHELL_SCRIPT = '.sh'
    ENERGY = '.ene'
    GEOMETRY = '.xyz'
    TRAJECTORY = '.t.xyz'
    ZMATRIX = '.zmat'
    VMATRIX = '.vmat'
    TORS = '.tors'
    RTORS = '.rtors'
    GRADIENT = '.grad'
    HESSIAN = '.hess'
    CUBIC_FC = '.cubic'
    QUARTIC_FC = '.quartic'
    HARMONIC_ZPVE = '.hzpve'
    ANHARMONIC_ZPVE = '.azpve'
    HARMONIC_FREQUENCIES = '.hfrq'
    ANHARMONIC_FREQUENCIES = '.afrq'
    ANHARMONICITY_MATRIX = '.xmat'
    VIBRO_ROT_MATRIX = '.vrmat'
    CENTRIF_DIST_CONSTS = '.qcd'
    LJ_EPSILON = '.eps'
    LJ_SIGMA = '.sig'
    EXTERNAL_SYMMETRY_FACTOR = '.esym'
    INTERNAL_SYMMETRY_FACTOR = '.isym'
    DIPOLE_MOMENT = '.dmom'
    POLARIZABILITY = '.polar'
    # Transformation files
    REACTION = '.r.yaml'
    # Instability Transformation files
    INSTAB = '.yaml'
    # Various VaReCoF files
    VRC_TST = '.tst'
    VRC_DIVSUR = '.divsur'
    VRC_MOLP = '.molpro'
    VRC_TML = '.tml'
    VRC_STRUCT = '.struct'
    VRC_POT = '.pot'
    VRC_FLUX = '.flux'
    JSON = '.json'


def information(file_name):
    """ adds information extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.INFORMATION)


def input_file(file_name):
    """ adds input file extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.INPUT_LOG)


def output_file(file_name):
    """ adds output file extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.OUTPUT_LOG)


def instability(file_name):
    """ adds instability extension, if missing
    """
    return _add_extension(file_name, Extension.INSTAB)


def run_script(file_name):
    """ adds run script extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.SHELL_SCRIPT)


def energy(file_name):
    """ adds energy extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.ENERGY)


def geometry(file_name):
    """ adds geometry extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.GEOMETRY)


def trajectory(file_name):
    """ adds trajectory extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.TRAJECTORY)


def zmatrix(file_name):
    """ adds zmatrix extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.ZMATRIX)


def vmatrix(file_name):
    """ adds variable zmatrix extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.VMATRIX)


def torsions(file_name):
    """ adds variable torsions extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.TORS)


def ring_torsions(file_name):
    """ adds variable torsions extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.RTORS)


def gradient(file_name):
    """ adds gradient extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.GRADIENT)


def hessian(file_name):
    """ adds hessian extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.HESSIAN)


def harmonic_zpve(file_name):
    """ adds harmonic zpve extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.HARMONIC_ZPVE)


def anharmonic_zpve(file_name):
    """ adds anharmonic zpve extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.ANHARMONIC_ZPVE)


def harmonic_frequencies(file_name):
    """ adds harmonic frequencies extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.HARMONIC_FREQUENCIES)


def anharmonic_frequencies(file_name):
    """ adds anharmonic frequencies extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.ANHARMONIC_FREQUENCIES)


def cubic_force_constants(file_name):
    """ adds cubic force constants extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.CUBIC_FC)


def quartic_force_constants(file_name):
    """ adds quartic force constants extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.QUARTIC_FC)


def anharmonicity_matrix(file_name):
    """ adds anharmonicity maxtrix extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.ANHARMONICITY_MATRIX)


def vibro_rot_alpha_matrix(file_name):
    """ adds vibro_rot_alpha maxtrix extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.VIBRO_ROT_MATRIX)


def quartic_centrifugal_dist_consts(file_name):
    """ adds quartic centrifugal distortion constants, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.CENTRIF_DIST_CONSTS)


def lennard_jones_epsilon(file_name):
    """ adds lennard-jones epsilon extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.LJ_EPSILON)


def lennard_jones_sigma(file_name):
    """ adds lennard-jones sigma extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.LJ_SIGMA)


def lennard_jones_input(file_name):
    """ adds lennard-jones input file extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.INPUT_LOG)


def lennard_jones_elstruct(file_name):
    """ adds lennard-jones sigma extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.TEMPLATE)


def external_symmetry_factor(file_name):
    """ adds external symmetry number extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.EXTERNAL_SYMMETRY_FACTOR)


def internal_symmetry_factor(file_name):
    """ adds internal symmetry number extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.INTERNAL_SYMMETRY_FACTOR)


def dipole_moment(file_name):
    """ adds dipole moment extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.DIPOLE_MOMENT)


def polarizability(file_name):
    """ adds dipole moment extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.POLARIZABILITY)


def reaction(file_name):
    """ adds reaction extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.REACTION)


def vrctst_tst(file_name):
    """ adds vrctst_tst extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.VRC_TST)


def vrctst_divsur(file_name):
    """ adds vrctst_divsur extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.VRC_DIVSUR)


def vrctst_molpro(file_name):
    """ adds vrctst_molpro extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.VRC_MOLP)


def vrctst_tml(file_name):
    """ adds vrctst_tml extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.VRC_TML)


def vrctst_struct(file_name):
    """ adds vrctst_struct extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.VRC_STRUCT)


def vrctst_pot(file_name):
    """ adds vrctst_pot extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.VRC_POT)


def vrctst_flux(file_name):
    """ adds vrctst_flux extension, if missing

    :param file_name: name of file
    :type file_name: str
    :returns: file with extension added
    :rtype: str
    """
    return _add_extension(file_name, Extension.VRC_FLUX)


def _add_extension(file_name, ext):
    if not str(file_name).endswith(ext):
        file_name = f'{file_name}{ext}'
    return file_name
