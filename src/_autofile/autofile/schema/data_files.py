""" DataFiles
"""
from autofile import model
import autofile.data_types
import autofile.info


def information(file_prefix, function=None):
    """ information DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :param function: optional information-generator function, for checking the
        function signature against the information object
    :type function: callable
    :return: instance of DataFile class
    :rtype: Datafile
    """
    def writer_(inf_obj):
        if function is not None:
            assert autofile.info.matches_function_signature(inf_obj, function)
        inf_str = autofile.data_types.swrite.information(inf_obj)
        return inf_str

    def reader_(inf_str):
        inf_obj = autofile.data_types.sread.information(inf_str)
        if function is not None:
            assert autofile.info.matches_function_signature(inf_obj, function)
        return inf_obj

    name = autofile.data_types.name.information(file_prefix)
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def locator(file_prefix, map_dct_, loc_keys):
    """ locator DataFile

    Specifiers are stored in information files according to `map_dct_` and read
    back out according to `loc_keys_`. The file may contain auxiliary
    information (such as SMILES along with InChI), but for the read to work it
    must contain each locator value.

    :param map_dct_: Maps on the locator list to the values stored in the
        information file, by key.
    :type map_dct_: dict[key: callable]
    :param loc_keys: Keys to the original locator values.
    :type loc_keys: tuple[str]
    """
    def writer_(locs):
        inf_dct = {key: map_(locs) for key, map_ in map_dct_.items()}
        inf_obj = autofile.info.object_(inf_dct)
        return autofile.data_types.swrite.information(inf_obj)

    def reader_(inf_str):
        inf_obj = autofile.data_types.sread.information(inf_str)
        inf_dct = dict(inf_obj)
        return list(map(inf_dct.__getitem__, loc_keys))

    name = autofile.data_types.name.information(file_prefix)
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def input_file(file_prefix):
    """ generate input file DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.input_file(file_prefix)
    return model.DataFile(name=name)


def output_file(file_prefix):
    """ generate output file DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.output_file(file_prefix)
    return model.DataFile(name=name)


def instability(file_prefix):
    """ Generate a data file for instabiliy
    """
    name = autofile.data_types.name.instability(file_prefix)
    writer_ = autofile.data_types.swrite.instability
    reader_ = autofile.data_types.sread.instability
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def energy(file_prefix):
    """ generate energy DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.energy(file_prefix)
    writer_ = autofile.data_types.swrite.energy
    reader_ = autofile.data_types.sread.energy
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def geometry(file_prefix):
    """ generate geometry DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.geometry(file_prefix)
    writer_ = autofile.data_types.swrite.geometry
    reader_ = autofile.data_types.sread.geometry
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def gradient(file_prefix):
    """ generate gradient DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.gradient(file_prefix)
    writer_ = autofile.data_types.swrite.gradient
    reader_ = autofile.data_types.sread.gradient
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def hessian(file_prefix):
    """ generate hessian DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.hessian(file_prefix)
    writer_ = autofile.data_types.swrite.hessian
    reader_ = autofile.data_types.sread.hessian
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def harmonic_frequencies(file_prefix):
    """ generate harmonic_frequencies DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.harmonic_frequencies(file_prefix)
    writer_ = autofile.data_types.swrite.harmonic_frequencies
    reader_ = autofile.data_types.sread.harmonic_frequencies
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def anharmonic_frequencies(file_prefix):
    """ generate anharmonic_frequencies DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.anharmonic_frequencies(file_prefix)
    writer_ = autofile.data_types.swrite.anharmonic_frequencies
    reader_ = autofile.data_types.sread.anharmonic_frequencies
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def anharmonic_zpve(file_prefix):
    """ generate anharmonic_zpve DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.anharmonic_zpve(file_prefix)
    writer_ = autofile.data_types.swrite.anharmonic_zpve
    reader_ = autofile.data_types.sread.anharmonic_zpve
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def cubic_force_constants(file_prefix):
    """ generate cubic_force_constants DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.cubic_force_constants(file_prefix)
    writer_ = autofile.data_types.swrite.cubic_force_constants
    reader_ = autofile.data_types.sread.cubic_force_constants
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def quartic_force_constants(file_prefix):
    """ generate quartic_force_constants DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.quartic_force_constants(file_prefix)
    writer_ = autofile.data_types.swrite.quartic_force_constants
    reader_ = autofile.data_types.sread.quartic_force_constants
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def anharmonicity_matrix(file_prefix):
    """ generate anharmonicity matrix DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.anharmonicity_matrix(file_prefix)
    writer_ = autofile.data_types.swrite.anharmonicity_matrix
    reader_ = autofile.data_types.sread.anharmonicity_matrix
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def vibro_rot_alpha_matrix(file_prefix):
    """ generate vibro_rot_alpha matrix DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.vibro_rot_alpha_matrix(file_prefix)
    writer_ = autofile.data_types.swrite.vibro_rot_alpha_matrix
    reader_ = autofile.data_types.sread.vibro_rot_alpha_matrix
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def quartic_centrifugal_dist_consts(file_prefix):
    """ generate vibro_rot_alpha matrix DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = (
        autofile.data_types.name.quartic_centrifugal_dist_consts(file_prefix))
    writer_ = autofile.data_types.swrite.quartic_centrifugal_dist_consts
    reader_ = autofile.data_types.sread.quartic_centrifugal_dist_consts
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def zmatrix(file_prefix):
    """ generate zmatrix DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.zmatrix(file_prefix)
    writer_ = autofile.data_types.swrite.zmatrix
    reader_ = autofile.data_types.sread.zmatrix
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def vmatrix(file_prefix):
    """ generate vmatrix DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """

    name = autofile.data_types.name.vmatrix(file_prefix)
    writer_ = autofile.data_types.swrite.vmatrix
    reader_ = autofile.data_types.sread.vmatrix
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def ring_torsions(file_prefix):
    """ generate ring torsions DataFile
    """
    name = autofile.data_types.name.ring_torsions(file_prefix)
    writer_ = autofile.data_types.swrite.ring_torsions
    reader_ = autofile.data_types.sread.ring_torsions
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def torsions(file_prefix):
    """ generate torsions DataFile
    """
    name = autofile.data_types.name.torsions(file_prefix)
    writer_ = autofile.data_types.swrite.torsions
    reader_ = autofile.data_types.sread.torsions
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def trajectory(file_prefix):
    """ generate trajectory DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.trajectory(file_prefix)
    writer_ = autofile.data_types.swrite.trajectory
    reader_ = autofile.data_types.sread.trajectory
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def reaction(file_prefix):
    """ generate reaction DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.reaction(file_prefix)
    writer_ = autofile.data_types.swrite.reaction
    reader_ = autofile.data_types.sread.reaction
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def lennard_jones_epsilon(file_prefix):
    """ generate lennard_jones_epsilon DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.lennard_jones_epsilon(file_prefix)
    writer_ = autofile.data_types.swrite.lennard_jones_epsilon
    reader_ = autofile.data_types.sread.lennard_jones_epsilon
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def lennard_jones_sigma(file_prefix):
    """ generate lennard_jones_sigma DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.lennard_jones_sigma(file_prefix)
    writer_ = autofile.data_types.swrite.lennard_jones_sigma
    reader_ = autofile.data_types.sread.lennard_jones_sigma
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def external_symmetry_number(file_prefix):
    """ generate external_symmetry_number DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.external_symmetry_factor(file_prefix)
    writer_ = autofile.data_types.swrite.external_symmetry_factor
    reader_ = autofile.data_types.sread.external_symmetry_factor
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def internal_symmetry_number(file_prefix):
    """ generate internal_symmetry_number DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.internal_symmetry_factor(file_prefix)
    writer_ = autofile.data_types.swrite.internal_symmetry_factor
    reader_ = autofile.data_types.sread.internal_symmetry_factor
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def lennard_jones_input(file_prefix):
    """ generate input file for the LJ params program

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.lennard_jones_input(file_prefix)
    return model.DataFile(name=name)


def lennard_jones_elstruct(file_prefix):
    """ generate elec struct template file for the LJ params program

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.lennard_jones_elstruct(file_prefix)
    return model.DataFile(name=name)


def dipole_moment(file_prefix):
    """ generate dipole_moment DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.dipole_moment(file_prefix)
    writer_ = autofile.data_types.swrite.dipole_moment
    reader_ = autofile.data_types.sread.dipole_moment
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def polarizability(file_prefix):
    """ generate polarizability DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.polarizability(file_prefix)
    writer_ = autofile.data_types.swrite.polarizability
    reader_ = autofile.data_types.sread.polarizability
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


#  vrctst
def vrctst_tst(file_prefix):
    """ generate vrcttst_tst DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.vrctst_tst(file_prefix)
    return model.DataFile(name=name)
    # return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def vrctst_divsur(file_prefix):
    """ generate vrctst_divsur DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.vrctst_divsur(file_prefix)
    return model.DataFile(name=name)
    # return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def vrctst_molpro(file_prefix):
    """ generate vrctst_molpro DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.vrctst_molpro(file_prefix)
    return model.DataFile(name=name)
    # return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def vrctst_tml(file_prefix):
    """ generate vrctst_tml DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.vrctst_tml(file_prefix)
    return model.DataFile(name=name)
    # return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def vrctst_struct(file_prefix):
    """ generate vrctst_struct DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.vrctst_struct(file_prefix)
    return model.DataFile(name=name)
    # return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def vrctst_pot(file_prefix):
    """ generate vrctst_pot DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.vrctst_pot(file_prefix)
    return model.DataFile(name=name)
    # return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def vrctst_flux(file_prefix):
    """ generate vrctst_flux DataFile

    :param file_prefix: path to file
    :type file_prefix: str
    :return: instance of DataFile class
    :rtype: Datafile
    """
    name = autofile.data_types.name.vrctst_flux(file_prefix)
    return model.DataFile(name=name)
    # return model.DataFile(name=name, writer_=writer_, reader_=reader_)
