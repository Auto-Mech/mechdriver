""" DataFiles
"""
import autofile.file
import autofile.info
import autofile.system.info
from autofile.system import model


def information(file_prefix, function=None):
    """ information DataFile

    :param function: optional information-generator function, for checking the
        function signature against the information object
    :type function: callable
    """
    def writer_(inf_obj):
        if function is not None:
            assert autofile.info.matches_function_signature(inf_obj, function)
        inf_str = autofile.file.write.information(inf_obj)
        return inf_str

    def reader_(inf_str):
        inf_obj = autofile.file.read.information(inf_str)
        if function is not None:
            assert autofile.info.matches_function_signature(inf_obj, function)
        return inf_obj

    name = autofile.file.name.information(file_prefix)
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
        return autofile.file.write.information(inf_obj)

    def reader_(inf_str):
        inf_obj = autofile.file.read.information(inf_str)
        inf_dct = dict(inf_obj)
        return list(map(inf_dct.__getitem__, loc_keys))

    name = autofile.file.name.information(file_prefix)
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def input_file(file_prefix):
    """ generate input file DataFile
    """
    name = autofile.file.name.input_file(file_prefix)
    return model.DataFile(name=name)


def output_file(file_prefix):
    """ generate output file DataFile
    """
    name = autofile.file.name.output_file(file_prefix)
    return model.DataFile(name=name)


def energy(file_prefix):
    """ generate energy DataFile
    """
    name = autofile.file.name.energy(file_prefix)
    writer_ = autofile.file.write.energy
    reader_ = autofile.file.read.energy
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def geometry(file_prefix):
    """ generate geometry DataFile
    """
    name = autofile.file.name.geometry(file_prefix)
    writer_ = autofile.file.write.geometry
    reader_ = autofile.file.read.geometry
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def gradient(file_prefix):
    """ generate gradient DataFile
    """
    name = autofile.file.name.gradient(file_prefix)
    writer_ = autofile.file.write.gradient
    reader_ = autofile.file.read.gradient
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def hessian(file_prefix):
    """ generate hessian DataFile
    """
    name = autofile.file.name.hessian(file_prefix)
    writer_ = autofile.file.write.hessian
    reader_ = autofile.file.read.hessian
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def harmonic_frequencies(file_prefix):
    """ generate harmonic_frequencies DataFile
    """
    name = autofile.file.name.harmonic_frequencies(file_prefix)
    writer_ = autofile.file.write.harmonic_frequencies
    reader_ = autofile.file.read.harmonic_frequencies
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def anharmonic_frequencies(file_prefix):
    """ generate anharmonic_frequencies DataFile
    """
    name = autofile.file.name.anharmonic_frequencies(file_prefix)
    writer_ = autofile.file.write.anharmonic_frequencies
    reader_ = autofile.file.read.anharmonic_frequencies
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def anharmonic_zpve(file_prefix):
    """ generate anharmonic_zpve DataFile
    """
    name = autofile.file.name.anharmonic_zpve(file_prefix)
    writer_ = autofile.file.write.anharmonic_zpve
    reader_ = autofile.file.read.anharmonic_zpve
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def anharmonicity_matrix(file_prefix):
    """ generate anharmonicity matrix DataFile
    """
    name = autofile.file.name.anharmonicity_matrix(file_prefix)
    writer_ = autofile.file.write.anharmonicity_matrix
    reader_ = autofile.file.read.anharmonicity_matrix
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def vibro_rot_alpha_matrix(file_prefix):
    """ generate vibro_rot_alpha matrix DataFile
    """
    name = autofile.file.name.vibro_rot_alpha_matrix(file_prefix)
    writer_ = autofile.file.write.vibro_rot_alpha_matrix
    reader_ = autofile.file.read.vibro_rot_alpha_matrix
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def quartic_centrifugal_dist_consts(file_prefix):
    """ generate vibro_rot_alpha matrix DataFile
    """
    name = autofile.file.name.quartic_centrifugal_dist_consts(file_prefix)
    writer_ = autofile.file.write.quartic_centrifugal_dist_consts
    reader_ = autofile.file.read.quartic_centrifugal_dist_consts
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def zmatrix(file_prefix):
    """ generate zmatrix DataFile
    """
    name = autofile.file.name.zmatrix(file_prefix)
    writer_ = autofile.file.write.zmatrix
    reader_ = autofile.file.read.zmatrix
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def vmatrix(file_prefix):
    """ generate vmatrix DataFile
    """
    name = autofile.file.name.vmatrix(file_prefix)
    writer_ = autofile.file.write.vmatrix
    reader_ = autofile.file.read.vmatrix
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def trajectory(file_prefix):
    """ generate trajectory DataFile
    """
    name = autofile.file.name.trajectory(file_prefix)
    writer_ = autofile.file.write.trajectory
    reader_ = _not_implemented
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def lennard_jones_epsilon(file_prefix):
    """ generate lennard_jones_epsilon DataFile
    """
    name = autofile.file.name.lennard_jones_epsilon(file_prefix)
    writer_ = autofile.file.write.lennard_jones_epsilon
    reader_ = autofile.file.read.lennard_jones_epsilon
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


def lennard_jones_sigma(file_prefix):
    """ generate lennard_jones_sigma DataFile
    """
    name = autofile.file.name.lennard_jones_sigma(file_prefix)
    writer_ = autofile.file.write.lennard_jones_sigma
    reader_ = autofile.file.read.lennard_jones_sigma
    return model.DataFile(name=name, writer_=writer_, reader_=reader_)


# helpers
def _not_implemented(*_args, **_kwargs):
    raise NotImplementedError
