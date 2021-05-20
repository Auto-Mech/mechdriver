""" Reads all of the MechDriver input files provided by the user
    at the start of the calculation.
"""

import os
import sys
import ioformat
import automol
from mechlib.amech_io import printer as ioprinter


# Names of input file strings (relative to the submission directory)
# INP_STR_DCT = {
RUN_INP = 'inp/run.dat'
CSV_INP = 'inp/species.csv'
DAT_INP = 'inp/species.dat'
MECH_INP = 'inp/mechanism.dat'
SORT_INP = 'inp/sort.dat'
MODEL_INP = 'inp/models.dat'
THY_INP = 'inp/theory.dat'


def read_amech_input(job_path):
    """ Reads all MechDriver input files provided by the user into strings.
        All whitespace and comment lines are stripped from the files.

        :param job_path: directory path where all input files exist
        :type job_path: str
        :rtype dict[str:str]
    """

    # Read required input strings
    run_str = ioformat.pathtools.read_file(
        job_path, RUN_INP, remove_comments='#', remove_whitespace=True)
    ioprinter.reading('run.dat...', newline=1)  # Add a Found <file> to msg

    thy_str = ioformat.pathtools.read_file(
        job_path, THY_INP, remove_comments='#', remove_whitespace=True)
    ioprinter.reading('theory.dat...', newline=1)

    mod_str = ioformat.pathtools.read_file(
        job_path, MODEL_INP, remove_comments='#', remove_whitespace=True)
    ioprinter.reading('model.dat...', newline=1)

    spc_str = ioformat.pathtools.read_file(job_path, CSV_INP)
    ioprinter.reading('species.csv...', newline=1)

    mech_str = ioformat.pathtools.read_file(
        job_path, MECH_INP, remove_comments='!', remove_whitespace=True)
    ioprinter.reading('mechanism.dat...', newline=1)

    # Read auxiliary input strings
    dat_str = ioformat.pathtools.read_file(
        job_path, DAT_INP, remove_comments='#')
    ioprinter.reading('species.dat...', newline=1)

    # Read structural and template files
    geo_dct = _geometry_dictionary(job_path)
    ioprinter.reading('geom.xyzs...', newline=1)

    act_dct = _active_space_dictionary(job_path)
    ioprinter.reading('active_space templates...', newline=1)

    inp_str_dct = {
        # required
        'run': run_str,
        'thy': thy_str,
        'mod': mod_str,
        'spc': spc_str,
        'mech': mech_str,
        # auxiliary
        'dat': dat_str,
        # structural/template
        'geo': geo_dct,
        'act': act_dct
    }

    # Assess if all required strings are present
    _check_input_avail(inp_str_dct)

    return inp_str_dct


def _check_input_avail(inp_str_dct):
    """ Checks if needed input files exist by seeing if inputs are
        none.

        TODO: Right now it ignores mechanism.dat, need to fix that
    """

    keys = ('run', 'thy', 'mod', 'spc')
    names = ('run.dat', 'theory.dat', 'models.dat', 'species.csv')

    inp_missing = False
    for key, name in zip(keys, names):
        if inp_str_dct[key] is None:
            print('No file {} found in inp directory'.format(name))
            inp_missing = True

    if inp_missing:
        sys.exit()


# formatters, dont know where to build this
def _geometry_dictionary(job_path):
    """ Reads any .xyz files provided by the user in the directory
        where other input files. The function extracts the mechanism name
        of the species/transition-state from the comment line of the .xyz file
        and indexes each .xyz string with this name.

        :param job_path: directory path where all input files exist
        :type job_path: str
        :rtype: dict[str: str]
    """

    geo_dct = {}
    for file_path, file_name in zip(_inp_file_paths(job_path)):
        if file_path.endswith('.xyz'):
            xyz_str = ioformat.pathtools.read_file(file_path, file_name)
            spc_name = automol.geom.comment_from_xyz_string(xyz_str)
            geo = automol.geom.from_xyz_string(xyz_str)
            if spc_name in geo_dct:
                print('Warning: Dupilicate xyz geometry for ', spc_name)
            geo_dct[spc_name] = geo

    return geo_dct


def _active_space_dictionary(job_path):
    """ Reads any active-space templates provided by the user in the directory
        where other input files. The function extracts the mechanism name
        from a top-line comment in the template file.

        Currently, only Molpro templates are supported.

        :param job_path: directory path where all input files exist
        :type job_path: str
        :rtype: tuple(str)
    """

    def _comment_name(aspace_str):
        """ read the species name from a comment line in template
            comment line in FIRST line of template of form
            ! {spc_name}
        """
        comm_line = aspace_str.splitlines()[0]
        comm_line.replace('!', '').strip()
        return comm_line

    aspace_dct = {}
    for file_path, file_name in zip(_inp_file_paths(job_path)):
        if file_path.endswith('.asp'):
            aspace_str = ioformat.pathtools.read_file(file_path, file_name)
            spc_name = _comment_name(aspace_str)
            if spc_name in aspace_dct:
                print('Warning: Dupilicate active space geometry for ',
                      spc_name)
            aspace_dct[spc_name] = aspace_str

    return aspace_dct


def _inp_file_paths(job_path):
    """ Utility function to get paths to all of the auxiliary input files.

        :param job_path: directory path where all input files exist
        :type job_path: str
        :rtype: tuple(str)
    """

    file_paths = ()
    file_names = ()

    geom_path = os.path.join(job_path, 'data')
    for dir_path, _, file_names in os.walk(geom_path):
        for file_name in file_names:
            file_paths += (dir_path,)
            file_names += (file_name,)

    return file_paths, file_names
