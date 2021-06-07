""" Reads all of the MechDriver input files and formats
    the data into usable Python objects that are used by
    all subsequent drivers in the workflow.
"""

import os
import ioformat
import automol
import autofile
from mechlib.amech_io import printer as ioprinter


# Names of input file strings
RUN_INP = 'inp/run.dat'
CSV_INP = 'inp/species.csv'
DAT_INP = 'inp/species.dat'
MECH_INP = 'inp/mechanism.dat'
SORT_INP = 'inp/sort.dat'
MODEL_INP = 'inp/models.dat'
THY_INP = 'inp/theory.dat'


def read_amech_input(job_path):
    """ Parse the run.dat file
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

    return {
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


# formatters, dont know where to build this
def _geometry_dictionary(job_path):
    """ read in dictionary of saved geometries
    """

    geo_dct = {}
    for path in _inp_file_paths(job_path):
        if file_path.endswith('.xyz'):
            xyz_str = ioformat.pathtools.read_file(file_path)
            spc_name = automol.geom.comment_from_xyz_string(xyz_str)
            geo = automol.geom.from_xyz_string(xyz_str)
            if spc_name in geo_dct:
                print('Warning: Dupilicate xyz geometry for ', spc_name)
            geo_dct[spc_name] = geo

    return geo_dct


def _active_space_dictionary(job_path):
    """ Read in files for active space calculations
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
    for path in _inp_file_paths(job_path):
        if file_path.endswith('.asp'):
            aspace_str = ioformat.pathtools.read_file(file_path)
            spc_name = _comment_name(aspace_str)
            if spc_name in aspace_dct:
                print('Warning: Dupilicate active space geometry for ', spc_name)
            aspace_dct[spc_name] = aspace_str

    return aspace_dct


def _inp_file_paths(job_path):
    """ Read the paths to all files 
    """

    file_paths = ()

    geom_path = os.path.join(job_path, 'data')
    for dir_path, _, file_names in os.walk(geom_path):
        for file_name in file_names:
            file_paths += (os.path.join(dir_path, file_name),)

    return file_paths 
