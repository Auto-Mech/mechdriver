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

    # Read the input string from the file
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

    dat_str = ioformat.pathtools.read_file(
        job_path, DAT_INP, remove_comments='#')
    ioprinter.reading('species.dat...', newline=1)

    geo_dct = _geometry_dictionary(job_path)
    ioprinter.reading('geom.xyzs...', newline=1)

    mech_str = ioformat.pathtools.read_file(
        job_path, MECH_INP, remove_comments='!', remove_whitespace=True)
    ioprinter.reading('mechanism.dat...', newline=1)

    return {
        'run': run_str,
        'thy': thy_str,
        'mod': mod_str,
        'spc': spc_str,
        'dat': dat_str,
        'geo': geo_dct,
        'mech': mech_str
    }


# formatters, dont know where to build this
def _geometry_dictionary(job_path):
    """ read in dictionary of saved geometries
    """
    geom_path = os.path.join(job_path, 'data')
    geom_dct = {}
    for dir_path, _, file_names in os.walk(geom_path):
        for file_name in file_names:
            file_path = os.path.join(dir_path, file_name)
            if file_path.endswith('.xyz'):
                xyz_str = autofile.io_.read_file(file_path)
                geo = automol.geom.from_xyz_string(xyz_str)
                ich = automol.geom.inchi(geo)
                if ich in geom_dct:
                    print('Warning: Dupilicate xyz geometry for ', ich)
                geom_dct[ich] = geo

    return geom_dct


def _active_space_dictionary(job_path):
    """ Read in files for active space calculations
    """
