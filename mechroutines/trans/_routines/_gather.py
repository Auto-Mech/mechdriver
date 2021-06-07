"""
routines for running OneDMin
"""

import os
import onedmin_io
from mechlib.amech_io import printer as ioprinter


def read_filesys(etrans_save_fs, etrans_locs):
    """ get the lj params thar are saved currently in the filesystem
    """
    _, _ = etrans_save_fs, etrans_locs
    sigmas, epsilons, geoms = [], [], []
    return sigmas, epsilons, geoms


def read_output(run_path):
    """ get the lj params from each run and average them together
    """

    all_sigmas, all_epsilons = [], []
    for jobdir in _jobdirs(run_path):

        # Read the output file strings
        lj_str = _output_str(jobdir, 'lj.out')
        geo_str = _output_str(jobdir, 'min_geoms.out')

        # Parse the sigma and epsilon values from the output
        sigmas, epsilons = onedmin_io.reader.lennard_jones(lj_str)
        if sigmas is not None and epsilons is not None:
            all_sigmas.extend(sigmas)
            all_epsilons.extend(epsilons)

        # Parse the geometries from the min geoms file
        geoms = geo_str.split()

    return sigmas, epsilons, geoms


def prog_version(run_path):
    """ read the program and version
    """
    for jobdir in _jobdirs(run_path):
        lj_str = _output_str(jobdir, 'lj.out')
        break
    version = onedmin_io.reader.program_version(lj_str)

    return version


def _jobdirs(run_path):
    """ Obtain a list of all the directory names where OneDMin
        jobs were run
    """
    return [os.path.join(run_path, directory)
            for directory in os.listdir(run_path)
            if 'build' not in directory and 'yaml' not in directory]


def _output_str(jobdir, output_name):
    """ Read the output file
    """

    output_file_name = os.path.join(jobdir, output_name)
    if os.path.exists(output_file_name):
        with open(output_file_name, 'r') as outfile:
            out_str = outfile.read()
    else:
        out_str = None

    return out_str


def print_lj_parms(sigmas, epsilons):
    """ Print the lj parameters out
    """
    if sigmas and epsilons:
        ioprinter.info_message(
            '{0:<14s}{1:<16s}'.format('\nSigma (Ang)', 'Epsilon (cm-1)'))
        for sig, eps in zip(sigmas, epsilons):
            ioprinter.info_message(
                '{0:<14.4f}{1:<16.4f}'.format(sig, eps))
