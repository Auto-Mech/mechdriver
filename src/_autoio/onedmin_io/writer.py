""" Writes files for OneDmin.
"""

import os
from random import randrange
from phydat import phycon
from ioformat import build_mako_str


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')


def input_file(nsamp, smin, smax, ranseed=None,
               target_xyz_name='target.xyz', bath_xyz_name='bath.xyz',
               spin_method=2):
    """ Writes the input file string for a single instance of OneDMin.

        :param nsamp: number of samples
        :type nsamp: int
        :param smin: minimum-allowed intermolecular distance (in bohr)
        :type smin: float
        :param smax: maximum-allowed intermolecular distance (in bohr)
        :type smax: float
        :param ranseed: 10-digit integer for orientation sampling
        :type ranseed: int
        :param: target_xyz_name: name of xyz file with target molecule
        :type target_xyz_name: str
        :param: bath_xyz_name: name of xyz file with bath molecule
        :type bath_xyz_name: str
        :param spin_method: parameter to control orientational spinning
        :type spin_method: int
        :rtype: str
    """

    if ranseed is not None:
        assert isinstance(ranseed, int)
    else:
        ranseed = randrange(int(1.0e7), int(9.9e7))

    assert isinstance(nsamp, int)
    assert target_xyz_name.endswith('.xyz')
    assert bath_xyz_name.endswith('.xyz')
    assert isinstance(smin, float)
    assert isinstance(smax, float)

    # Format smin and smax
    smin = f'{smin*phycon.BOHR2ANG:.3f}'
    smax = f'{smax*phycon.BOHR2ANG:.3f}'

    # Set the dictionary for the 1DMin input file
    inp_keys = {
        "ranseed": ranseed,
        "nsamples": nsamp,
        "target_xyz_name": target_xyz_name,
        "bath_xyz_name": bath_xyz_name,
        "smin": smin,
        "smax": smax,
        "spin_method": spin_method
    }

    return build_mako_str(
        template_file_name='onedmin_inp.mako',
        template_src_path=TEMPLATE_PATH,
        template_keys=inp_keys)


def submission_script(njobs, run_dir, exe_path):
    """ Writes a special BASH submission script for launching
        parallel instances of OneDMin.

        :param njobs: number of OneDMin processes to run
        :type njobs: int
        :param run_dir: directory to run each OneDMin process
        :type run_dir: str
        :param exe_path: full path to the OneDMin executable
        :type exe_path: str
        :rtype: str
    """

    # Write the string for running all of the job lines
    # job_lines = 'mkdir -p {0}/run1\n'.format(run_dir)
    job_lines = f'cd {run_dir}/run1\n'
    job_lines += '$ONEDMINEXE < input.dat > output.dat &\n'
    for i in range(njobs-1):
        job_lines += f'cd ../run{str(i+2)}\n'
        job_lines += '$ONEDMINEXE < input.dat > output.dat &\n'
    job_lines += 'wait\n'

    # Set the dictionary for the 1DMin input file
    exe_keys = {
        "exe_path": exe_path,
        "job_lines": job_lines
    }

    return build_mako_str(
        template_file_name='onedmin_sub.mako',
        template_src_path=TEMPLATE_PATH,
        template_keys=exe_keys)
