"""
routines for running OneDMin
"""

import os
import statistics
import py1dmin.interface


def lj_parameters(run_path):
    """ get the lj params from each run and average them together
    """

    avg_sigma = None
    avg_epsilon = None

    job_dirs = [os.path.join(run_path, directory)
                for directory in os.listdir(run_path)
                if 'build' not in directory and 'yaml' not in directory]
    sigmas, epsilons = [], []
    for job_dir in job_dirs:
        lj_file_name = os.path.join(job_dir, 'lj.out')
        if os.path.exists(lj_file_name):
            with open(lj_file_name, 'r') as lj_file:
                output_string = lj_file.read()
            sigs, epss = py1dmin.interface.reader.lennard_jones(
                output_string)
            if sigs is not None and epss is not None:
                for sig, eps in zip(sigs, epss):
                    sigmas.append(sig)
                    epsilons.append(eps)

    assert len(sigmas) == len(epsilons)
    if sigmas and epsilons:
        avg_sigma = statistics.mean(sigmas)
        avg_epsilon = statistics.mean(epsilons)
        print('{0:<14s}{1:<16s}'.format('\nSigma (Ang)', 'Epsilon (cm-1)'))
        for sig, eps in zip(sigmas, epsilons):
            print('{0:<14.4f}{1:<16.4f}'.format(sig, eps))
        print('\nAverage Sigma =', avg_sigma)
        print('Average Epsilon =', avg_epsilon)
        print('Number of values = ', len(sigmas))
    else:
        print('No Sigma and Epsilon Values obtained =(')

    return avg_sigma, avg_epsilon


def lj_well_geometries(run_path):
    """ get the miniumum geometries file from each run
    """

    geo_str = ''

    job_dirs = [os.path.join(run_path, directory)
                for directory in os.listdir(run_path)
                if 'build' not in directory and 'yaml' not in directory]
    for job_dir in job_dirs:
        geo_file_name = os.path.join(job_dir, 'lj.out')
        if os.path.exists(geo_file_name):
            with open(geo_file_name, 'r') as geo_file:
                geo_str += geo_file.read()

    return geo_str


def zero_energies(run_path):
    """ get the zero-energy
    """

    ene = None

    job_dirs = [os.path.join(run_path, directory)
                for directory in os.listdir(run_path)
                if 'build' not in directory and 'yaml' not in directory]
    for job_dir in job_dirs:
        ene_file_name = os.path.join(job_dir, 'zero.ene')
        if os.path.exists(ene_file_name):
            with open(ene_file_name, 'r') as ene_file:
                ene_str = ene_file.read()
            ene = float(ene_str)
            # use autofile above to read the string
            break

    return ene
