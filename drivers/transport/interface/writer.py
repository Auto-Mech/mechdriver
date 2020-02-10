"""
Executes the automation part of 1DMin
"""

import os
from mako.template import Template


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')


def onedmin_input(ranseed, nsamp,
                  target_xyz_name, bath_xyz_name, smin, smax):
    """ writes the 1dmin input file for each instance
    """

    # Set the dictionary for the 1DMin input file
    fill_vals = {
        "ranseed": ranseed,
        "nsamples": nsamp,
        "target_xyz_name": target_xyz_name,
        "bath_xyz_name": bath_xyz_name,
        "smin": smin,
        "smax": smax
    }

    # Set template name and path for the 1dmin input file
    template_file_name = '1dmin_inp.mako'
    template_file_path = os.path.join(TEMPLATE_PATH, template_file_name)

    # Build the 1dmin input string
    input_str = Template(filename=template_file_path).render(**fill_vals)

    return input_str


def submission_script(drive_path, run_path, njobs):
    """ launches the job
    """

    # Write the bottom of the string
    job_exe_lines = '# Run several onedmin.x instances\n'
    job_exe_lines += 'cd {0}/run1\n'.format(run_path)
    job_exe_lines += 'time $ONEDMINEXE < input.dat > output.dat &\n'
    for i in range(njobs-1):
        job_exe_lines += 'cd ../run{0}\n'.format(str(i+2))
        job_exe_lines += 'time $ONEDMINEXE < input.dat > output.dat &\n'
    job_exe_lines += 'wait\n'

    # Set the dictionary for the 1DMin input file
    fill_vals = {
        "job_exe_lines": job_exe_lines,
    }

    # Set template name and path for the 1dmin input file
    template_file_name = 'submit.mako'
    template_file_path = os.path.join(drive_path, template_file_name)

    # Build the 1dmin input string
    sub_str = Template(filename=template_file_path).render(**fill_vals)

    return sub_str
