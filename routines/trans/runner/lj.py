"""
input file writing routines for parallel OneDMin runs
"""

import os
import subprocess
import random
import moldr
import elstruct
import py1dmin.interface


def write_input(job_dir_path, nsamp,
                target_name='target.xyz', bath_name='bath,xyz',
                smin=2, smax=6):
    """ write the input file
    """

    ranseed = random.randrange(1E8, 1E9)
    inp_str = py1dmin.interface.writer.onedmin_input(
        ranseed, nsamp, target_name, bath_name, smin, smax)

    job_file_path = os.path.join(job_dir_path, 'input.dat')
    with open(job_file_path, 'w') as input_file:
        input_file.write(inp_str)


def write_xyz(job_dir_path, target_geo, bath_geo):
    """ write the target and bath xyz files
    """

    job_file_path = os.path.join(job_dir_path, 'target.xyz')
    with open(job_file_path, 'w') as xyz_file:
        xyz_file.write(target_geo)

    job_file_path = os.path.join(job_dir_path, 'bath.xyz')
    with open(job_file_path, 'w') as xyz_file:
        xyz_file.write(bath_geo)


def write_elstruct_inp(job_dir_path,
                       charge, mult, method, basis, thry_lvl,
                       prog, memory):
    """ writes the electronic structure input file
    """
    assert prog in ('g09', 'molpro')

    elstruct_inp_str = elstruct.writer.energy(
        geom='GEOMETRY',
        charge=charge,
        mult=mult,
        method=method,
        basis=basis,
        prog=prog,
        mol_options=('nosym', 'noorient', 'angstrom'),
        memory=memory,
        comment='SAMPLE GEOM',
        orb_restricted=moldr.util.orbital_restriction(
            ['', charge, mult], thry_lvl),
    )

    elstruct_inp_name = os.path.join(job_dir_path, 'qc.mol')
    with open(elstruct_inp_name, 'w') as elstruct_inp_file:
        elstruct_inp_file.write(elstruct_inp_str)


def write_elstruct_sub(job_dir_path, drive_path, prog):
    """ writes the elstruct submission string
    """

    if prog == 'molpro':
        sub_name = 'm.x'
    else:
        raise NotImplementedError

    # Read the original input script as a string
    elstruct_sub_name_in = os.path.join(drive_path, 'elstruct.x')
    with open(elstruct_sub_name_in, 'r') as elstruct_sub_file_in:
        in_sub_str = elstruct_sub_file_in.read()

    # Write the electronic structure sumbission script
    elstruct_sub_name = os.path.join(job_dir_path, sub_name)
    with open(elstruct_sub_name, 'w') as elstruct_sub_file:
        elstruct_sub_file.write(in_sub_str)

    # Make the auto1dmin.x and elstruct.x file executable
    subprocess.check_call(['chmod', '+x', elstruct_sub_name])


def submit_job(drive_path, run_path, njobs):
    """ submit the 1DMin job
    """

    submit_str = py1dmin.interface.writer.submission_script(
        drive_path, run_path, njobs)

    moldr.util.run_script(submit_str, run_path)
