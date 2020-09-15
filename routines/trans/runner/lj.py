"""
input file writing routines for parallel OneDMin runs
"""

import os
import stat
import subprocess
import random
import elstruct

from lib.amech_io.writer import write_files


def write():
    """ write all of the input files
    """

    # Write the 1DMin input file
    print('  Writing input files...')
    onedmin_str = _write_input(
        job_dir_path, NSAMPS,
        target_name='target.xyz', bath_name='bath.xyz',
        smin=SMIN, smax=SMAX)

    # Write the geometry files
    print('  Writing xyz files for target and bath....')
    xyz1_str, xyz2_str = _write_xyz(
        job_dir_path, TARGET_GEO, BATH_GEO)

    # Write the electronic structure template file
    print('  Writing electronic structure submission inp template...')
    elstruct_inp_str = _write_elstruct_inp(
        job_dir_path,
        chg, mlt, run_method, run_basis, theory_info,
        run_prog, run_memory)

    # Write the electronic structure sumbission script
    print('  Writing electronic structure submission script...')
    elstruct_sub_str = _write_elstruct_sub(
        job_dir_path, DRIVE_PATH, RUN_PROG)

    # Collate the input strings and write the remaining files
    input_strs = (
        onedmin_str,
        xyz1_str, xyz2_str,
        elstruct_inp_str, elstruct_sub_str)
    input_names = (
        'input.dat'
        'target.xyz', 'bath.xyz'
        'qc.mol', 'elstruct.x')
    inp = tuple(zip(input_strs, input_names))
    write_files(inp, vrc_path)


def _write_input(job_dir_path, nsamp,
                 target_name='target.xyz', bath_name='bath.xyz',
                 smin=2, smax=6):
    """ write the input file
    """

    ranseed = random.randrange(1E8, 1E9)
    inp_str = py1dmin.interface.writer.onedmin_input(
        ranseed, nsamp, target_name, bath_name, smin, smax)


def _write_elstruct_inp(job_dir_path,
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


def make_job_dirs(etrans_run_path, idx):
    """ write the files
    """
    job_dir_path = os.path.join(
        etrans_run_path, 'run{0}'.format(str(idx+1)))
    os.mkdir(job_dir_path)
    print('\n\nWriting files to'+job_dir_path)
