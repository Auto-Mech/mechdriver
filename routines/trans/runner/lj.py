"""
input file writing routines for parallel OneDMin runs
"""

import os
import random
import automol
import autofile
import elstruct
import onedmin_io
from lib.submission import qchem_params


def write_input(nsamp, smin=2.0, smax=6.0,
                target_name='target.xyz', bath_name='bath.xyz',
                spin_method=1):
    """ write the input file
    """

    ranseed = random.randrange(1E8, 1E9)

    inp_str = onedmin_io.writer.input_file(
        ranseed, nsamp, smin, smax,
        target_name, bath_name,
        spin_method=spin_method)

    return inp_str


def write_xyz(tgt_geo, bath_geo):
    """ Write the geom strings
    """

    tgt_str = automol.geom.string(tgt_geo)
    bath_str = automol.geom.string(bath_geo)

    return tgt_str, bath_str


def write_elstruct_inp(lj_info, mod_lj_thy_info):
    """ writes the electronic structure input file
    """

    # Unpack info objects
    [_, charge, mult] = lj_info
    [prog, method, basis, orb_lbl] = mod_lj_thy_info

    assert prog in ('gaussian09', 'gaussian16', 'molpro2015')

    # Get the job running options
    script_str, _, kwargs, _ = qchem_params(
        *mod_lj_thy_info[0:2])

    # Write the string
    elstruct_inp_str = elstruct.writer.energy(
        geom='GEOMETRY',
        charge=charge,
        mult=mult,
        method=method,
        basis=basis,
        prog=prog,
        mol_options=('nosym', 'noorient'),
        memory=kwargs['memory'],
        comment='SAMPLE GEOM',
        orb_type=orb_lbl
    )

    elstruct_sp_str = '\n'.join(script_str.splitlines()[1:])

    return elstruct_inp_str, elstruct_sp_str


def write_onedmin_sub(njobs, job_path, onedmin_path,
                      exe_name='onedmin-dd-molpro.x'):
    """ Write the OneDMin submission script
    """
    return onedmin_io.writer.submission_script(
        njobs, job_path, onedmin_path, exe_name=exe_name)


# Make the dirs
def build_rundir(etrans_run_path):
    """ build the run directory
    """

    # for idx in range(100):
    bld_locs = ['ONEDMIN', 0]
    bld_save_fs = autofile.fs.build(etrans_run_path)
    bld_save_fs[-1].create(bld_locs)
    run_path = bld_save_fs[-1].path(bld_locs)

    print('Build Path for OneDMin calculations')
    print(run_path)

    return run_path


def make_jobdir(etrans_run_path, idx):
    """ write the files
    """
    job_path = os.path.join(
        etrans_run_path, 'run{0}'.format(str(idx+1))
    )
    os.makedirs(job_path)

    return job_path
