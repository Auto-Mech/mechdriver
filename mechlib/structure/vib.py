"""
  Get frequencies
"""

import os
import random
import projrot_io
import autofile
from ioformat import run_script
from mechlib import filesys
from mechlib.submission import DEFAULT_SCRIPT_DCT


def projrot_freqs(geoms, hessians, run_path,
                  grads=((),), rotors_str='', coord_proj='cartesian',
                  script_str=DEFAULT_SCRIPT_DCT['projrot']):
    """ Get the projected frequencies from projrot code
        run path at thy later
    """

    # Set up the filesys
    bld_locs = ['PROJROT', 0]
    bld_save_fs = autofile.fs.build(run_path)
    bld_save_fs[-1].create(bld_locs)
    projrot_path = bld_save_fs[-1].path(bld_locs)

    projrot_path = os.path.join(projrot_path, str(random.randint(0,1234567)))
    if not os.path.exists(projrot_path):
        os.makedirs(projrot_path)

    # print('run path test for ProjRot:', run_path)
    # bld_fs, bld_locs = filesys.build.build_fs(
        # run_path, 'PROJROT', locs_idx=None)
    # bld_fs[-1].create(bld_locs)
    # projrot_path = bld_fs[-1].path(bld_locs)

    print('Run path for ProjRot:')
    print(projrot_path)

    # Write the ProjRot input file
    projrot_inp_str = projrot_io.writer.rpht_input(
        geoms, grads, hessians, rotors_str=rotors_str,
        coord_proj=coord_proj)
    proj_file_path = os.path.join(projrot_path, 'RPHt_input_data.dat')
    with open(proj_file_path, 'w') as proj_file:
        proj_file.write(projrot_inp_str)

    # Run ProjRot
    run_script(script_str, projrot_path)

    # Read vibrational frequencies from ProjRot output
    rtproj_file = os.path.join(projrot_path, 'RTproj_freq.dat')
    if os.path.exists(rtproj_file):
        with open(rtproj_file, 'r') as projfile:
            rtproj_str = projfile.read()
        rtproj_freqs, rt_imag_freq = projrot_io.reader.rpht_output(
            rtproj_str)
    else:
        rtproj_freqs, rt_imag_freq = [], []

    hrproj_file = os.path.join(projrot_path, 'hrproj_freq.dat')
    if os.path.exists(hrproj_file):
        with open(hrproj_file, 'r') as projfile:
            hrproj_str = projfile.read()
        hrproj_freqs, hr_imag_freq = projrot_io.reader.rpht_output(
            hrproj_str)
    else:
        hrproj_freqs, hr_imag_freq = [], []

    return rtproj_freqs, hrproj_freqs, rt_imag_freq, hr_imag_freq
