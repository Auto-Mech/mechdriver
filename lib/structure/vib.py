"""
  Get frequencies
"""

import os
import autofile
import projrot_io
from lib.submission import run_script
from lib.submission import DEFAULT_SCRIPT_DCT


def projrot_frequencies(geo, hess, thy_info, thy_run_fs,
                        script_str=DEFAULT_SCRIPT_DCT['projrot']):
    """ Get the projected frequencies from projrot code
    """

    # Write the string for the ProjRot input
    thy_run_fs[-1].create(thy_info[1:4])
    thy_run_path = thy_run_fs[-1].path(thy_info[1:4])

    coord_proj = 'cartesian'
    grad = ''
    rotors_str = ''
    projrot_inp_str = projrot_io.writer.rpht_input(
        geo, grad, hess, rotors_str=rotors_str,
        coord_proj=coord_proj)

    bld_locs = ['PROJROT', 0]
    bld_run_fs = autofile.fs.build(thy_run_path)
    bld_run_fs[-1].create(bld_locs)
    projrot_path = bld_run_fs[-1].path(bld_locs)

    proj_file_path = os.path.join(projrot_path, 'RPHt_input_data.dat')
    with open(proj_file_path, 'w') as proj_file:
        proj_file.write(projrot_inp_str)

    run_script(script_str, projrot_path)

    imag_freq = ''
    if os.path.exists(projrot_path+'/hrproj_freq.dat'):
        with open(projrot_path+'/hrproj_freq.dat', 'r') as projfile:
            hrproj_str = projfile.read()
        rthrproj_freqs, imag_freq = projrot_io.reader.rpht_output(
            hrproj_str)
        proj_freqs = rthrproj_freqs
    else:
        with open(projrot_path+'/RTproj_freq.dat', 'r') as projfile:
            rtproj_str = projfile.read()
        rtproj_freqs, imag_freq = projrot_io.reader.rpht_output(
            rtproj_str)
        proj_freqs = rtproj_freqs

    return proj_freqs, imag_freq
