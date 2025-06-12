"""
  Get frequencies
"""

import projrot_io.writer
from autorun._run import from_input_string


# Default names of input and output files
SCRIPT_NAME = 'run_projrot.sh'
INPUT_NAME = 'RPHt_input_data.dat'
OUTPUT_NAMES = (
        'RTproj_freq.dat',
        'hrproj_freq.dat',
        'RTproj_cd.dat',
        'hrproj_cd.dat',
        'imactint.dat'
)


# Specialized Runners
def pot_frequencies(script_str, geoms, grads, hessians, run_path, 
        rotors_str=''):
    """ Calculate the frequencies (need to replace)

        all input here are dictionaries (like potentials)
    """

    # Initialize hr freqs list
    hr_freqs = {}
    imag_freqs = ()
    for i, point in enumerate(geoms.keys()):
        _, proj_freqs, rt_imag_freqs, hr_imag_freqs = frequencies(
            script_str,
            run_path,
            [geoms[point]],
            [grads[point]],
            [hessians[point]],
            rotors_str=rotors_str)
        if i == 0:
            imag_freqs = rt_imag_freqs
        if len(imag_freqs) > 0 and len(hr_imag_freqs) < 1:
            print(imag_freqs, proj_freqs)
            print(f'Imaginary frequencies lost during projection, for pot#{i}')
            print('setting lowest frequency to imag')
            proj_freqs = proj_freqs[1:]
        hr_freqs[point] = proj_freqs

    return hr_freqs


def displacements(script_str, run_dir, geoms, grads, hessians,
                  rotors_str='', aux_dct=None):
    """ Calculate the vibrational frequencies for a single molecule
        for the RT-projections and RTHr-projections.
    """

    output_strs = direct(
        script_str, run_dir, geoms, grads, hessians,
        rotors_str=rotors_str, aux_dct=aux_dct)
    disp_str = output_strs[2]

    disps = projrot_io.reader.normal_coordinates(disp_str)

    return disp_str, disps


def frequencies(script_str, run_dir, geoms, grads, hessians,
                rotors_str='', aux_dct=None):
    """ Calculate the vibrational frequencies for a single molecule
        for the RT-projections and RTHr-projections.
    """

    output_strs = direct(
        script_str, run_dir, geoms, grads, hessians,
        rotors_str=rotors_str, aux_dct=aux_dct)

    rtproj_str = output_strs[0]
    if rtproj_str is not None:
        rtproj_freqs, rt_imag_freq = projrot_io.reader.rpht_output(
            rtproj_str)
    else:
        rtproj_freqs, rt_imag_freq = (), ()

    hrproj_str = output_strs[1]
    if hrproj_str is not None:
        hrproj_freqs, hr_imag_freq = projrot_io.reader.rpht_output(
            hrproj_str)
    else:
        hrproj_freqs, hr_imag_freq = (), ()
    return rtproj_freqs, hrproj_freqs, rt_imag_freq, hr_imag_freq


def small_curvature_tunneling(script_str, run_dir, geoms, grads, hessians,
                              rpath_coords, rpath_enes,  # sadpt_idx,
                              rotors_str=''):
    """ Determine the transmission
    """

    projrot_en_str = projrot_io.writer.rpht_path_coord_en(
        rpath_coords, rpath_enes,
        bnd1=(), bnd2=())
    aux_dct = {'RPHt_coord_en.dat': projrot_en_str}

    output_strs = direct(
        script_str, run_dir, geoms, grads, hessians,
        rotors_str=rotors_str, aux_dct=aux_dct,
        script_name=SCRIPT_NAME,
        input_name=INPUT_NAME,
        output_names=('imactint.dat',))
    # coord_proj=cart/int

    return output_strs[0]


# Generalized Runner
def direct(script_str, run_dir, geoms, grads, hessians,
           rotors_str='', aux_dct=None,
           script_name=SCRIPT_NAME,
           input_name=INPUT_NAME,
           output_names=OUTPUT_NAMES):
    """ Generates an input file for a ProjRot job, runs it directly, and
        obtains all of the possible output file strings
    """

    input_str = projrot_io.writer.rpht_input(
        geoms, grads, hessians, rotors_str=rotors_str)

    output_strs = from_input_string(
        script_str, run_dir, input_str,
        aux_dct=aux_dct,
        script_name=script_name,
        input_name=input_name,
        output_names=output_names)

    return output_strs
