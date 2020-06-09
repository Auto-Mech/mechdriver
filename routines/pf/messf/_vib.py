"""
  Handle vibrational data info
"""

import os
import projrot_io
import autofile
from lib import structure
from lib.phydat import phycon
from lib.submission import run_script
from lib.submission import DEFAULT_SCRIPT_DCT


def read_harmonic_freqs(geom, cnf_save_fs, cnf_save_locs, saddle=False):
    """ Read the harmonic frequencies
    """

    # Probably just read the freqs from the filesys

    # Do the freqs obtain for two species for fake and pst
    if cnf_save_locs is not None:

        # Obtain geom and freqs from filesys
        hess = cnf_save_fs[-1].file.hessian.read(cnf_save_locs)
        freqs = elstruct.util.harmonic_frequencies(
            geom, hess, project=False)

        # Modify freqs lst and get imaginary frequencies
        mode_start = 6
        if saddle:
            mode_start = mode_start + 1
            imag_freq = freqs[0]
        else:
            imag_freq = None

        # Grab the freqs from the lst, for cases of linear and nonlinear mol
        if automol.geom.is_linear(geom):
            mode_start = mode_start - 1
        freqs = freqs[mode_start:]
    else:
        print('ERROR: Reference geometry is missing for harmonic frequencies')

    return freqs, imag_freq


def read_anharmon_matrix(cnf_save_fs, cnf_save_locs):
    """ Read the anharmonicity matrix """
    if cnf_save_locs is not None:
        xmat = cnf_save_fs[-1].file.anharmonicity_matrix.read(
            cnf_save_locs)
    else:
        print('No anharm matrix')

    return xmat


def tors_projected_freqs_zpe(geo, hess, proj_rotors_str, tors_zpe, save_path,
                             has_tors=False, saddle=False):
    """ Get frequencies from one version of ProjRot
    """
    
    # Calculate ZPVES of the hindered rotors
    tors_zpe = calc_tors_freqs_zpe(
        geom, sym_factor, elec_levels,
        mess_hr_str, cnf_save_path)

    # Run ProjRot to get the frequencies v1
    rt_freqs1, rth_freqs1, rt_imag1, rth_imag1 = structure.vib.projrot_freqs(
        geo, hess, thy_info, thy_run_fs,
        grad=(), rotors_str=projrot_str, coord_proj='cartesian',
        script_str=DEFAULT_SCRIPT_DCT['projrot'])

    # Run ProjRot to get the frequencies v2
    projrot_script_str2 = (
        "#!/usr/bin/env bash\n"
        "RPHt2.exe >& /dev/null"
    )
    rt_freqs2, rth_freqs2, rt_imag2, rth_imag2 = structure.vib.projrot_freqs(
        geo, hess, thy_info, thy_run_fs,
        grad=(), rotors_str=projrot_str, coord_proj='cartesian',
        script_str=projrot_script_str2)

    # Set the correct frequency set if there are any torsions
    if has_tors:
        freqs1, freqs2 = rth_freqs1, rth_freqs2
        imag1, imag2 = rth_imag1, rth_imag2
    else:
        freqs1, freqs2 = rt_freqs1, rt_freqs2
        imag1, imag2 = rt_imag1, rt_imag2

    # Calculate harmonic ZPVE from all harmonic freqs, including torsionals
    harm_zpe = (sum(rt_freqs1) / 2.0) * phycon.WAVEN2KCAL

    # Calculate harmonic ZPVE from freqs where torsions have been projected out
    # Value from both projrot versions, which use different projection schemes
    harm_zpe_notors_1 = (sum(freqs1) / 2.0) * phycon.WAVEN2KCAL
    harm_zpe_notors_2 = (sum(freqs2) / 2.0) * phycon.WAVEN2KCAL

    # Calcuate the difference in the harmonic ZPVE from projecting out torsions
    harm_tors_zpe = harm_zpe - harm_zpe_notors_1
    harm_tors_zpe_2 = harm_zpe - harm_zpe_notors_2

    # Check to see which of the above ZPVEs match more closely with tors ZPVE
    # calculated directly by treating the torsions in MESS
    diff_tors_zpe = harm_tors_zpe - tors_zpe
    diff_tors_zpe_2 = harm_tors_zpe_2 - tors_zpe
    if diff_tors_zpe <= diff_tors_zpe_2:
        zpe = zpe_harm_no_tors + tors_zpe
        freqs = freqs1
        imag_freq = imag_freq1
    else:
        zpe = zpe_harm_no_tors_2 + tors_zpe
        freqs = freqs2
        imag_freq = imag_freq2

    # Check if there are significant differences caused by the rotor projection
    if abs(diff_tors_zpe) > 0.2 and abs(diff_tors_zpe_2) > 0.2:
        print('Warning: There is a difference of ',
              '{0:.2f} and {1:.2f}'.format(diff_tors_zpe, diff_tors_zpe_2),
              'kcal/mol between harmonic and hindered torsional ZPVEs')

    return freqs, imag_freq, zpe
