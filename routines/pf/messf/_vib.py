"""
  Handle vibrational data info
"""

import os
import projrot_io
import autofile
from lib.structure import tors as torsprep
from lib.structure import vib as vibprep
from lib.phydat import phycon
from lib.submission import run_script
from lib.submission import DEFAULT_SCRIPT_DCT


def read_harmonic_freqs(pf_filesystems, pf_levels):
    """ Read the harmonic frequencies
    """

    # Get the harmonic filesys information
    [cnf_fs, _, min_cnf_locs, _] = pf_filesystems['harm']

    # Do the freqs obtain for two species for fake and pst
    if min_cnf_locs is not None:

        # Obtain geom and freqs from filesys
        geo = cnf_fs[-1].file.hessian.read(min_cnf_locs)
        hess = cnf_fs[-1].file.hessian.read(min_cnf_locs)

        # Obtain the frequencies
        freqs, _, imag, _ = structure.vib.projrot_freqs(
            geo, hess, thy_info, thy_run_fs,
            grad=(), rotors_str='', coord_proj='cartesian',
            script_str=DEFAULT_SCRIPT_DCT['projrot'])

        # Calculate the zpve
        zpe = (sum(freqs) / 2.0) * phycon.WAVEN2KCAL

    else:
        print('ERROR: Reference geometry is missing for harmonic frequencies')

    return freqs, imag, zpe


def read_anharmon_matrix(pf_filesystems):
    """ Read the anharmonicity matrix """

    # Set up vpt2 level filesystem for rotational values
    [cnf_fs, cnf_path, min_cnf_locs, _h] = pf_filesystems['vpt2']

    if min_cnf_locs is not None:
        xmat = cnf_fs[-1].file.anharmonicity_matrix.read(
            min_cnf_locs)
        print('Anharm matrix at path')
        print(cnf_path)
    else:
        print('No anharm matrix at path')
        print(cnf_path)

    return xmat


def tors_projected_freqs_zpe(pf_filesystems, mess_hr_str, projrot_hr_str,
                             save_path, saddle=False):
    """ Get frequencies from one version of ProjRot
    """

    # Build the filesystems
    [harm_cnf_fs, _, harm_min_cnf_locs, _] = pf_filesystems['harm']
    [tors_cnf_fs, tors_path, tors_min_cnf_locs, _] = pf_filesystems['tors']

    # Read info from the filesystem that is needed
    harm_geo = harm_cnf_fs[-1].file.geometry.read(harm_min_cnf_locs)
    hess = harm_cnf_fs[-1].file.hessian.read(harm_min_cnf_locs)
    tors_geo = tors_cnf_fs[-1].file.geometry.read(tors_min_cnf_locs)

    # Calculate ZPVES of the hindered rotors
    tors_zpe = torsprep.calc_tors_freqs_zpe(tors_geo, mess_hr_str, tors_path)

    # Run ProjRot to get the frequencies v1
    rt_freqs1, rth_freqs1, _, rth_imag1 = vibprep.projrot_freqs(
        harm_geo, hess, thy_info, thy_run_fs,
        grad=(), rotors_str=projrot_hr_str, coord_proj='cartesian',
        script_str=DEFAULT_SCRIPT_DCT['projrot'])

    # Run ProjRot to get the frequencies v2
    projrot_script_str2 = (
        "#!/usr/bin/env bash\n"
        "RPHt2.exe >& /dev/null"
    )
    _, rth_freqs2, _, rth_imag2 = vibprep.projrot_freqs(
        harm_geo, hess, thy_info, thy_run_fs,
        grad=(), rotors_str=projrot_hr_str, coord_proj='cartesian',
        script_str=projrot_script_str2)

    # Calculate harmonic ZPVE from all harmonic freqs, including torsionals
    harm_zpe = (sum(rt_freqs1) / 2.0) * phycon.WAVEN2KCAL

    # Calculate harmonic ZPVE from freqs where torsions have been projected out
    # Value from both projrot versions, which use different projection schemes
    harm_zpe_notors_1 = (sum(rth_freqs1) / 2.0) * phycon.WAVEN2KCAL
    harm_zpe_notors_2 = (sum(rth_freqs2) / 2.0) * phycon.WAVEN2KCAL

    # Calcuate the difference in the harmonic ZPVE from projecting out torsions
    harm_tors_zpe = harm_zpe - harm_zpe_notors_1
    harm_tors_zpe_2 = harm_zpe - harm_zpe_notors_2

    # Check to see which of the above ZPVEs match more closely with tors ZPVE
    # calculated directly by treating the torsions in MESS
    diff_tors_zpe = harm_tors_zpe - tors_zpe
    diff_tors_zpe_2 = harm_tors_zpe_2 - tors_zpe
    if diff_tors_zpe <= diff_tors_zpe_2:
        freqs = rth_freqs1
        imag = rth_imag1
        zpe = harm_zpe_notors_1 + tors_zpe
    else:
        freqs = rth_freqs2
        imag = rth_imag2
        zpe = harm_zpe_notors_2 + tors_zpe

    # Check if there are significant differences caused by the rotor projection
    if abs(diff_tors_zpe) > 0.2 and abs(diff_tors_zpe_2) > 0.2:
        print('Warning: There is a difference of ',
              '{0:.2f} and {1:.2f}'.format(diff_tors_zpe, diff_tors_zpe_2),
              'kcal/mol between harmonic and hindered torsional ZPVEs')

    return freqs, imag, zpe
