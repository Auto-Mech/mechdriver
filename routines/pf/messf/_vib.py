"""
  Handle vibrational data info
"""

import autofile
from lib.structure import tors as torsprep
from lib.structure import vib as vibprep
from lib.phydat import phycon
from lib.submission import DEFAULT_SCRIPT_DCT


def read_harmonic_freqs(pf_filesystems, saddle=False):
    """ Read the harmonic frequencies
    """

    # Get the harmonic filesys information
    [cnf_fs, harm_path, min_cnf_locs, _, run_path] = pf_filesystems['harm']

    # probably should read freqs
    # Do the freqs obtain for two species for fake and pst
    if min_cnf_locs is not None:

        # Obtain geom and freqs from filesys
        geo = cnf_fs[-1].file.geometry.read(min_cnf_locs)
        hess = cnf_fs[-1].file.hessian.read(min_cnf_locs)
        hess_path = cnf_fs[-1].path(min_cnf_locs)
        print(' - Reading Hessian from path {}'.format(hess_path))

        # Build the run filesystem using locs
        run_fs = autofile.fs.conformer(harm_path)
        run_fs[-1].create(min_cnf_locs)
        run_path = run_fs[-1].path(min_cnf_locs)

        # Obtain the frequencies
        print('Calling ProjRot to diagonalize Hessian and get freqs...')
        freqs, _, imag_freqs, _ = vibprep.projrot_freqs(
            [geo], [hess], run_path,
            grads=[[]], rotors_str='', coord_proj='cartesian',
            script_str=DEFAULT_SCRIPT_DCT['projrot'])

        # Calculate the zpve
        print('Calculating ZPVE using harmonic frequencies...')
        zpe = (sum(freqs) / 2.0) * phycon.WAVEN2EH

        # Check imaginary frequencies and set freqs
        if saddle:
            if len(imag_freqs) > 1:
                print('WARNING: Saddle Point has more than',
                      'one imaginary frequency')
            imag = imag_freqs[0]
        else:
            imag = None

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
                             saddle=False):
    """ Get frequencies from one version of ProjRot
    """

    # Build the filesystems
    [harm_cnf_fs, harm_path, harm_min_cnf_locs, _, _] = pf_filesystems['harm']
    [tors_cnf_fs, tors_path, tors_min_cnf_locs, _, _] = pf_filesystems['tors']

    # Build the run filesystem using locs
    run_fs = autofile.fs.conformer(harm_path)
    run_fs[-1].create(harm_min_cnf_locs)
    run_path = run_fs[-1].path(harm_min_cnf_locs)

    # Read info from the filesystem that is needed
    harm_geo = harm_cnf_fs[-1].file.geometry.read(harm_min_cnf_locs)
    hess = harm_cnf_fs[-1].file.hessian.read(harm_min_cnf_locs)
    tors_geo = tors_cnf_fs[-1].file.geometry.read(tors_min_cnf_locs)
    hess_path = harm_cnf_fs[-1].path(harm_min_cnf_locs)
    print(' - Reading Hessian from path {}'.format(hess_path))

    # Calculate ZPVES of the hindered rotors
    print(' - Calculating the torsional ZPVES using MESS') 
    tors_zpe = torsprep.mess_tors_zpes(tors_geo, mess_hr_str, tors_path)
    tors_zpe *= phycon.KCAL2EH

    # Run ProjRot to get the frequencies v1
    rt_freqs1, rth_freqs1, _, rth_imag1 = vibprep.projrot_freqs(
        [harm_geo], [hess], run_path,
        grads=[[]], rotors_str=projrot_hr_str, coord_proj='cartesian',
        script_str=DEFAULT_SCRIPT_DCT['projrot'])

    # Run ProjRot to get the frequencies v2
    projrot_script_str2 = (
        "#!/usr/bin/env bash\n"
        "RPHt2.exe >& /dev/null"
    )
    _, rth_freqs2, _, rth_imag2 = vibprep.projrot_freqs(
        [harm_geo], [hess], run_path,
        grads=[[]], rotors_str=projrot_hr_str, coord_proj='cartesian',
        script_str=projrot_script_str2)

    # Calculate harmonic ZPVE from all harmonic freqs, including torsionals
    harm_zpe = (sum(rt_freqs1) / 2.0) * phycon.WAVEN2EH

    # Calculate harmonic ZPVE from freqs where torsions have been projected out
    # Value from both projrot versions, which use different projection schemes
    harm_zpe_notors_1 = (sum(rth_freqs1) / 2.0) * phycon.WAVEN2EH
    harm_zpe_notors_2 = (sum(rth_freqs2) / 2.0) * phycon.WAVEN2EH

    # Calcuate the difference in the harmonic ZPVE from projecting out torsions
    harm_tors_zpe = harm_zpe - harm_zpe_notors_1
    harm_tors_zpe_2 = harm_zpe - harm_zpe_notors_2

    # Check to see which of the above ZPVEs match more closely with tors ZPVE
    # calculated directly by treating the torsions in MESS
    diff_tors_zpe = harm_tors_zpe - tors_zpe
    diff_tors_zpe_2 = harm_tors_zpe_2 - tors_zpe
    if diff_tors_zpe <= diff_tors_zpe_2:
        freqs = rth_freqs1
        imag_freqs = rth_imag1
        zpe = harm_zpe_notors_1 + tors_zpe
    else:
        freqs = rth_freqs2
        imag_freqs = rth_imag2
        zpe = harm_zpe_notors_2 + tors_zpe

    # Check imaginary frequencies and set freqs
    if saddle:
        if len(imag_freqs) > 1:
            print('There is more than one imaginary frequency')
        imag = imag_freqs[0]
    else:
        imag = None

    # Check if there are significant differences caused by the rotor projection
    diff_tors_zpe *= phycon.EH2KCAL
    diff_tors_zpe_2 *= phycon.EH2KCAL
    if abs(diff_tors_zpe) > 0.2 and abs(diff_tors_zpe_2) > 0.2:
        print('Warning: There is a difference of ',
              '{0:.2f} and {1:.2f}'.format(diff_tors_zpe, diff_tors_zpe_2),
              'kcal/mol between harmonic and hindered torsional ZPVEs')

    return freqs, imag, zpe
