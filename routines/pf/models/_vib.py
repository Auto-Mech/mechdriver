"""
  Handle vibrational data info
"""

import numpy
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
    [cnf_fs, cnf_path, min_cnf_locs, _, _] = pf_filesystems['vpt2']

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
    [harm_cnf_fs, _, harm_min_locs, _, harm_run_fs] = pf_filesystems['harm']
    [tors_cnf_fs, _, tors_min_locs, _, tors_run_fs] = pf_filesystems['tors']

    # Build the run filesystem using locs
    harm_run_fs[-1].create(harm_min_locs)
    harm_run_path = harm_run_fs[-1].path(harm_min_locs)
    tors_run_fs[-1].create(tors_min_locs)
    tors_run_path = tors_run_fs[-1].path(tors_min_locs)

    # Read info from the filesystem that is needed
    harm_geo = harm_cnf_fs[-1].file.geometry.read(harm_min_locs)
    hess = harm_cnf_fs[-1].file.hessian.read(harm_min_locs)
    tors_geo = tors_cnf_fs[-1].file.geometry.read(tors_min_locs)
    hess_path = harm_cnf_fs[-1].path(harm_min_locs)
    print(' - Reading Hessian from path {}'.format(hess_path))

    # Read info for the hindered rotors
    print(' - Calculating the torsional ZPVES using MESS...')
    tors_zpes, tors_freqs = torsprep.mess_tors_zpes(
        tors_geo, mess_hr_str, tors_run_path)

    # Calculate ZPVES of the hindered rotors
    tors_zpe = sum(tors_zpes) if tors_zpes else 0.0
    tors_zpe *= phycon.KCAL2EH

    print(' - Calculating the RT and RT-rotor projected frequencies ProjRot')
    # Run ProjRot to get the frequencies v1
    rt_freqs1, rth_freqs1, rt_imag1, _ = vibprep.projrot_freqs(
        [harm_geo], [hess], harm_run_path,
        grads=[[]], rotors_str=projrot_hr_str, coord_proj='cartesian',
        script_str=DEFAULT_SCRIPT_DCT['projrot'])

    # Run ProjRot to get the frequencies v2
    projrot_script_str2 = (
        "#!/usr/bin/env bash\n"
        "RPHt2.exe >& /dev/null"
    )
    _, rth_freqs2, rt_imag2, _ = vibprep.projrot_freqs(
        [harm_geo], [hess], harm_run_path,
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
        imag_freqs = rt_imag1
        proj_zpe = harm_zpe_notors_1
    else:
        freqs = rth_freqs2
        imag_freqs = rt_imag2
        proj_zpe = harm_zpe_notors_2

    # Check imaginary frequencies and set freqs
    if saddle:
        if len(imag_freqs) > 1:
            print('There is more than one imaginary frequency')
        imag = max(imag_freqs)
    else:
        imag = None

    # Create a scaling factor for the frequencie
    log_rt_freq = [numpy.log(freq) for freq in rt_freqs1]
    log_rt_freq = sum(log_rt_freq)
    log_freq = [numpy.log(freq) for freq in freqs]
    log_freq = sum(log_freq)
    log_tors_freq = [numpy.log(freq) for freq in tors_freqs]
    log_tors_freq = sum(log_tors_freq)
    #unproj_prod = numpy.prod(rt_freqs1)
    #proj_prod = numpy.prod(freqs) * numpy.prod(tors_freqs)
    #print('proj_prod test:', unproj_prod, proj_prod)
    print('log_freq_tests:', log_rt_freq, log_freq, log_tors_freq)
    print('freq test:', freqs, tors_freqs, rt_freqs1)
    #scale_factor = unproj_prod / proj_prod
    scale_factor = numpy.exp(log_rt_freq - log_freq - log_tors_freq)
    print('scale fact test', scale_factor)

    # Check if there are significant differences caused by the rotor projection
    diff_tors_zpe *= phycon.EH2KCAL
    diff_tors_zpe_2 *= phycon.EH2KCAL
    if abs(diff_tors_zpe) > 0.2 and abs(diff_tors_zpe_2) > 0.2:
        print('Warning: There is a difference of ',
              '{0:.2f} and {1:.2f}'.format(diff_tors_zpe, diff_tors_zpe_2),
              'kcal/mol between harmonic and hindered torsional ZPVEs')

    return freqs, imag, tors_zpe, scale_factor


M3_COEFFS = {
    ('b2plypd3', 'cc-pvtz'): (1.066, 0.008045, 0.33),
    ('wb97xd', '6-31g*'): (1.657244, 0.56000691, 0.029624)
}


def scale_frequencies(freqs, tors_zpe,
                      chn_pf_levels, scale_method='3c'):
    """ Scale frequencies according to some method
        obtain a corrected zpe
    """

    # Scale the frequencies
    print('in scale freqs')
    thy_info = chn_pf_levels['harm'][1]
    thy_method = (thy_info[1], thy_info[2])
    print('thy_method', thy_method)
    if scale_method == '3c':
        cf1, cf2, cf3 = M3_COEFFS.get(thy_method, (1.0, 0.0, 0.0))
        scaled_freqs = []
        for freq in freqs:
            scale_factor = cf1 - (cf2 * freq**cf3)
            scaled_freqs.append(freq * scale_factor)
        print('coef2', cf1, cf2, cf3)
    else:
        scaled_freqs = freqs
    print('freqs', freqs)
    print('scaled_freqs', scaled_freqs)

    # Calculate the zpe
    freq_zpe = 0.0
    for freq, scfreq in zip(freqs, scaled_freqs):
        freq_zpe += ((freq / 2.0) - (1.0 / 8.0) * (scfreq - freq))
    freq_zpe *= phycon.WAVEN2EH
    tot_zpe = freq_zpe + tors_zpe

    return scaled_freqs, tot_zpe
