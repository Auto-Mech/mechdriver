"""
  Handle vibrational data info
"""

import os
import numpy
import autofile
import autorun
import automol.geom
import projrot_io
from phydat import phycon
from mechlib.amech_io import printer as ioprinter
# from mechlib.amech_io import job_path
from mechroutines.pf.models import typ
from mechroutines.pf.models import _tors as tors


def vib_analysis(spc_dct_i, pf_filesystems, chn_pf_models, pf_levels,
                 run_prefix, saddle=False):
    """ process to get freq
    """

    tors_strs = ['']

    rotors = tors.build_rotors(
        spc_dct_i, pf_filesystems, chn_pf_models, pf_levels)

    if typ.nonrigid_tors(chn_pf_models, rotors):
        tors_strs = tors.make_hr_strings(rotors)
        [_, hr_str, _, prot_str, _] = tors_strs

        freqs, imag, tors_zpe, pot_scalef = tors_projected_freqs_zpe(
            pf_filesystems, hr_str, prot_str, run_prefix, saddle=saddle)

        # Make final hindered rotor strings and get corrected tors zpe
        if typ.scale_1d(chn_pf_models):
            rotors = tors.scale_rotor_pots(rotors, scale_factor=pot_scalef)
            tors_strs = tors.make_hr_strings(rotors)
            [_, hr_str, _, prot_str, _] = tors_strs
            _, _, tors_zpe, _ = tors_projected_freqs_zpe(
                pf_filesystems, hr_str, prot_str, run_prefix, saddle=saddle)
            # Calculate current zpe assuming no freq scaling: tors+projfreq

        zpe = tors_zpe + (sum(freqs) / 2.0) * phycon.WAVEN2EH

        # For mdhrv model no freqs needed in MESS input, zero out freqs lst
        if 'mdhrv' in chn_pf_models['tors']:
            freqs = ()
    else:
        freqs, imag, zpe = read_harmonic_freqs(
            pf_filesystems, run_prefix, saddle=saddle)
        tors_zpe = 0.0

    return freqs, imag, zpe, tors_strs


def full_vib_analysis(spc_dct_i, pf_filesystems, chn_pf_models, pf_levels,
                      run_prefix, saddle=False):
    """ process to get freq
    """

    tors_strs = ['']
    freqs = []
    imag = []
    tors_zpe = 0.0
    scale_factor = 1.0
    tors_freqs = []
    rt_freqs1 = []

    rotors = tors.build_rotors(
        spc_dct_i, pf_filesystems, chn_pf_models, pf_levels)

    if typ.nonrigid_tors(chn_pf_models, rotors):
        tors_strs = tors.make_hr_strings(rotors)
        [_, hr_str, _, prot_str, _] = tors_strs

        freqs, imag, tors_zpe, pot_scalef = tors_projected_freqs_zpe(
            pf_filesystems, hr_str, prot_str, run_prefix, saddle=saddle)

        # Make final hindered rotor strings and get corrected tors zpe
        if typ.scale_1d(chn_pf_models):
            rotors = tors.scale_rotor_pots(rotors, scale_factor=pot_scalef)
            tors_strs = tors.make_hr_strings(rotors)
            [_, hr_str, _, prot_str, _] = tors_strs
            freqs, imag, tors_zpe, scale_factor, tors_freqs, rt_freqs1 = tors_projected_freqs(
                pf_filesystems, hr_str, prot_str, run_prefix, saddle=saddle)
            # Calculate current zpe assuming no freq scaling: tors+projfreq

        zpe = tors_zpe + (sum(freqs) / 2.0) * phycon.WAVEN2EH

        # For mdhrv model no freqs needed in MESS input, zero out freqs lst
        if 'mdhrv' in chn_pf_models['tors']:
            freqs = ()
    else:
        freqs, imag, zpe = read_harmonic_freqs(
            pf_filesystems, run_prefix, saddle=saddle)
        tors_zpe = 0.0

    return freqs, imag, tors_zpe, scale_factor, tors_freqs, rt_freqs1


def read_harmonic_freqs(pf_filesystems, run_prefix, saddle=False):
    """ Read the harmonic frequencies for the minimum
        energy conformer
    """
    # Get the harmonic filesys information
    [cnf_fs, harm_path, min_cnf_locs, _, _] = pf_filesystems['harm']
    freqs, imag, zpe = read_locs_harmonic_freqs(
        cnf_fs, harm_path, min_cnf_locs, run_prefix, saddle=saddle)
    return freqs, imag, zpe


def read_locs_harmonic_freqs(
        cnf_fs, harm_path, cnf_locs, run_prefix, saddle=False):
    """ Read the harmonic frequencies for a specific conformer
    """

    # probably should read freqs
    # Do the freqs obtain for two species for fake and pst
    if cnf_locs is not None:

        # Obtain geom and freqs from filesys
        geo = cnf_fs[-1].file.geometry.read(cnf_locs)
        hess = cnf_fs[-1].file.hessian.read(cnf_locs)
        ioprinter.reading('Hessian', cnf_fs[-1].path(cnf_locs))

        # Build the run filesystem using locs
        fml_str = automol.geom.formula_string(geo)
        vib_path = job_path(run_prefix, 'PROJROT', 'FREQ', fml_str)

        # Obtain the frequencies
        ioprinter.info_message(
            'Calling ProjRot to diagonalize Hessian and get freqs...')
        script_str = autorun.SCRIPT_DCT['projrot']
        freqs, _, imag_freqs, _ = autorun.projrot.frequencies(
            script_str, vib_path, [geo], [[]], [hess])

        # Calculate the zpve
        ioprinter.frequencies(freqs)
        zpe = (sum(freqs) / 2.0) * phycon.WAVEN2EH

        # Check imaginary frequencies and set freqs
        if saddle:
            if len(imag_freqs) > 1:
                ioprinter.warning_message(
                    'Saddle Point has more than',
                    'one imaginary frequency')
            imag = imag_freqs[0]
        else:
            imag = None

    else:
        ioprinter.error_message(
            'Reference geometry is missing for harmonic frequencies')

    return freqs, imag, zpe


def read_anharmon_matrix(pf_filesystems):
    """ Read the anharmonicity matrix """

    # Set up vpt2 level filesystem for rotational values
    [cnf_fs, cnf_path, min_cnf_locs, _, _] = pf_filesystems['vpt2']

    if cnf_path:
        xmat = cnf_fs[-1].file.anharmonicity_matrix.read(
            min_cnf_locs)
        ioprinter.reading('Anharm matrix', cnf_path)
    else:
        ioprinter.error_message('No anharm matrix at path', cnf_path)

    return xmat


def tors_projected_freqs_zpe(pf_filesystems, mess_hr_str, projrot_hr_str,
                             prefix, saddle=False, conf=None):
    """ Get frequencies from one version of ProjRot
    """
    ret = tors_projected_freqs(
        pf_filesystems, mess_hr_str, projrot_hr_str,
        prefix, saddle=saddle, conf=conf)
    freqs, imag, tors_zpe, scale_factor, _, _ = ret

    return freqs, imag, tors_zpe, scale_factor


def tors_projected_freqs(pf_filesystems, mess_hr_str, projrot_hr_str,
                         prefix, saddle=False, conf=None):
    """ Get frequencies from one version of ProjRot
    """
    run_prefix = pf_filesystems['run_prefix']

    # Build the filesystems
    [harm_cnf_fs, _, harm_min_locs, _, _] = pf_filesystems['harm']
    [tors_cnf_fs, _, tors_min_locs, _, _] = pf_filesystems['tors']
    if conf:
        harm_min_locs = conf[1]
        harm_cnf_fs = conf[2]

    # Read info from the filesystem that is needed
    harm_geo = harm_cnf_fs[-1].file.geometry.read(harm_min_locs)
    hess = harm_cnf_fs[-1].file.hessian.read(harm_min_locs)
    tors_geo = tors_cnf_fs[-1].file.geometry.read(tors_min_locs)
    ioprinter.reading('Hessian', harm_cnf_fs[-1].path(harm_min_locs))

    fml_str = automol.geom.formula_string(harm_geo)
    vib_path = job_path(
        run_prefix, 'PROJROT', 'PROJFREQ', fml_str, print_path=True)

    # Read info for the hindered rotors and calculate the ZPVE
    ioprinter.info_message(' - Calculating the torsional ZPVES using MESS...')

    # NEW autorun function for the frequencies
    mess_script_str = autorun.SCRIPT_DCT['messpf']
    projrot_script_str = autorun.SCRIPT_DCT['projrot']

    proj_freqs, proj_imags, proj_zpe, harm_freqs, tors_freqs = autorun.projected_frequencies(
        mess_script_str, projrot_script_str, vib_path,
        mess_hr_str, projrot_hr_str,
        tors_geo, harm_geo, hess)
    if saddle:
        proj_imag = proj_imags[0]
    else:
        proj_imag = []

    # NEW scale factor functions
    scale_factor = automol.prop.freq.rotor_scale_factor_from_harmonics(
        harm_freqs, proj_freqs, tors_freqs)

    return (proj_freqs, proj_imag, proj_zpe,
            scale_factor, tors_freqs, harm_freqs)


def scale_frequencies(freqs, tors_zpe,
                      chn_pf_levels, scale_method='3c'):
    """ Scale frequencies according to some method
        obtain a corrected zpe
    """

    thy_info = chn_pf_levels['harm'][1]
    method, basis = thy_info[1], thy_info[2]
    scaled_freqs, scaled_zpe = automol.prop.freq.anharm_by_scaling(
        freqs, method, basis, scale_method='c3')
    tot_zpe = scaled_zpe + tors_zpe

    return scaled_freqs, tot_zpe
