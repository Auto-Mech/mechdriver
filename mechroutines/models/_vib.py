"""
  Handle vibrational data info
"""

import os
import autorun
import automol.pot
import automol.geom
import autofile.fs
from phydat import phycon
from mechlib.amech_io import printer as ioprinter
from mechlib.amech_io._path import job_path
from mechroutines.models import typ
from mechroutines.models import _tors as tors


def full_vib_analysis(
        spc_dct_i, pf_filesystems, spc_mod_dct_i,
        run_prefix, zrxn=None):
    """ process to get freq
    """
    # Pack into big object to pass into functions and return
    tors_strs = ['']
    freqs = []
    imag = []
    tors_zpe = 0.0
    pot_scalef = 1.0
    tors_freqs = []
    harm_freqs = []

    rotors, mdhr_dct = tors.build_rotors(
        spc_dct_i, pf_filesystems, spc_mod_dct_i)
    # Squash the rotor potentials as necessary
    if rotors is not None:
        if typ.squash_tors_pot(spc_mod_dct_i):
            for rotor in rotors:
                for torsion in rotor:
                    torsion.pot = automol.pot.relax_scale(torsion.pot)
    if typ.nonrigid_tors(spc_mod_dct_i, rotors):

        # Build initial MESS+ProjRot HindRot strings; calc. projected freq info
        tors_strs = tors.make_hr_strings(rotors, mdhr_dct=mdhr_dct)
        [_, hr_str, _, prot_str, _] = tors_strs
        ret = tors_projected_freqs(
            pf_filesystems, hr_str, prot_str, run_prefix, zrxn=zrxn)

        if ret is not None:
            proj_freqs, harm_freqs, tors_freqs, imag, disps = ret

            # Make final hindered rotor strings and get corrected tors zpe
            if typ.scale_1d(spc_mod_dct_i):

                scaled_proj_freqs, _ = scale_frequencies(
                    proj_freqs, None, spc_mod_dct_i,
                    scale_method='c3_harm')
                scaled_harm_freqs, _ = scale_frequencies(
                    harm_freqs, None, spc_mod_dct_i,
                    scale_method='c3_harm')
                # print('scaling test:', scaled_harm_freqs,
                #       scaled_proj_freqs, tors_freqs)
                # print('tors string before scaling',tors_strs)
                pot_scalef = potential_scale_factor(
                    scaled_harm_freqs, scaled_proj_freqs, tors_freqs)
                # print('pot_scalef:', pot_scalef)
                # get the pot scale factor
                # print('before scaling')
                # for rotor in rotors:
                #     for _tors in rotor:
                #         print(_tors.pot)
                rotors = tors.scale_rotor_pots(rotors, scale_factor=pot_scalef)
                # print('after scaling')
                # for rotor in rotors:
                #     for _tors in rotor:
                #         print(_tors.pot)
                tors_strs = tors.make_hr_strings(rotors, mdhr_dct=mdhr_dct)
                [_, hr_str, _, prot_str, _] = tors_strs
                # print('tors string after scaling',tors_strs)
                tors_zpe = tors_projected_scaled_zpe(
                    pf_filesystems, hr_str, prot_str, run_prefix,
                    spc_mod_dct_i, zrxn=zrxn)
            freqs = ret[0]  # freqs equal to proj_freqs

        # For mdhrv model no freqs needed in MESS input, zero out freqs lst
        if 'mdhrv' in spc_mod_dct_i['tors']['mod']:
            freqs = ()
    else:
        ret = read_harmonic_freqs(
            pf_filesystems, run_prefix, zrxn=zrxn)
        freqs, imag, zpe, disps = ret

    if freqs:
        freqs, zpe = scale_frequencies(
            freqs, tors_zpe,
            spc_mod_dct_i, scale_method='c3')
        harm_freqs, _ = scale_frequencies(
            harm_freqs, None,
            spc_mod_dct_i, scale_method='c3')

        # harm_freqs, zpe = scale_frequencies(
        #     harm_freqs, 0,
        #    spc_mod_dct_i, scale_method='c3')

    return (freqs, imag, zpe, pot_scalef, tors_strs, tors_freqs,
            harm_freqs, disps)


# Read the Frequencies from the filesystem using the fs objs
def read_harmonic_freqs(pf_filesystems, run_prefix, zrxn=None):
    """ Read the harmonic frequencies for the minimum
        energy conformer
    """
    # Get the harmonic filesys information
    [cnf_fs, _, min_cnf_locs, _, _] = pf_filesystems['harm']
    return read_locs_harmonic_freqs(
        cnf_fs, min_cnf_locs, run_prefix, zrxn=zrxn)


def read_locs_harmonic_freqs(cnf_fs, cnf_locs, run_prefix, zrxn=None):
    """ Read the harmonic frequencies for a specific conformer
        Do the freqs obtain for two species for fake and pst?
    """

    if cnf_locs is not None:
        geo_exists = cnf_fs[-1].file.geometry.exists(cnf_locs)
        hess_exists = cnf_fs[-1].file.hessian.exists(cnf_locs)

        if not geo_exists:
            ioprinter.error_message(
                'No Reference geometry for harmonic frequencies at path',
                cnf_fs[-1].file.hessian.path(cnf_locs))
        if not hess_exists:
            ioprinter.error_message(
                'No Hessian available for harmonic frequencies at path',
                cnf_fs[-1].file.hessian.path(cnf_locs))
    else:
        geo_exists, hess_exists = False, False

    if geo_exists and hess_exists:
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

        # Obtain the displacements
        norm_coord_str, _ = autorun.projrot.displacements(
            script_str, vib_path, [geo], [[]], [hess])

        # Calculate the zpve
        ioprinter.frequencies(freqs)
        zpe = (sum(freqs) / 2.0) * phycon.WAVEN2EH

        # Check imaginary frequencies and set freqs
        if zrxn is not None:
            if len(imag_freqs) > 1:
                ioprinter.warning_message(
                    'Saddle Point has more than',
                    'one imaginary frequency')
            imag = imag_freqs[0]
        else:
            imag = None

        ret = (freqs, imag, zpe, norm_coord_str)
    else:
        ret = None

    return ret


def read_anharmon_matrix(pf_filesystems):
    """ Read a anharmonicity matrix from the SAVE filesystem for a
        species or transition state.
    """

    # Set up vpt2 level filesystem for rotational values
    [cnf_fs, cnf_path, min_cnf_locs, _, _] = pf_filesystems['vpt2']

    if cnf_path:
        xmat = cnf_fs[-1].file.anharmonicity_matrix.read(
            min_cnf_locs)
        ioprinter.reading('Anharm matrix', cnf_path)
    else:
        ioprinter.error_message('No anharm matrix at path', cnf_path)

    return xmat


def tors_projected_freqs(pf_filesystems, mess_hr_str, projrot_hr_str,
                         prefix, zrxn=None, conf=None):
    """ Get the projected frequencies from harmonic frequencies,
        which requires projrot run

        :param pf_filesytem: dictionary of the locations of various info
        :param runf_pfx: location to run projrot
    """
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
    vib_path = job_path(prefix, 'PROJROT', 'FREQ', fml_str, print_path=True)
    # print('proj test:', vib_path)
    # tors_path = job_path(run_pfx, 'MESS', 'TORS', fml_str, print_path=True)
    mess_script_str = autorun.SCRIPT_DCT['messpf']
    projrot_script_str = autorun.SCRIPT_DCT['projrot']
    dist_cutoff_dct1 = {('H', 'O'): 2.26767, ('H', 'C'): 2.26767}
    dist_cutoff_dct2 = {('H', 'O'): 2.83459, ('H', 'C'): 2.83459,
                        ('C', 'O'): 3.7807}
    proj_inf = autorun.projected_frequencies(
        mess_script_str, projrot_script_str, vib_path,
        mess_hr_str, projrot_hr_str,
        tors_geo, harm_geo, hess,
        dist_cutoff_dct1=dist_cutoff_dct1,
        dist_cutoff_dct2=dist_cutoff_dct2,
        saddle=(zrxn is not None))

    # Obtain the displacements
    disp_path = os.path.join(vib_path, 'DISP')
    harm_disps = autorun.projrot.displacements(
        projrot_script_str, disp_path, [harm_geo], [[]], [hess])

    proj_freqs, proj_imag, _, harm_freqs, tors_freqs = proj_inf

    return proj_freqs, harm_freqs, tors_freqs, proj_imag, harm_disps


def tors_projected_scaled_zpe(pf_filesystems, mess_hr_str, projrot_hr_str,
                             prefix, spc_mod_dct_i, zrxn=None, conf=None):
    """ Get frequencies from one version of ProjRot
    """
    ret = tors_projected_freqs(
        pf_filesystems, mess_hr_str, projrot_hr_str,
        prefix, zrxn=zrxn, conf=conf)
    _, _, tors_freqs, _, _ = ret
    scaled_tors_freqs, scaled_tors_zpe = scale_frequencies(
        tors_freqs, 0.0,
        spc_mod_dct_i, scale_method='c3')
    # tors_zpe = 0.5 * sum(tors_freqs) * phycon.WAVEN2EH

    return scaled_tors_zpe


def potential_scale_factor(harm_freqs, proj_freqs, tors_freqs):
    """ Get frequencies from one version of ProjRot
    """
    return automol.prop.freq.rotor_scale_factor_from_harmonics(
        harm_freqs, proj_freqs, tors_freqs)


def scale_frequencies(freqs, tors_zpe,
                      spc_mod_dct_i, scale_method='c3'):
    """ Empirically scale the harmonic vibrational frequencies and harmonic
        zero-point energy (ZPVE) of a species or transition to account
        for anharmonic effects. The final ZPVE value also includes the ZPVEs
        of internal rotations which are not scaled.

        Scaling factors determined by the electronic structure
        method used to calculate the frequencies and ZPVE as well as the
        requested scaling method.

        :param freqs: harmonic frequencies [cm-1]
        :type freqs: tuple(float)
        :param tors_zpe:

    """

    thy_info = spc_mod_dct_i['vib']['geolvl'][1][1]
    method, basis = thy_info[1], thy_info[2]
    if tors_zpe is None:
        ret = automol.prop.freq.scale_frequencies(
            freqs, method, basis, scale_method=scale_method)
        scaled_freqs = ret
        tot_zpe = None
    else:
        ret = automol.prop.freq.scale_frequencies_and_zpe(
            freqs, method, basis, scale_method=scale_method)
        scaled_freqs, scaled_zpe = ret
        tot_zpe = scaled_zpe + tors_zpe

    return scaled_freqs, tot_zpe


def ted(spc_dct_i, pf_filesystems, spc_mod_dct_i,
        run_prefix, zrxn=None):
    """ process to get freq
    """
    [cnf_fs, _, min_cnf_locs, _, _] = pf_filesystems['harm']

    if min_cnf_locs is not None:
        geo_exists = cnf_fs[-1].file.geometry.exists(min_cnf_locs)
        hess_exists = cnf_fs[-1].file.hessian.exists(min_cnf_locs)

        if not geo_exists:
            ioprinter.error_message(
                'No Reference geometry for harmonic frequencies at path',
                cnf_fs[-1].file.hessian.path(min_cnf_locs))
        if not hess_exists:
            ioprinter.error_message(
                'No Hessian available for harmonic frequencies at path',
                cnf_fs[-1].file.hessian.path(min_cnf_locs))
    else:
        geo_exists, hess_exists = False, False

    if geo_exists and hess_exists:

        geo = cnf_fs[-1].file.geometry.read(min_cnf_locs)
        hess = cnf_fs[-1].file.hessian.read(min_cnf_locs)

        geo_path = cnf_fs[-1].path(min_cnf_locs)
        zma_fs = autofile.fs.zmatrix(geo_path)
        zma = zma_fs[-1].file.zmatrix.read((0,))

        fml_str = automol.geom.formula_string(geo)
        ted_path = job_path(run_prefix, 'INTDER', 'TED', fml_str)

        ioprinter.running(
            f'Reaction Coordinate Check with INTDER at {ted_path}')
        script_str = autorun.SCRIPT_DCT['intder']
        _ = autorun.intder.reaction_coordinate_check_idxs(
            script_str, ted_path, geo, zma, hess, zrxn)
