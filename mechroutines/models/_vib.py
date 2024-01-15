"""
  Handle vibrational data info
"""

import os
from copy import deepcopy
import autorun
import automol.pot
import automol.geom
import autofile.fs
from phydat import phycon
from mechlib.reaction import _util as rxn_util
from mechlib.amech_io import printer as ioprinter
from mechlib.reaction import _util as rxn_util
from mechlib.amech_io._path import job_path
from mechroutines.models import typ
from mechroutines.models import _tors as tors
from mechroutines.models import _rot as rot
from copy import deepcopy
import numpy
import itertools


def full_vib_analysis(
        spc_dct_i, pf_filesystems, spc_mod_dct_i,
        run_prefix, zrxn=None):
    """ process to get freq
    """
    # Pack into big object to pass into functions and return
    tors_strs = ['', '', '', '', '']
    imag = []           # harmonic imaginary frequency
    pot_scalef = 1.0    # scaling applied to torsional potentials
    unproj_hfreqs = []  # all vibrational and torsional harmonic freqs
    proj_hfreqs = []    # only vibrational  harmonic frequencies
    proj_ffreqs = []    # all vibrational and torsional vpt2/scaled freqs
    unproj_ffreqs = []  # only vibrational vpt2/scaled frequencies
    tors_freqs = []     # torsional oscillator frequencies
    total_zpe = None    # vibrational + torsional zero point
    tors_zpe = 0.0      # torsional zero point
    disps = []
    xmat = None
    rovib_mat = None
    rot_dists = None

    rotors, mdhr_dct, zma_locs = tors.build_rotors(
        spc_dct_i, pf_filesystems, spc_mod_dct_i)

    if typ.nonrigid_tors(spc_mod_dct_i, rotors):
        # Build initial MESS+ProjRot HindRot strings; calc. projected freq info
        tors_strs = tors.make_hr_strings(rotors, mdhr_dct=mdhr_dct)
        [_, hr_str, _, prot_str, _] = tors_strs
        ret = tors_projected_freqs(
            pf_filesystems, hr_str, prot_str, run_prefix,
            zrxn=zrxn, zma_locs=zma_locs)

        if ret is not None:
            proj_hfreqs, unproj_hfreqs, tors_freqs, imag, disps = ret

            # Make final hindered rotor strings and get corrected tors zpe
            if typ.scale_1d(spc_mod_dct_i):

                # scale harmonic from low-level to high-level estimation
                scaled_proj_hfreqs, _ = scale_frequencies(
                   proj_hfreqs, None, spc_mod_dct_i,
                   scale_method='c3_harm')
                scaled_unproj_hfreqs, _ = scale_frequencies(
                    unproj_hfreqs, None, spc_mod_dct_i,
                    scale_method='c3_harm')

                # use harmonic frequencies to determine torsional scaling
                pot_scalef = potential_scale_factor(
                    scaled_unproj_hfreqs, scaled_proj_hfreqs, tors_freqs)
                rotors = tors.scale_rotor_pots(
                    rotors, scale_factor=pot_scalef,
                    scale_override=None)
                tors_strs = tors.make_hr_strings(rotors, mdhr_dct=mdhr_dct)

                # The lines commented below are to use torsional oscillator
                # frequencies from projrot for the ZPVE
                # [_, hr_str, _, prot_str, _] = tors_strs
                # tors_zpe = tors_projected_scaled_zpe(
                #     pf_filesystems, hr_str, prot_str, run_prefix,
                #     spc_mod_dct_i, zrxn=zrxn, zma_locs=zma_locs)

        # For mdhrv model no freqs needed in MESS input, zero out freqs lst
        if 'mdhrv' in spc_mod_dct_i['tors']['mod']:
            proj_ffreqs = ()
            unproj_ffreqs = ()
    else:
        ret = read_harmonic_freqs(
            pf_filesystems, run_prefix, zrxn=zrxn)
        if ret is not None:
            proj_hfreqs, imag, _, disps = ret
            unproj_hfreqs = proj_hfreqs

    # get fundamental/scaled vibrational frequecies and zpe
    if proj_hfreqs:
        if typ.anharm_vib(spc_mod_dct_i):
            ret = fund_frequencies(
                proj_hfreqs, unproj_hfreqs,
                pf_filesystems)
            proj_ffreqs, unproj_ffreqs, total_zpe, mats = ret
            xmat, rovib_mat, rot_dists = mats

        else:
            proj_ffreqs, total_zpe = scale_frequencies(
                proj_hfreqs, tors_zpe,
                spc_mod_dct_i, scale_method='c3')

            unproj_ffreqs, total_zpe = scale_frequencies(
                unproj_hfreqs, 0,
                spc_mod_dct_i, scale_method='c3')

            # Note order of previous functions determines whether
            # torsional zpe is based on unprojected frequencies
            # or on projected frequecies + torsional zpe
            # but lines 79--82 need to be uncommented as well
    vib_anal_dct = {
        'fund_proj_RTimag': unproj_ffreqs,
        'harm_proj_RTimag': unproj_hfreqs,
        'fund_proj_RTimagTors': proj_ffreqs,
        'harm_proj_RTimagTors': proj_hfreqs,
        'harm_tors': tors_freqs,
        'harm_imag': imag,
        'rotors': rotors,
        'anharm_zpe': total_zpe,
        'x_mat': xmat,
        'rovib_mat': rovib_mat,
        'rot_dists': rot_dists,
        'disps': disps,
        'pot_scale_fact': pot_scalef,
        'mess_tors_strs': tors_strs}
    return vib_anal_dct


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
    xmat = None

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
                         prefix, zrxn=None, conf=None, zma_locs=None):
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
    harm_geo, hess, tors_geo = _morph(
        harm_geo, hess, tors_geo,
        zrxn, pf_filesystems, zma_locs)
    fml_str = automol.geom.formula_string(harm_geo)
    vib_path = job_path(prefix, 'PROJROT', 'FREQ', fml_str, print_path=True)
    # print('proj test:', vib_path)
    # tors_path = job_path(run_pfx, 'MESS', 'TORS', fml_str, print_path=True)
    mess_script_str = autorun.SCRIPT_DCT['messpf']
    projrot_script_str = autorun.SCRIPT_DCT['projrot']
    dist_cutoff_dct1 = {('H', 'O'): 2.26767, ('H', 'C'): 2.26767}
    # dist_cutoff_dct2 = {('H', 'O'): 2.83459, ('H', 'C'): 2.83459,
    dist_cutoff_dct2 = {('H', 'O'): 2.83459, ('H', 'C'): 3.023,
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


def tors_projected_scaled_zpe(
        pf_filesystems, mess_hr_str, projrot_hr_str,
        prefix, spc_mod_dct_i, zrxn=None, conf=None,
        zma_locs=None):
    """ Get frequencies from one version of ProjRot
    """
    ret = tors_projected_freqs(
        pf_filesystems, mess_hr_str, projrot_hr_str,
        prefix, zrxn=zrxn, conf=conf, zma_locs=zma_locs)
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
    if spc_mod_dct_i['vib']['scale']=='off':
        scale_method = 'no_scale'
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


def _morph(hess_geo, hess, tors_geo, zrxn, pf_filesystems, zma_locs):
    ret = hess_geo, hess, tors_geo
    if zma_locs not in [(0,), [0]]:
        [cnf_fs, _, min_cnf_locs, _, _] = pf_filesystems['tors']
        zma_fs = autofile.fs.zmatrix(cnf_fs[-1].path(min_cnf_locs))
        print('zma path', zma_fs[0].path(), zma_locs)
        zma = zma_fs[-1].file.zmatrix.read(zma_locs)
        zma_gra = automol.zmat.graph(zma)
        hess_gra = automol.geom.graph(hess_geo)
        tors_gra = automol.geom.graph(tors_geo)
        hess_iso_dct = rxn_util.zmatrix_conversion_keys(hess_gra, zma_gra)
        tors_iso_dct = rxn_util.zmatrix_conversion_keys(tors_gra, zma_gra)
        zma_hess_geo = automol.geom.reorder(hess_geo, hess_iso_dct)
        zma_tors_geo = automol.geom.reorder(tors_geo, tors_iso_dct)
        zma_hess = _reorder_hessian(hess, hess_iso_dct)
        ret = zma_hess_geo, zma_hess, zma_tors_geo
    return ret


def _reorder_hessian(hess, idx_dct):
    """ Reorder the atoms of a hessian using
        the mapping of an input dictionary.

        :param hess: tuple of tuples of the hessian
        :param idx_dct: The new order of the atoms, by index
        :type idx_dct: dict
        :rtype: tuple of tuples for hessiane
    """

    update_hess = deepcopy(hess)
    for orig_i, new_i in idx_dct.items():
        for orig_j, new_j in idx_dct.items():
            update_hess[new_i][new_j] = hess[orig_i][orig_j]
            update_hess[new_i + 1][new_j + 1] = hess[orig_i + 1][orig_j + 1]
            update_hess[new_i + 2][new_j + 2] = hess[orig_i + 2][orig_j + 2]

    return update_hess


def predict_hind_modes(proj_freqs, unproj_freqs):
    """
    determine which modes are projected out by minimizing overlap
    between projected and unprojected frequencies
    INPUTS:
    :param proj_freqs:  frequencies after projection
    :param unproj_freqs:  unprojected frequencies
    :rtype proj_modes: list of idxes of modes
    """
    best_overlap = 100000
    best_combo = None
    num_high_modes = sum(i > 1000 for i in unproj_freqs)
    red_proj_freqs = proj_freqs[:-num_high_modes]
    red_unproj_freqs = unproj_freqs[:-num_high_modes]
    for combo in itertools.permutations(red_unproj_freqs, len(red_proj_freqs)):
        curr_overlap = 0
        for mode, freq in enumerate(combo):
            curr_overlap += abs(freq - red_proj_freqs[mode])
        if curr_overlap < best_overlap:
            best_overlap = curr_overlap
            best_combo = combo
    if best_combo is None:
        proj_modes = []
    else:
        nonproj_modes = [red_unproj_freqs.index(freq) for freq in best_combo]
        proj_modes = list(set(range(len(red_unproj_freqs))) - set(nonproj_modes))
        # next statement checks if procedure failed bc of duplicate freqs
        if len(proj_modes) > len(red_unproj_freqs) - len(red_proj_freqs):
            for idx, val in enumerate(nonproj_modes):
                if idx > 0:
                    if val == nonproj_modes[idx - 1]:
                        if abs(best_combo[idx] - best_combo[idx - 1]) < .1:
                            nonproj_modes[idx] = val + 1
        
            proj_modes = list(set(range(len(red_unproj_freqs))) - set(nonproj_modes))
    return proj_modes


def remove_modes_from_mat(mat, modes, dim=1):
    """
    Removes specified modes from anharmonic constant matrix or vibrot matrix
    INPUTS:
    :param xmat: anharmonic constant matrix
    :param modes: the modes to delete from the matrix
    :xmat  - anharmonic constant matrix with columns and rows deleted
    """
    modes.sort()
    mat = numpy.array(mat)
    for index in modes[::-1]:
        mat = numpy.delete(mat, index, 0)
        if dim > 1:
            mat = numpy.delete(mat, index, 1)
    mat = tuple([tuple(row) for row in mat])
    return mat


def zero_out_modes(mat, modes):
    """
    """
    mat = numpy.array(mat)
    for mode in modes:
        mat[mode, :] = numpy.zeros(len(mat))
        mat[:, mode] = numpy.zeros(len(mat))
    mat = tuple([tuple(row) for row in mat])
    return mat


def compute_lambda(freqs, xmat, temp=298.15):
    """
    Finds perturbation parameters, lambda, for all modes
    INPUT:
    xmat - anharmonic x constant matrix
    freqs - harmonic frequencies
    temp    - temperature (in K)
    OUTPUT:
    lambda_mat  - matrix of lambda_ij values
    """
    def _lambda_ij(w_i, w_j, x_ij, temp):
        return (
            0.69503476*float(temp)*float(x_ij)) / (float(w_i)*float(w_j))

    length = min(len(freqs), len(xmat))
    lambda_mat = numpy.zeros((length, length))
    for i in range(length):
        for j in range(length):
            lambda_mat[i][j] = _lambda_ij(
                freqs[i], freqs[j], xmat[i][j], temp)
    return lambda_mat


def evaluate_lambda(lambda_mat, freqs, thresh=0.01):
    """
    Determines which modes cause bad lambda values
    """
    bad_modes = []
    for i, _ in enumerate(lambda_mat):
        if abs(lambda_mat[i][i]) > thresh:
            print(
               'lambda value of {:.4f} is greater than thresh of {:.4f}. Mode v_{:d}={:.2f} is flagged'.format(
                   lambda_mat[i][i], thresh, i, freqs[i]))
            bad_modes.append(i)

    for _ in range(len(lambda_mat)):
        largest = 0
        index = None
        for i, lambda_row in enumerate(lambda_mat):
            if i in bad_modes:
                continue
            for j, lam in enumerate(lambda_row):
                if j in bad_modes:
                    continue
                if abs(lam) < thresh or abs(lam) < largest:
                    continue
                largest = abs(lam)
                lambda_mat[i][j] = 0
                index = j
                if abs(lambda_mat[i][i]) > abs(lambda_mat[j][j]):
                    index = i
        if index is not None:
            print(
               'lambda value of {:.4f} is greater than thresh of {:.4f}. Mode v_{:d}={:.2f} is flagged'.format(
                   largest, thresh, index, freqs[index]))
            bad_modes.append(index)
    for i in range(len(freqs)):
        for j in range(i):
            if abs(freqs[i] - freqs[j]) < .05:
                if i in bad_modes and not j in bad_modes:
                    bad_modes.append(j)
                    print('Mode v_{:d}={:.2f} is flagged because it is degenerate to flagged mode'.format(j, freqs[j]))
                if j in bad_modes and not i in bad_modes:
                    bad_modes.append(i)
                    print('Mode v_{:d}={:.2f} is flagged because it is degenerate to flagged mode'.format(i, freqs[i]))
    return list(set(bad_modes))


def anharm_freqs(freqs, xmat):
    """
    Uses anharmonic frequency matrix and harmonic frequencies to compute VPT2 anharmonic frequencies
    INPUT:
    freqs   - harmonic frequencies
    xmat    - anharmonic constant matrix
    OUTPUT:
    anharms - VPT2 anharmonic frequencies
    """
    anharms = numpy.zeros(len(freqs))
    for i, freq in enumerate(freqs):
        anharms[i] = freq
        anharms[i] += 2. * xmat[i][i]
        tmp = 0
        for j in range(len(freqs)):
            if j != i:
                tmp += xmat[i][j]
        anharms[i] += 1./2 * tmp
    return tuple(anharms)


def anharm_zpe(freqs, xmat):
    xmat = numpy.array(xmat)
    freqs = numpy.array(freqs)
    zpe = 1/2. * numpy.sum(freqs)
    zpe += 1/4. * numpy.sum(numpy.triu(xmat, k=0))
    return zpe * phycon.WAVEN2EH


def fund_frequencies(
        proj_hfreqs, unproj_hfreqs, pf_filesystems):
    """
    """
    # read in mats
    xmat = read_anharmon_matrix(pf_filesystems)
    rovib_mat, rot_dists = rot.read_rotational_values(pf_filesystems)

    # zero out values that blow up vpt2, often umbrella modes
    lambda_mat = compute_lambda(unproj_hfreqs, xmat)
    high_lambda_modes = evaluate_lambda(lambda_mat, unproj_hfreqs)
    xmat = zero_out_modes(xmat, high_lambda_modes)

    # get unprojected fundamental freqs and zpve
    unproj_ffreqs = anharm_freqs(unproj_hfreqs, xmat)
    total_zpe = anharm_zpe(unproj_hfreqs, xmat)
    print('total (unprojected) anharmonic zpve', total_zpe)
    [cnf_fs, cnf_path, min_cnf_locs, _, _] = pf_filesystems['vpt2']
    if cnf_path:
        print(
            'anharm zpve with G0 from elec. output',
            cnf_fs[-1].file.anharmonic_zpve.read(
                min_cnf_locs))
    # delete out projected hindered rotor modes from
    # xmat and rovib mat to getprojected fund. freqs
    tors_proj_modes = predict_hind_modes(proj_hfreqs, unproj_hfreqs)
    xmat = remove_modes_from_mat(xmat, tors_proj_modes, dim=2)
    rovib_mat = remove_modes_from_mat(rovib_mat, tors_proj_modes)
    proj_ffreqs = anharm_freqs(proj_hfreqs, xmat)

    return (
        proj_ffreqs, unproj_ffreqs, total_zpe,
        (xmat, rovib_mat, rot_dists,))
