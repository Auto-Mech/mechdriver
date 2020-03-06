"""
  Build the species string based on the model
"""

import elstruct
import automol
from routines.es.scan import hr_prep
from routines.pf.messf import _tors as tors
from routines.pf.messf import _vib as vib
from routines.pf.messf import _tau as tau
# from routines.pf.messf import _vpt2 as vpt2
from lib.phydat import phycon


def vib_harm_tors_rigid(spc_info, harm_min_cnf_locs, harm_cnf_save_fs,
                        saddle=False):
    """ Build the species string for a model of
        harmonic vibrational frequencies and rigid torsions
    """
    # Do the freqs obtain for two species for fake and pst
    if harm_min_cnf_locs is not None:
        # Obtain geom and freqs from filesys
        harm_geo = harm_cnf_save_fs[-1].file.geometry.read(harm_min_cnf_locs)
        hess = harm_cnf_save_fs[-1].file.hessian.read(harm_min_cnf_locs)
        freqs = elstruct.util.harmonic_frequencies(
            harm_geo, hess, project=False)
        # Modify freqs lst and get imaginary frequencies
        mode_start = 6
        if saddle:
            mode_start = mode_start + 1
            imag_freq = freqs[0]
        else:
            imag_freq = ''
        if automol.geom.is_linear(harm_geo):
            mode_start = mode_start - 1
        freqs = freqs[mode_start:]
    else:
        print('ERROR: Reference geometry is missing for harmonic frequencies ',
              'for species {}'.format(spc_info[0]))
        raise ValueError

    return harm_geo, freqs, imag_freq


def vib_harm_tors_1dhr(harm_min_cnf_locs, harm_cnf_save_fs,
                       tors_min_cnf_locs, tors_cnf_save_fs,
                       tors_save_path, tors_cnf_save_path,
                       spc_dct_i, spc_info,
                       frm_bnd_key, brk_bnd_key,
                       sym_factor, elec_levels,
                       saddle=False):
    """ Build the species string for a model: Harm, 1DHR
    """
    if harm_min_cnf_locs is not None:
        harm_geo = harm_cnf_save_fs[-1].file.geometry.read(
            harm_min_cnf_locs)
        min_ene = harm_cnf_save_fs[-1].file.energy.read(
            harm_min_cnf_locs)
        hess = harm_cnf_save_fs[-1].file.hessian.read(
            harm_min_cnf_locs)
        freqs = elstruct.util.harmonic_frequencies(
            harm_geo, hess, project=False)

        if tors_min_cnf_locs is not None:

            # Get geometry for the torsional minimum
            zma = tors_cnf_save_fs[-1].file.zmatrix.read(
                tors_min_cnf_locs)
            tors_geo = tors_cnf_save_fs[-1].file.geometry.read(
                tors_min_cnf_locs)

            # Set torsional stuff
            tors_names = tors.get_tors_names(
                spc_dct_i, tors_cnf_save_fs, saddle=saddle)
            tors_grids = tors.get_tors_grids(
                spc_dct_i, zma, tors_names, frm_bnd_key, brk_bnd_key)
            tors_sym_nums = tors.get_tors_sym_nums(
                spc_dct_i, tors_min_cnf_locs, tors_cnf_save_fs,
                frm_bnd_key, brk_bnd_key, saddle=False)

            # Set ts bond
            ts_bnd = None
            if saddle:
                dist_name = spc_dct_i['dist_info'][0]
                ts_bnd = automol.zmatrix.bond_idxs(zma, dist_name)

            # Write strings containing rotor info for MESS and ProjRot
            hind_rot_str, proj_rotors_str = tors.write_1dhr_tors_mess_strings(
                harm_geo, spc_info, spc_dct_i, ts_bnd, zma,
                tors_names, tors_grids, tors_sym_nums,
                tors_cnf_save_path, min_ene,
                saddle=False, hind_rot_geo=None)

            # Calculate ZPVES of the hindered rotors
            if saddle and tors_names is not None:
                tors_zpe = tors.calc_tors_freqs_zpe(
                    tors_geo, sym_factor, elec_levels,
                    hind_rot_str, tors_save_path)
            else:
                tors_zpe = 0.0

            # Run one vers ProjRot to proj freqs for that version
            freqs1, imag_freq1, zpe_harm_no_tors = vib.projrot_freqs_1(
                tors_geo, hess,
                proj_rotors_str,
                tors_save_path, pot=True, saddle=False)

            # Now run the other version of ProjRot
            pfreqs2 = vib.projrot_freqs_2(
                tors_save_path, pot=True, saddle=saddle)
            [freqs2, imag_freq2,
             zpe_harm_no_tors_2, harm_zpe] = pfreqs2

            # Determine freqs and imag_freqs
            freqs, imag_freq, zpe = vib.determine_freqs_zpe(
                freqs1, freqs2, imag_freq1, imag_freq2,
                zpe_harm_no_tors, zpe_harm_no_tors_2,
                harm_zpe, tors_zpe)
    else:
        print('ERROR: Reference geometry is missing for harmonic frequencies',
              ' for species {}'.format(spc_info[0]))
        tors_geo, freqs, imag_freq, hind_rot_str = (), (), 0.0, ''
        raise ValueError

    return tors_geo, tuple(freqs), imag_freq, hind_rot_str, zpe


def vib_harm_tors_mdhr(harm_min_cnf_locs, harm_cnf_save_fs,
                       tors_min_cnf_locs, tors_cnf_save_fs,
                       tors_save_path, tors_cnf_save_path,
                       spc_dct_i, spc_info,
                       frm_bnd_key, brk_bnd_key,
                       sym_factor, elec_levels,
                       tors_mod,
                       saddle=False):
    """ Build the species string for a model: Harm, 1DHR
    """

    # Set up some info with the torsions
    ndim_tors, _ = tors_mod
    if 'hind_def' in spc_dct_i:
        run_tors_names = spc_dct_i['hind_def']
    else:
        run_tors_names = ()
    if 'hind_inc' in spc_dct_i:
        scan_increment = spc_dct_i['hind_inc'] * phycon.DEG2RAD
    else:
        scan_increment = 30. * phycon.DEG2RAD

    if harm_min_cnf_locs is not None:
        harm_geo = harm_cnf_save_fs[-1].file.geometry.read(
            harm_min_cnf_locs)
        min_ene = harm_cnf_save_fs[-1].file.energy.read(
            harm_min_cnf_locs)
        hess = harm_cnf_save_fs[-1].file.hessian.read(
            harm_min_cnf_locs)
        freqs = elstruct.util.harmonic_frequencies(
            harm_geo, hess, project=False)

        if tors_min_cnf_locs is not None:

            # Get geometry for the torsional minimum
            zma = tors_cnf_save_fs[-1].file.zmatrix.read(
                tors_min_cnf_locs)
            tors_geo = tors_cnf_save_fs[-1].file.geometry.read(
                tors_min_cnf_locs)

            # Get the hr prep stuff
            tors_name_grps, tors_grid_grps = hr_prep(
                zma, tors_geo, run_tors_names=run_tors_names,
                scan_increment=scan_increment, ndim_tors=ndim_tors,
                saddle=saddle,
                frm_bnd_key=(), brk_bnd_key=())
            # frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key)

            # Set torsional stuff
            tors_sym_nums = tors.get_tors_sym_nums(
                spc_dct_i, tors_min_cnf_locs, tors_cnf_save_fs,
                frm_bnd_key, brk_bnd_key, saddle=False)

            # Set ts bond
            ts_bnd = None
            if saddle:
                dist_name = spc_dct_i['dist_info'][0]
                ts_bnd = automol.zmatrix.bond_idxs(zma, dist_name)

            # Write strings containing rotor info for MESS and ProjRot
            hind_rot_str, proj_rotors_str, mdhr_dat_str = tors.write_mdhr_tors_mess_strings(
                harm_geo, spc_info, 1.0, spc_dct_i, ts_bnd, zma,
                tors_name_grps, tors_grid_grps, tors_sym_nums,
                tors_cnf_save_path, min_ene,
                saddle=False, hind_rot_geo=None)

            # Calculate ZPVES of the hindered rotors
            if saddle:  # and tors_names is not None:
                tors_zpe = tors.calc_tors_freqs_zpe(
                    tors_geo, sym_factor, elec_levels,
                    hind_rot_str, tors_save_path)
            else:
                tors_zpe = 0.0

            # Run one vers ProjRot to proj freqs for that version
            freqs1, imag_freq1, zpe_harm_no_tors = vib.projrot_freqs_1(
                tors_geo, hess,
                proj_rotors_str,
                tors_save_path, pot=True, saddle=False)

            # Now run the other version of ProjRot
            pfreqs2 = vib.projrot_freqs_2(
                tors_save_path, pot=True, saddle=saddle)
            [freqs2, imag_freq2,
             zpe_harm_no_tors_2, harm_zpe] = pfreqs2

            # Determine freqs and imag_freqs
            freqs, imag_freq, zpe = vib.determine_freqs_zpe(
                freqs1, freqs2, imag_freq1, imag_freq2,
                zpe_harm_no_tors, zpe_harm_no_tors_2,
                harm_zpe, tors_zpe)
    else:
        print('ERROR: Reference geometry is missing for harmonic frequencies',
              ' for species {}'.format(spc_info[0]))
        tors_geo, freqs, imag_freq, hind_rot_str = (), (), 0.0, ''
        raise ValueError

    return tors_geo, freqs, imag_freq, hind_rot_str, zpe, mdhr_dat_str


def vib_harm_tors_tau(harm_min_cnf_locs, harm_cnf_save_fs,
                      tors_min_cnf_locs, tors_cnf_save_fs,
                      tors_save_path, tors_cnf_save_path,
                      spc_dct_i, spc_info,
                      frm_bnd_key, brk_bnd_key,
                      sym_factor, elec_levels,
                      save_prefix,
                      tau_dat_file_name,
                      saddle=False):
    """ Build the species string for a model: Harm, 1DHR
    """
    if harm_min_cnf_locs is not None:
        harm_geo = harm_cnf_save_fs[-1].file.geometry.read(
            harm_min_cnf_locs)
        min_ene = harm_cnf_save_fs[-1].file.energy.read(
            harm_min_cnf_locs)
        hess = harm_cnf_save_fs[-1].file.hessian.read(
            harm_min_cnf_locs)
        freqs = elstruct.util.harmonic_frequencies(
            harm_geo, hess, project=False)

        if tors_min_cnf_locs is not None:

            # Get geometry for the torsional minimum
            zma = tors_cnf_save_fs[-1].file.zmatrix.read(
                tors_min_cnf_locs)
            tors_geo = tors_cnf_save_fs[-1].file.geometry.read(
                tors_min_cnf_locs)

            # Set torsional stuff
            tors_names = tors.get_tors_names(
                spc_dct_i, tors_cnf_save_fs, saddle=saddle)
            tors_grids = tors.get_tors_grids(
                spc_dct_i, zma, tors_names, frm_bnd_key, brk_bnd_key)
            tors_sym_nums = tors.get_tors_sym_nums(
                spc_dct_i, tors_min_cnf_locs, tors_cnf_save_fs,
                frm_bnd_key, brk_bnd_key, saddle=False)

            # Set ts bond
            ts_bnd = None
            if saddle:
                dist_name = spc_dct_i['dist_info'][0]
                ts_bnd = automol.zmatrix.bond_idxs(zma, dist_name)

            # Write strings containing rotor info for ProjRot should split
            mess_str, proj_rotors_str = tors.write_1dhr_tors_mess_strings(
                harm_geo, spc_info, spc_dct_i, ts_bnd, zma,
                tors_names, tors_grids, tors_sym_nums,
                tors_cnf_save_path, min_ene,
                saddle=False, hind_rot_geo=None)

            # Calculate ZPVES of the hindered rotors
            if saddle and tors_names is not None:
                tors_zpe = tors.calc_tors_freqs_zpe(
                    tors_geo, sym_factor, elec_levels,
                    mess_str, tors_save_path)
            else:
                tors_zpe = 0.0

            # Run one vers ProjRot to proj freqs for that version
            freqs1, imag_freq1, zpe_harm_no_tors = vib.projrot_freqs_1(
                tors_geo, hess,
                proj_rotors_str,
                tors_save_path, pot=True, saddle=False)

            # Now run the other version of ProjRot
            pfreqs2 = vib.projrot_freqs_2(
                tors_save_path, pot=True, saddle=saddle)
            [freqs2, imag_freq2,
             zpe_harm_no_tors_2, harm_zpe] = pfreqs2

            # Determine freqs and imag_freqs
            freqs, imag_freq, zpe = vib.determine_freqs_zpe(
                freqs1, freqs2, imag_freq1, imag_freq2,
                zpe_harm_no_tors, zpe_harm_no_tors_2,
                harm_zpe, tors_zpe)

            # Write the tau Monte Carlo Core section
            monte_carlo_str = tau.write_monte_carlo_mess_strings(
                tors_min_cnf_locs, tors_cnf_save_fs,
                spc_dct_i,
                frm_bnd_key, brk_bnd_key,
                elec_levels,
                tau_dat_file_name,
                freqs=freqs,
                saddle=saddle)

            # Write the file containing the tau enes, geos, grads, hess
            spc_formula = automol.inchi.formula(spc_dct_i['ich'])
            tau_dat_str = tau.write_tau_data_str(
                spc_formula, save_prefix, gradient=False, hessian=False)

    else:
        print('ERROR: Reference geometry is missing for harmonic frequencies',
              ' for species {}'.format(spc_info[0]))
        tors_geo, freqs, imag_freq, hind_rot_str = (), (), 0.0, ''
        raise ValueError

    return tors_geo, freqs, imag_freq, hind_rot_str, zpe, sym_factor, monte_carlo_str, tau_dat_str


def vib_tau_tors_tau(tors_min_cnf_locs, tors_cnf_save_fs,
                     spc_dct_i,
                     frm_bnd_key, brk_bnd_key,
                     sym_factor, elec_levels,
                     save_prefix,
                     tau_dat_file_name,
                     saddle=False):
    """ Build the species string for a model: Harm, 1DHR
    """

    # Write the Monte Carlo Core section
    monte_carlo_str = tau.write_monte_carlo_mess_strings(
        tors_min_cnf_locs, tors_cnf_save_fs,
        spc_dct_i,
        frm_bnd_key, brk_bnd_key,
        elec_levels,
        tau_dat_file_name,
        freqs=(),
        saddle=saddle)

    # Write the file containing the tau enes, geos, grads, hess
    spc_formula = automol.inchi.formula(spc_dct_i['ich'])
    tau_dat_str = tau.write_tau_data_str(
        spc_formula, save_prefix, gradient=True, hessian=True)

    return monte_carlo_str, tau_dat_str


def vib_vpt2_tors_rigid(harm_min_cnf_locs, harm_cnf_save_fs,
                        vpt2_min_cnf_locs, vpt2_cnf_save_fs,
                        saddle=False):
    """ Build the species info for a model: VPT2, Rigid
    """
    if harm_min_cnf_locs is not None:
        harm_geo = harm_cnf_save_fs[-1].file.geometry.read(
            harm_min_cnf_locs)
        hess = harm_cnf_save_fs[-1].file.hessian.read(
            harm_min_cnf_locs)
        freqs = elstruct.util.harmonic_frequencies(
            harm_geo, hess, project=False)
        # Modify freqs lst and get imaginary frequencies
        mode_start = 6
        if saddle:
            mode_start = mode_start + 1
            imag_freq = freqs[0]
        else:
            imag_freq = ''
        if automol.geom.is_linear(harm_geo):
            mode_start = mode_start - 1
        freqs = freqs[mode_start:]
        harm_zpve = sum(freqs)/2.0

    if vpt2_min_cnf_locs is not None:
        xmat = vpt2_cnf_save_fs[-1].file.anharmonicity_matrix.read(
            vpt2_min_cnf_locs)
        vpt2_anharm_zpve = vpt2_cnf_save_fs[-1].file.anharmonic_zpve.read(
            vpt2_min_cnf_locs)
        hess = harm_cnf_save_fs[-1].file.hessian.read(
            harm_min_cnf_locs)
        freqs = elstruct.util.harmonic_frequencies(
            harm_geo, hess, project=False)
        # Modify freqs lst and get imaginary frequencies
        mode_start = 6
        if saddle:
            mode_start = mode_start + 1
            imag_freq = freqs[0]
        else:
            imag_freq = ''
        if automol.geom.is_linear(harm_geo):
            mode_start = mode_start - 1
        freqs = freqs[mode_start:]
        vpt2_harm_zpve = sum(freqs)/2.0

    zpve = harm_zpve + (vpt2_anharm_zpve - vpt2_harm_zpve)

    return harm_geo, freqs, imag_freq, xmat, zpve


def vib_vpt2_tors_1dhr(harm_min_cnf_locs, harm_cnf_save_fs,
                       vpt2_min_cnf_locs, vpt2_cnf_save_fs,
                       saddle=False):
    """ Build the species info for a model: VPT2, Rigid
    """
    # Get the unprojected frequencies
    if harm_min_confs:
        harm_geo = harm_cnf_save_fs[-1].file.geometry.read(
            harm_min_cnf_locs)
        hess = harm_cnf_save_fs[-1].file.hessian.read(
            harm_min_cnf_locs)
        unproj_freqs = elstruct.util.harmonic_frequencies(
            harm_geo, hess, project=False)

    # Get the projected frequencies
    proj_freqs = 10
    
    # Determine what modes were projected out
    modes = vpt2.find_hinfreqs(proj, unproj, order)
    xmat_mod = vpt2.remove_modes(xmat, modes)
   
    # Determine the anharmonic freqs
    anharm_freqs = anharm_freq(freqs, xmat)

    return harm_geo, anharm_freqs, imag_freq, xmat_mod, zpve
