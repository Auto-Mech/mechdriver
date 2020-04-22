"""
  Functions handling hindered rotor model calculations
"""

import os
import numpy
from scipy.interpolate import interp1d
import automol
import mess_io
import projrot_io
import autofile

# New libs
from lib.phydat import phycon
from lib.runner import script
from lib.struct import torsprep


# Function to deal with setting up all of the torsion info since it is a pain
def make_tors_info(spc_dct_i, cnf_save_fs, cnf_save_locs, tors_model,
                   saddle=False, frm_bnd_key=(), brk_bnd_key=()):
    """ get tors stuff
    """

    # Set up some info with the torsions
    if 'hind_def' in spc_dct_i:
        run_tors_names = spc_dct_i['hind_def']
    else:
        run_tors_names = ()
    if 'hind_inc' in spc_dct_i:
        scan_increment = spc_dct_i['hind_inc'] * phycon.DEG2RAD
    else:
        scan_increment = 30. * phycon.DEG2RAD
    if 'tors_names' in spc_dct_i:
        dct_tors_names = spc_dct_i['tors_names']
        run_tors_names = [[name] for name in dct_tors_names]
    else:
        run_tors_names = ()

    tors_name_grps, tors_grid_grps, tors_sym_grps = (), (), ()
    if cnf_save_locs is not None:

        # Get geometry for the torsional minimum
        zma = cnf_save_fs[-1].file.zmatrix.read(cnf_save_locs)
        geom = cnf_save_fs[-1].file.geometry.read(cnf_save_locs)

        # Get the hr prep stuff
        tors_name_grps, tors_grid_grps, tors_sym_grps = torsprep.hr_prep(
            zma, geom, run_tors_names=run_tors_names,
            scan_increment=scan_increment, tors_model=tors_model,
            saddle=saddle,
            frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key)

    return tors_name_grps, tors_grid_grps, tors_sym_grps


# Functions to write MESS strings
def make_hr_strings(spc_dct_i, tors_name_grps, tors_grid_grps, tors_sym_grps,
                    tors_model, cnf_save_path, saddle=False):
    """ Procedure for building the MESS strings
    """

    # Get the torsion info

    # Set ts bond
    ts_bnd = None
    if saddle:
        dist_name = spc_dct_i['dist_info'][0]

    # Get minimum-energy
    min_ene = 0.0

    # Write strings containing rotor info for MESS and ProjRot
    if tors_model in ('1dhr', '1dhrf'):
        mess_str, projrot_str = _make_1dhr_tors_strs(
            geom, spc_info, spc_dct_i, ts_bnd, zma,
            tors_names, tors_grids, tors_sym_nums,
            cnf_save_path, min_ene,
            saddle=False, hind_rot_geo=None)
        mdhr_dat_str_lst= []

    elif tors_model in ('mdhr', 'mdhrv'):
        mess_str, projrot_str, mdhr_dat_str_lst = _make_mdhr_tors_strs(
            geom, spc_info, sym_num, spc_dct_i,
            ts_bnd, zma,
            tors_name_grps, tors_grid_grps, tors_sym_nums,
            cnf_save_path, min_ene,
            saddle=False, hind_rot_geo=None,
            vib_adiabatic=False)

    return hind_rot_str, proj_rotors_str, mdhr_dat_str_lst


def _make_1dhr_tors_strs(harm_geo, spc_info, rxn_class, ts_bnd, zma,
                         tors_names, tors_grids, tors_sym_nums,
                         cnf_save_path, min_ene,
                         saddle=False, hind_rot_geo=None):
    """ Gather the 1DHR torsional data and gather them into a MESS file
    """

    # Initialize empty strings
    mess_hr_str, projrot_hr_str = '', ''

    # Loop over the torsions
    tors_info = zip(tors_names, tors_grids, tors_sym_nums)
    for tors_name, tors_grid, tors_sym in tors_info:

        # Read the hindered rotor potential
        pot, _ = _read_hr_pot(
            spc_info, [tors_name], tors_grid,
            cnf_save_path, min_ene)

        # Build potential lst from only successful calculations
        pot = _hrpot_spline_fitter(pot)

        # Get the HR groups and axis for the rotor
        group, axis, atm_key = torsprep.set_groups_ini(
            zma, tors_name, ts_bnd, saddle)
        if saddle:
            group, axis, pot = torsprep.check_saddle_groups(
                zma, rxn_class, group, axis,
                pot, ts_bnd, tors_sym)
        group = list(numpy.add(group, 1))
        axis = list(numpy.add(axis, 1))
        if (atm_key+1) != axis[1]:
            axis.reverse()

        # Check for dummy transformations
        remdummy = torsprep.check_dummy_trans(zma)

        # Write the MESS and ProjRot strings for the rotor
        hrgeo = harm_geo if hind_rot_geo else None
        mess_hr_str += mess_io.writer.rotor_hindered(
            group, axis, tors_sym, pot, remdummy=remdummy, geom=hrgeo)
        projrot_hr_str += projrot_io.writer.rotors(
            axis, group, remdummy=remdummy)

    return mess_hr_str, projrot_hr_str


def _make_mdhr_tors_strs(geom, spc_info, sym_num, spc_dct_i,
                         ts_bnd, zma,
                         tors_name_grps, tors_grid_grps, tors_sym_nums,
                         tors_cnf_save_path, min_ene,
                         saddle=False, hind_rot_geo=None,
                         vib_adiabatic=False):
    """ Gather the MDHR torsional data and gather them into a MESS file
    """

    # Loop over the torsion groups and get the int rot strings and potentials
    mess_hr_str, projrot_hr_str = '', ''
    mdhr_dat_str_lst = []

    # Loop over the torsions
    tors_info = zip(tors_name_grps, tors_grid_grps)
    for idx, (tors_names, tors_grids) in enumerate(tors_info):

        # Read the hindered rotor potential and add to master list
        hr_pot, hr_freqs = _read_hr_pot(
            spc_info, tors_names, tors_grids,
            tors_cnf_save_path, min_ene,
            read_freqs=vib_adiabatic)

        # Write the MDHR potential file for each rotor set; append to lst
        mdhr_dat_str_lst.append(mess_io.mdhr_data(hr_pot, hr_freqs))

        # Check for dummy transformations
        remdummy = torsprep.check_dummy_trans(zma)

        # Loop over the rotors in the group and write the internal rotor strs
        for tors_name, _ in zip(tors_names, tors_grids):

            # Set pot to empty list (may need fix)
            pot = ()

            # Get the HR groups and axis for the rotor
            group, axis, atm_key = torsprep.set_groups_ini(
                zma, tors_name, ts_bnd, saddle)
            if saddle:
                group, axis, pot = torsprep.check_saddle_groups(
                    zma, spc_dct_i, group, axis,
                    pot, ts_bnd, tors_sym_nums[idx])
            group = list(numpy.add(group, 1))
            axis = list(numpy.add(axis, 1))
            if (atm_key+1) != axis[1]:
                axis.reverse()

            # Write the MESS and ProjRot strings for the rotor
            mess_hr_str += mess_io.writer.mol_data.rotor_internal(
                group, axis, tors_sym_nums[idx],
                rotor_id='', remdummy=remdummy,
                mass_exp_size=5, pot_exp_size=5,
                hmin=13, hmax=101,
                grid_size=100)
            projrot_hr_str += projrot_io.writer.rotors(
                axis, group, remdummy=remdummy)

    return mess_hr_str, projrot_hr_str, mdhr_dat_str_lst


def make_flux_str(tors_min_cnf_locs, tors_cnf_save_fs,
                  tors_names, tors_sym_nums):
    """ Write out the input string for tau samling
    """
    # Loop over the torsions to get the flux strings
    flux_mode_str = ''
    if tors_min_cnf_locs is not None:

        # Get geometry for the torsional minimum
        zma = tors_cnf_save_fs[-1].file.zmatrix.read(
            tors_min_cnf_locs)
        name_matrix = automol.zmatrix.name_matrix(zma)
        key_matrix = automol.zmatrix.key_matrix(zma)

        # Write the MESS flux strings for each of the modes
        for tors_name, tors_sym in zip(tors_names, tors_sym_nums):

            # Get the idxs for the atoms used to define the dihedrals
            # Move code at some point to automol
            mode_idxs = [0, 0, 0, 0]
            for i, name_row in enumerate(name_matrix):
                if tors_name in name_row:
                    mode_idxs[0] = i
                    mode_idxs[1], mode_idxs[2], mode_idxs[3] = key_matrix[i]
                    break
            mode_idxs = [idx+1 for idx in mode_idxs]

            # Determine the flux span using the symmetry number
            mode_span = 360.0 / tors_sym

            # Write flux string for mode_ to overall flux string
            flux_mode_str += mess_io.writer.fluxional_mode(
                mode_idxs, span=mode_span)

    return flux_mode_str


# Use the torsions and vibrational frequenices to calculate the zpve
def calc_tors_zpves(geom, sym_factor, elec_levels, hess,
                    mess_hr_str, projrot_hr_str,
                    cnf_save_path, pot=True, saddle=False):
    """ proc to calculate freqs and zpves
    """
    # Calculate ZPVES of the hindered rotors
    if saddle:  # and tors_names is not None:
        tors_zpe = calc_tors_freqs_zpe(
            geom, sym_factor, elec_levels,
            mess_hr_str, cnf_save_path)
    else:
        tors_zpe = 0.0

    # Run one vers ProjRot to proj freqs for that version
    freqs1, imag_freq1, zpe_harm_no_tors = vib.projrot_freqs_1(
        geom, hess, projrot_hr_str,
        cnf_save_path, pot=True, saddle=False)

    # Now run the other version of ProjRot
    pfreqs2 = vib.projrot_freqs_2(
        cnf_save_path, pot=True, saddle=saddle)
    [freqs2, imag_freq2,
     zpe_harm_no_tors_2, harm_zpe] = pfreqs2

    # Determine freqs and imag_freqs
    freqs, imag_freq, zpe = vib.determine_freqs_zpe(
        freqs1, freqs2, imag_freq1, imag_freq2,
        zpe_harm_no_tors, zpe_harm_no_tors_2,
        harm_zpe, tors_zpe)

    return freqs, imag_freq, zpe


def calc_tors_freqs_zpe(tors_geo, sym_factor, elec_levels,
                        hind_rot_str, tors_save_path):
    """ Calculate the frequencies and ZPVES of the hindered rotors
        create a messpf input and run messpf to get tors_freqs and tors_zpes
    """
    dummy_freqs = [1000.]
    dummy_zpe = 0.0
    core = mess_io.writer.core_rigidrotor(tors_geo, sym_factor)
    spc_str = mess_io.writer.molecule(
        core, dummy_freqs, elec_levels,
        hind_rot=hind_rot_str,
        )
    temp_step = 100.
    ntemps = 5
    zpe_str = '{0:<8.2f}\n'.format(dummy_zpe)
    zpe_str = ' ZeroEnergy[kcal/mol] ' + zpe_str
    zpe_str += 'End\n'
    global_pf_str = mess_io.writer.global_pf(
        [], temp_step, ntemps, rel_temp_inc=0.001,
        atom_dist_min=0.6)
    spc_head_str = 'Species ' + ' Tmp'
    pf_inp_str = '\n'.join(
        [global_pf_str, spc_head_str,
         spc_str, zpe_str])

    bld_locs = ['PF', 0]
    bld_save_fs = autofile.fs.build(tors_save_path)
    bld_save_fs[-1].create(bld_locs)
    pf_path = bld_save_fs[-1].path(bld_locs)

    # run messpf
    with open(os.path.join(pf_path, 'pf.inp'), 'w') as pf_file:
        pf_file.write(pf_inp_str)
    pf_script_str = ("#!/usr/bin/env bash\n"
                     "export OMP_NUM_THREADS=10\n"
                     "messpf pf.inp pf.out >> stdout.log &> stderr.log")

    script.run_script(pf_script_str, pf_path)

    with open(os.path.join(pf_path, 'pf.log'), 'r') as mess_file:
        output_string = mess_file.read()

    # Read the freqs and zpes
    # tors_freqs = mess_io.reader.tors.freqs(output_string)
    tors_zpes = mess_io.reader.tors.zpves(output_string)

    # Calculate the torsional zpe
    tors_zpe = sum(tors_zpes) if tors_zpes else 0.0
    # tors_zpe_cor = 0.0
    # tors_zpe = 0.0
    # for (tors_freq, tors_1dhr_zpe) in zip(tors_freqs, tors_zpes):
    #     tors_zpe_cor += tors_1dhr_zpe - tors_freq*phycon.WAVEN2KCAL/2
    #     tors_zpe += tors_1dhr_zpe

    return tors_zpe


# Functions to obtain values of the HR potentials from the filesystem
def _read_hr_pot(tors_names, tors_grids, cnf_save_path, min_ene,
                 read_freqs=False):
    """ Get the potential for a hindered rotor
    """

    # Build template pot lst and freqs list into a list-of-lists if ndim > 1
    if len(tors_names) == 1:
        dims = (len(tors_grids),)
    elif len(tors_names) == 2:
        dims = (len(tors_grids[0]), len(tors_grids[1]))
    elif len(tors_names) == 3:
        dims = (len(tors_grids[0]), len(tors_grids[1]), len(tors_grids[2]))
    elif len(tors_names) == 4:
        dims = (len(tors_grids[0]), len(tors_grids[1]),
                len(tors_grids[2]), len(tors_grids[3]))
    pot = numpy.zeros(dims).tolist()
    if read_freqs:
        freqs = numpy.zeros(dims).tolist()
    else:
        freqs = []

    # Read the energies from the filesystem
    scn_fs = autofile.fs.scan(cnf_save_path)
    if len(tors_names) == 1:
        for i, grid_val_i in enumerate(tors_grids):
            locs = [tors_names, [grid_val_i]]
            if scn_fs[-1].exists(locs):
                ene = scn_fs[-1].file.energy.read(locs)
                pot[i] = (ene - min_ene) * phycon.EH2KCAL
            else:
                pot[i] = 10.0
            if read_freqs:
                freqs[i] = scn_fs[-1].file.harmonic_frequencies.read(locs)
    elif len(tors_names) == 2:
        for i, grid_val_i in enumerate(tors_grids[0]):
            for j, grid_val_j in enumerate(tors_grids[1]):
                locs = [tors_names, [grid_val_i, grid_val_j]]
                if scn_fs[-1].exists(locs):
                    ene = scn_fs[-1].file.energy.read(locs)
                    pot[i][j] = (ene - min_ene) * phycon.EH2KCAL
                else:
                    pot[i][j] = 10.0
                if read_freqs:
                    freqs = scn_fs[-1].file.harmonic_frequencies.read(locs)
                    freqs[i][j] = freqs
    elif len(tors_names) == 3:
        for i, grid_val_i in enumerate(tors_grids[0]):
            for j, grid_val_j in enumerate(tors_grids[1]):
                for k, grid_val_k in enumerate(tors_grids[2]):
                    locs = [tors_names, [grid_val_i, grid_val_j, grid_val_k]]
                    if scn_fs[-1].exists(locs):
                        ene = scn_fs[-1].file.energy.read(locs)
                        pot[i][j][k] = (ene - min_ene) * phycon.EH2KCAL
                    else:
                        pot[i][j][k] = 10.0
                    if read_freqs:
                        freqs = scn_fs[-1].file.harmonic_frequencies.read(locs)
                        freqs[i][j][k] = freqs
    elif len(tors_names) == 4:
        for i, grid_val_i in enumerate(tors_grids[0]):
            for j, grid_val_j in enumerate(tors_grids[1]):
                for k, grid_val_k in enumerate(tors_grids[2]):
                    for lm, grid_val_l in enumerate(tors_grids[3]):
                        locs = [tors_names,
                                [grid_val_i, grid_val_j,
                                 grid_val_k, grid_val_l]]
                        if scn_fs[-1].exists(locs):
                            ene = scn_fs[-1].file.energy.read(locs)
                            pot[i][j][k][lm] = (ene - min_ene) * phycon.EH2KCAL
                        else:
                            pot[i][j][k][lm] = 10.0
                        if read_freqs:
                            freqs = scn_fs[-1].file.harmonic_frequencies.read(
                                locs)
                            freqs[i][j][k][lm] = freqs

    return pot, freqs


def _hrpot_spline_fitter(pot, thresh=-0.05):
    """ Get a physical hindered rotor potential via a series of spline fits
    """

    # Initialize a variable for the size of the potential
    lpot = len(pot)+1

    # Build a potential list from only successful calculations
    idx_success = []
    pot_success = [0.0]
    for idx in range(lpot):
        if pot[idx] < 600.:
            idx_success.append(idx)
            pot_success.append(pot[idx])
    idx_success.append(lpot)
    pot_success.append(pot[0])

    # Build a new potential list using a spline fit of the HR potential
    pot_spl = interp1d(
        numpy.array(idx_success), numpy.array(pot_success), kind='cubic')
    for idx in range(lpot):
        pot[idx] = pot_spl(idx)

    # Do second spline fit of only positive values if any negative values found
    if any(val < thresh for val in pot):
        print('Found pot vals below',
              ' {0} kcal. Refit w/ positives'.format(thresh))
        print('Potential before spline:', pot)
        x_pos = numpy.array([i for i in range(lpot)
                             if pot[i] >= thresh])
        y_pos = numpy.array([pot[i] for i in range(lpot)
                             if pot[i] >= thresh])
        pos_pot_spl = interp1d(x_pos, y_pos, kind='cubic')
        pot_pos_fit = []
        for idx in range(lpot):
            pot_pos_fit.append(pos_pot_spl(idx))

        print('Potential after spline:', pot_pos_fit)
        # Perform second check to see if negative potentials have been fixed
        if any(val < thresh for val in pot_pos_fit):
            print('Found values below {0} kcal again.'.format(thresh),
                  ' Trying linear interp of positive vals')
            neg_idxs = [i for i in range(lpot) if pot_pos_fit[i] < thresh]
            clean_pot = []
            for i in range(lpot):
                if i in neg_idxs:
                    # Find the indices for positive vals around negative value
                    idx_0 = i - 1
                    while idx_0 in neg_idxs:
                        idx_0 = idx_0 - 1
                    for j in range(i, lpot):
                        if pot_pos_fit[j] >= thresh:
                            idx_1 = j
                            break
                    # Get a new value for this point on the potential by
                    # doing a linear interp of positives
                    interp_val = (
                        pot_pos_fit[idx_0] * (1.0-((i-idx_0)/(idx_1-idx_0))) +
                        pot_pos_fit[idx_1] * ((i-idx_0)/(idx_1-idx_0))
                    )
                    clean_pot.append(interp_val)
                else:
                    clean_pot.append(pot[i])
            final_potential = clean_pot.copy()

        else:
            final_potential = pot_pos_fit.copy()

    else:
        final_potential = pot.copy()

    final_potential = final_potential[:-1]

    return final_potential
