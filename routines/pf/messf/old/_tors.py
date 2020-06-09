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
from lib import structure
from lib.struct import torsprep
from lib.phydat import phycon
from lib.submission import run_script
from lib.submission import DEFAULT_SCRIPT_DCT


# MESS strings
def write_1dhr_tors_mess_strings(harm_geo, spc_info, spc_dct_i, ts_bnd, zma,
                                 tors_names, tors_grids, tors_sym_nums,
                                 tors_cnf_save_path, min_ene,
                                 saddle=False, hind_rot_geo=None,
                                 frz_tors=False):
    """ Gather the 1DHR torsional data and gather them into a MESS file
    """

    # Build constraint dct
    if frz_tors:
        constraint_dct = structure.tors.build_constraint_dct(zma, tors_names)
    else:
        constraint_dct = None

    # Loop over the torsions
    hind_rot_str = ""
    proj_rotors_str = ""
    tors_info = zip(tors_names, tors_grids, tors_sym_nums)
    new_tors_sym_nums = []
    for tors_name_lst, tors_grid_lst, tors_sym in tors_info:

        # Grab zero elment because of formatting
        tors_name = tors_name_lst[0]
        tors_grid = tors_grid_lst[0]

        # Read the hindered rotor potential
        pot, _ = read_hr_pot(
            spc_info, [tors_name], tors_grid,
            tors_cnf_save_path, min_ene,
            saddle=saddle, read_freqs=False,
            frz_tors=frz_tors, constraint_dct=constraint_dct)

        # Build potential lst from only successful calculations
        pot = hrpot_spline_fitter(pot)
        # Get the HR groups and axis for the rotor
        group, axis, atm_key = set_groups_ini(
            zma, tors_name, ts_bnd, saddle)
        if saddle:
            group, axis, pot, sym_num = check_saddle_groups(
                zma, spc_dct_i, group, axis,
                pot, ts_bnd, tors_sym)
        else:
            sym_num = tors_sym
        group = list(numpy.add(group, 1))
        axis = list(numpy.add(axis, 1))
        if (atm_key+1) != axis[1]:
            axis.reverse()

        # Check for dummy transformations
        remdummy = check_dummy_trans(zma)

        # Write the MESS and ProjRot strings for the rotor
        hrgeo = harm_geo if hind_rot_geo else None
        hind_rot_str += mess_io.writer.rotor_hindered(
            group, axis, sym_num, pot,
            remdummy=remdummy, geom=hrgeo, use_quantum_weight=True)
        proj_rotors_str += projrot_io.writer.rotors(
            axis, group, remdummy=remdummy)
        new_tors_sym_nums.append(sym_num)

    return hind_rot_str, proj_rotors_str, new_tors_sym_nums


def write_mdhr_tors_mess_strings(geom, spc_info, sym_num, spc_dct_i,
                                 ts_bnd, zma,
                                 tors_name_grps, tors_grid_grps, tors_sym_nums,
                                 tors_cnf_save_path, min_ene,
                                 saddle=False, hind_rot_geo=None,
                                 vib_adiabatic=False):
    """ Gather the MDHR torsional data and gather them into a MESS file
    """

    # Loop over the torsion groups and get the int rot strings and potentials
    rotor_internal_str = ''
    proj_rotors_str = ''
    mdhr_dat_str_lst = []
    tors_idx = 0
    for tors_names, tors_grids in zip(tors_name_grps, tors_grid_grps):

        # Read the hindered rotor potential and add to master list
        vib_adiabatic=True
        hr_pot, hr_freqs = read_hr_pot(
            spc_info, tors_names, tors_grids,
            tors_cnf_save_path, min_ene,
            saddle=saddle, read_freqs=vib_adiabatic)

        # Write the MDHR potential file for each rotor set
        mdhr_dat_str = write_mdhr_dat_file(hr_pot, hr_freqs)

        # Check for dummy transformations
        remdummy = check_dummy_trans(zma)

        # Loop over the rotors in the group and write the internal rotor strs
        for tors_name, tors_grid in zip(tors_names, tors_grids):

            # Set pot to empty list (may need fix)
            pot = ()

            # Get the HR groups and axis for the rotor
            group, axis, atm_key = set_groups_ini(
                zma, tors_name, ts_bnd, saddle)
            if saddle:
                group, axis, pot, sym_num = check_saddle_groups(
                    zma, spc_dct_i, group, axis,
                    pot, ts_bnd, tors_sym_nums[tors_idx])
            else:
                sym_num = tors_sym_nums[tors_idx]
            group = list(numpy.add(group, 1))
            axis = list(numpy.add(axis, 1))
            if (atm_key+1) != axis[1]:
                axis.reverse()

            # Write the MESS and ProjRot strings for the rotor
            rotor_internal_str += mess_io.writer.mol_data.rotor_internal(
                group, axis, tors_sym_nums[tors_idx],
                rotor_id='', remdummy=remdummy,
                mass_exp_size=5, pot_exp_size=5,
                hmin=13, hmax=101,
                grid_size=100)
            proj_rotors_str += projrot_io.writer.rotors(
                axis, group, remdummy=remdummy)

            # Increment tors idx to keep track of the sym number
            tors_idx +=1

    return rotor_internal_str, proj_rotors_str, mdhr_dat_str


# Calculating certain quantities on the torsions
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
    pf_script_str = DEFAULT_SCRIPT_DCT['messpf']

    run_script(pf_script_str, pf_path)

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

    # print('tors_zpe test:', tors_zpe)
    return tors_zpe

