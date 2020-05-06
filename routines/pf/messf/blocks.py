""" drivers
"""

import os
import numpy
import autofile
import mess_io
from lib import filesys
from lib.amech_io import cleaner
from lib.phydat import phycon
from routines.pf.messf import models as pfmodels
from routines.pf.messf import _vib as vib
from routines.pf.messf import _tors as tors
from routines.pf.messf import _sym as sym
from routines.pf.messf import _fake as fake
from routines.pf.messf import _util as messfutil


def species_block(spc, spc_dct_i, spc_info, spc_model,
                  pf_levels, save_prefix):
    """ prepare the species input for pfpf
    """
    print('spc name in block', spc)
    print('path', save_prefix)
    # Unpack the models and levels
    [_, _, harm_level, vpt2_level, sym_level, tors_level] = pf_levels
    tors_model, vib_model, sym_model = spc_model

    # Lazy set for tors mod
    if tors_model == '1dhrf':
        tors_mod = ('1dhr', False)

    # Set theory filesystem used throughout
    thy_save_fs = autofile.fs.theory(save_prefix)

    # Set boolean to account for rad-rad reaction (not supported by vtst)
    rad_rad_ts = False
    if 'ts_' in spc:
        if spc_dct_i['rad_rad']:
            rad_rad_ts = True

    # Set the filesystem objects for various species models
    harmfs = set_model_filesys(
        thy_save_fs, spc_info, harm_level, saddle=('ts_' in spc))
    [harm_cnf_save_fs, _,
     harm_min_cnf_locs, harm_save_path] = harmfs
    if sym_level:
        symfs = set_model_filesys(
            thy_save_fs, spc_info, sym_level, saddle=('ts_' in spc))
        [sym_cnf_save_fs, _,
         sym_min_cnf_locs, _] = symfs
    if tors_level and not rad_rad_ts:
        torsfs = set_model_filesys(
            thy_save_fs, spc_info, tors_level[0], saddle=('ts_' in spc))
        [tors_cnf_save_fs, tors_cnf_save_path,
         tors_min_cnf_locs, tors_save_path] = torsfs
    if vpt2_level:
        vpt2fs = set_model_filesys(
            thy_save_fs, spc_info, vpt2_level, saddle=('ts_' in spc))
        [vpt2_cnf_save_fs, vpt2_cnf_save_path,
         vpt2_min_cnf_locs, vpt2_save_path] = vpt2fs

    # Check if any torsions to set the model

    # Set additional info for a saddle point
    saddle = False
    dist_names = []
    tors_names = []
    if 'ts_' in spc:
        saddle = True
        tors_names = spc_dct_i['amech_ts_tors_names']
        mig = 'migration' in spc_dct_i['class']
        elm = 'elimination' in spc_dct_i['class']
        if mig or elm:
            dist_names.append(spc_dct_i['dist_info'][0])
            dist_names.append(spc_dct_i['dist_info'][3])

    no_tors = not bool(tors.get_tors_names(spc_dct_i, tors_cnf_save_fs, saddle=saddle))
    print('no_tors test:', no_tors, spc_dct_i)
    # Set TS information
    frm_bnd_key, brk_bnd_key = messfutil.get_bnd_keys(spc_dct_i, saddle)

    # Initialize electronic energy levels
    elec_levels = messfutil.ini_elec_levels(spc_dct_i, spc_info)

    # Determine the species symmetry factor using the given model
    sym_factor = sym.symmetry_factor(
        sym_model, spc_dct_i, spc_info, dist_names,
        saddle, frm_bnd_key, brk_bnd_key, tors_names,
        tors_cnf_save_fs, tors_min_cnf_locs,
        sym_cnf_save_fs, sym_min_cnf_locs)

    # Build the species string and get the imaginary frequency
    # Initialize various auxiliary MESS data strings and names
    dat_str_dct = {}
    tau_dat_str = ''
    mdhr_dat_str = ''
    sct_dat_str = ''
    tau_dat_file_name = '{}_tau.dat'.format(spc[0])
    mdhr_dat_file_name = '{}_mdhr.dat'.format(spc[0])
    sct_dat_file_name = '{}_sct.dat'.format(spc[0])
    # Initialize various variables and strings
    symf = sym_factor
    imag = 0.0
    freqs = []
    hr_str = ""
    xmat = ()

    if messfutil.is_atom(harm_min_cnf_locs, harm_cnf_save_fs):
        mass = messfutil.atom_mass(harm_min_cnf_locs, harm_cnf_save_fs)
        spc_str = mess_io.writer.atom(
            mass, elec_levels)
    else:
        if (vib_model == 'harm' and tors_model == 'rigid') or rad_rad_ts:
            geo, freqs, imag = pfmodels.vib_harm_tors_rigid(
                spc_info, harm_min_cnf_locs, harm_cnf_save_fs, saddle=saddle)
            symf = sym_factor
        elif vib_model == 'harm' and (tors_model == '1dhr' or tors_model == '1dhrf'):
            frz_tors = True if tors_model == '1dhrf' else False
            if no_tors:
                geo, freqs, imag = pfmodels.vib_harm_tors_rigid(
                    spc_info, harm_min_cnf_locs, harm_cnf_save_fs, saddle=saddle)
                hr_str = ""
                symf = sym_factor
            else:
                geo, freqs, imag, hr_str, _ = pfmodels.vib_harm_tors_1dhr(
                    harm_min_cnf_locs, harm_cnf_save_fs,
                    tors_min_cnf_locs, tors_cnf_save_fs,
                    tors_save_path, tors_cnf_save_path,
                    spc_dct_i, spc_info,
                    frm_bnd_key, brk_bnd_key,
                    sym_factor, elec_levels,
                    saddle=saddle,
                    frz_tors=frz_tors)
                sym_nums = tors.get_tors_sym_nums(
                    spc_dct_i, tors_min_cnf_locs, tors_cnf_save_fs,
                    frm_bnd_key, brk_bnd_key, saddle=saddle)
                    # frm_bnd_key, brk_bnd_key, saddle=False)
                symf = sym_factor
                for num in sym_nums:
                    symf /= num
        elif vib_model == 'harm' and tors_model == 'mdhr':
            if no_tors:
                geo, freqs, imag = pfmodels.vib_harm_tors_rigid(
                    spc_info, harm_min_cnf_locs, harm_cnf_save_fs, saddle=saddle)
                hr_str = ""
                symf = sym_factor
            else:
                geo, freqs, imag, core_hr_str, _, mdhr_dat_str = pfmodels.vib_harm_tors_mdhr(
                    harm_min_cnf_locs, harm_cnf_save_fs,
                    tors_min_cnf_locs, tors_cnf_save_fs,
                    tors_save_path, tors_cnf_save_path,
                    spc_dct_i, spc_info,
                    frm_bnd_key, brk_bnd_key,
                    sym_factor, elec_levels,
                    tors_mod=tors_mod,
                    saddle=saddle)
                sym_nums = tors.get_tors_sym_nums(
                    spc_dct_i, tors_min_cnf_locs, tors_cnf_save_fs,
                    frm_bnd_key, brk_bnd_key, saddle=saddle)
                    #frm_bnd_key, brk_bnd_key, saddle=False)
                symf = sym_factor
                for num in sym_nums:
                    symf /= num
                mdhr_str = mess_io.writer.mol_data.core_multirotor(
                    geo, sym_factor, mdhr_dat_file_name, core_hr_str,
                    interp_emax=100, quant_lvl_emax=9)  # , forceq=False)
        elif vib_model == 'harm' and tors_model == 'tau':
            if no_tors:
                geo, freqs, imag = pfmodels.vib_harm_tors_rigid(
                    spc_info, harm_min_cnf_locs, harm_cnf_save_fs, saddle=saddle)
                hr_str = ""
                symf = sym_factor
            else:
                _, _, _, _, _, _, mc_str, tau_dat_str = pfmodels.vib_harm_tors_tau(
                    harm_min_cnf_locs, harm_cnf_save_fs,
                    tors_min_cnf_locs, tors_cnf_save_fs,
                    tors_save_path, tors_cnf_save_path,
                    spc_dct_i, spc_info,
                    frm_bnd_key, brk_bnd_key,
                    sym_factor, elec_levels,
                    tau_dat_file_name,
                    hind_rot_geo=False,
                    saddle=saddle)
                    #saddle=False)
                print('HARM and TAU combination is not yet implemented')
        elif vib_model == 'tau' and tors_model == 'tau':
            mc_str, tau_dat_str = pfmodels.vib_tau_tors_tau(
                tors_min_cnf_locs, tors_cnf_save_fs,
                spc_dct_i,
                frm_bnd_key, brk_bnd_key,
                sym_factor, elec_levels,
                tors_save_path,  # save_prefix,
                tau_dat_file_name,
                saddle=saddle)
                #saddle=False)
        elif vib_model == 'vpt2' and tors_model == 'rigid':
            geo, freqs, imag, xmat, _ = pfmodels.vib_vpt2_tors_rigid(
                harm_min_cnf_locs, harm_cnf_save_fs,
                vpt2_min_cnf_locs, vpt2_cnf_save_fs,
                saddle=saddle)
        elif vib_model == 'vpt2' and tors_model == '1dhr':
            print('VPT2 and 1DHR combination is not yet implemented')
        elif vib_model == 'vpt2' and tors_model == 'tau':
            print('VPT2 and TAU combination is not yet implemented')

        # Write the species string for the molecule
        if tors_model == 'tau':
            spc_str = mc_str
        else:
            if tors_model == 'mdhr':
                core = mdhr_str
            else:
                core = mess_io.writer.core_rigidrotor(geo, symf)
            spc_str = mess_io.writer.molecule(
                core, freqs, elec_levels,
                hind_rot=hr_str, xmat=xmat)

        # Combine various data strings into a dct
        dat_str_dct = {
            'tau': (tau_dat_str, tau_dat_file_name),
            'mdhr': (mdhr_dat_str, mdhr_dat_file_name),
            'sct': (sct_dat_str, sct_dat_file_name)
        }

    return spc_str, dat_str_dct, imag


def vtst_with_no_saddle_block(
        ts_dct, ts_label, reac_label, prod_label,
        spc_ene, projrot_script_str, multi_info):
    """ prepare the pf input string for a variational TS that does not have
    a saddle point. Do it by calling the species block for each grid point
    in the scan file system
    """

    ts_info = ['', ts_dct['chg'], ts_dct['mul']]
    multi_level = filesys.inf.modify_orb_restrict(ts_info, multi_info)

    rxn_save_path = ts_dct['rxn_fs'][3]
    thy_save_fs = autofile.fs.theory(rxn_save_path)
    thy_save_fs[-1].create(multi_level[1:4])
    thy_save_path = thy_save_fs[-1].path(multi_level[1:4])
    scn_save_fs = autofile.fs.scan(thy_save_path)

    # Read the scan save filesystem to get the molecular info
    sym_factor = 1.
    irc_pt_strs = []
    proj_rotors_str = ''
    pot = []

    elec_levels = [[0., ts_dct['mul']]]
    grid = ts_dct['grid']
    grid = numpy.append(grid[0], grid[1])
    dist_name = ts_dct['dist_info'][0]

    # Read the infinite separation energy
    inf_locs = [[dist_name], [1000.]]
    inf_sep_ene = scn_save_fs[-1].file.energy.read(inf_locs)

    grid[::-1].sort()
    for idx, grid_val in enumerate(grid):
        # Set the filesystem locators for each grid point
        locs = [[dist_name], [grid_val]]

        # Get geometry, energy, and vibrational freqs
        if not scn_save_fs[-1].file.geometry.exists(locs):
            continue
        else:
            geom = scn_save_fs[-1].file.geometry.read(locs)
        if not scn_save_fs[-1].file.energy.exists(locs):
            continue
        else:
            ene = scn_save_fs[-1].file.energy.read(locs)
        if not scn_save_fs[-1].file.hessian.exists(locs):
            continue
        else:
            hess = scn_save_fs[-1].file.hessian.read(locs)
            scn_save_path = scn_save_fs[-1].path(locs)
            freqs, _, _ = pfmodels.projrot_freqs_1(
                geom, hess, pot,
                proj_rotors_str, projrot_script_str,
                scn_save_path, saddle=True)

        # Calculate standard harmonic ZPVE
        zpe = sum(freqs)*phycon.WAVEN2KCAL/2.

        # Calcuate infinite separation ZPVE
        # Assumes the ZPVE = ZPVE(1st grid pt) as an approximation
        if idx == 0:
            rct_zpe = zpe

        # Calculate the reference energies
        ene_rel = (ene - inf_sep_ene) * phycon.EH2KCAL
        zpe_rel = zpe - rct_zpe
        eref_abs = ene_rel + zpe_rel + spc_ene

        # Iniialize the header of the string
        irc_pt_str = '!----------------------------------------------- \n'
        irc_pt_str += '! IRC Point {0}\n'.format(str(idx+1))

        # Write the MESS string for the molecule section for each irc point
        core = mess_io.writer.mol_data.core_rigidrotor(
            geom, sym_factor, interp_emax='')
        irc_pt_str += mess_io.writer.species.molecule(
            core, freqs, elec_levels,
            hind_rot='', xmat=None, rovib_coups='', rot_dists='')

        # Append the zero energy for the molecule
        irc_pt_str += ('    ZeroEnergy[kcal/mol]      ',
                       '{0:<8.2f}\n'.format(eref_abs))
        if grid_val != grid[-1]:
            irc_pt_str += 'End \n'

        # Append string to list
        irc_pt_strs.append(irc_pt_str)

    # Write the MESS string for the entire variational section
    variational_str = mess_io.writer.rxnchan.ts_variational(
        ts_label, reac_label, prod_label, irc_pt_strs)

    return variational_str


def vtst_saddle_block(ts_dct, ene_thy_level, geo_thy_level,
                      ts_label, reac_label, prod_label, first_ground_ene):
    """ prepare the pf input string for a variational TS where there is a
        saddle point on the MEP.
        In this case, there is limited torsional information.
    """
    [geo_thy_level, ene_thy_level, _, _, _, _] = pf_levels

    irc_idxs = ts_dct['irc_idxs']
    ts_info = ['', ts_dct['chg'], ts_dct['mul']]

    # Build the TS scan file system
    scn_save_fs, _, _, _ = irc.ts_scn_fs(
        ts_dct, ts_info, geo_thy_level)

    # Set the distance name for the reaction coordinate
    dist_name = 'RC'

    # Loop over all the points of the irc and build MESS strings
    full_irc_str = []
    for idx in irc_idxs:

        # Set the filesystem locators for each grid point
        locs = [[dist_name], [idx]]
        print(scn_save_fs[-1].path(locs))

        # Get geometry, energy, vibrational freqs, and zpe
        if scn_save_fs[-1].file.geometry.exists(locs):
            geom = scn_save_fs[-1].file.geometry.read(locs)
        else:
            print('no geom')
            continue
        if scn_save_fs[-1].file.energy.exists(locs):
            if ene_thy_level == geo_thy_level:
                ene = scn_save_fs[-1].file.energy.read(locs)
                print('ene', ene)
            else:
                scn_save_path = scn_save_fs[-1].path(locs)
                sp_save_fs = autofile.fs.single_point(scn_save_path)
                sp_level = filesys.inf.modify_orb_restrict(ts_info, ene_thy_level)
                if sp_save_fs[-1].file.energy.exists(sp_level[1:4]):
                    ene = sp_save_fs[-1].file.energy.read(sp_level[1:4])
                    print('ene-high', ene)
                else:
                    print('no energy')
                    continue
        else:
            print('no energy')
            continue
        if scn_save_fs[-1].file.hessian.exists(locs):
            proj_rotors_str = ''
            hess = scn_save_fs[-1].file.hessian.read(locs)
            scn_save_path = scn_save_fs[-1].path(locs)
            freqs, _, _ = vib.projrot_freqs_1(
                geom, hess,
                proj_rotors_str,
                scn_save_path, pot=False, saddle=True)
            zpe = sum(freqs)*phycon.WAVEN2KCAL/2.
        else:
            print('no hessian')
            continue

        # Calculate the relative energy
        erel = ene * phycon.EH2KCAL + zpe - first_ground_ene

        # Set the symmetry factor and electronic levels
        sym_factor = 1.0
        elec_levels = [[0., ts_dct['mul']]]

        # Iniialize the header of the string
        irc_pt_str = '!-----------------------------------------------\n'
        irc_pt_str += '! IRC Point {0}\n'.format(str(int(idx)))

        # Write the molecule section for each irc point
        core = mess_io.writer.mol_data.core_rigidrotor(
            geom, sym_factor, interp_emax='')
        irc_pt_str += mess_io.writer.species.molecule(
            core, freqs, elec_levels,
            hind_rot='', xmat=None, rovib_coups='', rot_dists='')

        # Append the zero point energy string for the point
        irc_pt_str += ('  ZeroEnergy[kcal/mol]      '
                       '{0:<8.2f}'.format(erel))

        # Append string to list
        full_irc_str.append(irc_pt_str)

    # Write the MESS string for the variational sections
    variational_str = mess_io.writer.rxnchan.ts_variational(
        ts_label, reac_label, prod_label, full_irc_str)

    return variational_str


def pst_block(spc_dct_i, spc_dct_j, spc_model, pf_levels,
              spc_save_fs, pst_params=(1.0, 6)):
    """ prepare a Phase Space Theory species block
    """

    # Unpack the models and levels
    [_, _, harm_level, _, sym_level, tors_level] = pf_levels
    tors_model, vib_model, sym_model = spc_model

    # prepare the four sets of file systems
    spc_info_i = (spc_dct_i['ich'], spc_dct_i['chg'], spc_dct_i['mul'])
    spc_info_j = (spc_dct_j['ich'], spc_dct_j['chg'], spc_dct_j['mul'])
    spc_save_fs[-1].create(spc_info_i)
    spc_save_fs[-1].create(spc_info_j)
    save_path_i = spc_save_fs[-1].path(spc_info_i)
    save_path_j = spc_save_fs[-1].path(spc_info_j)

    # Set theory filesystem used throughout
    thy_save_fs_i = autofile.fs.theory(save_path_i)
    thy_save_fs_j = autofile.fs.theory(save_path_j)

    # Set the filesystem objects for the two species
    harmfs_i = set_model_filesys(
        thy_save_fs_i, spc_info_i, harm_level, saddle=False)
    harmfs_j = set_model_filesys(
        thy_save_fs_j, spc_info_j, harm_level, saddle=False)
    [harm_cnf_save_fs_i, _,
     harm_min_cnf_locs_i, _] = harmfs_i
    [harm_cnf_save_fs_j, _,
     harm_min_cnf_locs_j, _] = harmfs_j

    if sym_level:
        symfs_i = set_model_filesys(
            thy_save_fs_i, spc_info_i, sym_level, saddle=False)
        symfs_j = set_model_filesys(
            thy_save_fs_j, spc_info_j, sym_level, saddle=False)
        [sym_cnf_save_fs_i, _,
         sym_min_cnf_locs_i, _] = symfs_i
        [sym_cnf_save_fs_j, _,
         sym_min_cnf_locs_j, _] = symfs_j

    if tors_level:
        torsfs_i = set_model_filesys(
            thy_save_fs_i, spc_info_i, tors_level[0], saddle=False)
        torsfs_j = set_model_filesys(
            thy_save_fs_j, spc_info_j, tors_level[0], saddle=False)
        [tors_cnf_save_fs_i, tors_cnf_save_path_i,
         tors_min_cnf_locs_i, tors_save_path_i] = torsfs_i
        [tors_cnf_save_fs_j, tors_cnf_save_path_j,
         tors_min_cnf_locs_j, tors_save_path_j] = torsfs_j

    # if vpt2_level:
    #     vpt2fs_i = set_model_filesys(
    #         thy_save_fs_i, spc_info_i, vpt2_level, saddle=False)
    #     vpt2fs_j = set_model_filesys(
    #         thy_save_fs_j, spc_info_j, vpt2_level, saddle=False)
    #     [vpt2_cnf_save_fs_i, vpt2_cnf_save_path_i,
    #      vpt2_min_cnf_locs_i, vpt2_save_path_i] = vpt2fs_i
    #     [vpt2_cnf_save_fs_j, vpt2_cnf_save_path_j,
    #      vpt2_min_cnf_locs_j, vpt2_save_path_j] = vpt2fs_j

    # are there any torsions
    no_tors_i = not bool(tors.get_tors_names(spc_dct_i, tors_cnf_save_fs_i, saddle=False))
    no_tors_j = not bool(tors.get_tors_names(spc_dct_j, tors_cnf_save_fs_j, saddle=False))

    # Get the combined electronic energy levels
    elec_levels = messfutil.combine_elec_levels(spc_dct_i, spc_dct_j)

    # Determine the species symmetry factor using the given model
    dist_names = []
    tors_names = []
    frm_bnd_key = []
    brk_bnd_key = []
    saddle = False
    sym_factor_i = sym.symmetry_factor(
        sym_model, spc_dct_i, spc_info_i, dist_names,
        saddle, frm_bnd_key, brk_bnd_key, tors_names,
        tors_cnf_save_fs_i, tors_min_cnf_locs_i,
        sym_cnf_save_fs_i, sym_min_cnf_locs_i)
    sym_factor_j = sym.symmetry_factor(
        sym_model, spc_dct_j, spc_info_j, dist_names,
        saddle, frm_bnd_key, brk_bnd_key, tors_names,
        tors_cnf_save_fs_j, tors_min_cnf_locs_j,
        sym_cnf_save_fs_j, sym_min_cnf_locs_j)
    sym_factor = sym_factor_i * sym_factor_j

    # Get the stoichiometry
    stoich = messfutil.get_stoich(
        harm_min_cnf_locs_i, harm_min_cnf_locs_j,
        harm_cnf_save_fs_i, harm_cnf_save_fs_j)

    spc_str = ''
    if vib_model == 'harm' and tors_model == 'rigid':
        if messfutil.is_atom(harm_min_cnf_locs_i, harm_cnf_save_fs_i):
            freqs_i = ()
        else:
            geo_i, freqs_i, _ = pfmodels.vib_harm_tors_rigid(
                spc_info_i, harm_min_cnf_locs_i,
                harm_cnf_save_fs_i, saddle=False)
        if messfutil.is_atom(harm_min_cnf_locs_j, harm_cnf_save_fs_j):
            freqs_j = ()
        else:
            geo_j, freqs_j, _ = pfmodels.vib_harm_tors_rigid(
                spc_info_j, harm_min_cnf_locs_j,
                harm_cnf_save_fs_j, saddle=False)
        if harm_min_cnf_locs_i is not None:
            geo_i = harm_cnf_save_fs_i[-1].file.geometry.read(
                harm_min_cnf_locs_i)
        if harm_min_cnf_locs_j is not None:
            geo_j = harm_cnf_save_fs_j[-1].file.geometry.read(
                harm_min_cnf_locs_j)
        freqs = freqs_i + freqs_j
        hind_rot_str = ""

    if vib_model == 'harm' and (tors_model == '1dhr' or tors_model == '1dhrf'):
        if messfutil.is_atom(harm_min_cnf_locs_i, harm_cnf_save_fs_i):
            geo_i = harm_cnf_save_fs_i[-1].file.geometry.read(
                harm_min_cnf_locs_i)
            freqs_i = []
            hr_str_i = ''
            symf_i = sym_factor_i
        elif no_tors_i:
            geo_i = harm_cnf_save_fs_i[-1].file.geometry.read(
                harm_min_cnf_locs_i)
            _, freqs_i, _ = pfmodels.vib_harm_tors_rigid(
                spc_info_i, harm_min_cnf_locs_i,
                harm_cnf_save_fs_i, saddle=False)
            hr_str_i = ''
            symf_i = sym_factor_i
        else:
            frz_tors = True if tors_model == '1dhrf' else False
            geo_i, freqs_i, _, hr_str_i, _ = pfmodels.vib_harm_tors_1dhr(
                harm_min_cnf_locs_i, harm_cnf_save_fs_i,
                tors_min_cnf_locs_i, tors_cnf_save_fs_i,
                tors_save_path_i, tors_cnf_save_path_i,
                spc_dct_i, spc_info_i,
                frm_bnd_key, brk_bnd_key,
                sym_factor_i, elec_levels,
                saddle=False,
                frz_tors = frz_tors)
            sym_nums_i = tors.get_tors_sym_nums(
                spc_dct_i, tors_min_cnf_locs_i, tors_cnf_save_fs_i,
                frm_bnd_key, brk_bnd_key, saddle=False)
            symf_i = sym_factor_i
            for num in sym_nums_i:
                symf_i /= num
        if messfutil.is_atom(harm_min_cnf_locs_j, harm_cnf_save_fs_j):
            geo_j = harm_cnf_save_fs_j[-1].file.geometry.read(
                harm_min_cnf_locs_j)
            freqs_j = []
            hr_str_j = ''
            symf_j = sym_factor_j
        elif no_tors_j:
            geo_j = harm_cnf_save_fs_j[-1].file.geometry.read(
                harm_min_cnf_locs_j)
            hr_str_j = ''
            symf_j = sym_factor_j
            _, freqs_j, _ = pfmodels.vib_harm_tors_rigid(
                spc_info_j, harm_min_cnf_locs_j,
                harm_cnf_save_fs_j, saddle=False)
        else:
            frz_tors = True if tors_model == '1dhrf' else False
            geo_j, freqs_j, _, hr_str_j, _ = pfmodels.vib_harm_tors_1dhr(
                harm_min_cnf_locs_j, harm_cnf_save_fs_j,
                tors_min_cnf_locs_j, tors_cnf_save_fs_j,
                tors_save_path_j, tors_cnf_save_path_j,
                spc_dct_j, spc_info_j,
                frm_bnd_key, brk_bnd_key,
                sym_factor_j, elec_levels,
                saddle=False,
                frz_tors = frz_tors)
            sym_nums_j = tors.get_tors_sym_nums(
                spc_dct_j, tors_min_cnf_locs_j, tors_cnf_save_fs_j,
                frm_bnd_key, brk_bnd_key, saddle=False)
            symf_j = sym_factor_j
            for num in sym_nums_j:
                symf_j /= num
        freqs = list(freqs_i) + list(freqs_j)
        hind_rot_str = hr_str_i + hr_str_j
        sym_factor = symf_i * symf_j

    # Write the MESS input strings
    core = mess_io.writer.core_phasespace(
        geo_i, geo_j, sym_factor, stoich,
        pot_prefactor=pst_params[0], pot_power_exp=pst_params[1])
    spc_str = mess_io.writer.molecule(
        core, freqs, elec_levels,
        hind_rot=hind_rot_str)

    return spc_str


def fake_species_block(
        spc_dct_i, spc_dct_j, spc_info_i, spc_info_j,
        spc_model, pf_levels,
        save_prefix_i, save_prefix_j):
    """ prepare a fake species block corresponding to the
        van der Waals well between two fragments
    """
    [_, _, harm_level, _, sym_level, tors_level] = pf_levels
    tors_model, vib_model, sym_model = spc_model

    # prepare the four sets of file systems
    har_levelp_i = filesys.inf.modify_orb_restrict(
        spc_info_i, harm_level)
    har_levelp_j = filesys.inf.modify_orb_restrict(
        spc_info_j, harm_level)

    # Set theory filesystem used throughout
    thy_save_fs_i = autofile.fs.theory(save_prefix_i)
    thy_save_fs_j = autofile.fs.theory(save_prefix_j)

    # Set the filesystem objects for the two species
    harmfs_i = set_model_filesys(
        thy_save_fs_i, spc_info_i, harm_level, saddle=False)
    harmfs_j = set_model_filesys(
        thy_save_fs_j, spc_info_j, harm_level, saddle=False)
    [harm_cnf_save_fs_i, _,
     harm_min_cnf_locs_i, _] = harmfs_i
    [harm_cnf_save_fs_j, _,
     harm_min_cnf_locs_j, _] = harmfs_j

    if sym_level:
        symfs_i = set_model_filesys(
            thy_save_fs_i, spc_info_i, sym_level, saddle=False)
        symfs_j = set_model_filesys(
            thy_save_fs_j, spc_info_j, sym_level, saddle=False)
        [sym_cnf_save_fs_i, _,
         sym_min_cnf_locs_i, _] = symfs_i
        [sym_cnf_save_fs_j, _,
         sym_min_cnf_locs_j, _] = symfs_j

    if tors_level:
        torsfs_i = set_model_filesys(
            thy_save_fs_i, spc_info_i, tors_level[0], saddle=False)
        torsfs_j = set_model_filesys(
            thy_save_fs_j, spc_info_j, tors_level[0], saddle=False)
        [tors_cnf_save_fs_i, tors_cnf_save_path_i,
         tors_min_cnf_locs_i, tors_save_path_i] = torsfs_i
        [tors_cnf_save_fs_j, tors_cnf_save_path_j,
         tors_min_cnf_locs_j, tors_save_path_j] = torsfs_j

    spc_str = ''

    # are there any torsion
    no_tors_i = not bool(tors.get_tors_names(spc_dct_i, tors_cnf_save_fs_i, saddle=False))
    no_tors_j = not bool(tors.get_tors_names(spc_dct_j, tors_cnf_save_fs_j, saddle=False))

    # Get the combined electronic energy levels
    elec_levels = messfutil.combine_elec_levels(spc_dct_i, spc_dct_j)

    # Determine the species symmetry factor using the given model
    dist_names = []
    tors_names = []
    frm_bnd_key = []
    brk_bnd_key = []
    saddle=False
    sym_factor_i = sym.symmetry_factor(
        sym_model, spc_dct_i, spc_info_i, dist_names,
        saddle, frm_bnd_key, brk_bnd_key, tors_names,
        tors_cnf_save_fs_i, tors_min_cnf_locs_i,
        sym_cnf_save_fs_i, sym_min_cnf_locs_i)
    sym_factor_j = sym.symmetry_factor(
        sym_model, spc_dct_j, spc_info_j, dist_names,
        saddle, frm_bnd_key, brk_bnd_key, tors_names,
        tors_cnf_save_fs_j, tors_min_cnf_locs_j,
        sym_cnf_save_fs_j, sym_min_cnf_locs_j)
    sym_factor = sym_factor_i * sym_factor_j

    # Get the freqs
    freqs = fake.set_fake_freqs(
        harm_min_cnf_locs_i, harm_min_cnf_locs_j,
        harm_cnf_save_fs_i, harm_cnf_save_fs_j)
    geo = fake.combine_geos_in_fake_well(
        harm_min_cnf_locs_i, harm_min_cnf_locs_j,
        harm_cnf_save_fs_i, harm_cnf_save_fs_j)

    if vib_model == 'harm' and tors_model == 'rigid':
        if messfutil.is_atom(harm_min_cnf_locs_i, harm_cnf_save_fs_i):
            freqs_i = ()
        else:
            _, freqs_i, _ = pfmodels.vib_harm_tors_rigid(
                spc_info_i, harm_min_cnf_locs_i,
                harm_cnf_save_fs_i, saddle=False)
        if messfutil.is_atom(harm_min_cnf_locs_j, harm_cnf_save_fs_j):
            freqs_j = ()
        else:
            _, freqs_j, _ = pfmodels.vib_harm_tors_rigid(
                spc_info_j, harm_min_cnf_locs_j,
                harm_cnf_save_fs_j, saddle=False)
        freqs = freqs + freqs_i + freqs_j
        hind_rot_str = ""

    if vib_model == 'harm' and (tors_model == '1dhr' or tors_model == '1dhrf'):
        if messfutil.is_atom(harm_min_cnf_locs_i, harm_cnf_save_fs_i):
            freqs_i = ()
            hr_str_i = ''
            symf_i = sym_factor_i
        elif no_tors_i:
            hr_str_i = ''
            symf_i = sym_factor_i
            _, freqs_i, _ = pfmodels.vib_harm_tors_rigid(
                spc_info_i, harm_min_cnf_locs_i,
                harm_cnf_save_fs_i, saddle=False)
        else:
            frz_tors = True if tors_model == '1dhrf' else False
            _, freqs_i, _, hr_str_i, _ = pfmodels.vib_harm_tors_1dhr(
                harm_min_cnf_locs_i, harm_cnf_save_fs_i,
                tors_min_cnf_locs_i, tors_cnf_save_fs_i,
                tors_save_path_i, tors_cnf_save_path_i,
                spc_dct_i, spc_info_i,
                frm_bnd_key, brk_bnd_key,
                sym_factor_i, elec_levels,
                saddle=False,
                frz_tors=frz_tors)
            sym_nums_i = tors.get_tors_sym_nums(
                spc_dct_i, tors_min_cnf_locs_i, tors_cnf_save_fs_i,
                frm_bnd_key, brk_bnd_key, saddle=False)
            symf_i = sym_factor_i
            for num_i in sym_nums_i:
                symf_i /= num_i

        if messfutil.is_atom(harm_min_cnf_locs_j, harm_cnf_save_fs_j):
            freqs_j = ()
            hr_str_j = ''
            symf_j = sym_factor_j
        elif no_tors_j:
            hr_str_j = ''
            symf_j = sym_factor_j
            _, freqs_j, _ = pfmodels.vib_harm_tors_rigid(
                spc_info_j, harm_min_cnf_locs_j,
                harm_cnf_save_fs_j, saddle=False)
        else:
            frz_tors = True if tors_model == '1dhrf' else False
            _, freqs_j, _, hr_str_j, _ = pfmodels.vib_harm_tors_1dhr(
                harm_min_cnf_locs_j, harm_cnf_save_fs_j,
                tors_min_cnf_locs_j, tors_cnf_save_fs_j,
                tors_save_path_j, tors_cnf_save_path_j,
                spc_dct_j, spc_info_j,
                frm_bnd_key, brk_bnd_key,
                sym_factor_j, elec_levels,
                saddle=False,
                frz_tors=frz_tors)
            sym_nums_j = tors.get_tors_sym_nums(
                spc_dct_j, tors_min_cnf_locs_j, tors_cnf_save_fs_j,
                frm_bnd_key, brk_bnd_key, saddle=False)
            symf_j = sym_factor_j
            for num_j in sym_nums_j:
                symf_j /= num_j
        freqs = freqs + freqs_i + freqs_j
        hind_rot_str = hr_str_i + hr_str_j
        sym_factor = symf_i * symf_j

    core = mess_io.writer.core_rigidrotor(geo, sym_factor)
    spc_str = mess_io.writer.molecule(
        core, freqs, elec_levels,
        hind_rot=hind_rot_str)

    return spc_str


# Write blocks for MESSPF
def get_pf_header(temps):
    """ prepare partition function header string
    """
    global_pf_str = mess_io.writer.global_pf(
        temps, rel_temp_inc=0.001, atom_dist_min=0.6)
    return global_pf_str


def get_pf_input(spc, spc_str, global_pf_str, zpe_str):
    """ prepare the full pf input string for running pfpf
    """

    # create a pfpf input file
    spc_head_str = 'Species ' + '\n' + spc
    pf_inp_str = '\n'.join(
        [global_pf_str,
         spc_head_str,
         spc_str,
         zpe_str,
         '\n'])
    print(pf_inp_str)
    return pf_inp_str


def write_pf_input(pf_inp_str, data_str_dct, pf_path):
    """ write the pf.inp file
    """
    # Write the data file
    for data in data_str_dct.values():
        string, name = data
        if string:
            data_file_path = os.path.join(pf_path, name)
            with open(data_file_path, 'w') as data_file:
                data_file.write(string)

    # Write the MESSPF input file
    pf_inp_str = cleaner.remove_whitepsace(pf_inp_str)
    with open(os.path.join(pf_path, 'pf.inp'), 'w') as pf_file:
        pf_file.write(pf_inp_str)



####################
# Helper functions #
####################

def set_model_filesys(thy_save_fs, spc_info, level, saddle=False):
    """ Gets filesystem objects for torsional calculations
    """
    # Set the level for the model
    levelp = filesys.inf.modify_orb_restrict(spc_info, level)

    # Get the save fileystem path
    save_path = thy_save_fs[-1].path(levelp[1:4])
    if saddle:
        save_fs = autofile.fs.ts(save_path)
        save_fs[0].create()
        save_path = save_fs[0].path()

    # Get the fs object and the locs
    cnf_save_fs = autofile.fs.conformer(save_path)
    min_cnf_locs = filesys.mincnf.min_energy_conformer_locators(cnf_save_fs)

    # Get the save path for the conformers
    if min_cnf_locs:
        cnf_save_path = cnf_save_fs[-1].path(min_cnf_locs)
    else:
        cnf_save_path = ''

    return cnf_save_fs, cnf_save_path, min_cnf_locs, save_path
