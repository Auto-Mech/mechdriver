"""
Find a TS from the grid as well as associated vdW wells
"""

import automol
import elstruct
from lib.reaction import grid as rxngrid
from lib.reaction import ts as lts
from lib.phydat import phycon
from lib.runner import driver
from lib.filesystem import minc as fsmin
from routines.es import conformer
from routines.es import variational
from routines.es import scan


def find_ts(spc_dct, ts_dct, ts_zma, ts_info,
            mod_ini_thy_info, mod_thy_info,
            multi_opt_info, multi_sp_info,
            thy_save_fs,
            cnf_run_fs, cnf_save_fs,
            scn_run_fs, scn_save_fs,
            run_fs,
            ts_save_fs, ts_save_path,
            run_prefix, save_prefix,
            opt_script_str,
            overwrite, vrc_dct,
            rad_rad_ts='vtst',
            **opt_kwargs):
    """ find the ts geometry
    """

    # Set various TS information using the dictionary
    typ = ts_dct['class']
    dist_info = ts_dct['dist_info']
    grid = ts_dct['grid']
    rad_rad = bool('radical radical' in typ)
    low_spin = bool('low' in typ)

    # Unpack the dist info
    dist_name, _, update_guess, brk_name = dist_info

    # Get TS from filesys or get it from some procedure
    min_cnf_locs = fsmin.min_energy_conformer_locators(cnf_save_fs)
    if min_cnf_locs and not overwrite:
        # ts_class, ts_original_zma, ts_tors_names, ts_dist_info
        # geo, zma, final_dist = check_filesys_for_ts(
        #     ts_dct, ts_zma, cnf_save_fs, overwrite,
        #     typ, dist_info, dist_name, bkp_ts_class_data)
        zma = cnf_save_fs[-1].file.zmatrix.read(min_cnf_locs)

        # Add an angle check which is added to spc dct for TS (crap code...)
        vals = automol.zmatrix.values(zma)
        final_dist = vals[dist_name]
        dist_info[1] = final_dist
        angle = lts.check_angle(
            ts_dct['zma'],
            ts_dct['dist_info'],
            ts_dct['class'])
        ts_dct['dist_info'][1] = final_dist
        ts_dct['dist_info'].append(angle)

        # Setting to true for now
        ts_found = True
    else:
        if rad_rad and low_spin and 'elimination' not in ts_dct['class']:
            print('Running Scan for Barrierless TS:')
            ts_found = find_barrierless_transition_state(
                ts_info, ts_zma, ts_dct, spc_dct,
                grid,
                dist_name,
                rxn_run_path, rxn_save_path,
                rad_rad_ts,
                mod_ini_thy_info, mod_thy_info,
                multi_opt_info, multi_sp_info,
                run_prefix, save_prefix,
                scn_run_fs, scn_save_fs,
                overwrite, vrc_dct,
                update_guess, **opt_kwargs)
            # Have code analyze the path and switch to a sadpt finder if needed
        else:
            print('Running Scan for Fixed TS:')
            max_zma = run_sadpt_scan(
                typ, grid, dist_name, brk_name, ts_zma, ts_info,
                mod_thy_info,
                scn_run_fs, scn_save_fs, opt_script_str,
                overwrite, update_guess, **opt_kwargs)
            ts_found = find_sadpt_transition_state(
                opt_script_str,
                run_fs,
                max_zma,
                ts_info,
                mod_thy_info,
                overwrite,
                ts_save_path,
                ts_save_fs,
                dist_name,
                dist_info,
                thy_save_fs,
                cnf_run_fs,
                cnf_save_fs,
                **opt_kwargs)

    return ts_found


def find_barrierless_transition_state(ts_info, ts_zma, ts_dct, spc_dct,
                                      grid,
                                      dist_name,
                                      rad_rad_ts,
                                      ini_thy_info, thy_info,
                                      mod_multi_opt_info, mod_multi_sp_info,
                                      run_prefix, save_prefix,
                                      scn_run_fs, scn_save_fs,
                                      overwrite, vrc_dct,
                                      update_guess, **opt_kwargs):
    """ Run TS finder for barrierless reactions
    """

    ts_formula = automol.geom.formula(automol.zmatrix.geometry(ts_zma))
    [grid1, grid2] = grid

    # Run PST, VTST, VRC-TST based on RAD_RAD_TS model
    if rad_rad_ts.lower() == 'pst':
        ts_found = True
        print('Phase Space Theory Used, No ES calculations are needed')
    elif rad_rad_ts.lower() == 'vtst':
        print('Beginning Calculations for VTST Treatments')
        ts_found = variational.vtst.run_vtst_scan(
            ts_zma, ts_formula, ts_info, ts_dct, spc_dct,
            high_mul, grid1, grid2, dist_name,
            multi_level, multi_sp_info,
            multi_info, ini_thy_info, thy_info,
            run_prefix, save_prefix, scn_run_fs, scn_save_fs,
            overwrite, update_guess, **opt_kwargs)
        if ts_found:
            print('Scans for VTST succeeded')
        else:
            print('Scans for VTST failed')
    elif rad_rad_ts.lower() == 'vrctst':
        print('Beginning Calculations for VRC-TST Treatments')
        ts_found = variational.vrctst.calc_vrctst_flux(
            ts_zma, ts_formula, ts_info, ts_dct, spc_dct,
            ts_dct['high_mul'], grid1, grid2, dist_name,
            multi_level, mod_multi_opt_info, mod_multi_sp_info,
            mod_ini_thy_info, mod_thy_info,
            thy_run_path, thy_save_path,
            overwrite, update_guess,
            run_prefix, save_prefix,
            vrc_dct,
            corr_pot=True)
        if ts_found:
            print('VaReCoF run successful and flux file was obtained')
        else:
            print('VaReCoF run failed')

    return ts_found


def run_sadpt_scan(typ, grid, dist_name, brk_name, ts_zma, ts_info, ref_level,
                   scn_run_fs, scn_save_fs, opt_script_str,
                   overwrite, update_guess, **opt_kwargs):
    """ saddle point scan code
    """
    if 'elimination' in typ:
        grid1, grid2 = grid
        grid_dct = {dist_name: grid1, brk_name: grid2}
    else:
        grid_dct = {dist_name: grid}
    scan.run_scan(
        zma=ts_zma,
        spc_info=ts_info,
        thy_level=ref_level,
        grid_dct=grid_dct,
        scn_run_fs=scn_run_fs,
        scn_save_fs=scn_save_fs,
        script_str=opt_script_str,
        saddle=False,
        overwrite=overwrite,
        update_guess=update_guess,
        reverse_sweep=False,
        fix_failures=False,
        **opt_kwargs,
        )
    if 'elimination' in typ:
        coo_names = [dist_name, brk_name]
    else:
        coo_names = [dist_name]
    scan.save_scan(
        scn_run_fs=scn_run_fs,
        scn_save_fs=scn_save_fs,
        coo_names=coo_names,
        )

    # Find the structure at the maximum on the grid opt scan
    if 'elimination' in typ:
        max_zma, max_ene = rxngrid.find_max_2d(
            grid1, grid2, dist_name, brk_name, scn_save_fs)
    else:
        max_zma, max_ene = rxngrid.find_max_1d(
            typ, grid, ts_zma, dist_name, scn_save_fs)
    print('geometry for maximum along scan:', max_zma)
    print('energy for maximum along scan:', max_ene)

    return max_zma


def find_sadpt_transition_state(
        opt_script_str,
        run_fs,
        max_zma,
        ts_info,
        ref_level,
        overwrite,
        ts_save_path,
        ts_save_fs,
        dist_name,
        dist_info,
        thy_save_fs,
        cnf_run_fs,
        cnf_save_fs,
        **opt_kwargs):
    """ Optimize the transition state structure obtained from the grid search
    """

    # Run the transition state optimization
    driver.run_job(
        job='optimization',
        script_str=opt_script_str,
        run_fs=run_fs,
        geom=max_zma,
        spc_info=ts_info,
        thy_level=ref_level,
        saddle=True,
        overwrite=overwrite,
        **opt_kwargs,
        )

    # Read the contents of the optimization
    opt_ret = driver.read_job(
        job='optimization',
        run_fs=run_fs,
    )

    if opt_ret is not None:

        # Set ts found
        ts_found = True

        # If successful, Read the geom and energy from the optimization
        inf_obj, _, out_str = opt_ret
        prog = inf_obj.prog
        method = inf_obj.method
        ene = elstruct.reader.energy(prog, method, out_str)
        geo = elstruct.reader.opt_geometry(prog, out_str)
        zma = elstruct.reader.opt_zmatrix(prog, out_str)

        # Save the information into the filesystem
        print(" - Saving...")
        print(" - Save path: {}".format(ts_save_path))

        ts_save_fs[0].file.energy.write(ene)
        ts_save_fs[0].file.geometry.write(geo)
        ts_save_fs[0].file.zmatrix.write(zma)

        # Run single conformer to get intitial conformer in filesystem
        vals = automol.zmatrix.values(zma)
        final_dist = vals[dist_name]
        dist_info[1] = final_dist
        conformer.conformer_sampling(
            spc_info=ts_info,
            thy_info=ref_level,
            thy_save_fs=ts_save_fs,
            cnf_run_fs=cnf_run_fs,
            cnf_save_fs=cnf_save_fs,
            script_str=opt_script_str,
            overwrite=overwrite,
            nsamp_par=[False, 0, 0, 0, 0, 1],
            saddle=True,
            dist_info=dist_info,
            two_stage=True,
            **opt_kwargs
        )
    else:
        ts_found = False

    return ts_found


# HELPER FUNCTIONS FOR THE MAIN FINDER FUNCTIONS
def check_ts_zma(zma, ts_zma):
    """ Check to see if zma in filesystem matches guess ts zma
        check to see if rxn class for already found ts is of expected class
        do this by comparing names
    """
    chk_bkp = False
    if automol.zmatrix.names(zma) == automol.zmatrix.names(ts_zma):
        if not automol.zmatrix.almost_equal(zma, ts_zma, 4e-1, True):
            if 'babs1' in automol.zmatrix.names(ts_zma):
                babs1 = 170. * phycon.DEG2RAD
                if automol.zmatrix.values(ts_zma)['babs1'] == babs1:
                    babs1 = 85. * phycon.DEG2RAD
                ts_zma = automol.zmatrix.set_valuess(
                    ts_zma, {'babs1': babs1})
                if not automol.zmatrix.almost_equal(zma, ts_zma, 4e-1):
                    chk_bkp = True
            else:
                chk_bkp = True
    else:
        chk_bkp = True

    return chk_bkp
# def check_filesys_for_ts(ts_dct, ts_zma, cnf_save_fs, overwrite,
#                          typ, dist_info, dist_name, bkp_ts_class_data):
#     """ Check if TS is in filesystem and matches original guess
#     """
#     update_dct = {}
#
#     # Check if TS is in filesystem and check if there is a match
#     min_cnf_locs = fsmin.min_energy_conformer_locators(cnf_save_fs)
#     if min_cnf_locs and not overwrite:
#
#         print('Found TS at {}'.format(cnf_save_fs[0].path()))
#
#         # Check if TS matches original guess
#         zma = cnf_save_fs[-1].file.zmatrix.read(min_cnf_locs)
#         chk_bkp = check_ts_zma(zma, ts_zma)
#
#         # Check if TS matches original guess from back reaction
#         if chk_bkp and bkp_ts_class_data:
#            [bkp_typ, bkp_ts_zma, _, _, bkp_tors_names, _] = bkp_ts_class_data
#             is_bkp = check_ts_zma(zma, bkp_ts_zma)
#
#         # Set information in ts_dct as needed
#         update_dct['class'] = ts_dct['class'] if not is_bkp else bkp_typ
#         update_dct['zma'] = ts_dct['zma'] if not is_bkp else bkp_ts_zma
#        update_dct['tors_names'] = ts_dct['zma'] if not is_bkp else bkp_ts_zma
#         # ts_dct['original_zma'] = ts_zma
#         if is_bkp:
#             print('updating reaction class to {}'.format(bkp_typ))
#             update_dct['class'] = ts_dct['class'] if not is_bkp else bkp_typ
#             update_dct['zma'] = ts_dct['zma'] if not is_bkp else bkp_ts_zma
#             # ts_dct['class'] = bkp_typ
#             # ts_dct['original_zma'] = bkp_ts_zma
#             # ts_dct['tors_names'] = bkp_tors_names
#             # if not is_typ or not is
#         else:
#             print("TS may not be original type or backup type")
#             print("Some part of the z-matrices have changed")
#
#         print('class test:', ts_dct['class'])
#         vals = automol.zmatrix.values(zma)
#         final_dist = vals[dist_name]
#         dist_info[1] = final_dist
#
#         # Add an angle check which is added to spc dct for TS
#         angle = lts.check_angle(
#             ts_dct['original_zma'],
#             ts_dct['dist_info'],
#             ts_dct['class'])
#         ts_dct['dist_info'][1] = final_dist
#         ts_dct['dist_info'].append(angle)
#
#     return ts_dct
# SOME SECOND ATTEMPT REACTION BASED ON REACTION TYPES
# def aa
#     """
#     """
#     elif 'addition' in typ and bkp_ts_class_data and attempt > 2:
#         # Try to find addition rxn TS in reverse direction
#         bkp_typ, bkp_ts_zma, bkp_dist_name, bkp_grid, bkp_tors_names,
#         bkp_update_guess = bkp_ts_class_data
#         print('TS find failed. Attempting to find with new',
#               'reaction class: {}'.format(bkp_typ))
#         bkp_dist_info = [bkp_dist_name, 0., bkp_update_guess]
#         ts_dct['class'] = bkp_typ
#         ts_dct['original_zma'] = bkp_ts_zma
#         ts_dct['dist_info'] = bkp_dist_info
#         ts_dct['tors_names'] = bkp_tors_names
#         attempt += 1
#         geo, zma, final_dist = find_ts(
#             spc_dct, ts_dct, ts_info, bkp_ts_zma, bkp_typ, bkp_dist_info,
#             bkp_grid, None, ini_thy_info, thy_info, run_prefix,
#             save_prefix, rxn_run_path, rxn_save_path, overwrite=True,
#             attempt=attempt)
#     elif ('addition ' in typ or 'abstraction' in typ) and attempt < 3:
#         # Try to find addition rxn TS in reverse direction with addn tricks
#         babs1 = 170. * phycon.DEG2RAD
#         if automol.zmatrix.values(ts_zma)['babs1'] == babs1:
#             babs1 = 85. * phycon.DEG2RAD
#         print('TS find failed. Attempting to find with '
#               'new angle of attack: {:.1f}'.format(babs1))
#         ts_zma = automol.zmatrix.set_values(ts_zma, {'babs1': babs1})
#         ts_dct['original_zma'] = ts_zma
#         attempt += 1
#         geo, zma, final_dist = find_ts(
#             spc_dct, ts_dct, ts_info, ts_zma, typ, dist_info, grid,
#             bkp_ts_class_data, ini_thy_info, thy_info, run_prefix,
#             save_prefix, rxn_run_path, rxn_save_path, overwrite=True,
#             attempt=attempt)
#     elif 'beta scission' in typ and bkp_ts_class_data and attempt < 2:
#         # Run reverse for a beta scission reaction
#         [bkp_typ, bkp_ts_zma, bkp_dist_name, bkp_grid,
#          bkp_tors_names, bkp_update_guess] = bkp_ts_class_data
#         print('TS find failed. Attempting to find with'
#               'new reaction class: {}'.format(bkp_typ))
#         bkp_dist_info = [bkp_dist_name, 0., bkp_update_guess]
#         ts_dct['class'] = bkp_typ
#         ts_dct['original_zma'] = bkp_ts_zma
#         ts_dct['dist_info'] = bkp_dist_info
#         ts_dct['tors_names'] = bkp_tors_names
#         attempt += 1
#         geo, zma, final_dist = find_ts(
#             spc_dct, ts_dct, ts_info, bkp_ts_zma, bkp_typ, bkp_dist_info,
#             bkp_grid, None, ini_thy_info, thy_info, run_prefix,
#             save_prefix, rxn_run_path, rxn_save_path, overwrite=True,
#             attempt=attempt)
