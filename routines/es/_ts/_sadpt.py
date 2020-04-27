""" Functions for sadpt
"""

import automol
import elstruct
from routines.es import conformer
from routines.es import _scan as scan
from runners import es as es_runner
from lib import filesys
from lib.reaction import grid as rxngrid


def check_filesys_for_guess(ini_thy_save_path, es_keyword_dct):
    """ Check if the filesystem for any TS structures at the input
        level of theory
    """
    guess_zmas = []

    # Check and see if a zma is found from the filesystem
    ini_cnf_save_fs, ini_cnf_save_locs = filesys.build.cnf_fs_from_prefix(
        ini_thy_save_path, cnf='min')
    if ini_cnf_save_locs:
        if ini_cnf_save_fs[-1].file.zmatrix.exists(ini_cnf_save_locs):
            geo_path = ini_cnf_save_fs[-1].path(ini_cnf_save_locs)
            print(' - Z-Matrix found.')
            print(' - Reading Z-Matrix from path {}'.format(geo_path))
            guess_zmas.append(
                ini_cnf_save_fs[-1].file.zmatrix.read(
                    ini_cnf_save_locs))

    return guess_zmas


def scan_for_guess(typ, grid, dist_name, brk_name,
                   ts_zma, ts_info, ref_level,
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
        max_zma = rxngrid.find_max_2d(
            grid1, grid2, dist_name, brk_name, scn_save_fs)
        guess_zmas = [max_zma]
    else:
        guess_zmas = rxngrid.find_max_1d(
            typ, grid, ts_zma, dist_name, scn_save_fs)

    return guess_zmas


def optimize_transition_state(
        guess_zmas, ts_info, mod_thy_info,
        cnf_run_fs, cnf_save_fs,
        ts_save_fs, ts_save_path, run_fs,
        dist_name, dist_info,
        opt_script_str, overwrite, **opt_kwargs):
    """ Optimize the transition state structure obtained from the grid search
    """

    if len(guess_zmas) == 1:
        print('\nThere is 1 guess Z-Matrix',
              'to attempt to find saddle point.')
    else:
        print('\nThere are {} guess Z-Matrices'.format(len(guess_zmas)),
             'to attempt to find saddle point.')

    # Loop over all the guess zmas to find a TS
    ts_found = False
    for idx, zma in enumerate(guess_zmas):
        print('\nAttempting optimization of',
              'guess Z-Matrix {}...'.format(idx+1))

        # Run the transition state optimization
        es_runner.run_job(
            job='optimization',
            script_str=opt_script_str,
            run_fs=run_fs,
            geom=zma,
            spc_info=ts_info,
            thy_level=mod_thy_info,
            saddle=True,
            overwrite=overwrite,
            **opt_kwargs,
            )

        # Read the contents of the optimization
        opt_ret = es_runner.read_job(
            job='optimization',
            run_fs=run_fs,
        )

        if opt_ret is not None:

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

            # Save this structure as first conformer
            # cnf_save_fs[-1].create(locs)
            # cnf_save_fs[-1].file.geometry_info.write(inf_obj, locs)
            # cnf_save_fs[-1].file.geometry_input.write(inp_str, locs)
            # cnf_save_fs[-1].file.energy.write(ene, locs)
            # cnf_save_fs[-1].file.geometry.write(geo, locs)
            # cnf_save_fs[-1].file.zmatrix.write(zma, locs)

            # Run single conformer to get intitial conformer in filesystem
            vals = automol.zmatrix.values(zma)
            final_dist = vals[dist_name]
            dist_info[1] = final_dist
            conformer.single_conformer(
                zma=zma,
                spc_info=ts_info,
                thy_info=mod_thy_info,
                thy_save_fs=ts_save_fs,
                cnf_run_fs=cnf_run_fs,
                cnf_save_fs=cnf_save_fs,
                overwrite=overwrite,
                saddle=True,
                dist_info=dist_info
            )

            # Exit the for loop if a TS has been found
            ts_found = True
            break

    return ts_found


# HELPER FUNCTIONS FOR THE MAIN FINDER FUNCTIONS
# def check_ts_zma(zma, ts_zma):
#     """ Check to see if zma in filesystem matches guess ts zma
#         check to see if rxn class for already found ts is of expected class
#         do this by comparing names
#     """
#     chk_bkp = False
#     if automol.zmatrix.names(zma) == automol.zmatrix.names(ts_zma):
#         if not automol.zmatrix.almost_equal(zma, ts_zma, 4e-1, True):
#             if 'babs1' in automol.zmatrix.names(ts_zma):
#                 babs1 = 170. * phycon.DEG2RAD
#                 if automol.zmatrix.values(ts_zma)['babs1'] == babs1:
#                     babs1 = 85. * phycon.DEG2RAD
#                 ts_zma = automol.zmatrix.set_valuess(
#                     ts_zma, {'babs1': babs1})
#                 if not automol.zmatrix.almost_equal(zma, ts_zma, 4e-1):
#                     chk_bkp = True
#             else:
#                 chk_bkp = True
#     else:
#         chk_bkp = True
#
#     return chk_bkp
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
