"""
Find a TS from the grid as well as associated vdW wells
"""

import automol
import elstruct
from lib.reaction import grid as rxngrid
from lib.phydat import phycon
from lib.runner import driver
from routines.es import conformer
from routines.es import scan
from routines.es.ts import _sadpt as sadpt
from routines.es.ts import _vtst as vtst
from routines.es.ts import _vrctst as vrctst


# Main TS finder functions
def sadpt_transition_state(
    opt_script_str,
    run_fs,
    guess_zmas,
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
    **opt_kwargs):
    """ aaa """

    # Initialize empty list of guess zmas
    guess_zmas = []

    # Check and see if a zma is found from the filesystem
    ini_cnf_save_fs, ini_cnf_save_locs = fbuild.cnf_fs_from_prefix(
        ini_thy_save_path, cnf='min')
    if ini_cnf_save_locs:
        if ini_cnf_save_fs[-1].file.zmatrix.exists(ini_cnf_save_locs):
            print('Z-Matrix calculated at {} found'.format(
                es_keyword_dct['inplvl']))
            geo_path = ini_cnf_save_fs[-1].path(ini_cnf_save_locs)
            print('Reading Z-Matrix from path {}'.format(geo_path))
            guess_zma = ini_cnf_save_fs[-1].file.zmatrix.read(
                ini_cnf_save_locs)
            guess_zmas.append(guess_zma)
    
    # If no guess zma, run a TS searching algorithm
    if not guess_zmas:
        print('No Z-Matrix in filesys for {} level'.format(
            es_keyword_dct['inplvl']))
        scn_run_fs = autofile.fs.scan(thy_run_path)
        scn_save_fs = autofile.fs.scan(thy_save_path)
        print('Running scan to generate guess Z-Matrix for opt...')
        guess_zma = ts.find.run_sadpt_scan(
            typ, grid, dist_name, brk_name, ts_dct['zma'], ts_info,
            mod_thy_info,
            scn_run_fs, scn_save_fs, opt_script_str,
            overwrite, update_guess, **opt_kwargs)
        guess_zmas.append(guess_zma)

    # Optimize the saddle point
    print('Optimiziing Guess Z-Matrix from scan or filesys...')
    ts.find.find_sadpt_transition_state(
        opt_script_str,
        run_fs,
        guess_zmas,
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


def barrierless_transition_state(ts_info, ts_zma, ts_dct, spc_dct,
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
        ts_found = vtst.run_vtst_scan(
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
        ts_found = vrctst.calc_vrctst_flux(
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
