"""
Find a TS from the grid as well as associated vdW wells
"""

import automol
import autofile
from routines.es._routines import _sadpt as sadpt
from routines.es._routines import _vtst as vtst
from routines.es._routines import _vrctst as vrctst


# Main TS finder functions
def sadpt_transition_state(
        ini_zma, ts_info, mod_thy_info,
        thy_run_path, thy_save_path,
        ini_thy_save_path,
        cnf_run_fs, cnf_save_fs,
        ts_save_fs, ts_save_path, run_fs,
        typ, grid, update_guess,
        dist_name, dist_info, brk_name,
        opt_script_str, overwrite,
        es_keyword_dct, **opt_kwargs):
    """ Find a sadddle point
    """

    # Check filesystem for input level of theory
    print('\nSearching save filesys for guess Z-Matrix calculated',
          'at {} level...'.format(es_keyword_dct['inplvl']))
    guess_zmas = sadpt.check_filesys_for_guess(ini_thy_save_path)

    # If no guess zma, run a TS searching algorithm
    if not guess_zmas:

        print(' - No Z-Matrix in found in save filesys.')
        print('\nRunning scan to generate guess Z-Matrix for opt...')
        scn_run_fs = autofile.fs.scan(thy_run_path)
        scn_save_fs = autofile.fs.scan(thy_save_path)
        guess_zmas = sadpt.scan_for_guess(
            typ, grid, dist_name, brk_name, ini_zma, ts_info,
            mod_thy_info,
            scn_run_fs, scn_save_fs, opt_script_str,
            overwrite, update_guess, **opt_kwargs)

    # Optimize the saddle point
    print('\nOptimizing guess Z-Matrix obtained from scan or filesys...')
    sadpt.optimize_transition_state(
        guess_zmas, ts_info, mod_thy_info,
        cnf_run_fs, cnf_save_fs,
        ts_save_fs, ts_save_path, run_fs,
        dist_name, dist_info,
        opt_script_str, overwrite, **opt_kwargs)


def barrierless_transition_state(
        ts_info, ts_zma, ts_dct, spc_dct,
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
    switch = False

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

    return ts_found, switch
