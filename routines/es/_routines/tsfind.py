"""
Find a TS from the grid as well as associated vdW wells
"""

import automol
import autofile
from routines.es._routines import _sadpt as sadpt
from routines.es._routines import _vtst as vtst
# from routines.es._routines import _vrctst as vrctst
from routines.es._routines import _wfn as wfn


# Main TS finder functions
def sadpt_transition_state(
        ini_zma, ts_info, mod_thy_info,
        thy_run_path, thy_save_path,
        ini_ts_save_path,
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
    guess_zmas = sadpt.check_filesys_for_guess(ini_ts_save_path)

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
        mod_ini_thy_info, mod_thy_info,
        mod_var_scn_thy_info,
        mod_var_sp1_thy_info, mod_var_sp2_thy_info,
        hs_var_scn_thy_info,
        hs_var_sp1_thy_info,
        hs_var_sp2_thy_info,
        run_prefix, save_prefix,
        scn_run_fs, scn_save_fs,
        overwrite, vrc_dct,
        update_guess, **opt_kwargs):
    """ Run TS finder for barrierless reactions
    """
    switch = False

    # Set information from the transition state
    ts_formula = automol.geom.formula(automol.zmatrix.geometry(ts_zma))
    [grid1, grid2] = grid

    # Get info from the reactants
    # rct_zmas = ts_dct['rct_zmas']
    rcts = ts_dct['reacs']
    high_mul = ts_dct['high_mul']
    spc_1_info = [spc_dct[rcts[0]]['ich'],
                  spc_dct[rcts[0]]['chg'],
                  spc_dct[rcts[0]]['mul']]
    spc_2_info = [spc_dct[rcts[1]]['ich'],
                  spc_dct[rcts[1]]['chg'],
                  spc_dct[rcts[1]]['mul']]

    # Set the active space
    num_act_orb, num_act_elc = wfn.active_space(
        ts_dct, spc_dct, ts_dct['high_mul'])

    # Run PST, VTST, VRC-TST based on RAD_RAD_TS model
    if rad_rad_ts.lower() == 'pst':
        ts_found = True
        print('Phase Space Theory Used, No ES calculations are needed')
    elif rad_rad_ts.lower() == 'vtst':
        ts_found = True
        print('Beginning Calculations for VTST Treatments')
        vtst.run_scan(
            ts_zma, ts_info, ts_formula, high_mul,
            spc_1_info, spc_2_info,
            grid1, grid2, dist_name,
            num_act_orb, num_act_elc,
            mod_var_scn_thy_info,
            mod_var_sp1_thy_info, mod_var_sp2_thy_info,
            hs_var_scn_thy_info,
            hs_var_sp1_thy_info,
            hs_var_sp2_thy_info,
            mod_ini_thy_info, mod_thy_info,
            scn_run_fs, scn_save_fs,
            run_prefix, save_prefix,
            overwrite, update_guess,
            **opt_kwargs)
        # if ts_found:
        #     print('Scans for VTST succeeded')
        # else:
        #     print('Scans for VTST failed')
    # elif rad_rad_ts.lower() == 'vrctst':
    #     print('Beginning Calculations for VRC-TST Treatments')
    #     ts_found = vrctst.calc_vrctst_flux(
    #         ts_zma, ts_formula, ts_info, ts_dct, spc_dct,
    #         ts_dct['high_mul'], grid1, grid2, dist_name,
    #         multi_level, mod_multi_opt_info, mod_multi_sp_info,
    #         mod_ini_thy_info, mod_thy_info,
    #         thy_run_path, thy_save_path,
    #         overwrite, update_guess,
    #         run_prefix, save_prefix,
    #         vrc_dct,
    #         corr_pot=True)
    #     if ts_found:
    #         print('VaReCoF run successful and flux file was obtained')
    #     else:
    #         print('VaReCoF run failed')

    return ts_found, switch
