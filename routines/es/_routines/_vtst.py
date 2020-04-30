""" Run and Read the scans from VTST calculations
"""

from routines.es._routines import _mrscan as mrscan
from routines.es._routines import _scan as scan


def run_scan(ts_zma, ts_info, ts_formula, high_mul,
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
             **opt_kwargs):
    """ Run the scan for VTST calculations
    """

    mrscan.run_multiref_rscan(
        ts_zma=ts_zma,
        ts_info=ts_info,
        ts_formula=ts_formula,
        high_mul=high_mul,
        grid1=grid1,
        grid2=grid2,
        dist_name=dist_name,
        num_act_orb=num_act_orb,
        num_act_elc=num_act_elc,
        multi_level=mod_var_scn_thy_info,
        scn_run_fs=scn_run_fs,
        scn_save_fs=scn_save_fs,
        overwrite=overwrite,
        update_guess=update_guess,
        **opt_kwargs
    )

    scan.save_scan(
        scn_run_fs=scn_run_fs,
        scn_save_fs=scn_save_fs,
        coo_names=[dist_name],
    )

    # Calculate and store the infinite separation energy
    locs = [[dist_name], [grid1[0]]]
    print('ts zma locs', locs)
    ts_zma = scn_save_fs[-1].file.zmatrix.read(locs)

    # set up all the file systems for the TS
    # start with the geo and reference theory info
    geo_run_path = scn_run_fs[-1].path(locs)
    geo_save_path = scn_save_fs[-1].path(locs)
    geo = scn_save_fs[-1].file.geometry.read(locs)

    inf_sep_ene = mrscan.infinite_separation_energy(
        spc_1_info, spc_2_info, ts_info, high_mul, ts_zma,
        mod_var_scn_thy_info,
        mod_var_sp1_thy_info, mod_var_sp2_thy_info,
        hs_var_scn_thy_info,
        hs_var_sp1_thy_info,
        hs_var_sp2_thy_info,
        mod_ini_thy_info,
        geo, geo_run_path, geo_save_path,
        run_prefix, save_prefix,
        num_act_orb, num_act_elc)
    inf_locs = [[dist_name], [1000.]]
    scn_save_fs[-1].create(inf_locs)
    scn_save_fs[-1].file.energy.write(inf_sep_ene, inf_locs)
    # geo = automol.zmatrix.geometry(ts_zma)
    # zma = ts_zma
    # final_dist = grid1[0]
