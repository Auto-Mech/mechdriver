""" Run and Read the scans from VTST calculations
"""

from lib import filesys
from routines.es.runner import par as runpar
from routines.es._routines import sp
from routines.es._routines import _varscan as varscan
from routines.es._routines import _scan as scan


def radrad_scan(ts_zma, ts_info, ts_formula, high_mul,
                spc_1_info, spc_2_info,
                grid1, grid2, dist_name,
                num_act_orb, num_act_elc,
                mod_var_scn_thy_info,
                mod_var_sp1_thy_info,
                mod_var_sp2_thy_info,
                hs_var_scn_thy_info,
                hs_var_sp1_thy_info,
                hs_var_sp2_thy_info,
                mod_ini_thy_info,
                scn_run_fs, scn_save_fs,
                run_prefix, save_prefix,
                overwrite, update_guess,
                **opt_kwargs):
    """ Run the scan for VTST calculations
    """

    varscan.run_multiref_rscan(
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
        thy_info=mod_var_scn_thy_info
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

    inf_sep_ene = varscan.infinite_separation_energy(
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

    # Run the Hessians


def molrad_scan(ts_zma, ts_info,
                spc1_info, spc2_info,
                grid1, grid2, dist_name,
                thy_info, ini_thy_info,
                var_sp1_thy_info,
                scn_run_fs, scn_save_fs,
                run_prefix, save_prefix,
                overwrite, update_guess,
                **opt_kwargs):
    """ Run the scan for VTST calculations
    """

    # Modify the theory
    mod_ts_thy_info = filesys.inf.modify_orb_restrict(ts_info, thy_info)
    mod_ts_ini_thy_info = filesys.inf.modify_orb_restrict(ts_info, ini_thy_info)
    if var_sp1_thy_info is not None:
        mod_ts_sp_thy_info = filesys.inf.modify_orb_restrict(
            ts_info, var_sp1_thy_info)
        inf_thy_info = var_sp1_thy_info
    else:
        mod_ts_sp_thy_info = mod_ts_thy_info
        inf_thy_info = thy_info

    # Set script
    _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
        *mod_ts_thy_info[0:2])

    # Setup and run the first part of the scan to shorte
    varscan.run_two_way_scan(
        ts_zma, ts_info, mod_ts_thy_info,
        grid1, grid2, dist_name,
        scn_run_fs, scn_save_fs,
        opt_script_str, overwrite,
        update_guess=update_guess,
        reverse_sweep=False,
        fix_failures=True,
        saddle=False,
        constraint_dct=None,
        **opt_kwargs
    )

    # Infinite seperation energy calculation
    print('\nCalculating infinite separation energy...')
    inf_sep_ene = varscan.molrad_inf_sep_ene(
        spc1_info, spc2_info,
        run_prefix, save_prefix,
        inf_thy_info, ini_thy_info,
        overwrite)

    inf_locs = [[dist_name], [1000.]]
    scn_save_fs[-1].create(inf_locs)
    scn_save_fs[-1].file.energy.write(inf_sep_ene, inf_locs)

    # Run the hessians
    print('\nRunning Hessians and energies...')
    script_str, _, kwargs, _ = runpar.run_qchem_par(*mod_ts_sp_thy_info[0:2])
    scn_locs = filesys.build.scn_locs_from_fs(
        scn_save_fs, [dist_name], constraint_dct=None)
    print('\n Running Hessians...')
    script_str, _, kwargs, _ = runpar.run_qchem_par(*mod_ts_thy_info[0:2])
    for locs in scn_locs:
        if locs != inf_locs:
            geo_run_path = scn_run_fs[-1].path(locs)
            geo_save_path = scn_save_fs[-1].path(locs)
            zma, geo = filesys.inf.get_zma_geo(scn_save_fs, locs)
            scn_run_fs[-1].create(locs)
            zma, geo = filesys.inf.get_zma_geo(scn_save_fs, locs)
            sp.run_hessian(zma, geo, ts_info, mod_ts_thy_info,
                          scn_save_fs, geo_run_path, geo_save_path, locs,
                          script_str, overwrite, **kwargs)
    print('\n Running Energies...')
    script_str, _, kwargs, _ = runpar.run_qchem_par(*mod_ts_sp_thy_info[0:2])
    for locs in scn_locs:
        if locs != inf_locs:
            geo_run_path = scn_run_fs[-1].path(locs)
            geo_save_path = scn_save_fs[-1].path(locs)
            zma, geo = filesys.inf.get_zma_geo(scn_save_fs, locs)
            scn_run_fs[-1].create(locs)
            zma, geo = filesys.inf.get_zma_geo(scn_save_fs, locs)
            sp.run_energy(zma, geo, ts_info, mod_ts_sp_thy_info,
                          scn_save_fs, geo_run_path, geo_save_path, locs,
                          script_str, overwrite, **kwargs)


