""" Run and Read the scans from VTST calculations
"""

from lib import filesys
from routines.es import runner as es_runner
from routines.es._routines import sp
from routines.es._routines import _varscan as varscan
from routines.es._routines import _scan as scan


def radrad_scan(ts_zma, ts_info, ts_formula, high_mul,
                grid1, grid2, coord_name,
                mod_var_scn_thy_info,
                mod_var_sp1_thy_info,
                mod_var_sp2_thy_info,
                hs_var_scn_thy_info,
                hs_var_sp1_thy_info,
                hs_var_sp2_thy_info,
                mod_ini_thy_info,
                vscnlvl_thy_save_fs,
                scn_run_fs, scn_save_fs,
                rcts_cnf_fs,
                run_prefix, save_prefix,
                overwrite, update_guess):
    """ Run the scan for VTST calculations
    """

    varscan.multiref_rscan(
        ts_zma=ts_zma,
        ts_info=ts_info,
        ts_formula=ts_formula,
        high_mul=high_mul,
        grid1=grid1,
        grid2=grid2,
        coord_name=coord_name,
        mod_var_scn_thy_info=mod_var_scn_thy_info,
        vscnlvl_thy_save_fs=vscnlvl_thy_save_fs,
        scn_run_fs=scn_run_fs,
        scn_save_fs=scn_save_fs,
        overwrite=overwrite,
        update_guess=update_guess,
    )

    # Calculate and store the infinite separation energy
    locs = [[coord_name], [grid1[0]]]
    print('ts zma locs', locs)
    ts_zma = scn_save_fs[-1].file.zmatrix.read(locs)

    # set up all the file systems for the TS
    # start with the geo and reference theory info
    geo_run_path = scn_run_fs[-1].path(locs)
    geo_save_path = scn_save_fs[-1].path(locs)
    geo = scn_save_fs[-1].file.geometry.read(locs)

    inf_sep_ene = varscan.radrad_inf_sep_ene(
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

    inf_locs = [[coord_name], [1000.]]
    scn_save_fs[-1].create(inf_locs)
    scn_save_fs[-1].file.energy.write(inf_sep_ene, inf_locs)

    # Run the Hessians
    # _vtst_hess_ene(ts_info, mod_thy_info, mod_vsp1_thy_info,
    #                scn_save_fs, scn_run_fs, scn_locs, inf_locs,
    #                overwrite)


def molrad_scan(ts_zma, ts_info,
                grid1, grid2, coord_name,
                mod_thy_info, mod_ini_thy_info,
                mod_vsp1_thy_info,
                thy_save_fs,
                scn_run_fs, scn_save_fs,
                rcts_cnf_fs,
                run_prefix, save_prefix,
                overwrite, update_guess, retryfail):
    """ Run the scan for VTST calculations
    """

    if mod_vsp1_thy_info is not None:
        inf_thy_info = mod_vsp1_thy_info
    else:
        inf_thy_info = mod_thy_info

    # Set script
    _, opt_script_str, _, opt_kwargs = es_runner.qchem_params(
        *mod_thy_info[0:2])

    # Setup and run the first part of the scan to shorte
    varscan.run_two_way_scan(
        ts_zma, ts_info, mod_thy_info,
        grid1, grid2, coord_name,
        thy_save_fs,
        scn_run_fs, scn_save_fs,
        opt_script_str, overwrite,
        update_guess=update_guess,
        reverse_sweep=False,
        saddle=False,
        constraint_dct=None,
        retryfail=retryfail,
        **opt_kwargs
    )

    # Infinite seperation energy calculation
    print('\nCalculating infinite separation energy...')
    inf_sep_ene = varscan.molrad_inf_sep_ene(
        rcts_cnf_fs,
        run_prefix, save_prefix,
        inf_thy_info, mod_ini_thy_info,
        overwrite)

    inf_locs = [[coord_name], [1000.]]
    scn_save_fs[-1].create(inf_locs)
    scn_save_fs[-1].file.energy.write(inf_sep_ene, inf_locs)

    print('\nRunning Hessians and energies...')
    scn_locs = filesys.build.scn_locs_from_fs(
        scn_save_fs, [coord_name], constraint_dct=None)

    _vtst_hess_ene(ts_info, mod_thy_info, mod_vsp1_thy_info,
                   scn_save_fs, scn_run_fs, scn_locs, inf_locs,
                   overwrite)


def _vtst_hess_ene(ts_info, mod_thy_info, mod_vsp1_thy_info,
                   scn_save_fs, scn_run_fs, scn_locs, inf_locs,
                   overwrite):
    """ VTST Hessians and Energies
    """

    print('\n Running Hessians...')
    script_str, _, kwargs, _ = es_runner.qchem_params(
        *mod_thy_info[0:2])
    for locs in scn_locs:
        if locs != inf_locs:
            geo_run_path = scn_run_fs[-1].path(locs)
            geo_save_path = scn_save_fs[-1].path(locs)
            scn_run_fs[-1].create(locs)
            zma, geo = filesys.inf.cnf_fs_zma_geo(scn_save_fs, locs)
            sp.run_hessian(zma, geo, ts_info, mod_thy_info,
                           scn_save_fs, geo_run_path, geo_save_path, locs,
                           script_str, overwrite, **kwargs)

    print('\n Running Energies...')
    script_str, _, kwargs, _ = es_runner.qchem_params(
        *mod_vsp1_thy_info[0:2])
    for locs in scn_locs:
        if locs != inf_locs:
            geo_run_path = scn_run_fs[-1].path(locs)
            geo_save_path = scn_save_fs[-1].path(locs)
            scn_run_fs[-1].create(locs)
            zma, geo = filesys.inf.cnf_fs_zma_geo(scn_save_fs, locs)
            sp.run_energy(zma, geo, ts_info, mod_vsp1_thy_info,
                          scn_save_fs, geo_run_path, geo_save_path, locs,
                          script_str, overwrite, **kwargs)
