""" Run and Read the scans from VTST calculations
"""

import automol
import autofile
from lib import filesys
from routines.es._routines import sp
from routines.es._routines import _wfn as wfn
from routines.es._routines import _scan as scan
from lib.submission import qchem_params


def radrad_scan(ts_zma, ts_info, hs_info,
                ts_formula, high_mul, active_space,
                rct_info, rct_ichs, rcts_cnf_fs, rcts_gra,
                grid1, grid2, coord_name, frm_bnd_keys,
                mod_var_scn_thy_info,
                mod_var_sp1_thy_info,  # Need an unmodifie
                var_sp1_thy_info,
                var_sp2_thy_info,
                hs_var_sp1_thy_info,
                hs_var_sp2_thy_info,
                mod_thy_info,
                vscnlvl_thy_save_fs,
                vscnlvl_ts_save_fs,
                scn_run_fs, scn_save_fs,
                overwrite, update_guess,
                constraint_dct=None,
                zma_locs=(0,)):
    """ Run the scan for VTST calculations
    """

    # Set up the casscf options
    ref_zma = automol.zmatrix.set_values(ts_zma, {coord_name: grid1[0]})
    cas_kwargs = wfn.build_wfn(ref_zma, ts_info, ts_formula, high_mul,
                               rct_ichs, rct_info,
                               active_space, mod_var_scn_thy_info)

    # Run the scan along the reaction coordinate
    scan.multiref_rscan(
        ts_zma=ts_zma,
        ts_info=ts_info,
        grid1=grid1,
        grid2=grid2,
        coord_name=coord_name,
        mod_var_scn_thy_info=mod_var_scn_thy_info,
        vscnlvl_thy_save_fs=vscnlvl_thy_save_fs,
        scn_run_fs=scn_run_fs,
        scn_save_fs=scn_save_fs,
        overwrite=overwrite,
        update_guess=update_guess,
        constraint_dct=constraint_dct,
        **cas_kwargs
    )

    # Calculate and store the infinite separation energy
    far_locs = [[coord_name], [grid1[0]]]
    ts_zma = scn_save_fs[-1].file.zmatrix.read(far_locs)
    geo = scn_save_fs[-1].file.geometry.read(far_locs)

    geo_run_path = scn_run_fs[-1].path(far_locs)
    geo_save_path = scn_save_fs[-1].path(far_locs)

    inf_sep_ene = scan.radrad_inf_sep_ene(
        hs_info, ts_zma,
        rct_info, rcts_cnf_fs,
        var_sp1_thy_info, var_sp2_thy_info,
        hs_var_sp1_thy_info, hs_var_sp2_thy_info,
        geo, geo_run_path, geo_save_path,
        scn_save_fs, far_locs,
        overwrite=overwrite,
        **cas_kwargs)

    inf_locs = [[coord_name], [1000.]]
    scn_save_fs[-1].create(inf_locs)
    scn_save_fs[-1].file.energy.write(inf_sep_ene, inf_locs)

    # Save the vmatrix for use in reading
    _save_traj(ts_zma, frm_bnd_keys, rcts_gra, 
               vscnlvl_ts_save_fs, zma_locs=zma_locs)

    print('\nRunning Hessians and energies...')
    scn_locs = filesys.build.scn_locs_from_fs(
        scn_save_fs, [coord_name], constraint_dct=None)

    _vtst_hess_ene(ts_info, mod_var_scn_thy_info, mod_var_sp1_thy_info,
                   scn_save_fs, scn_run_fs, scn_locs, inf_locs,
                   overwrite, **cas_kwargs)


def molrad_scan(ts_zma, ts_info,
                rct_info, rcts_cnf_fs, rcts_gra,
                grid1, grid2, coord_name, frm_bnd_keys,
                thy_info, vsp1_thy_info,
                thy_save_fs,
                ts_save_fs,
                scn_run_fs, scn_save_fs,
                overwrite, update_guess, retryfail,
                zma_locs=(0,)):
    """ Run the scan for VTST calculations
    """

    # Set the thy info objects appropriately
    if vsp1_thy_info is not None:
        inf_thy_info = vsp1_thy_info
    else:
        inf_thy_info = thy_info
    mod_thy_info = filesys.inf.modify_orb_restrict(ts_info, thy_info)
    mod_vsp1_thy_info = filesys.inf.modify_orb_restrict(ts_info, vsp1_thy_info)

    # Set script
    _, opt_script_str, _, opt_kwargs = qchem_params(
        *mod_thy_info[0:2])

    # Setup and run the first part of the scan to shorte
    scan.run_two_way_scan(
        ts_zma, ts_info, mod_thy_info,
        grid1, grid2, coord_name,
        thy_save_fs,
        scn_run_fs, scn_save_fs,
        opt_script_str, overwrite,
        update_guess=update_guess,
        reverse_sweep=False,
        saddle=False,   # opts along scan are min, not sadpt opts
        constraint_dct=None,
        retryfail=retryfail,
        **opt_kwargs
    )

    # Infinite seperation energy calculation
    print('\nCalculating infinite separation energy...')
    print('inf_thy_info', inf_thy_info)
    inf_sep_ene = scan.molrad_inf_sep_ene(
        rct_info, rcts_cnf_fs,
        inf_thy_info, overwrite)

    inf_locs = [[coord_name], [1000.]]
    scn_save_fs[-1].create(inf_locs)
    scn_save_fs[-1].file.energy.write(inf_sep_ene, inf_locs)

    # Save the vmatrix for use in reading
    _save_traj(ts_zma, frm_bnd_keys, rcts_gra,
               ts_save_fs, zma_locs=zma_locs)

    print('\nRunning Hessians and energies...')
    cas_kwargs = {}
    scn_locs = filesys.build.scn_locs_from_fs(
        scn_save_fs, [coord_name], constraint_dct=None)

    _vtst_hess_ene(ts_info, mod_thy_info, mod_vsp1_thy_info,
                   scn_save_fs, scn_run_fs, scn_locs, inf_locs,
                   overwrite, **cas_kwargs)


def _vtst_hess_ene(ts_info, mod_thy_info, mod_vsp1_thy_info,
                   scn_save_fs, scn_run_fs, scn_locs, inf_locs,
                   overwrite, **cas_kwargs):
    """ VTST Hessians and Energies
    """

    print('\n Running Hessians and Gradients...')
    hess_script_str, _, hess_kwargs, _ = qchem_params(
        *mod_thy_info[0:2])
    hess_kwargs.update(cas_kwargs)
    for x, y in hess_kwargs.items():
        print(x)
        print(y)
    for locs in scn_locs:
        if locs != inf_locs:
            geo_run_path = scn_run_fs[-1].path(locs)
            geo_save_path = scn_save_fs[-1].path(locs)
            scn_run_fs[-1].create(locs)
            zma, geo = filesys.inf.cnf_fs_zma_geo(scn_save_fs, locs)
            sp.run_hessian(zma, geo, ts_info, mod_thy_info,
                           scn_save_fs, geo_run_path, geo_save_path, locs,
                           hess_script_str, overwrite, **hess_kwargs)
            sp.run_gradient(zma, geo, ts_info, mod_thy_info,
                            scn_save_fs, geo_run_path, geo_save_path, locs,
                            hess_script_str, overwrite, **hess_kwargs)

    print('\n Running Energies...')
    script_str, _, ene_kwargs, _ = qchem_params(
        *mod_vsp1_thy_info[0:2])
    ene_kwargs.update(cas_kwargs)
    for locs in scn_locs:
        if locs != inf_locs:
            geo_run_path = scn_run_fs[-1].path(locs)
            geo_save_path = scn_save_fs[-1].path(locs)
            scn_run_fs[-1].create(locs)
            zma, geo = filesys.inf.cnf_fs_zma_geo(scn_save_fs, locs)
            sp.run_energy(zma, geo, ts_info, mod_vsp1_thy_info,
                          scn_save_fs, geo_run_path, geo_save_path, locs,
                          script_str, overwrite, **ene_kwargs)


def _save_traj(ts_zma, frm_bnd_keys, rcts_gra, ts_save_fs, zma_locs=(0,)):
    """ save trajectory and zma stuff
    """

    print('\nSaving the V-Matrix into the filesystem...')
    ts_fs, _ = ts_save_fs
    ts_path = ts_fs[-1].path()
    zma_fs = autofile.fs.manager(ts_path, 'ZMATRIX')
    zma_fs[-1].create(zma_locs)
    zma_fs[-1].file.vmatrix.write(automol.zmatrix.var_(ts_zma), zma_locs)

    print('\nSaving the trajectory into the filesystem...')
    tra = (frozenset({frm_bnd_keys}),
           frozenset({}))
    zma_fs[-1].file.transformation.write(tra, zma_locs)
    zma_fs[-1].file.reactant_graph.write(rcts_gra, zma_locs)
