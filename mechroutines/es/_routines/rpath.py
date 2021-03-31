"""
IRC calcs
"""

import elstruct
import autofile
from mechroutines.es import runner as es_runner
from mechroutines.es.runner import qchem_params


# Sadpt scan
def reac_coord_scan(ts_zma, ts_info, zrxn, method_dct,
                    scn_run_fs, scn_save_fs,
                    es_keyword_dct):
    """ saddle point scan code
    """

    # Build grid and names appropriate for reaction type
    scan_inf = automol.reac.build_scan_info(zrxn, ts_zma)
    coord_names, constraint_dct, coord_grids, update_guess = scan_inf

    # Set up script string and kwargs
    mod_thy_info = tinfo.modify_orb_label(tinfo.from_dct(method_dct), ts_info)
    script_str, kwargs = qchem_params(
        method_dct, job=elstruct.Job.OPTIMIZATION)

    es_runner.scan.execute_scan(
        zma=ts_zma,
        spc_info=ts_info,
        mod_thy_info=mod_thy_info,
        thy_save_fs=(),
        coord_names=coord_names,
        coord_grids=coord_grids,
        scn_run_fs=scn_run_fs,
        scn_save_fs=scn_save_fs,
        scn_typ='relaxed',
        script_str=script_str,
        overwrite=es_keyword_dct['overwrite'],
        update_guess=update_guess,
        reverse_sweep=False,
        saddle=False,
        constraint_dct=constraint_dct,
        retryfail=False,
        chkstab=False,
        **kwargs,
        )

    # Find the structure at the maximum on the grid opt scan
    # if 'elimination' in rxn_typ:
    #     [grid1, grid2] = coord_grids
    #     max_zma = rxngrid.find_max_2d(
    #         grid1, grid2, frm_name, brk_name, scn_save_fs,
    #         mod_thy_info, constraint_dct)
    #     guess_zmas = [max_zma]
    # else:
    guess_zmas = rxngrid.find_max_1d(
        zrxn.class_, coord_grids[0], ts_zma, coord_names[0], scn_save_fs,
        mod_thy_info, constraint_dct)

    return guess_zmas



# Intrinsic Reaction Coordinates
def irc_scan(zma, ts_info, coord_name,
             mod_ini_thy_info, ini_method_dct,
             ini_scn_save_fs, geo_run_path,
             overwrite):
    """ Run the IRC
    """

    # Set up run filesys
    run_fs = autofile.fs.run(geo_run_path)

    # Run and Read the IRC in the forward and reverse direction
    for irc_job in (elstruct.Job.IRCF, elstruct.Job.IRCR):

        # Set up the script
        script_str, kwargs = qchem_params(
            ini_method_dct, job=irc_job)

        run_irc(
            zma,
            irc_job,
            coord_name,
            run_fs,
            ini_scn_save_fs,
            ts_info,
            mod_ini_thy_info,
            overwrite,
            script_str,
            **kwargs
        )
        save_irc(
            irc_job,
            coord_name,
            run_fs,
            ini_scn_save_fs,
            mod_ini_thy_info
        )


def run_irc(zma, irc_job, coord_name, run_fs, ini_scn_save_fs,
            ts_info, mod_ini_thy_info, overwrite,
            opt_script_str, **opt_kwargs):
    """ Run the irc job
    """

    # Maybe check for positive coords
    if not _irc_ran(ini_scn_save_fs, coord_name, irc_job):
        print('No IRC calculation in save filesystem')
        opt_success, _ = es_runner.read_job(
            job=irc_job,
            run_fs=run_fs,
        )
        need_irc = bool(not opt_success)
    else:
        print('Found IRC directory at {}'.format(
            ini_scn_save_fs[1].path([coord_name])))
        need_irc = False

    if need_irc:
        print('Running IRC calculation')
        es_runner.run_job(
            job=irc_job,
            script_str=opt_script_str,
            run_fs=run_fs,
            geo=zma,
            spc_info=ts_info,
            thy_info=mod_ini_thy_info,
            overwrite=overwrite,
            **opt_kwargs,
        )
    else:
        print('Skipping IRC calculation')


def save_irc(irc_job, coord_name,
             run_fs, ini_scn_save_fs, mod_ini_thy_info):
    """ Read IRC output and store data in filesystem
    """

    opt_success, opt_ret = es_runner.read_job(
        job=irc_job,
        run_fs=run_fs,
    )
    if opt_success is not None:

        # Read the IRC output file
        inf_obj, inp_str, out_str = opt_ret
        prog = inf_obj.prog
        geos, gras, hessians = elstruct.reader.irc_points(prog, out_str)
        coord_vals, enes = elstruct.reader.irc_path(prog, out_str)

        # Write the data for each geom along IRC to the filesystem
        save_path = ini_scn_save_fs[1].path([coord_name])
        print(" - Saving...")
        print(" - Save path: {}".format(save_path))
        locs_lst = []
        for idx, val in enumerate(coord_vals):

            # Set locs idx; for reverse, ignore SadPt and flip idx to negative
            locs_idx = idx
            if irc_job == elstruct.Job.IRCR:
                if locs_idx == 0:
                    continue
                # val *= -1  coord vals negative from elstruct fxn for g09

            # Scale the coordinates so rounding to .2f number is non-zero
            locs = [coord_name, [val*100.0]]
            locs_lst.append(locs)

            # Save files
            ini_scn_save_fs[-1].create(locs)
            ini_scn_save_fs[-1].file.energy.write(enes[idx], locs)
            ini_scn_save_fs[-1].file.geometry.write(geos[idx], locs)
            ini_scn_save_fs[-1].file.geometry_input.write(inp_str, locs)
            ini_scn_save_fs[-1].file.geometry_info.write(inf_obj, locs)
            if gras:
                ini_scn_save_fs[-1].file.gradient.write(gras[idx], locs)
                ini_scn_save_fs[-1].file.gradient_info.write(inf_obj, locs)
            if hessians:
                ini_scn_save_fs[-1].file.hessian.write(hessians[idx], locs)
                ini_scn_save_fs[-1].file.hessian_info.write(inf_obj, locs)

            scn_save_path = ini_scn_save_fs[-1].path(locs)
            sp_save_fs = autofile.fs.single_point(scn_save_path)
            sp_save_fs[-1].create(mod_ini_thy_info[1:4])
            sp_save_fs[-1].file.input.write(inp_str, mod_ini_thy_info[1:4])
            sp_save_fs[-1].file.info.write(inf_obj, mod_ini_thy_info[1:4])
            sp_save_fs[-1].file.energy.write(enes[idx], mod_ini_thy_info[1:4])
        
        # Build the trajectory file
        if locs_lst:
            es_runner.write_traj(
                coord_name, ini_scn_save_fs, mod_ini_thy_info, locs_lst
            )


def _irc_ran(ini_scn_save_fs, coord_name, irc_job):
    """ See if coords are available
    """

    coords = ini_scn_save_fs[-1].existing([coord_name])
    if irc_job == elstruct.Job.IRCF:
        ran_coords = [coord[1][0] for coord in coords if coord[1][0] > 0.0]
    else:
        ran_coords = [coord[1][0] for coord in coords if coord[1][0] < 0.0]

    return bool(ran_coords)


# Infinite Separation Energy 
def inf_sep_ene(ts_dct, thy_inf_dct, savefs_dct, runfs_dct, es_keyword_dct):
    """ complete scan calcs
    """

    # Get info from the reactants
    ts_info = thy_inf_dct['ts_info']
    ini_zma = ts_dct['zma']
    rct_info = thy_inf_dct['rct_info']
    overwrite = es_keyword_dct['overwrite']
    update_guess = False  # check

    # Grid
    [grid1, grid2] = grid

    # Get thy_inf_dct stuff
    mod_thy_info = thy_inf_dct['mod_runlvl']
    mod_var_scn_thy_info = thy_inf_dct['mod_var_scnlvl']
    mod_var_sp1_thy_info = thy_inf_dct['mod_var_splvl1']
    var_sp1_thy_info = thy_inf_dct['var_splvl2']
    var_sp2_thy_info = thy_inf_dct['var_splvl2']
    hs_var_sp1_thy_info = thy_inf_dct['hs_var_splvl1']
    hs_var_sp2_thy_info = thy_inf_dct['hs_var_splvl2']

    # Get the filesys stuff
    var_scn_save_fs = savefs_dct['vscnlvl_scn_fs']
    var_scn_run_fs = runfs_dct['vscnlvl_scn_fs']
    rcts_cnf_fs = savefs_dct['rcts_cnf_fs']
    vscnlvl_thy_save_fs = savefs_dct['vscnlvl_thy_fs']
    vscnlvl_ts_save_fs = savefs_dct['vscnlvl_ts_fs']

    radrad = False
    if radrad:
        high_mul = ts_dct['high_mult']
        hs_info = info_dct['hs_info']
        rct_ichs = [spc_dct[rct]['inchi'] for rct in ts_dct['reacs']]
        ts_formula = automol.geom.formula(automol.zmatrix.geometry(ini_zma))
        active_space = ts_dct['active_space']
        pot_thresh = es_keyword_dct['pot_thresh']

        scan.radrad_inf_sep_ene(
            hs_info, ts_zma,
            rct_info, rcts_cnf_fs,
            var_sp1_thy_info, var_sp2_thy_info,
            hs_var_sp1_thy_info, hs_var_sp2_thy_info,
            geo, geo_run_path, geo_save_path,
            scn_save_fs, far_locs,
            overwrite=overwrite,
            **cas_kwargs)
    else:
        scan.molrad_inf_sep_ene(
            rct_info, rcts_cnf_fs,
            inf_thy_info, overwrite)

    # Probably just move into the tasks from splitting the tasks initially
    # _vtst_hess_ene(ts_info, coord_name,
    #                mod_var_scn_thy_info, mod_var_sp1_thy_info,
    #                scn_save_fs, scn_run_fs,
    #                overwrite, **cas_kwargs)


def radrad_inf_sep_ene(hs_info, ref_zma,
                       rct_info, rcts_cnf_fs,
                       var_sp1_thy_info, var_sp2_thy_info,
                       hs_var_sp1_thy_info, hs_var_sp2_thy_info,
                       geo, geo_run_path, geo_save_path,
                       scn_save_fs, inf_locs,
                       overwrite=False,
                       **cas_kwargs):
    """ Obtain the infinite separation energy from the multireference energy
        at a given reference point, the high-spin low-spin splitting at that
        reference point, and the high level energy for the high spin state
        at the reference geometry and for the fragments
        scn = thy for optimizations
        sp1 = low-spin single points
        sp2 = high-spin single points for inf sep

        inf = spc0 + spc1 - hs_sr_e + hs_mr_ene
    """

    # Initialize infinite sep energy
    inf_sep_ene = -1.0e12

    # Calculate the energies for the two cases
    hs_inf = (hs_var_sp2_thy_info, hs_var_sp1_thy_info)
    for idx, thy_info in enumerate(hs_inf):

        if idx == 0:
            ioprinter.info_message(
                " - Running high-spin single reference energy ...")
        else:
            ioprinter.info_message(
                " - Running high-spin multi reference energy ...")

        # Calculate the single point energy
        script_str, _, kwargs, _ = qchem_params(*thy_info[0:2])
        cas_kwargs.update(kwargs)

        sp.run_energy(ref_zma, geo, hs_info, thy_info,
                      scn_save_fs, geo_run_path, geo_save_path, inf_locs,
                      script_str, overwrite, highspin=True, **cas_kwargs)

        # Read the energty from the filesystem
        hs_save_fs = autofile.fs.high_spin(geo_save_path)
        if not hs_save_fs[-1].file.energy.exists(thy_info[1:4]):
            ioprinter.error_message(
                'High-spin energy job failed: ',
                'energy is needed to evaluate infinite separation energy')
            ene = None
        else:
            ioprinter.info_message(
                " - Reading high-spin energy from filesystem...")
            ene = hs_save_fs[-1].file.energy.read(thy_info[1:4])

        if idx == 0:
            hs_sr_ene = ene
        else:
            hs_mr_ene = ene

    # Get the single reference energy for each of the reactant configurations
    for thy_inf in (var_sp1_thy_info, var_sp2_thy_info):
        sp_script_str, _, kwargs, _ = qchem_params(
            *var_sp2_thy_info[0:2])
        reac_ene = reac_sep_ene(
            rct_info, rcts_cnf_fs,
            thy_inf, overwrite, sp_script_str, **kwargs)

    # Calculate the infinite seperation energy
    all_enes = (reac_ene, hs_sr_ene, hs_mr_ene)
    if all(ene is not None for ene in all_enes):
        inf_sep_ene = reac_ene - hs_sr_ene + hs_mr_ene
    else:
        inf_sep_ene = None

    return inf_sep_ene


def molrad_inf_sep_ene(rct_info, rcts_cnf_fs,
                       mod_thy_info, overwrite):
    """ Calculate the inf spe ene for a mol-rad ene
    """
    script_str, kwargs = qchem_params(*mod_thy_info[0:2])
    return reac_sep_ene(
        rct_info, rcts_cnf_fs,
        mod_thy_info, overwrite, script_str, **kwargs)


def reac_sep_ene(rct_info, rcts_cnf_fs, thy_info,
                 overwrite, sp_script_str, **kwargs):
    """ Calculate the sum of two reactants
    """

    # get the single reference energy for each of the reactant configurations
    spc_enes = []
    for (run_fs, save_fs, mlocs, mpath), inf in zip(rcts_cnf_fs, rct_info):

        # Set the modified thy info
        mod_thy_info = filesys.inf.modify_orb_restrict(inf, thy_info)

        # Build filesys
        ioprinter.debug_message('locs test', mlocs)
        zma_fs = autofile.fs.zmatrix(mpath)

        # Read the geometry and set paths
        zma = zma_fs[-1].file.zmatrix.read([0])
        geo = save_fs[-1].file.geometry.read(mlocs)
        geo_run_path = run_fs[-1].path(mlocs)
        geo_save_path = save_fs[-1].path(mlocs)

        # Build the single point filesys objects
        sp_save_fs = autofile.fs.single_point(mpath)

        # Calculate the save single point energy
        sp.run_energy(zma, geo, inf, mod_thy_info,
                      save_fs, geo_run_path, geo_save_path,
                      mlocs, sp_script_str, overwrite,
                      **kwargs)
        exists = sp_save_fs[-1].file.energy.exists(mod_thy_info[1:4])
        if not exists:
            ioprinter.warning_message('No ene found')
            ene = None
        else:
            ioprinter.info_message('Reading energy')
            ene = sp_save_fs[-1].file.energy.read(mod_thy_info[1:4])

        # Append ene to list
        spc_enes.append(ene)

    # Analyze the energies in the list
    inf_ene = 0.0
    for ene, inf in zip(spc_enes, rct_info):
        if ene is not None:
            inf_ene += ene
        else:
            ioprinter.error_message(
                'Single reference energy job fails',
                'for {}: '.format(inf),
                'Energy needed to evaluate infinite separation energy')
            inf_ene = None

    return inf_ene
