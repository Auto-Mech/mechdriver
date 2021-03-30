""" es_runners for coordinate scans
"""

import autofile
from mechroutines.es._routines import sp
from mechroutines.es.runner import qchem_params
from mechlib import filesys
from mechlib.amech_io import printer as ioprinter


# CALCULATE INFINITE SEPARATION ENERGY #
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
