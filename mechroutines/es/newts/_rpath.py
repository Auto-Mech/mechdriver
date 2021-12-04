""" Calculations along reaction coordinates x
"""

import autofile
import elstruct
import automol
# from mechanalyzer.inf import rxn as rinfo
from mechanalyzer.inf import thy as tinfo
from mechlib.reaction import grid as rxngrid
from mechlib.amech_io import printer as ioprinter
from mechroutines.es import runner as es_runner
from mechroutines.es._routines import sp
from mechroutines.es.runner import qchem_params


# Scans along coordinate(s)
def internal_coordinates_scan(ts_zma, zrxn,
                              ts_info, rxn_class,
                              method_dct, mref_params,
                              scn_run_fs, scn_save_fs,
                              es_keyword_dct,
                              find_max=True):
    """ Scan along the internal coordinates that correspond to the
        reaction coordinate. Additional constraints will be used as needed.

        Both the internal coordinate and constrained coordinates are set
        according to reaction class.
    """

    # Determine if scan should be for variational reaction class

    # Build grid and names appropriate for reaction type
    var = (automol.par.is_radrad(rxn_class) and
           automol.par.is_low_spin(rxn_class))
    scan_inf = automol.reac.build_scan_info(zrxn, ts_zma, var=var)
    coord_names, constraint_dct, coord_grids, update_guess = scan_inf

    # Set up script string and kwargs
    mod_thy_info = tinfo.modify_orb_label(tinfo.from_dct(method_dct), ts_info)
    script_str, kwargs = qchem_params(
        method_dct, job=elstruct.Job.OPTIMIZATION)
    kwargs.update(mref_params)

    es_runner.scan.execute_scan(
        zma=ts_zma,
        spc_info=ts_info,
        mod_thy_info=mod_thy_info,
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
        **kwargs,
    )

    if find_max:
        include_endpts = not mref_params
        max_zmas = rxngrid.grid_maximum_zmatrices(
            zrxn.class_, ts_zma, coord_grids, coord_names, scn_save_fs,
            mod_thy_info, constraint_dct, include_endpts=include_endpts)
    else:
        max_zmas = None

    return max_zmas


# Infinite Separation Energy
def inf_sep_ene(ts_dct, thy_inf_dct, thy_method_dct, mref_params,
                savefs_dct, runfs_dct, es_keyword_dct):
    """ Determine the total electronic energy of two reacting species that
        at infinite separation.
    """

    # Get info from the reactants
    rct_info = thy_inf_dct['rct_info']
    overwrite = es_keyword_dct['overwrite']

    # Get thy_inf_dct stuff
    var_sp1_thy_info = thy_inf_dct['var_splvl2']
    var_sp2_thy_info = thy_inf_dct['var_splvl2']
    hs_var_sp1_thy_info = thy_inf_dct['hs_var_splvl1']
    hs_var_sp2_thy_info = thy_inf_dct['hs_var_splvl2']
    vscn_method_dct = thy_method_dct['var_scnlvl']
    var_sp1_method_dct = thy_method_dct['var_splvl1']
    var_sp2_method_dct = thy_method_dct['var_splvl2']

    print('v1 dct', var_sp1_method_dct)
    print('v2 dct', var_sp2_method_dct)

    # Get the filesys stuff
    rcts_cnf_fs = savefs_dct['rcts_cnf']
    vscnlvl_scn_run_fs = runfs_dct['vscnlvl_scn']
    vscnlvl_scn_save_fs = savefs_dct['vscnlvl_scn']
    run_prefix = runfs_dct['prefix']

    # Get kwargs for the calculation
    rxn_class = ts_dct['class']
    var = (automol.par.is_radrad(rxn_class) and
           automol.par.is_low_spin(rxn_class))
    if var:
        ts_zma, zrxn = ts_dct['zma'], ts_dct['zrxn']
        # rxn_info = ts_dct['rxn_info']
        # aspace = ts_dct['active']

        # ts_info = rinfo.ts_info(rxn_info)
        hs_ts_info = thy_inf_dct['hs_var_scnlvl']
        # mod_thy_info = thy_inf_dct['mod_vscnlvl_thy_info']

        # Build grid and names appropriate for reaction type
        names, _, grids, _ = automol.reac.build_scan_info(zrxn, ts_zma)
        far_locs = [[names[0]], [grids[0]]]

        # cas_kwargs = es_runner.multireference_calculation_parameters(
        #     ts_zma, ts_info, hs_ts_info, rxn_info,
        #     aspace, mod_thy_info)
        cas_kwargs = mref_params['var_scnlvl']

        _inf_sep_ene = _multiref_inf_sep_ene(
            hs_ts_info, ts_zma,
            rct_info, rcts_cnf_fs, run_prefix,
            var_sp1_thy_info, var_sp2_thy_info,
            hs_var_sp1_thy_info, hs_var_sp2_thy_info,
            var_sp1_method_dct, var_sp2_method_dct,
            vscnlvl_scn_run_fs, vscnlvl_scn_save_fs, far_locs,
            overwrite=overwrite,
            **cas_kwargs)
    else:
        vscn_thy_info = thy_inf_dct['var_scnlvl']
        _inf_sep_ene = _singleref_inf_sep_ene(
            rct_info, vscn_thy_info,
            rcts_cnf_fs, run_prefix,
            vscn_method_dct, overwrite)

    if _inf_sep_ene is not None:
        ioprinter.energy(_inf_sep_ene)
    else:
        ioprinter.warning_message('Infinite separation ene not computed')

    return _inf_sep_ene


def _multiref_inf_sep_ene(hs_info, ref_zma,
                          rct_info, rcts_cnf_fs, run_prefix,
                          var_sp1_thy_info, var_sp2_thy_info,
                          hs_var_sp1_thy_info, hs_var_sp2_thy_info,
                          var_sp1_method_dct, var_sp2_method_dct,
                          scn_run_fs, scn_save_fs, inf_locs,
                          overwrite=False,
                          **cas_kwargs):
    """ Obtain the electronic energy for a set of reactants at infinite
        separation for a multi-reference electronic structure method.

        Since multireference methods are not size-consistent, we cannot
        sum the electronic energies of the two reactants from individual
        calculations. One could determine this energy by calculating the
        it where the two reactants are in the same input but the intermolecular
        distance is set arbitrarily large; unfortunately, this leads to
        convergence issues.

        To resulve this, the following approach is used:
        At a given reference point, the high-spin low-spin splitting at that
        reference point, and the high level energy for the high spin state
        at the reference geometry and for the fragments
        scn = thy for optimizations
        sp1 = low-spin single points
        sp2 = high-spin single points for inf sep

        inf = spc0 + spc1 - hs_sr_e + hs_mr_ene
    """

    # Set groups for loops
    hs_thy_infs = (hs_var_sp2_thy_info, hs_var_sp1_thy_info)
    thy_infs = (var_sp2_thy_info, var_sp1_thy_info)
    method_dcts = (var_sp2_method_dct, var_sp1_method_dct)

    # Calculate the energies for the two cases
    for idx, (meth_dct, thy_inf) in enumerate(zip(method_dcts, hs_thy_infs)):

        if idx == 0:
            ioprinter.info_message(
                " - Running high-spin single reference energy ...")
        else:
            ioprinter.info_message(
                " - Running high-spin multi reference energy ...")

        # Calculate the single point energy
        print('meth dct test', meth_dct)
        script_str, kwargs = qchem_params(meth_dct)
        cas_kwargs.update(kwargs)

        geo = automol.zmat.geometry(ref_zma)
        sp.run_energy(ref_zma, geo, hs_info, thy_inf,
                      scn_run_fs, scn_save_fs, inf_locs, run_prefix,
                      script_str, overwrite, highspin=True, **cas_kwargs)

        # Read the energty from the filesystem
        geo_save_path = scn_save_fs[-1].path(inf_locs)
        hs_save_fs = autofile.fs.high_spin(geo_save_path)
        if not hs_save_fs[-1].file.energy.exists(thy_inf[1:4]):
            ioprinter.error_message(
                'High-spin energy job failed: ',
                'energy is needed to evaluate infinite separation energy')
            ene = None
        else:
            ioprinter.info_message(
                " - Reading high-spin energy from filesystem...")
            ene = hs_save_fs[-1].file.energy.read(thy_inf[1:4])

        if idx == 0:
            hs_sr_ene = ene
        else:
            hs_mr_ene = ene

    # Get the single reference energy for each of the reactant configurations
    for method_dct, thy_inf in zip(method_dcts, thy_infs):
        sp_script_str, kwargs = qchem_params(method_dct)
        reac_ene = _reac_sep_ene(
            rct_info, thy_inf,
            rcts_cnf_fs, run_prefix,
            overwrite, sp_script_str, **kwargs)

    # Calculate the infinite seperation energy
    all_enes = (reac_ene, hs_sr_ene, hs_mr_ene)
    if all(ene is not None for ene in all_enes):
        _inf_sep_ene = reac_ene - hs_sr_ene + hs_mr_ene
    else:
        _inf_sep_ene = None

    return _inf_sep_ene


def _singleref_inf_sep_ene(rct_info, thy_info,
                           rcts_cnf_fs, run_prefix,
                           method_dct, overwrite):
    """ Obtain the electronic energy for a set of reactants at infinite
        separation for a single-reference electronic structure method.

        Since single-reference methods are size-consistent, this simply
        amounts to summing the electronic energies of the reactants
        calculated independently.

        Function will attempt to read the energy from the SAVE filesystem
        and calculate it if it does not currently exist.
    """
    script_str, kwargs = qchem_params(method_dct)
    return _reac_sep_ene(
        rct_info, thy_info,
        rcts_cnf_fs, run_prefix,
        overwrite, script_str, **kwargs)


def _reac_sep_ene(rct_info, thy_info, rcts_cnf_fs, run_prefix,
                  overwrite, sp_script_str, **kwargs):
    """ Determine the sum of electronic energies of two reactants specified
        at the level of theory described in the theory info object. Will
        calculate the energy if it is not currently in the SAVE filesystem.
    """

    # get the single reference energy for each of the reactant configurations
    spc_enes = []
    for (run_fs, save_fs, mlocs, mpath), inf in zip(rcts_cnf_fs, rct_info):

        # Set the modified thy info
        mod_thy_info = tinfo.modify_orb_label(thy_info, inf)

        # Build filesys
        zma_fs = autofile.fs.zmatrix(mpath)

        # Read the geometry and set paths
        zma = zma_fs[-1].file.zmatrix.read([0])
        geo = save_fs[-1].file.geometry.read(mlocs)

        # Build the single point filesys objects
        sp_save_fs = autofile.fs.single_point(mpath)

        # Calculate the save single point energy
        sp.run_energy(zma, geo, inf, mod_thy_info,
                      run_fs, save_fs, mlocs, run_prefix,
                      sp_script_str, overwrite,
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
                f'for {inf}: ',
                'Energy needed to evaluate infinite separation energy')
            inf_ene = None

    return inf_ene
