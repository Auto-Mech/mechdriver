""" Calculations along reaction coordinates x
"""

import autofile
import elstruct
import automol
from mechanalyzer.inf import rxn as rinfo
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
                              cscn_run_fs, cscn_save_fs,
                              es_keyword_dct,
                              find_max=True):
    """ Scan along the internal coordinates that correspond to the
        reaction coordinate. Additional constraints will be used as needed.

        Both the internal coordinate and constrained coordinates are set
        according to reaction class.
    """

    # Determine if scan should be for variational reaction class

    # Build grid and names appropriate for reaction type
    var = (automol.ReactionInfo.is_radical_radical(rxn_class) and
           automol.ReactionInfo.is_low_spin(rxn_class))
    scan_inf = automol.reac.build_scan_info(zrxn, ts_zma, var=var)
    coord_names, constraint_dct, coord_grids, update_guess = scan_inf

    # Set the filesystem
    if constraint_dct is None:
        _scn_run_fs, _scn_save_fs = scn_run_fs, scn_save_fs
    else:
        _scn_run_fs, _scn_save_fs = cscn_run_fs, cscn_save_fs
        # these lines are to help ring forming scissions, coudl be a bad idea
        for key, val in constraint_dct.items():
            if abs(val - 1.57) < .05:
                val = val - .25
                constraint_dct[key] = val
    # Set up script string and kwargs
    mod_thy_info = tinfo.modify_orb_label(tinfo.from_dct(method_dct), ts_info)
    script_str, kwargs = qchem_params(
        method_dct, job=elstruct.Job.OPTIMIZATION,
        geo=automol.zmat.geometry(ts_zma), spc_info=ts_info)
    kwargs.update(mref_params)

    es_runner.scan.execute_scan(
        zma=ts_zma,
        spc_info=ts_info,
        mod_thy_info=mod_thy_info,
        coord_names=coord_names,
        coord_grids=coord_grids,
        scn_run_fs=_scn_run_fs,
        scn_save_fs=_scn_save_fs,
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
    scan_finished_ = es_runner.scan.scan_finished(
        coord_names=coord_names,
        coord_grids=coord_grids,
        scn_save_fs=_scn_save_fs,
        constraint_dct=constraint_dct,
    )
    if scan_finished_ and find_max:
        include_endpts = not mref_params
        max_zmas = rxngrid.grid_maximum_zmatrices(
            automol.reac.class_(zrxn), ts_zma, coord_grids, coord_names, _scn_save_fs,
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

    ioprinter.info_message('Calculating the Infinite Separation Energy...')
    ioprinter.info_message('')

    # Get info from the reactants
    rct_info = thy_inf_dct['rct_info']
    overwrite = es_keyword_dct['overwrite']

    # Get thy_inf_dct stuff
    thy_info = thy_inf_dct['runlvl']
    var_scn_thy_info = thy_inf_dct['var_scnlvl']
    var_sp1_thy_info = thy_inf_dct['var_splvl1']
    var_sp2_thy_info = thy_inf_dct['var_splvl2']
    hs_var_scn_thy_info = thy_inf_dct['hs_var_scnlvl']
    hs_var_sp1_thy_info = thy_inf_dct['hs_var_splvl1']
    hs_var_sp2_thy_info = thy_inf_dct['hs_var_splvl2']
    vscn_method_dct = thy_method_dct['var_scnlvl']
    var_sp1_method_dct = thy_method_dct['var_splvl1']
    var_sp2_method_dct = thy_method_dct['var_splvl2']

    # Get the filesys stuff
    rcts_cnf_fs = savefs_dct['rcts_cnf']
    vscnlvl_scn_run_fs = runfs_dct['vscnlvl_scn']
    vscnlvl_scn_save_fs = savefs_dct['vscnlvl_scn']
    run_prefix = runfs_dct['prefix']

    # Get kwargs for the calculation
    rxn_class = ts_dct['class']
    var = (automol.ReactionInfo.is_radical_radical(rxn_class) and
           automol.ReactionInfo.is_low_spin(rxn_class))
    if var:
        ts_zma, zrxn = ts_dct['zma'], ts_dct['zrxn']
        rxn_info = ts_dct['rxn_info']

        ts_info = rinfo.ts_info(rxn_info)
        high_mult = rinfo.ts_mult(rxn_info, rxn_mul='high')
        hs_ts_info = (ts_info[0], ts_info[1], high_mult)

        # Build grid and names appropriate for reaction type
        names, _, grids, _ = automol.reac.build_scan_info(
            zrxn, ts_zma, var=var)
        inf_locs = (names, (grids[1][-1],))

        cas_kwargs = mref_params['var_scnlvl']

        _inf_sep_ene = _multiref_inf_sep_ene(
            hs_ts_info, ts_zma,
            rct_info, rcts_cnf_fs, run_prefix,
            thy_info,
            var_scn_thy_info,
            var_sp1_thy_info, var_sp2_thy_info,
            hs_var_scn_thy_info,
            hs_var_sp1_thy_info, hs_var_sp2_thy_info,
            vscn_method_dct,
            var_sp1_method_dct, var_sp2_method_dct,
            vscnlvl_scn_run_fs, vscnlvl_scn_save_fs, inf_locs,
            overwrite=overwrite,
            **cas_kwargs)
    else:
        _inf_sep_ene = _singleref_inf_sep_ene(
            rct_info, rcts_cnf_fs,
            vscn_method_dct, thy_info,
            run_prefix, overwrite)

    if _inf_sep_ene is not None:
        ioprinter.energy(_inf_sep_ene)
    else:
        ioprinter.warning_message('Infinite separation ene not computed')

    return _inf_sep_ene


def _multiref_inf_sep_ene(hs_info, ref_zma,
                          rct_info, rcts_cnf_fs, run_prefix,
                          thy_info,
                          var_scn_thy_info,
                          var_sp1_thy_info, var_sp2_thy_info,
                          hs_var_scn_thy_info,
                          hs_var_sp1_thy_info, hs_var_sp2_thy_info,
                          var_scn_method_dct,
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
    # hs_thy_infs = (hs_var_sp2_thy_info, hs_var_sp1_thy_info)
    # # thy_infs = (var_sp2_thy_info, var_sp1_thy_info)
    # method_dcts = (var_sp2_method_dct, var_sp1_method_dct)
    hs_thy_infs = (hs_var_sp2_thy_info, hs_var_scn_thy_info)
    # thy_infs = (var_sp2_thy_info, var_sp1_thy_info)
    method_dcts = (var_sp2_method_dct, var_scn_method_dct)

    # Calculate the energies for the two cases
    for idx, (meth_dct, thy_inf) in enumerate(zip(method_dcts, hs_thy_infs)):

        if idx == 0:
            ioprinter.info_message(
                " - Running high-spin single reference energy ...")
        else:
            ioprinter.info_message(
                " - Running high-spin multi reference energy ...")

        ioprinter.info_message(
            ' - Method:', tinfo.string(var_scn_thy_info, thy_inf))

        geo = scn_save_fs[-1].file.geometry.read(inf_locs)
        zma = scn_save_fs[-1].file.zmatrix.read(inf_locs)

        # Calculate the single point energy
        script_str, kwargs = qchem_params(
            meth_dct, geo=geo, spc_info=hs_info)
        cas_kwargs.update(kwargs)

        # geo = automol.zmat.geometry(ref_zma)
        sp.run_energy(zma, geo, hs_info, thy_inf,
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
            ene = hs_save_fs[-1].file.energy.read(thy_inf[1:4])

        if idx == 0:
            hs_sr_ene = ene
        else:
            hs_mr_ene = ene

    # Get the single reference energy for each of the reactant configurations
    ioprinter.info_message('')
    ioprinter.info_message(
        'Running single-point calculations for reactants (need DFT)...')
    ioprinter.info_message('Method:', tinfo.string(thy_info, var_sp2_thy_info))
    reac_ene = _reac_sep_ene(
        rct_info, rcts_cnf_fs,
        var_sp2_method_dct, var_sp2_thy_info,
        run_prefix, overwrite)

    # Calculate the infinite seperation energy
    all_enes = (reac_ene, hs_sr_ene, hs_mr_ene)
    if all(ene is not None for ene in all_enes):
        _inf_sep_ene = reac_ene - hs_sr_ene + hs_mr_ene
        ioprinter.info_message('Infinite Separation Energy [au]: '
                               f'{_inf_sep_ene}')
    else:
        _inf_sep_ene = None

    print('inf ene components')
    print('reac', reac_ene)
    print('hs sr', hs_sr_ene)
    print('hs mr', hs_mr_ene)

    return _inf_sep_ene


def _singleref_inf_sep_ene(rct_info, rcts_cnf_fs,
                           method_dct, sp_thy_info,
                           run_prefix, overwrite):
    """ Obtain the electronic energy for a set of reactants at infinite
        separation for a single-reference electronic structure method.

        Since single-reference methods are size-consistent, this simply
        amounts to summing the electronic energies of the reactants
        calculated independently.

        Function will attempt to read the energy from the SAVE filesystem
        and calculate it if it does not currently exist.
    """
    return _reac_sep_ene(
        rct_info, rcts_cnf_fs,
        method_dct, sp_thy_info,
        run_prefix, overwrite)


def _reac_sep_ene(rct_info, rcts_cnf_fs, method_dct, sp_thy_info,
                  run_prefix, overwrite):
    """ Determine the sum of electronic energies of two reactants specified
        at the level of theory described in the theory info object. Will
        calculate the energy if it is not currently in the SAVE filesystem.
    """

    # get the single reference energy for each of the reactant configurations
    spc_enes = []
    for (run_fs, save_fs, mlocs, mpath), inf in zip(rcts_cnf_fs, rct_info):

        # Set the modified thy info
        mod_sp_thy_info = tinfo.modify_orb_label(sp_thy_info, inf)

        # Build filesys
        zma_fs = autofile.fs.zmatrix(mpath)

        # Read the geometry and set paths
        zma = zma_fs[-1].file.zmatrix.read([0])
        geo = save_fs[-1].file.geometry.read(mlocs)

        # Set up the run script
        sp_script_str, kwargs = qchem_params(
            method_dct, geo=geo, spc_info=inf)

        # Build the single point filesys objects
        sp_save_fs = autofile.fs.single_point(mpath)

        # Calculate the save single point energy
        sp.run_energy(zma, geo, inf, mod_sp_thy_info,
                      run_fs, save_fs, mlocs, run_prefix,
                      sp_script_str, overwrite,
                      **kwargs)
        exists = sp_save_fs[-1].file.energy.exists(mod_sp_thy_info[1:4])
        if not exists:
            ioprinter.warning_message('No ene found')
            ene = None
        else:
            ene = sp_save_fs[-1].file.energy.read(mod_sp_thy_info[1:4])

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

    if inf_ene is not None:
        ioprinter.info_message(f'Reactant Energy [au]: {inf_ene}')

    return inf_ene
