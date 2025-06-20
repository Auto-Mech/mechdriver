""" Set up objects used for TS calculations
"""

import elstruct.par
from mechanalyzer.inf import rxn as rinfo
from mechanalyzer.inf import thy as tinfo
from mechlib.filesys import build_fs, rcts_cnf_fs
from mechlib import filesys
from mechroutines.es.runner import multireference_calculation_parameters


def thy_dcts(tsname, ts_dct, thy_dct, es_keyword_dct,
             run_prefix, save_prefix):
    """ set the theory
    """

    # Get transition state info
    ts_zma = ts_dct['zma']
    aspace = ts_dct['active']

    # Set reaction info
    rxn_info = ts_dct['canon_rxn_info']
    ts_info = rinfo.ts_info(rxn_info)
    rct_info = rinfo.rgt_info(rxn_info, 'reacs')
    rxn_info_sort = rinfo.sort(rxn_info)

    # Set the high-sping ts info
    high_mult = rinfo.ts_mult(rxn_info, rxn_mul='high')
    hs_info = (ts_info[0], ts_info[1], high_mult)

    # Set the filesystem locs
    ts_locs = (int(tsname.split('_')[-1]),)
    zma_locs = (ts_dct.get('zma_idx', 0),)

    # Initialize the theory objects
    ini_thy_info, mod_ini_thy_info = None, None
    thy_info, mod_thy_info = None, None
    vscnlvl_thy_info, mod_vscnlvl_thy_info = None, None
    vsp1lvl_thy_info, mod_vsp1lvl_thy_info = None, None
    vsp2lvl_thy_info, mod_vsp2lvl_thy_info = None, None
    hs_vscnlvl_thy_info = None
    hs_vsp1lvl_thy_info = None
    hs_vsp2lvl_thy_info = None
    hs_thy_info = None

    # Method dicts
    ini_method_dct = None
    run_method_dct = None
    vscn_method_dct = None
    vsp1_method_dct = None
    vsp2_method_dct = None

    # Multireference Params
    inplvl_mref_params = {}
    runlvl_mref_params = {}
    vscn_mref_params = {}
    vsp1_mref_params = {}
    vsp2_mref_params = {}

    # Initialize the necessary run filesystem
    runlvl_scn_run_fs = None
    runlvl_cscn_run_fs = None
    vscnlvl_ts_run_fs = None
    vscnlvl_scn_run_fs = None
    vscnlvl_cscn_run_fs = None
    vrctst_run_fs = None

    # Initialize the necessary save filesystem
    inplvl_cnf_save_fs = None
    runlvl_scn_save_fs = None   # above cnf filesys, for search scans
    runlvl_cscn_save_fs = None   # above cnf filesys, for search scans
    runlvl_cnf_save_fs = None
    vscnlvl_scn_save_fs = None
    vscnlvl_cscn_save_fs = None
    vrctst_save_fs = None

    if es_keyword_dct.get('inplvl') is not None:

        ini_method_dct = thy_dct.get(es_keyword_dct['inplvl'])
        ini_thy_info = tinfo.from_dct(ini_method_dct)
        mod_ini_thy_info = tinfo.modify_orb_label(
            ini_thy_info, ts_info)

        # Conformer Layer for saddle point structures
        _, inplvl_cnf_save_fs = build_fs(
            run_prefix, save_prefix, 'CONFORMER',
            rxn_locs=rxn_info_sort, ts_locs=ts_locs,
            thy_locs=mod_ini_thy_info[1:])

        inplvl_loc_info = filesys.mincnf.min_energy_conformer_locators(
            inplvl_cnf_save_fs, mod_ini_thy_info)
        inplvl_min_cnf_locs, _ = inplvl_loc_info
        inplvl_cnf_save_tuple = (inplvl_cnf_save_fs, inplvl_min_cnf_locs)

        if elstruct.par.Method.is_multiref(ini_thy_info[1]):
            inplvl_mref_params = multireference_calculation_parameters(
                ts_zma, ts_info, hs_info,
                aspace, mod_ini_thy_info, rxn_info=rxn_info)

    if es_keyword_dct.get('runlvl') is not None:

        run_method_dct = thy_dct.get(es_keyword_dct['runlvl'])
        thy_info = tinfo.from_dct(run_method_dct)
        mod_thy_info = tinfo.modify_orb_label(
            thy_info, ts_info)
        hs_thy_info = tinfo.modify_orb_label(
            thy_info, hs_info)

        # Conformer Layer for saddle point structures
        runlvl_cnf_run_fs, runlvl_cnf_save_fs = build_fs(
            run_prefix, save_prefix, 'CONFORMER',
            rxn_locs=rxn_info_sort, ts_locs=ts_locs,
            thy_locs=mod_thy_info[1:])
        
        runlvl_loc_info = filesys.mincnf.min_energy_conformer_locators(
            runlvl_cnf_save_fs, mod_thy_info)
        runlvl_min_cnf_locs, _ = runlvl_loc_info
        runlvl_cnf_save_tuple = (runlvl_cnf_save_fs, runlvl_min_cnf_locs)

        # TS/IDX/Z Layer used for coordinate scans (perhaps for guesses)
        _, runlvl_z_save_fs = build_fs(
            run_prefix, save_prefix, 'ZMATRIX',
            rxn_locs=rxn_info_sort, ts_locs=ts_locs,
            thy_locs=mod_thy_info[1:])

        runlvl_scn_run_fs, runlvl_scn_save_fs = build_fs(
            run_prefix, save_prefix, 'SCAN',
            rxn_locs=rxn_info_sort, ts_locs=ts_locs,
            thy_locs=mod_thy_info[1:], zma_locs=zma_locs)

        runlvl_cscn_run_fs, runlvl_cscn_save_fs = build_fs(
            run_prefix, save_prefix, 'CSCAN',
            rxn_locs=rxn_info_sort, ts_locs=ts_locs,
            thy_locs=mod_thy_info[1:], zma_locs=zma_locs)

        if elstruct.par.Method.is_multiref(thy_info[1]):
            runlvl_mref_params = multireference_calculation_parameters(
                ts_zma, ts_info, hs_info,
                aspace, mod_thy_info, rxn_info=rxn_info)

    if es_keyword_dct.get('var_scnlvl') is not None:

        vscn_method_dct = thy_dct.get(es_keyword_dct['var_scnlvl'])
        vscnlvl_thy_info = tinfo.from_dct(vscn_method_dct)
        mod_vscnlvl_thy_info = tinfo.modify_orb_label(
            vscnlvl_thy_info, ts_info)
        hs_vscnlvl_thy_info = tinfo.modify_orb_label(
            vscnlvl_thy_info, hs_info)

        vscnlvl_scn_run_fs, vscnlvl_scn_save_fs = build_fs(
            run_prefix, save_prefix, 'SCAN',
            rxn_locs=rxn_info_sort, ts_locs=ts_locs,
            thy_locs=mod_vscnlvl_thy_info[1:], zma_locs=zma_locs)

        vscnlvl_cscn_run_fs, vscnlvl_cscn_save_fs = build_fs(
            run_prefix, save_prefix, 'CSCAN',
            rxn_locs=rxn_info_sort, ts_locs=ts_locs,
            thy_locs=mod_vscnlvl_thy_info[1:], zma_locs=zma_locs)

        vrctst_run_fs, vrctst_save_fs = build_fs(
            run_prefix, save_prefix, 'VRCTST',
            rxn_locs=rxn_info_sort, ts_locs=ts_locs,
            thy_locs=mod_vscnlvl_thy_info[1:])

        if elstruct.par.Method.is_multiref(vscnlvl_thy_info[1]):
            vscn_mref_params = multireference_calculation_parameters(
                ts_zma, ts_info, hs_info,
                aspace, mod_vscnlvl_thy_info, rxn_info=rxn_info)

    if es_keyword_dct.get('var_splvl1') is not None:

        vsp1_method_dct = thy_dct.get(es_keyword_dct['var_splvl1'])
        vsp1lvl_thy_info = tinfo.from_dct(vsp1_method_dct)
        mod_vsp1lvl_thy_info = tinfo.modify_orb_label(
            vsp1lvl_thy_info, ts_info)
        hs_vsp1lvl_thy_info = tinfo.modify_orb_label(
            vsp1lvl_thy_info, hs_info)

        if elstruct.par.Method.is_multiref(vsp1lvl_thy_info[1]):
            vsp1_mref_params = multireference_calculation_parameters(
                ts_zma, ts_info, hs_info,
                aspace, mod_vsp1lvl_thy_info, rxn_info=rxn_info)

    if es_keyword_dct.get('var_splvl2') is not None:

        vsp2_method_dct = thy_dct.get(es_keyword_dct['var_splvl2'])
        vsp2lvl_thy_info = tinfo.from_dct(vsp2_method_dct)
        mod_vsp2lvl_thy_info = tinfo.modify_orb_label(
            vsp2lvl_thy_info, ts_info)
        hs_vsp2lvl_thy_info = tinfo.modify_orb_label(
            vsp2lvl_thy_info, hs_info)

        if elstruct.par.Method.is_multiref(vsp2lvl_thy_info[1]):
            vsp2_mref_params = multireference_calculation_parameters(
                ts_zma, ts_info, hs_info,
                aspace, mod_vsp2lvl_thy_info, rxn_info=rxn_info)

    # Get the conformer filesys for the reactants
    _rcts_cnf_fs = rcts_cnf_fs(rct_info, ini_thy_info, run_prefix, save_prefix)
    if _rcts_cnf_fs is None:
        _rcts_cnf_fs = rcts_cnf_fs(rct_info, thy_info, run_prefix, save_prefix)

    thy_inf_dct = {
        'inplvl': ini_thy_info,
        'runlvl': thy_info,
        'var_scnlvl': vscnlvl_thy_info,
        'var_splvl1': vsp1lvl_thy_info,
        'var_splvl2': vsp2lvl_thy_info,
        'mod_inplvl': mod_ini_thy_info,
        'mod_runlvl': mod_thy_info,
        'mod_var_scnlvl': mod_vscnlvl_thy_info,
        'mod_var_splvl1': mod_vsp1lvl_thy_info,
        'mod_var_splvl2': mod_vsp2lvl_thy_info,
        'hs_var_scnlvl': hs_vscnlvl_thy_info,
        'hs_var_splvl1': hs_vsp1lvl_thy_info,
        'hs_var_splvl2': hs_vsp2lvl_thy_info,
        'hs_runlvl': hs_thy_info,
        'rct_info': rct_info
    }

    thy_method_dct = {
        'inplvl': ini_method_dct,
        'runlvl': run_method_dct,
        'var_scnlvl': vscn_method_dct,
        'var_splvl1': vsp1_method_dct,
        'var_splvl2': vsp2_method_dct,
    }

    mref_params_dct = {
        'inplvl': inplvl_mref_params,
        'runlvl': runlvl_mref_params,
        'var_scnlvl': vscn_mref_params,
        'var_splvl1': vsp1_mref_params,
        'var_splvl2': vsp2_mref_params,
    }

    runfs_dct = {
        'runlvl_scn': runlvl_scn_run_fs,
        'runlvl_cscn': runlvl_cscn_run_fs,
        'runlvl_cnf': runlvl_cnf_run_fs,
        'vscnlvl_ts': vscnlvl_ts_run_fs,
        'vscnlvl_scn': vscnlvl_scn_run_fs,
        'vscnlvl_cscn': vscnlvl_cscn_run_fs,
        'vrctst': vrctst_run_fs,
        'prefix': run_prefix
    }

    savefs_dct = {
        'runlvl_scn': runlvl_scn_save_fs,
        'runlvl_cscn': runlvl_cscn_save_fs,
        'runlvl_ts_zma': runlvl_z_save_fs,
        'inplvl_cnf_tuple': inplvl_cnf_save_tuple,
        'runlvl_cnf_tuple': runlvl_cnf_save_tuple,
        'vscnlvl_scn': vscnlvl_scn_save_fs,
        'vscnlvl_cscn': vscnlvl_cscn_save_fs,
        'vrctst': vrctst_save_fs,
        'rcts_cnf': _rcts_cnf_fs,
        'prefix': save_prefix
    }

    return (thy_inf_dct, thy_method_dct,
            mref_params_dct,
            runfs_dct, savefs_dct)
