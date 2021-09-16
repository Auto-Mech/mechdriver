"""
  TS Finding algorithms

    2-TS
    (1) mult. sadpts for TS (comes into task from ID)
    (2) ring puckering (ID after initial sadpt found)
    (3) different stereo in migration (ID after sadpt. I think)

"""

import autofile
from mechanalyzer.inf import rxn as rinfo
from mechanalyzer.inf import thy as tinfo
from mechlib.filesys import build_fs, rcts_cnf_fs
from mechlib import filesys


def set_thy_inf_dcts(tsname, ts_dct, thy_dct, es_keyword_dct,
                     run_prefix, save_prefix):
    """ set the theory
    """

    rxn_info = ts_dct['rxn_info']
    ts_info = rinfo.ts_info(rxn_info)
    rct_info = rinfo.rgt_info(rxn_info, 'reacs')
    rxn_info = rinfo.sort(rxn_info)

    ts_locs = (int(tsname.split('_')[-1]),)

    high_mult = rinfo.ts_mult(rxn_info, rxn_mul='high')

    # Set the hs info
    hs_info = (ts_info[0], ts_info[1], high_mult)

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

    # Initialize the necessary run filesystem
    runlvl_scn_run_fs = None
    runlvl_cscn_run_fs = None
    vscnlvl_ts_run_fs = None
    vscnlvl_scn_run_fs = None
    vscnlvl_cscn_run_fs = None
    vrctst_run_fs = None

    # Initialize the necessary save filesystem
    ini_zma_save_fs = None
    runlvl_scn_save_fs = None   # above cnf filesys, for search scans
    runlvl_cscn_save_fs = None   # above cnf filesys, for search scans
    runlvl_cnf_save_fs = None
    vscnlvl_scn_save_fs = None
    vscnlvl_cscn_save_fs = None
    vrctst_save_fs = None

    if es_keyword_dct.get('inplvl', None) is not None:

        ini_method_dct = thy_dct.get(es_keyword_dct['inplvl'])
        ini_thy_info = tinfo.from_dct(ini_method_dct)
        mod_ini_thy_info = tinfo.modify_orb_label(
            ini_thy_info, ts_info)

        _, ini_cnf_save_fs = build_fs(
            run_prefix, save_prefix, 'CONFORMER',
            rxn_locs=rxn_info, ts_locs=ts_locs,
            thy_locs=mod_ini_thy_info[1:])

        ini_loc_info = filesys.mincnf.min_energy_conformer_locators(
            ini_cnf_save_fs, mod_ini_thy_info)
        ini_min_cnf_locs, ini_path = ini_loc_info

        if ini_path:
            ini_zma_save_fs = autofile.fs.zmatrix(
                ini_cnf_save_fs[-1].path(ini_min_cnf_locs))

    if es_keyword_dct.get('runlvl', None) is not None:

        method_dct = thy_dct.get(es_keyword_dct['runlvl'])
        thy_info = tinfo.from_dct(method_dct)
        mod_thy_info = tinfo.modify_orb_label(
            thy_info, ts_info)
        hs_thy_info = tinfo.modify_orb_label(
            thy_info, hs_info)

        runlvl_cnf_run_fs, runlvl_cnf_save_fs = build_fs(
            run_prefix, save_prefix, 'CONFORMER',
            rxn_locs=rxn_info, ts_locs=ts_locs,
            thy_locs=mod_thy_info[1:])

        _, runlvl_ts_zma_save_fs = build_fs(
            run_prefix, save_prefix, 'ZMATRIX',
            rxn_locs=rxn_info, ts_locs=ts_locs,
            thy_locs=mod_thy_info[1:])

        runlvl_loc_info = filesys.mincnf.min_energy_conformer_locators(
            runlvl_cnf_save_fs, mod_thy_info)
        runlvl_min_cnf_locs, _ = runlvl_loc_info
        runlvl_cnf_save_fs = (runlvl_cnf_save_fs, runlvl_min_cnf_locs)

        runlvl_scn_run_fs, runlvl_scn_save_fs = build_fs(
            run_prefix, save_prefix, 'SCAN',
            rxn_locs=rxn_info, ts_locs=ts_locs,
            thy_locs=mod_thy_info[1:], zma_locs=(0,))

        runlvl_cscn_run_fs, runlvl_cscn_save_fs = build_fs(
            run_prefix, save_prefix, 'CSCAN',
            rxn_locs=rxn_info, ts_locs=ts_locs,
            thy_locs=mod_thy_info[1:], zma_locs=(0,))

    if es_keyword_dct.get('var_scnlvl', None) is not None:

        method_dct = thy_dct.get(es_keyword_dct['var_scnlvl'])
        vscnlvl_thy_info = tinfo.from_dct(method_dct)
        mod_vscnlvl_thy_info = tinfo.modify_orb_label(
            vscnlvl_thy_info, ts_info)
        hs_vscnlvl_thy_info = tinfo.modify_orb_label(
            vscnlvl_thy_info, hs_info)

        vscnlvl_scn_run_fs, vscnlvl_scn_save_fs = build_fs(
            run_prefix, save_prefix, 'SCAN',
            rxn_locs=rxn_info, ts_locs=ts_locs,
            thy_locs=mod_thy_info[1:], zma_locs=(0,))

        vscnlvl_cscn_run_fs, vscnlvl_cscn_save_fs = build_fs(
            run_prefix, save_prefix, 'CSCAN',
            rxn_locs=rxn_info, ts_locs=ts_locs,
            thy_locs=mod_thy_info[1:], zma_locs=(0,))

        vrctst_run_fs, vrctst_save_fs = build_fs(
            run_prefix, save_prefix, 'VRCTST',
            rxn_locs=rxn_info, ts_locs=ts_locs,
            thy_locs=mod_thy_info[1:])

        if es_keyword_dct.get('var_splvl1', None) is not None:

            method_dct = thy_dct.get(es_keyword_dct['var_scnlvl1'])
            vsp1lvl_thy_info = tinfo.from_dct(method_dct)
            mod_vsp1lvl_thy_info = tinfo.modify_orb_label(
                vsp1lvl_thy_info, ts_info)
            hs_vsp1lvl_thy_info = tinfo.modify_orb_label(
                vsp1lvl_thy_info, hs_info)

        if es_keyword_dct.get('var_splvl2', None) is not None:

            method_dct = thy_dct.get(es_keyword_dct['var_scnlvl2'])
            vsp2lvl_thy_info = tinfo.from_dct(method_dct)
            mod_vsp2lvl_thy_info = tinfo.modify_orb_label(
                vsp2lvl_thy_info, ts_info)
            hs_vsp2lvl_thy_info = tinfo.modify_orb_label(
                vsp2lvl_thy_info, hs_info)

    # Get the conformer filesys for the reactants
    _rcts_cnf_fs = rcts_cnf_fs(
        rct_info, thy_dct, es_keyword_dct, run_prefix, save_prefix)

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
        'hs_runlvl': hs_thy_info
    }

    runfs_dct = {
        'runlvl_scn_fs': runlvl_scn_run_fs,
        'runlvl_cscn_fs': runlvl_cscn_run_fs,
        'runlvl_cnf_fs': runlvl_cnf_run_fs,
        'vscnlvl_ts_fs': vscnlvl_ts_run_fs,
        'vscnlvl_scn_fs': vscnlvl_scn_run_fs,
        'vscnlvl_cscn_fs': vscnlvl_cscn_run_fs,
        'vrctst_fs': vrctst_run_fs,
    }

    savefs_dct = {
        'inilvl_zma_fs': ini_zma_save_fs,
        'runlvl_scn_fs': runlvl_scn_save_fs,
        'runlvl_cscn_fs': runlvl_cscn_save_fs,
        'runlvl_cnf_fs': runlvl_cnf_save_fs,
        'vscnlvl_scn_fs': vscnlvl_scn_save_fs,
        'vscnlvl_cscn_fs': vscnlvl_cscn_save_fs,
        'vrctst_fs': vrctst_save_fs,
        'rcts_cnf_fs': _rcts_cnf_fs,
        'runlvl_ts_zma_fs': runlvl_ts_zma_save_fs
    }

    return thy_inf_dct, runfs_dct, savefs_dct


