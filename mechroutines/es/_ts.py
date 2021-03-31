"""
  TS Finding algorithms

    2-TS
    (1) mult. sadpts for TS (comes into task from ID)
    (2) ring puckering (ID after initial sadpt found)
    (3) different stereo in migration (ID after sadpt. I think)

"""

import automol
import autofile
from mechanalyzer.inf import rxn as rinfo
from mechanalyzer.inf import thy as tinfo
from mechroutines.es._routines import _sadpt as sadpt
# from mechroutines.es._routines import _vrctst as vrctst
# from mechroutines.es._routines import _vtst as vtst
from mechlib.filesys import build_fs
from mechlib import filesys


def findts(tsk, spc_dct, tsname, thy_dct, es_keyword_dct,
           run_prefix, save_prefix):
    """ New run function
    """

    method_dct = thy_dct.get(es_keyword_dct['runlvl'])
    # ini_method_dct = thy_dct.get(es_keyword_dct['inplvl'])

    # Set the TS searching algorithm to use: (1) Check dct, (2) Set by Class
    search_thy_inf_dct = 'sadpt'
    # search_thy_inf_dct = _ts_finder_match(tsk, spc_dct[tsname])

    # Build necessary objects
    thy_inf_dct, runfs_dct, savefs_dct = _set_thy_inf_dcts(
        tsname, spc_dct[tsname], thy_dct, es_keyword_dct,
        run_prefix, save_prefix)

    # Find the transition state
    if search_thy_inf_dct == 'sadpt':
        run_sadpt(spc_dct, tsname, method_dct, es_keyword_dct,
                  thy_inf_dct, runfs_dct, savefs_dct)
    # elif search_thy_inf_dct == 'vtst':
    #     run_vtst(spc_dct, tsname, es_keyword_dct,
    #              thy_inf_dct, runfs_dct, savefs_dct)
    elif search_thy_inf_dct == 'vrctst':
        run_vrctst(spc_dct, tsname, es_keyword_dct,
                   thy_inf_dct, runfs_dct, savefs_dct)
    elif search_thy_inf_dct is None:
        print('No TS search algorithm was specified or able to determined')


def run_sadpt(spc_dct, tsname, method_dct, es_keyword_dct,
              thy_inf_dct, runfs_dct, savefs_dct):
    """ find a transition state
    """

    # Get objects for the calculations
    ts_dct = spc_dct[tsname]

    # Find the TS
    cnf_info = savefs_dct['runlvl_cnf_fs']
    cnf_save_fs, cnf_save_locs = cnf_info
    overwrite = es_keyword_dct['overwrite']
    if not cnf_save_locs[0]:
        print('No transition state found in filesys',
              'at {} level...'.format(es_keyword_dct['runlvl']),
              'Proceeding to find it...')
        _run = True
    elif overwrite:
        print('User specified to overwrite energy with new run...')
        _run = True
    else:
        print('TS found and saved previously in ',
              cnf_save_fs[-1].path(cnf_save_locs))
        _run = False

    if _run:
        # split below in guess, scan
        guess_zmas = sadpt.generate_guess_structure(
            ts_dct, method_dct, es_keyword_dct,
            runfs_dct, savefs_dct)
        sadpt.obtain_saddle_point(
            guess_zmas, ts_dct, method_dct,
            runfs_dct, savefs_dct, es_keyword_dct)

        # Generate a second sadpt


def run_vrctst(spc_dct, tsname, es_keyword_dct,
               thy_inf_dct, runfs_dct, savefs_dct,
               info_dct, grid, run_prefix, save_prefix):
    """ find a transition state
    """

    # Get dct for specific species task is run for
    ts_dct = spc_dct[tsname]

    # Get info from the reactants
    high_mul = ts_dct['high_mult']
    rct_zmas = ts_dct['rct_zmas']
    ts_info = info_dct['ts_info']
    hs_info = info_dct['hs_info']
    rct_info = info_dct['rct_info']
    rct_ichs = [spc_dct[rct]['inchi'] for rct in ts_dct['reacs']]

    # Set information from the transition state
    high_mul = ts_dct['high_mult']
    ini_zma = ts_dct['zma']
    frm_bnd_keys = ts_dct['frm_bnd_keys']
    ts_formula = automol.geom.formula(automol.zmatrix.geometry(ini_zma))
    active_space = ts_dct['active_space']

    # Get reaction coordinates
    frm_name = automol.zmatrix.bond_key_from_idxs(ini_zma, frm_bnd_keys)

    # Get es options
    overwrite = es_keyword_dct['overwrite']
    update_guess = False  # check

    # Grid
    [grid1, grid2] = grid

    # Get thy_inf_dct stuff
    mod_var_scn_thy_info = thy_inf_dct['var_scnlvl']
    mod_var_sp1_thy_info = thy_inf_dct['var_splvl1']
    mod_var_sp2_thy_info = thy_inf_dct['var_splvl2']
    hs_var_sp1_thy_info = thy_inf_dct['hs_var_splvl1']
    hs_var_sp2_thy_info = thy_inf_dct['hs_var_splvl2']

    # Get the filesys stuff
    vscnlvl_scn_save_fs = savefs_dct['vscnlvl_scn_fs']
    vscnlvl_scn_run_fs = runfs_dct['vscnlvl_scn_fs']
    vscnlvl_cscn_save_fs = savefs_dct['vscnlvl_cscn_fs']
    vscnlvl_cscn_run_fs = runfs_dct['vscnlvl_cscn_fs']
    rcts_cnf_fs = savefs_dct['rcts_cnf_fs']
    vscnlvl_thy_save_fs = savefs_dct['vscnlvl_thy_fs']
    vscnlvl_ts_save_fs = savefs_dct['vscnlvl_ts_fs']
    vscnlvl_ts_run_fs = runfs_dct['vscnlvl_ts_fs']

    print('Beginning Calculations for VRC-TST Treatments')
    vrctst.calc_vrctst_flux(
        ini_zma, ts_info, hs_info,
        ts_formula, high_mul, active_space,
        rct_info, rct_ichs, rct_zmas, rcts_cnf_fs,
        grid1, grid2, frm_name,
        mod_var_scn_thy_info,
        mod_var_sp1_thy_info, mod_var_sp2_thy_info,
        hs_var_sp1_thy_info, hs_var_sp2_thy_info,
        vscnlvl_thy_save_fs,
        vscnlvl_ts_save_fs,
        vscnlvl_ts_run_fs,
        vscnlvl_scn_run_fs, vscnlvl_scn_save_fs,
        vscnlvl_cscn_run_fs, vscnlvl_cscn_save_fs,
        run_prefix, save_prefix,
        overwrite, update_guess)


# SET THE SEARCHING ALGORITHM
def _ts_finder_match(tsk, ts_dct):
    """ Determine the algorithm that should be used for a given transition
        state by looking at the requested user input or determining
    """

    print('Determining if TS search algorithm matches the task requested')

    # Set search algorithm to one specified by the user, if specified
    if 'ts_search' in ts_dct:
        ini_thy_inf_dct = [ts_dct['ts_search']]
        print('Running search algorithm according to {},'.format(
            ts_dct['ts_search']),
              'as requested by the user')
    else:
        ini_thy_inf_dct = None
        print('No search algorithm requested')
    print()

    # ID search algorithm if user did not specify one (wrong)
    if ini_thy_inf_dct is None:
        if _nobarrier(ts_dct):
            ini_thy_inf_dct = ['radrad_vtst', 'vrctst']
            print('Reaction is low-spin, radical-radical addition/abstraction')
            print('Assuming reaction is barrierless...')
            print('Finding a transition state according to either vtst or '
                  'vrctst, depending on the current task')
        else:
            ini_thy_inf_dct = ['sadpt']
            print('Assuming reaction has saddle point on potential surface...')
            print('Use species.dat to specify VTST search for mol-rad rxn...')
            print('Finding the geometry of the saddle point...')

    # Print message if no algorithm found
    if ini_thy_inf_dct is None:
        print('No TS search algorithm was specified or able to determined')

    # Set return for ts searching algorithm if there is one
    if tsk in ini_thy_inf_dct:
        print('Search algorithm matches task')
        search_thy_inf_dct = tsk
    else:
        print('Algorithm does not match task')
        search_thy_inf_dct = None

    # Refine ret vtst thy_inf_dct if that is what is being used
    if search_thy_inf_dct == 'vtst':
        search_thy_inf_dct = 'radrad_vtst' if _radrad(ts_dct) else 'molrad_vtst'

    return search_thy_inf_dct


# CHECKS FOR TYPE OF TRANSITION STATE
def _nobarrier(ts_dct):
    """ Determine if reaction is barrierless
    """
    radrad = _radrad(ts_dct)
    low_spin = bool('low' in ts_dct['class'])
    return radrad and low_spin


def _radrad(ts_dct):
    return bool('radical radical' in ts_dct['class'])


# SET OPTIONS FOR THE TRANSITION STATE
def _set_thy_inf_dcts(tsname, ts_dct, thy_dct, es_keyword_dct,
                      run_prefix, save_prefix,
                      zma_locs=(0,)):
    """ set the theory
    """

    rxn_info = ts_dct['rxn_info']
    ts_info = rinfo.ts_info(rxn_info)
    rct_info = rinfo.rgt_info(rxn_info, 'reacs')
    rxn_info = rinfo.sort(rxn_info)

    ts_locs = ()
    # ts_locs = (int(tsname.split('_')[-1]),)

    # high_mult = rinfo.ts_high_mult(rxn_info)
    high_mult = 2

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
    runlvl_ts_run_fs = None
    runlvl_scn_run_fs = None
    runlvl_cscn_run_fs = None
    vscnlvl_ts_run_fs = None
    vscnlvl_scn_run_fs = None
    vscnlvl_cscn_run_fs = None
    vrctst_run_fs = None

    # Initialize the necessary save filesystem
    ini_zma_save_fs = None
    runlvl_ts_save_fs = None
    runlvl_scn_save_fs = None   # above cnf filesys, for search scans
    runlvl_cscn_save_fs = None   # above cnf filesys, for search scans
    runlvl_cnf_save_fs = None
    vscnlvl_thy_save_fs = None
    vscnlvl_ts_save_fs = None
    vscnlvl_scn_save_fs = None
    vscnlvl_cscn_save_fs = None
    vrctst_save_fs = None

    if es_keyword_dct.get('inplvl', None) is not None:

        ini_method_dct = thy_dct.get(es_keyword_dct['inplvl'])
        ini_thy_info = tinfo.from_dct(ini_method_dct)
        mod_ini_thy_info = tinfo.modify_orb_label(
            ini_thy_info, ts_info)

        ini_cnf_run_fs, ini_cnf_save_fs = build_fs(
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
            vsplvl1_thy_info = tinfo.from_dct(method_dct)
            mod_vsplvl1_thy_info = tinfo.modify_orb_label(
                vsplvl1_thy_info, ts_info)
            hs_vsplvl1_thy_info = tinfo.modify_orb_label(
                vsplvl1_thy_info, hs_info)

        if es_keyword_dct.get('var_splvl2', None) is not None:

            method_dct = thy_dct.get(es_keyword_dct['var_scnlvl2'])
            vsplvl2_thy_info = tinfo.from_dct(method_dct)
            mod_vsplvl2_thy_info = tinfo.modify_orb_label(
                vsplvl2_thy_info, ts_info)
            hs_vsplvl2_thy_info = tinfo.modify_orb_label(
                vsplvl2_thy_info, hs_info)

    # Get the conformer filesys for the reactants
    rcts_cnf_fs = _reac_cnf_fs(
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
        'rcts_cnf_fs': rcts_cnf_fs
    }

    return thy_inf_dct, runfs_dct, savefs_dct
