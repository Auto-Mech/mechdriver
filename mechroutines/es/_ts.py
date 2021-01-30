"""
  TS Finding algorithms
"""

import automol
import autofile
from mechanalyzer.inf import spc as sinfo
from mechroutines.es._routines import _sadpt as sadpt
from mechroutines.es._routines import _vrctst as vrctst
from mechroutines.es._routines import _vtst as vtst
from mechlib import filesys
from mechlib.structure import tors as torsprep


def findts(tsk, spc_dct, tsname, thy_dct, es_keyword_dct,
           run_prefix, save_prefix):
    """ New run function
    """

    # Set the TS searching algorithm to use: (1) Check dct, (2) Set by Class
    search_method = _ts_finder_match(tsk, spc_dct[tsname])

    # Build necessary objects
    method_dct, runfs_dct, savefs_dct = _set_methods(
        spc_dct[tsname], thy_dct, es_keyword_dct,
        run_prefix, save_prefix)

    # Find the transition state
    if search_method == 'sadpt':
        run_sadpt(spc_dct, tsname, es_keyword_dct,
                  method_dct, runfs_dct, savefs_dct)
    elif search_method == 'molrad_vtst':
        run_vtst(spc_dct, tsname, es_keyword_dct,
                 method_dct, runfs_dct, savefs_dct)
    elif search_method == 'radrad_vtst':
        run_vrctst(spc_dct, tsname, es_keyword_dct,
                   method_dct, runfs_dct, savefs_dct)
    elif search_method is None:
        print('No TS search algorithm was specified or able to determined')


def run_sadpt(spc_dct, tsname, es_keyword_dct,
              method_dct, runfs_dct, savefs_dct):
    """ find a transition state
    """

    # Get objects for the calculations
    ts_dct = spc_dct[tsname]
    mod_thy_info = method_dct['runlvl']

    # Find the TS
    cnf_save_fs, cnf_save_locs = savefs_dct['runlvl_cnf_fs']
    overwrite = es_keyword_dct['overwrite']
    if not cnf_save_locs:
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
        guess_zmas = sadpt.generate_guess_structure(
            ts_dct, mod_thy_info,
            runfs_dct, savefs_dct, es_keyword_dct)
        sadpt.obtain_saddle_point(
            guess_zmas, ts_dct, mod_thy_info,
            runfs_dct, savefs_dct, es_keyword_dct)


def run_vtst(spc_dct, tsname, es_keyword_dct,
             method_dct, runfs_dct, savefs_dct,
             info_dct, grid):
    """ find a transition state
    """

    # Get dct for specific species task is run for
    ts_dct = spc_dct[tsname]

    # Get info from the reactants
    ts_info = info_dct['ts_info']
    rct_info = info_dct['rct_info']
    rcts_gra = ts_dct['rcts_gra']
    if radrad:
        high_mul = ts_dct['high_mult']
        hs_info = info_dct['hs_info']
        rct_ichs = [spc_dct[rct]['inchi'] for rct in ts_dct['reacs']]

    # Set information from the transition state
    ini_zma = ts_dct['zma']
    frm_bnd_keys = ts_dct['frm_bnd_keys']
    if radrad:
        ts_formula = automol.geom.formula(automol.zmatrix.geometry(ini_zma))
        active_space = ts_dct['active_space']

    # Get reaction coordinates
    frm_name = automol.zmatrix.bond_key_from_idxs(ini_zma, frm_bnd_keys)

    # Get es options
    overwrite = es_keyword_dct['overwrite']
    update_guess = False  # check
    if radrad:
        pot_thresh = es_keyword_dct['pot_thresh']

    # Grid
    print('newts class', ts_dct['class'])
    print('newts grid', grid)
    [grid1, grid2] = grid

    # Get method stuff
    mod_thy_info = method_dct['mod_runlvl']
    mod_var_scn_thy_info = method_dct['mod_var_scnlvl']
    mod_var_sp1_thy_info = method_dct['mod_var_splvl1']
    var_sp1_thy_info = method_dct['var_splvl2']
    var_sp2_thy_info = method_dct['var_splvl2']
    hs_var_sp1_thy_info = method_dct['hs_var_splvl1']
    hs_var_sp2_thy_info = method_dct['hs_var_splvl2']

    # Get the filesys stuff
    var_scn_save_fs = savefs_dct['vscnlvl_scn_fs']
    var_scn_run_fs = runfs_dct['vscnlvl_scn_fs']
    rcts_cnf_fs = savefs_dct['rcts_cnf_fs']
    vscnlvl_thy_save_fs = savefs_dct['vscnlvl_thy_fs']
    vscnlvl_ts_save_fs = savefs_dct['vscnlvl_ts_fs']

    # Check for the saddle point
    if not cnf_save_locs:
        print('No energy found in save filesys. Running energy...')
        _run = True
    elif overwrite:
        print('User specified to overwrite transition state with new run...')
        _run_scan = True
    else:
        _run = False

    # Find the TS (check the path)
    if not cnf_save_locs:
        print('No energy found in save filesys. Running energy...')
        _run = True
    elif overwrite:
        print('User specified to overwrite transition state with new run...')
        _run_scan = True
    else:
        _run = False

    # Run the scan if needed
    if _run:
        if radrad:
            vtst.radrad_scan(
                ts_zma, ts_info,
                ts_formula, high_mul, active_space,
                rct_info, rct_ichs,
                grid1, grid2, coord_name,
                mod_var_scn_thy_info,
                scn_run_fs, scn_save_fs,
                overwrite, update_guess,
                constraint_dct=None,
                zma_locs=(0,))
        else:
            vtst.molrad_scan(ts_zma, ts_info,
                             rct_info, rcts_cnf_fs, rcts_gra,
                             grid1, grid2, coord_name, frm_bnd_keys,
                             thy_info, vsp1_thy_info,
                             thy_save_fs,
                             ts_save_fs,
                             scn_run_fs, scn_save_fs,
                             overwrite, update_guess, retryfail,
                             zma_locs=(0,))

        sadpt_zma = rxngrid.vtst_max(
            list(grid1)+list(grid2), coord_name, scn_save_fs,
            mod_var_scn_thy_info, constraint_dct,
            ethresh=pot_thresh)

    if sadpt_zma is None:

        if radrad:
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

        _save_traj(ts_zma, frm_bnd_keys, rcts_gra,
                   vscnlvl_ts_save_fs, zma_locs=zma_locs)

        _vtst_hess_ene(ts_info, coord_name,
                       mod_var_scn_thy_info, mod_var_sp1_thy_info,
                       scn_save_fs, scn_run_fs,
                       overwrite, **cas_kwargs)
    else:
        obtain_saddle_point()


def run_vrctst(spc_dct, tsname, es_keyword_dct,
               method_dct, runfs_dct, savefs_dct,
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

    # Get method stuff
    mod_var_scn_thy_info = method_dct['var_scnlvl']
    mod_var_sp1_thy_info = method_dct['var_splvl1']
    mod_var_sp2_thy_info = method_dct['var_splvl2']
    hs_var_sp1_thy_info = method_dct['hs_var_splvl1']
    hs_var_sp2_thy_info = method_dct['hs_var_splvl2']

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
        ini_method = [ts_dct['ts_search']]
        print('Running search algorithm according to {},'.format(
            ts_dct['ts_search']),
              'as requested by the user')
    else:
        ini_method = None
        print('No search algorithm requested')
    print()

    # ID search algorithm if user did not specify one (wrong)
    if ini_method is None:
        if _nobarrier(ts_dct):
            ini_method = ['radrad_vtst', 'vrctst']
            print('Reaction is low-spin, radical-radical addition/abstraction')
            print('Assuming reaction is barrierless...')
            print('Finding a transition state according to either vtst or '
                  'vrctst, depending on the current task')
        else:
            ini_method = ['sadpt']
            print('Assuming reaction has saddle point on potential surface...')
            print('Use species.dat to specify VTST search for mol-rad rxn...')
            print('Finding the geometry of the saddle point...')

    # Print message if no algorithm found
    if ini_method is None:
        print('No TS search algorithm was specified or able to determined')

    # Set return for ts searching algorithm if there is one
    if ini_method == 'pst':
        print('Phase Space Theory Used, No ES calculations are needed')
    if tsk in ini_method:
        print('Search algorithm matches task')
        search_method = tsk
    else:
        print('Algorithm does not match task')
        search_method = None

    # Refine ret vtst method if that is what is being used
    if search_method == 'vtst':
        search_method = 'radrad_vtst' if _radrad(ts_dct) else 'molrad_vtst'

    return search_method


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
def _set_grid(ts_search, ts_dct):
    """ Set the TS grid
    """

    if ts_search in ('molrad_vtst', 'radrad_vtst', 'vrctst'):
        if 'rad' in ts_dct['class']:
            grid = ts_dct['grid']
        else:
            grid = ts_dct['var_grid']
    else:
        grid = ts_dct['grid']

    return grid


def _set_methods(ts_dct, thy_dct, es_keyword_dct, info_dct,
                 run_prefix, save_prefix,
                 zma_locs=(0,)):
    """ set the theory
    """

    ts_info = info_dct['ts_info']
    rxn_info = info_dct['rxn_info']
    rct_info = info_dct['rct_info']

    # Get the name
    # ini_zma = ts_dct['zma']
    # frm_bnd_keys = ts_dct['frm_bnd_keys']
    # brk_bnd_keys = ts_dct['brk_bnd_keys']
    # frm_name = automol.zmatrix.bond_key_from_idxs(ini_zma, frm_bnd_keys)
    # brk_name = automol.zmatrix.bond_key_from_idxs(ini_zma, brk_bnd_keys)

    # Set the hs info
    hs_info = (ts_info[0], ts_info[1], ts_dct['high_mult'])

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
    vscnlvl_ts_run_fs = None
    vscnlvl_scn_run_fs = None
    vscnlvl_cscn_run_fs = None
    vrctst_run_fs = None

    # Initialize the necessary save filesystem
    ini_zma_save_fs = None
    runlvl_ts_save_fs = None
    runlvl_scn_save_fs = None   # above cnf filesys, for search scans
    runlvl_cnf_save_fs = None
    vscnlvl_thy_save_fs = None
    vscnlvl_ts_save_fs = None
    vscnlvl_scn_save_fs = None
    vscnlvl_cscn_save_fs = None
    vrctst_save_fs = None

    if es_keyword_dct.get('inplvl', None) is not None:

        ini_thy_info = filesys.inf.get_es_info(
            es_keyword_dct['inplvl'], thy_dct)
        mod_ini_thy_info = filesys.inf.modify_orb_restrict(
            ts_info, ini_thy_info)

        ini_thy_save_fs = filesys.build.rxn_thy_fs_from_root(
            save_prefix, rxn_info, mod_ini_thy_info)
        ini_ts_save_fs = filesys.build.ts_fs_from_thy(
            ini_thy_save_fs[1])
        ini_cnf_save_fs = autofile.fs.conformer(ini_ts_save_fs[1])
        ini_min_cnf_locs, _ = filesys.mincnf.min_energy_conformer_locators(
            ini_cnf_save_fs, mod_ini_thy_info)

        if ini_min_cnf_locs:
            ini_zma_save_fs = autofile.fs.manager(
                ini_cnf_save_fs[-1].path(ini_min_cnf_locs), 'ZMATRIX')

    if es_keyword_dct.get('runlvl', None) is not None:

        thy_info = filesys.inf.get_es_info(
            es_keyword_dct['runlvl'], thy_dct)
        mod_thy_info = filesys.inf.modify_orb_restrict(
            ts_info, thy_info)
        hs_thy_info = filesys.inf.modify_orb_restrict(
            hs_info, thy_info)

        runlvl_thy_run_fs = filesys.build.rxn_thy_fs_from_root(
            run_prefix, rxn_info, mod_thy_info)
        runlvl_thy_save_fs = filesys.build.rxn_thy_fs_from_root(
            save_prefix, rxn_info, mod_thy_info)

        runlvl_ts_save_fs = filesys.build.ts_fs_from_thy(
            runlvl_thy_save_fs[1])
        runlvl_ts_run_fs = filesys.build.ts_fs_from_thy(
            runlvl_thy_run_fs[1])

        runlvl_fs = autofile.fs.conformer(ini_ts_save_fs[1])
        ini_min_cnf_locs, _ = filesys.mincnf.min_energy_conformer_locators(
            runlvl_fs, mod_thy_info)
        runlvl_cnf_save_fs = (runlvl_fs, ini_min_cnf_locs)

        _, runlvl_zma_run_path = filesys.build.zma_fs_from_prefix(
            runlvl_thy_run_fs[1], zma_idxs=zma_locs)
        _, runlvl_zma_save_path = filesys.build.zma_fs_from_prefix(
            runlvl_thy_save_fs[1], zma_idxs=zma_locs)
        runlvl_scn_run_fs = filesys.build.scn_fs_from_cnf(
            runlvl_zma_run_path, constraint_dct=None)
        runlvl_scn_save_fs = filesys.build.scn_fs_from_cnf(
            runlvl_zma_save_path, constraint_dct=None)

    if es_keyword_dct.get('var_scnlvl', None) is not None:

        vscnlvl_thy_info = filesys.inf.get_es_info(
            es_keyword_dct['var_scnlvl'], thy_dct)
        mod_vscnlvl_thy_info = filesys.inf.modify_orb_restrict(
            ts_info, vscnlvl_thy_info)
        hs_vscnlvl_thy_info = filesys.inf.modify_orb_restrict(
            hs_info, vscnlvl_thy_info)

        vscnlvl_thy_run_fs = filesys.build.rxn_thy_fs_from_root(
            run_prefix, rxn_info, mod_vscnlvl_thy_info)
        vscnlvl_thy_save_fs = filesys.build.rxn_thy_fs_from_root(
            save_prefix, rxn_info, mod_vscnlvl_thy_info)

        vscnlvl_ts_save_fs = filesys.build.ts_fs_from_thy(
            vscnlvl_thy_save_fs[1])
        vscnlvl_ts_run_fs = filesys.build.ts_fs_from_thy(
            vscnlvl_thy_run_fs[1])

        # put the scan filesys
        _, vscnlvl_zma_run_path = filesys.build.zma_fs_from_prefix(
            vscnlvl_thy_run_fs[1], zma_idxs=zma_locs)
        _, vscnlvl_zma_save_path = filesys.build.zma_fs_from_prefix(
            vscnlvl_thy_save_fs[1], zma_idxs=zma_locs)

        vscnlvl_scn_run_fs = autofile.fs.scan(vscnlvl_zma_run_path)
        vscnlvl_scn_save_fs = autofile.fs.scan(vscnlvl_zma_save_path)
        vscnlvl_cscn_run_fs = autofile.fs.cscan(vscnlvl_zma_run_path)
        vscnlvl_cscn_save_fs = autofile.fs.cscan(vscnlvl_zma_save_path)

        vrctst_save_fs = filesys.build.vrc_fs_from_thy(
            vscnlvl_ts_save_fs[1])
        vrctst_run_fs = filesys.build.vrc_fs_from_thy(
            vscnlvl_ts_run_fs[1])

        if es_keyword_dct.get('var_splvl1', None) is not None:

            vsp1lvl_thy_info = filesys.inf.get_es_info(
                es_keyword_dct['var_splvl1'], thy_dct)
            mod_vsp1lvl_thy_info = filesys.inf.modify_orb_restrict(
                ts_info, vsp1lvl_thy_info)
            hs_vsp1lvl_thy_info = filesys.inf.modify_orb_restrict(
                hs_info, vsp1lvl_thy_info)

            # vscnlvl_scn_run_fs = filesys.build.scn_fs_from_cnf(
            #     vscnlvl_zma_run_path, constraint_dct=None)
            # vscnlvl_scn_save_fs = filesys.build.scn_fs_from_cnf(
            #    vscnlvl_zma_save_path, constraint_dct=None)
            #     var_sp1_thy_run_fs = filesys.build.rxn_thy_fs_from_root(
            #         run_prefix, rxn_info, mod_var_sp1_thy_info)
            #     var_sp1_thy_save_fs = filesys.build.rxn_thy_fs_from_root(
            #         save_prefix, rxn_info, mod_var_sp1_thy_info)

        if es_keyword_dct.get('var_splvl2', None) is not None:

            vsp2lvl_thy_info = filesys.inf.get_es_info(
                es_keyword_dct['var_splvl2'], thy_dct)
            mod_vsp2lvl_thy_info = filesys.inf.modify_orb_restrict(
                ts_info, vsp2lvl_thy_info)
            hs_vsp2lvl_thy_info = filesys.inf.modify_orb_restrict(
                hs_info, vsp2lvl_thy_info)

            #     var_sp2_thy_run_fs = filesys.build.rxn_thy_fs_from_root(
            #         run_prefix, rxn_info, mod_var_sp2_thy_info)
            #     var_sp2_thy_save_fs = filesys.build.rxn_thy_fs_from_root(
            #         save_prefix, rxn_info, mod_var_sp2_thy_info)

    # Get the conformer filesys for the reactants
    reac_cnf_fs = _reac_cnf_fs(
        rct_info, thy_dct, es_keyword_dct, run_prefix, save_prefix)

    # Build the dictionaries for the return
    method_dct = {
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
        'runlvl_ts_fs': runlvl_ts_run_fs,
        'runlvl_scn_fs': runlvl_scn_run_fs,
        'vscnlvl_ts_fs': vscnlvl_ts_run_fs,
        'vscnlvl_scn_fs': vscnlvl_scn_run_fs,
        'vscnlvl_cscn_fs': vscnlvl_cscn_run_fs,
        'vrctst_fs': vrctst_run_fs,
    }

    savefs_dct = {
        'inilvl_zma_fs': ini_zma_save_fs,
        'inilvl_ts_fs': ini_ts_save_fs,
        'runlvl_thy_fs': runlvl_thy_save_fs,
        'runlvl_ts_fs': runlvl_ts_save_fs,
        'runlvl_scn_fs': runlvl_scn_save_fs,
        'runlvl_cnf_fs': runlvl_cnf_save_fs,
        'vscnlvl_thy_fs': vscnlvl_thy_save_fs,
        'vscnlvl_ts_fs': vscnlvl_ts_save_fs,
        'vscnlvl_scn_fs': vscnlvl_scn_save_fs,
        'vscnlvl_cscn_fs': vscnlvl_cscn_save_fs,
        'vrctst_fs': vrctst_save_fs,
        'rcts_cnf_fs': reac_cnf_fs
    }

    return method_dct, runfs_dct, savefs_dct


def _reac_cnf_fs(rct_info, thy_dct, es_keyword_dct, run_prefix, save_prefix):
    """ set reactant method stuff
    """

    ini_thy_info = filesys.inf.get_es_info(
        es_keyword_dct['inplvl'], thy_dct)

    rct_cnf_fs = ()

    for rinfo in rct_info:

        mod_ini_thy_info = filesys.inf.modify_orb_restrict(
            rinfo, ini_thy_info)

        # Build filesys for ini thy info
        _, ini_thy_run_path = filesys.build.spc_thy_fs_from_root(
            run_prefix, rinfo, mod_ini_thy_info)
        _, ini_thy_save_path = filesys.build.spc_thy_fs_from_root(
            save_prefix, rinfo, mod_ini_thy_info)

        # Build conformer filesys
        ini_cnf_run_fs = autofile.fs.conformer(ini_thy_run_path)
        ini_cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)
        ini_loc_info = filesys.mincnf.min_energy_conformer_locators(
            ini_cnf_save_fs, mod_ini_thy_info)
        ini_min_cnf_locs, ini_cnf_path = ini_loc_info

        ini_cnf_run_fs[-1].create(ini_min_cnf_locs)

        rct_cnf_fs += ((ini_cnf_run_fs, ini_cnf_save_fs,
                        ini_min_cnf_locs, ini_cnf_path),)

    return rct_cnf_fs