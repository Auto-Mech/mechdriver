"""
  TS Finding algorithms
"""

import automol
import autofile
from routines.es._routines import _sadpt as sadpt
from routines.es._routines import _vrctst as vrctst
from routines.es._routines import _vtst as vtst
from routines.es._routines import _wfn as wfn
from routines.es import runner as es_runner
from lib import filesys


def run(tsk, spc_dct, tsname, thy_dct, es_keyword_dct,
        run_prefix, save_prefix):
    """ New run function
    """

    # Set the TS searching algorithm to use: (1) Check dct, (2) Set by Class
    search_method = _ts_finder_match(tsk, spc_dct[tsname], tsname)

    # Build necessary objects
    method_dct, runfs_dct, savefs_dct = _set_methods(
        thy_dct, es_keyword_dct)
    info_dct = _set_info(spc_dct, tsname)
    grid = _set_grid(search_method, spc_dct[tsname])

    # Find the transition state
    if search_method == 'sadpt':
        run_sadpt(spc_dct, tsname, es_keyword_dct,
                  method_dct, runfs_dct, savefs_dct,
                  info_dct, grid, run_prefix, save_prefix)
    elif search_method == 'molrad_vtst':
        run_molrad_vtst()
    elif search_method == 'radrad_vtst':
        run_radrad_vtst()
    elif search_method == 'vrctst':
        run_vrctst()
    elif search_method is None:
        print('No TS search algorithm was specified or able to determined')


def run_sadpt(spc_dct, tsname, es_keyword_dct,
              method_dct, runfs_dct, savefs_dct,
              info_dct, grid, run_prefix, save_prefix):
    """ find a transition state
    """

    # Get es options
    overwrite = es_keyword_dct['overwrite']
    update_guess = False  # check

    # Get dct for specific species task is run for
    ts_dct = spc_dct[tsname]

    # Build inf objects for the rxn and ts
    ts_info = info_dct['ts_info']

    # Set various TS information using the dictionary
    ini_zma = ts_dct['zma']
    typ = ts_dct['class']
    frm_bnd_keys = ts_dct['frm_bnd_keys']
    brk_bnd_keys = ts_dct['brk_bnd_keys']
    rcts_gra = ts_dct['rcts_gra']

    # Get reaction coordinates
    frm_name = automol.zmatrix.bond_key_from_idxs(ini_zma, frm_bnd_keys)
    brk_name = automol.zmatrix.bond_key_from_idxs(ini_zma, brk_bnd_keys)

    # Get method stuff
    mod_ini_thy_info = method_dct['inplvl']
    mod_thy_info = method_dct['runlvl']

    # Get filesys stuff
    _, ts_run_path = runfs_dct['ts_fs']
    scn_run_fs, _ = runfs_dct['scn_fs']

    thy_save_fs, _ = savefs_dct['thy_fs']
    cnf_save_fs, cnf_save_locs = savefs_dct['cnf_fs']
    _, ini_ts_save_path = savefs_dct['ini_ts_fs']
    ts_save_fs, ts_save_path = savefs_dct['ts_fs']
    scn_save_fs, _ = savefs_dct['scn_fs']

    run_fs = autofile.fs.run(ts_run_path)

    # Find the TS
    if cnf_save_locs and not overwrite:

        print('TS found and saved previously in ',
              cnf_save_fs[-1].path(cnf_save_locs))
    else:

        print('No transition state found in filesys',
              'at {} level...'.format(es_keyword_dct['runlvl']),
              'Proceeding to find it...')
        script_str, opt_script_str, _, opt_kwargs = es_runner.qchem_params(
            *mod_thy_info[0:2])
        sadpt_transition_state(
            ini_zma, ts_info,
            mod_ini_thy_info, mod_thy_info,
            thy_save_fs,
            ini_ts_save_path,
            cnf_save_fs,
            scn_save_fs, scn_run_fs,
            ts_save_fs, ts_save_path, run_fs,
            typ, grid, update_guess,
            frm_name, brk_name,
            frm_bnd_keys, brk_bnd_keys, rcts_gra,
            opt_script_str, script_str, overwrite,
            es_keyword_dct, **opt_kwargs
        )


# SADPT FINDER FUNCTIONS
def sadpt_transition_state(
        ini_zma, ts_info,
        mod_ini_thy_info, mod_thy_info,
        thy_save_fs,
        ini_ts_save_path,
        cnf_save_fs,
        scn_save_fs, scn_run_fs,
        ts_save_fs, ts_save_path, run_fs,
        typ, grid, update_guess,
        dist_name, brk_name,
        frm_bnd_keys, brk_bnd_keys, rcts_gra,
        opt_script_str, script_str, overwrite,
        es_keyword_dct, **opt_kwargs):
    """ Find a sadddle point
    """

    # Check filesystem for input level of theory
    print('\nSearching save filesys for guess Z-Matrix calculated',
          'at {} level...'.format(es_keyword_dct['inplvl']))
    guess_zmas = sadpt.check_filesys_for_guess(
        ini_ts_save_path, mod_ini_thy_info)

    # If no guess zma, run a TS searching algorithm
    if not guess_zmas:
        print(' - No Z-Matrix in found in save filesys.')
        print('\nRunning scan to generate guess Z-Matrix for opt...')
        guess_zmas = sadpt.scan_for_guess(
            typ, grid, dist_name, brk_name, ini_zma, ts_info,
            mod_thy_info, thy_save_fs,
            scn_run_fs, scn_save_fs, opt_script_str,
            overwrite, update_guess, **opt_kwargs)

    # Optimize the saddle point
    print('\nOptimizing guess Z-Matrix obtained from scan or filesys...')
    opt_ret = sadpt.optimize_saddle_point(
        guess_zmas, ts_info, mod_thy_info,
        run_fs, opt_script_str, overwrite, **opt_kwargs)

    # Calculate the Hessian for the optimized structure
    if opt_ret is not None:
        print('\nCalculating Hessian for the optimized geometry...')
        hess_ret, freqs, imags = sadpt.saddle_point_hessian(
            opt_ret, ts_info, mod_thy_info,
            run_fs, script_str, overwrite, **opt_kwargs)

        # Assess saddle point, save it if viable
        print('Assessing the saddle point...')
        saddle = sadpt.saddle_point_checker(imags)
        if saddle:
            sadpt.save_saddle_point(
                opt_ret, hess_ret, freqs, imags,
                mod_thy_info,
                cnf_save_fs,
                ts_save_fs, ts_save_path,
                frm_bnd_keys, brk_bnd_keys, rcts_gra,
                zma_locs=[0])
    else:
        print('\n TS optimization failed. No geom to check and save.')


def run_molrad_vtst(spc_dct, tsname, es_keyword_dct,
                    method_dct, runfs_dct, savefs_dct,
                    info_dct, grid, run_prefix, save_prefix):
    """ find a transition state
    """

    # Get dct for specific species task is run for
    ts_dct = spc_dct[tsname]

    # Build inf objects for the rxn and ts
    ts_info = info_dct['ts_info']
    rct1_info = info_dct['rct1_info']
    rct2_info = info_dct['rct2_info']

    # Set various TS information using the dictionary
    ini_zma = ts_dct['zma']
    frm_bnd_keys = ts_dct['frm_bnd_keys']

    # Get reaction coordinates
    frm_name = automol.zmatrix.bond_key_from_idxs(ini_zma, frm_bnd_keys)

    # Get es options
    overwrite = es_keyword_dct['overwrite']
    update_guess = False  # check

    # Make grid
    [grid1, grid2] = grid

    # Get method stuff
    mod_ini_thy_info = method_dct['inplvl']
    mod_thy_info = method_dct['runlvl']
    mod_var_sp1_thy_info = method_dct['var_splvl1']

    # Get filesys stuff
    scn_save_fs, _ = savefs_dct['scn_fs']
    scn_run_fs, _ = runfs_dct['scn_fs']

    # Run single reference mol-rad VTST Search
    vtst.molrad_scan(
        ini_zma, ts_info,
        rct1_info, rct2_info,
        grid1, grid2, frm_name,
        mod_thy_info, mod_ini_thy_info,
        mod_var_sp1_thy_info,
        scn_run_fs, scn_save_fs,
        run_prefix, save_prefix,
        overwrite, update_guess
    )


def run_radrad_vtst(spc_dct, tsname, es_keyword_dct,
                    method_dct, runfs_dct, savefs_dct,
                    info_dct, grid, run_prefix, save_prefix):
    """ find a transition state
    """

    switch = False
    s = vrc_dct

    # Set information from the transition state
    ts_formula = automol.geom.formula(automol.zmatrix.geometry(ts_zma))
    [grid1, grid2] = grid

    # Get info from the reactants
    # rct_zmas = ts_dct['rct_zmas']
    rcts = ts_dct['reacs']
    high_mul = ts_dct['high_mult']
    rct1_info = info_dct['rct1_info']
    rct2_info = info_dct['rct2_info']

    # Set the active space
    num_act_orb, num_act_elc = wfn.active_space(
        ts_dct, spc_dct, ts_dct['high_mult'])

    vtst.radrad_scan(
        ts_zma, ts_info, ts_formula, high_mul,
        rct1_info, rct2_info,
        grid1, grid2, dist_name,
        num_act_orb, num_act_elc,
        mod_var_scn_thy_info,
        mod_var_sp1_thy_info, mod_var_sp2_thy_info,
        hs_var_scn_thy_info,
        hs_var_sp1_thy_info,
        hs_var_sp2_thy_info,
        mod_ini_thy_info, mod_thy_info,
        scn_run_fs, scn_save_fs,
        run_prefix, save_prefix,
        overwrite, update_guess,
        **opt_kwargs)


def run_vrctst(spc_dct, tsname, es_keyword_dct,
               method_dct, runfs_dct, savefs_dct,
               info_dct, grid, run_prefix, save_prefix):
    """ find a transition state
    """

    switch = False
    _ = vrc_dct

    # Set information from the transition state
    ts_formula = automol.geom.formula(automol.zmatrix.geometry(ts_zma))
    [grid1, grid2] = grid

    # Get info from the reactants
    high_mul = ts_dct['high_mult']
    rct1_info = info_dct['rct1_info']
    rct2_info = info_dct['rct2_info']

    # Set the active space
    num_act_orb, num_act_elc = wfn.active_space(
        ts_dct, spc_dct, ts_dct['high_mult'])

    print('Beginning Calculations for VRC-TST Treatments')
    vrctst.calc_vrctst_flux(
        ts_zma, ts_formula, ts_info, ts_dct, spc_dct,
        high_mul, grid1, grid2, dist_name,
        mod_var_scn_thy_info,
        mod_var_scn_thy_info,
        mod_var_sp1_thy_info,
        mod_ini_thy_info, mod_thy_info,
        thy_run_path, thy_save_path,
        overwrite, update_guess,
        run_prefix, save_prefix,
        vrc_dct,
        corr_pot=True)


# SET THE SEARCHING ALGORITHM
def _ts_finder_match(tsk, ts_dct, tsname):
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
        # Add the tsname
        print('requested by the user for {}'.format(tsname))
    else:
        ini_method = None
        print('No search algorithm requested')
    print()

    # ID search algorithm if user did not specify one (wrong)
    if ini_method is None:
        if _nobarrier(ts_dct):
            ini_method = ['vtst', 'vrctst']
            print('Reaction is low-spin, radical-radical addition/abstraction')
            print('Assuming reaction is barrierless...')
            print('Finding a transition state according to either vtst or'
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
        search_method = ini_method
    elif rad_rad_ts.lower() == 'pst':
        print('Phase Space Theory Used, No ES calculations are needed')
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
    print('cla', ts_dct['class'])
    radrad = _radrad(ts_dct)
    low_spin = bool('low' in ts_dct['class'])
    return radrad and low_spin


def _radrad(ts_dct):
    return bool('radical radical' in ts_dct['class'])


# SET OPTIONS FOR THE TRANSITION STATE
def _set_grid(ts_search, ts_dct):
    """ Set the TS grid
    """

    if ts_search in ('vtst', 'vrctst'):
        if 'rad' in ts_dct['class']:
            grid = ts_dct['grid']
        else:
            grid = ts_dct['var_grid']
    else:
        grid = ts_dct['grid']

    return grid


def _set_info(spc_dct, tsname):
    """ Build info objects
    """

    # Get needed objs from ts dict
    ts_dct = spc_dct[tsname]
    chg = ts_dct['charge']
    mult, high_mult = ts_dct['mult'], ts_dct['high_mult']
    reacs, prods = ts_dct['reacs'], ts_dct['prods']

    # Build dct holding all the info objects
    info_dct = {
        'ts_info': ('', chg, mult),
        'hs_info': ('', chg, high_mult),
        'rxn_info': filesys.inf.rxn_info(reacs, prods, spc_dct),
        'rct1_info': filesys.get_spc_info(reacs[0]),
        'rct2_info': filesys.get_spc_info(reacs[1])
    }

    return info_dct


def _set_methods(ts_dct, thy_dct, es_keyword_dct,
                 ts_info, rxn_info,
                 run_prefix, save_prefix):
    """ set the theory
    """

    # Set the hs info
    hs_info = (ts_info[0], ts_info[1], ts_dct['high_mult'])

    # Build all of the theory objects
    mod_ini_thy_info = None
    mod_thy_info = None
    mod_var_scn_thy_info = None
    mod_var_sp1_thy_info = None
    mod_var_sp2_thy_info = None
    hs_var_scn_thy_info = None
    hs_var_sp1_thy_info = None
    hs_var_sp2_thy_info = None

    if 'inplvl' in es_keyword_dct:
        ini_thy_info = filesys.inf.get_es_info(
            es_keyword_dct['inplvl'], thy_dct)
        mod_ini_thy_info = filesys.inf.modify_orb_restrict(
            ts_info, ini_thy_info)

    if 'runlvl' in es_keyword_dct:
        thy_info = filesys.inf.get_es_info(
            es_keyword_dct['runlvl'], thy_dct)
        mod_thy_info = filesys.inf.modify_orb_restrict(
            ts_info, thy_info)

    if 'var_scnlvl' in es_keyword_dct:
        var_scn_thy_info = filesys.inf.get_es_info(
            es_keyword_dct['var_scnlvl'], thy_dct)
        mod_var_scn_thy_info = filesys.inf.modify_orb_restrict(
            ts_info, var_scn_thy_info)
        hs_var_scn_thy_info = filesys.inf.modify_orb_restrict(
            hs_info, var_scn_thy_info)

    if 'var_splvl1' in es_keyword_dct:
        var_sp1_thy_info = filesys.inf.get_es_info(
            es_keyword_dct['var_splvl1'], thy_dct)
        mod_var_sp1_thy_info = filesys.inf.modify_orb_restrict(
            ts_info, var_sp1_thy_info)
        hs_var_sp1_thy_info = filesys.inf.modify_orb_restrict(
            hs_info, var_sp1_thy_info)

    if 'var_splvl2' in es_keyword_dct:
        var_sp2_thy_info = filesys.inf.get_es_info(
            es_keyword_dct['var_splvl2'], thy_dct)
        mod_var_sp2_thy_info = filesys.inf.modify_orb_restrict(
            ts_info, var_sp2_thy_info)
        hs_var_sp2_thy_info = filesys.inf.modify_orb_restrict(
            hs_info, var_sp2_thy_info)

    # Build the filesys objects for the ini thy lvl
    # Only need save to see if guess zmat can be found
    ini_thy_save_fs = filesys.build.rxn_thy_fs_from_root(
        save_prefix, rxn_info, mod_ini_thy_info)

    ini_ts_save_fs = filesys.build.ts_fs_from_thy(
        ini_thy_save_fs[1])

    # Build the filesys objects for the thy lvl
    thy_run_fs = filesys.build.rxn_thy_fs_from_root(
        run_prefix, rxn_info, mod_thy_info)
    thy_save_fs = filesys.build.rxn_thy_fs_from_root(
        save_prefix, rxn_info, mod_thy_info)

    # Build the filesys objects for variational treatments
    if mod_var_scn_thy_info is not None:
        var_scn_thy_run_fs = filesys.build.rxn_thy_fs_from_root(
            run_prefix, rxn_info, mod_var_scn_thy_info)
        var_scn_thy_save_fs = filesys.build.rxn_thy_fs_from_root(
            save_prefix, rxn_info, mod_var_scn_thy_info)

    # Set up TS filesystem objects for sadpts
    ts_save_fs = filesys.build.ts_fs_from_thy(thy_save_fs[1])
    ts_run_fs = filesys.build.ts_fs_from_thy(thy_run_fs[1])
    cnf_save_fs = filesys.build.cnf_fs_from_prefix(
        ts_save_fs[1], mod_thy_info, cnf='min')

    # Set up TS filesystem objects for multiref objects
    mref_ts_save_fs = filesys.build.ts_fs_from_thy(var_scn_thy_save_fs[1])
    mref_ts_run_fs = filesys.build.ts_fs_from_thy(var_scn_thy_run_fs[1])

    vrctst_save_fs = filesys.build.vrc_fs_from_thy(mref_ts_save_fs[1])
    vrctst_run_fs = filesys.build.vrc_fs_from_thy(mref_ts_run_fs[1])

    # SP filesys may be better to build in function
    # if mod_var_sp1_thy_info is not None:
    #     var_sp1_thy_run_fs = filesys.build.rxn_thy_fs_from_root(
    #         run_prefix, rxn_info, mod_var_sp1_thy_info)
    #     var_sp2_thy_save_fs = filesys.build.rxn_thy_fs_from_root(
    #         save_prefix, rxn_info, mod_var_sp1_thy_info)

    # if mod_var_sp2_thy_info is not None:
    #     var_sp2_thy_run_fs = filesys.build.rxn_thy_fs_from_root(
    #         run_prefix, rxn_info, mod_var_sp2_thy_info)
    #     var_sp2_thy_save_fs = filesys.build.rxn_thy_fs_from_root(
    #         save_prefix, rxn_info, mod_var_sp2_thy_info)

    # Set up the scan filesys for single and multiref
    _, zma_run_path = filesys.build.zma_fs_from_prefix(
        thy_run_path, zma_idxs=[0])
    _, zma_save_path = filesys.build.zma_fs_from_prefix(
        thy_save_path, zma_idxs=[0])
    sadpt_scn_run_fs = filesys.build.scn_fs_from_cnf(
        zma_run_path, constraint_dct=None)
    sadpt_scn_save_fs = filesys.build.scn_fs_from_cnf(
        zma_save_path, constraint_dct=None)

    # fix
    _, zma_run_path = filesys.build.zma_fs_from_prefix(
        thy_run_path, zma_idxs=[0])
    _, zma_save_path = filesys.build.zma_fs_from_prefix(
        thy_save_path, zma_idxs=[0])
    var_scn_run_fs = autofile.fs.scan(var_scn_thy_run_fs[1])
    var_scn_save_fs = autofile.fs.scan(var_scn_thy_save_fs[1])

    # Build the dictionaries for the return
    method_dct = {
        'inplvl': mod_ini_thy_info,
        'runlvl': mod_thy_info,
        'var_scnlvl': mod_var_scn_thy_info,
        'hs_var_scnlvl': hs_var_scn_thy_info,
        'var_splvl1': mod_var_sp1_thy_info,
        'hs_var_splvl1': hs_var_sp1_thy_info,
        'var_splvl2': mod_var_sp2_thy_info,
        'hs_var_splvl2': hs_var_sp2_thy_info
    }

    runfs_dct = {
        'thy_fs': thy_run_fs,
        'var_scn_thy_fs': var_scn_thy_run_fs,
        'ts_fs': ts_run_fs,
        'mref_ts_fs': mref_ts_run_fs,
        'vrctst_fs': vrctst_run_fs
    }

    savefs_dct = {
        'ini_thy_fs': ini_thy_save_fs,
        'ini_ts_fs': ini_ts_save_fs,
        'thy_fs': thy_save_fs,
        'var_scn_thy_fs': var_scn_thy_save_fs,
        'ts_fs': ts_save_fs,
        'cnf_save_fs': cnf_save_fs,
        'mref_ts_fs': mref_ts_save_fs,
        'vrctst_fs': vrctst_save_fs
    }

    return method_dct, runfs_dct, savefs_dct
