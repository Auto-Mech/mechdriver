""" TS Finding algorithms

    Newish approach to eventually integrate
    - KBM
"""

import automol.par
# from mechanalyzer.inf import thy as tinfo
from mechanalyzer.inf import rxn as rinfo
from mechlib.amech_io import printer as ioprinter
from mechroutines.es.runner import qchem_params
from mechroutines.es.ts import _sadpt as sadpt
from mechroutines.es.ts import _irc as irc
from mechroutines.es.ts import _vrctst as vrctst
from mechroutines.es.ts import _rpath as rpath
from mechroutines.es.ts._util import thy_dcts
from mechroutines.es.ts._fs import rpath_fs


ES_TSKS = {}


# Main callable function
def findts(tsk, spc_dct, tsname, thy_dct, es_keyword_dct,
           run_prefix, save_prefix):
    """ run TS task
    """

    # Build necessary objects for transition states
    thy_inf_dct, thy_method_dct, mref_dct, runfs_dct, savefs_dct = thy_dcts(
        tsname, spc_dct[tsname], thy_dct, es_keyword_dct,
        run_prefix, save_prefix)

    # Run the task
    if 'find' in tsk:
        _findts(spc_dct, tsname,
                thy_inf_dct, thy_method_dct, mref_dct,
                es_keyword_dct,
                runfs_dct, savefs_dct)
    elif 'rpath' in tsk:
        job = tsk.split('_', 1)[1]
        rpath_tsk(job, spc_dct, tsname,
                  thy_inf_dct, thy_method_dct, mref_dct,
                  es_keyword_dct,
                  runfs_dct, savefs_dct)


# Task Functions for finding transition states
def _findts(spc_dct, tsname,
            thy_inf_dct, thy_method_dct, mref_dct,
            es_keyword_dct,
            runfs_dct, savefs_dct):
    """ Launches one of several TS finding algorithms based on
        the request of the user, or if missing, aspects about the
        reactants and spin-state of the reaction.

        :param spc_dct: mechanism species dictionary
        :type spc_dct: dict[str: dict[str: obj]]
        :param tsname: mechanism name of transition state
        :type tsname: str
        :param thy_dct: parameters of all user-defined elec. struct levels
        :type thy_dct: dict[str: dict[str:obj]]
        :param es_keyword_dct: keyword-values for electronic structure task
        :type es_keyword_dct: dict[str: obj]
        :param run_prefix: root-path to the run-filesystem
        :type run_prefix: str
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
    """

    # # Build necessary objects
    # thy_inf_dct, thy_method_dct, mref_dct, runfs_dct, savefs_dct = thy_dcts(
    #     tsname, spc_dct[tsname], thy_dct, es_keyword_dct,
    #     run_prefix, save_prefix)

    # Determine the TS finding algorithm to use
    search_method = _ts_search_method(spc_dct[tsname])

    # Calculate all required transition state information
    if search_method == 'sadpt':
        success = run_sadpt(
            spc_dct, tsname,
            thy_inf_dct, thy_method_dct, mref_dct,
            es_keyword_dct,
            runfs_dct, savefs_dct)
    elif search_method == 'pst':
        success = run_pst(
            spc_dct, tsname, savefs_dct)
    elif search_method == 'rpvtst':
        success = run_rpvtst(
            spc_dct, tsname,
            thy_inf_dct, thy_method_dct, mref_dct,
            es_keyword_dct,
            runfs_dct, savefs_dct)
    elif search_method == 'vrctst':
        success = run_vrctst(
            spc_dct, tsname,
            thy_inf_dct, thy_method_dct, mref_dct,
            es_keyword_dct,
            runfs_dct, savefs_dct)

    return success


def run_sadpt(spc_dct, tsname,
              thy_inf_dct, thy_method_dct, mref_dct,
              es_keyword_dct,
              runfs_dct, savefs_dct):
    """ Find the saddle-point for a reaction
    """

    # Check filesystem for existing zmatrixes
    run_zma, ini_zma = sadpt.read_existing_saddle_points(
        spc_dct, tsname, savefs_dct)

    if sadpt.search_required(run_zma, es_keyword_dct):
        success = sadpt.search(ini_zma, spc_dct, tsname,
                               thy_inf_dct, thy_method_dct, mref_dct,
                               es_keyword_dct,
                               runfs_dct, savefs_dct)
    else:
        success = True

    return success


def run_rpvtst(spc_dct, tsname,
               thy_inf_dct, thy_method_dct, mref_dct,
               es_keyword_dct,
               runfs_dct, savefs_dct):
    """ generate a reaction path

        need some way of putting in a different level of theory,
        maybe switch out the run_sadpt, for its compontnets.
        sadpt.search could take a level as input
    """

    # Try and first locate a saddle point
    print('First attempting to locate a saddle point')
    success = run_sadpt(spc_dct, tsname,
                        thy_inf_dct, thy_method_dct, mref_dct,
                        es_keyword_dct,
                        runfs_dct, savefs_dct)
    print('rpvtst (sadpt) success', success)

    if success:
        print('Sadpoint located. Will launch IRC from there')
    else:
        print('No saddle-point. Will launch IRC from max of potential')

    # Run rpath task to run a scan along the IRC
    # Internal logic should allow it to determine if there is sadpt or not
    # based on what happened in above function call
    job = 'scan'
    es_keyword_dct.update({'rxncoord': 'irc'})
    rpath_tsk(job, spc_dct, tsname,
              thy_inf_dct, thy_method_dct, mref_dct, es_keyword_dct,
              runfs_dct, savefs_dct)

    return success


def run_pst(spc_dct, tsname, savefs_dct):
    """ Saves all of the information required to model a transition
        state with Phase Space Theory into the SAVE filesystem.

        Currently, only a Z-Matrix and a automol.reac.Reaction object
        are saved, as all other required information to build a MESS
        input for the PST transition state are generated with kTPDriver.
        :rtype: bool
    """

    # Obtain necessary objects from various dictionaries
    zma_save_fs = savefs_dct['runlvl_ts_zma']

    ts_dct = spc_dct[tsname]
    zrxn, zma = ts_dct['zrxn'], ts_dct['zma']
    zma_locs = (ts_dct.get('zma_idx', 0),)

    # Save a Z-Matrix and Reaction object if missing
    zma_path = zma_save_fs[-1].path(zma_locs)
    if (
        not zma_save_fs[-1].file.reaction.exists(zma_locs) and
        not zma_save_fs[-1].file.zmatrix.exists(zma_locs)
    ):
        print('Saving info required for Phase Space Theory ',
              f'at path {zma_path}')

        zma_save_fs[-1].create(zma_locs)
        zma_save_fs[-1].file.reaction.write(zrxn, zma_locs)
        zma_save_fs[-1].file.zmatrix.write(zma, zma_locs)
    else:
        print('All info required for Phase Space Theory currently saved '
              f'at path {zma_path}')

    return True  # As long as the if-block passes, code should be successful


def run_vrctst(spc_dct, tsname,
               thy_inf_dct, thy_method_dct, mref_dct,
               es_keyword_dct,
               runfs_dct, savefs_dct):
    """ generate a reaction path
    """

    success = vrctst.calc_vrctst_flux(
        spc_dct[tsname],
        thy_inf_dct, thy_method_dct, mref_dct,
        es_keyword_dct,
        runfs_dct, savefs_dct)

    return success


# def run_isc(spc_dct, tsname, method_dct, es_keyword_dct,
#             runfs_dct, savefs_dct):
#     """ Search for an ISC TS by looking for
#         the minimum on the crossing seam
#     """
#     ts_dct = spc_dct[tsname]
#     # Get all of the structural information
#     ts_geo = automol.zmatrix.geometry(ts_zma['zma'])
#     ts_chg = rinfo.ts_chg(ts_dct['rxn_info'])
#     ts_mult = rinfo.ts_mult(ts_dct['rxn_info'])
#     # Set quantum chemistry parameters
#     thy_info = method_dct['runlvl']
#     prog, method, basis, orb_label = thy_info
#     if elstruct.Method.is_multireference(method):
#         cas_kwargs = multireference_calculation_parameters(
#             ref_zma, ts_info, ts_formula, high_mul,
#             rct_ichs, rct_info,
#             aspace, mod_var_scn_thy_info)
#     else:
#         cas_kwargs = None
#     _, kwargs = qchem_params(
#         method_dct, job=elstruct.Job.OPTIMIZATION)
#     # Get the run directory
#     runlvl_cnf_run_fs, _ = runfs_dct['runlvl_cnf_fs']
#     run_fs = autofile.fs.run()
#     run_path = run_fs[-1].path(['NST'])
#     # Perform the ISC TS search
#     msx_geo, hessians, flux_str = autorun.isc_flux(
#         run_dir, prog, ts_geo, ts_chg, mults,
#         method, basis, orb_label, ini_kwargs)
#     print('')
#     print('geo')
#     if msx_geo is not None:
#         print(automol.geom.string(msx_geo))
#     print('hessians')
#     if hessians is not None:
#         for hess in hessians:
#             print(hess)
#     print('flux')
#     if flux_str is not None:
#         print(flux_str)


# Task Functions for calculations involving reaction paths
def rpath_tsk(job, spc_dct, spc_name,
              thy_inf_dct, thy_method_dct, mref_dct, es_keyword_dct,
              runfs_dct, savefs_dct):
    """ run a scan over the specified torsional coordinate

        :param job:
        :type job:
        :param spc_dct:
        :type spc_dct:
        :param spc_name:
        :type spc_name:
        :param thy_dct:
        :type thy_dct:
        :param es_keyword_dct: keyword-value pairs for electronic structure tsk
        :type es_keyword_dct: dict[str:str]
        :param run_prefix: root-path to the run-filesystem
        :type run_prefix: str
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
    """

    # Get dct for specific species task is run for
    ts_dct = spc_dct[spc_name]
    ts_info = rinfo.ts_info(ts_dct['rxn_info'])

    # # Build various thy and filesystem objects
    # thy_inf_dct, thy_method_dct, mref_dct, runfs_dct, savefs_dct = thy_dcts(
    #     spc_name, spc_dct[spc_name], thy_dct, es_keyword_dct,
    #     run_prefix, save_prefix)

    mod_thy_info = thy_inf_dct['mod_runlvl']
    mod_ini_thy_info = thy_inf_dct['mod_inplvl']
    method_dct = thy_method_dct['runlvl']
    ini_method_dct = thy_method_dct['inplvl']
    mref_params = mref_dct['runlvl']

    # Build filesys objects
    scn_alg, scn_fs, cscn_fs, cnf_fs, cnf_locs = rpath_fs(
        ts_dct, spc_name,
        mod_ini_thy_info,
        es_keyword_dct,
        runfs_dct['prefix'], savefs_dct['prefix'])
    print('alg set test', scn_alg)

    # Run job
    if job == 'scan':
        if 'irc' in scn_alg:
            zmas = irc.launch_point_zmatrices(
                ts_dct, mod_thy_info, scn_alg,
                scn_fs, cnf_fs, cnf_locs)
            for zma in zmas:
                success = irc.execute_irc(
                    zma, ts_info,
                    mod_ini_thy_info, ini_method_dct,
                    scn_fs[0], scn_fs[1],
                    es_keyword_dct)
                if success:
                    break
        else:
            zma, zrxn, rclass, = ts_dct['zma'], ts_dct['zrxn'], ts_dct['class']
            _ = rpath.internal_coordinates_scan(
                zma, zrxn,
                ts_info, rclass,
                method_dct, mref_params,
                scn_fs[0], scn_fs[1],
                cscn_fs[0], cscn_fs[1],
                es_keyword_dct)
    elif job in ('energy', 'grad', 'hess'):
        # Run along the scan and calculate desired quantities
        ini_scn_run_fs, ini_scn_save_fs = scn_fs
        for locs in ini_scn_save_fs[-1].existing():
            geo = ini_scn_save_fs[-1].file.geometry.read(locs)
            script_str, kwargs = qchem_params(
                geo=geo, spc_info=ts_info)
            ini_scn_run_fs[-1].create(locs)
            ES_TSKS[job](
                None, geo, ts_info, mod_thy_info,
                ini_scn_run_fs, ini_scn_save_fs, locs,
                script_str, es_keyword_dct['overwrite'],
                **kwargs)
            ioprinter.obj('vspace')
    elif job == 'infene':
        rpath.inf_sep_ene(
            ts_dct, thy_inf_dct, mref_dct,
            savefs_dct, runfs_dct, es_keyword_dct)


# Helper functions
def _ts_search_method(ts_dct):
    """ Determine the algorithm that should be used for finding the
        transition state of a reaction.

        First the ts_dct is checked to see if an algorithm has been
        specified by the user from an input file. If an algorithm
        has not been requested, then one is automatically chosen based
        on (1) spin-state of reaction being low or high, and (2)
        whether the reactants consist of two-or-more radicals.

        :param ts_dct: species dictionary for reaction transition state
        :type ts_dct: dict[str:obj]
        :rtype: str
    """

    # Set search algorithm to one specified by the user, if specified
    _search_method = ts_dct.get('ts_search')
    if _search_method is not None:
        print(('User requested the use of TS finding algorithm ' +
               _search_method +
               '. Using this algorithm for the search'))
    else:
        print('No TS finding algorithm requested by the user.')

        # ID search algorithm if user did not specify one (wrong)
        if automol.par.isc(ts_dct['class']):
            _search_method = 'isc'
            print('Reactant and Product spins differ...')
            print('Using intersystem crossing search sceme')
        elif automol.par.has_nobarrier(ts_dct['class']):
            _search_method = 'pst'
            print()
            print('Reaction is low-spin, radical-radical reaction')
            print('Assuming reaction is barrierless...')
            print('Generating information required for phase space theory...')
        else:
            _search_method = 'sadpt'
            print()
            print('Reaction is either mol-rad or high-spin radical-radical.')
            print('Assuming reaction has saddle point on potential surface...')
            print('Finding the geometry of the saddle point...')

    return _search_method
