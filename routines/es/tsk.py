""" eletronic structure routines modules
"""

# electronic structure routines
import sys
import autofile
import automol
from routines.es import conformer
from routines.es import geom
from routines.es import scan
from routines.es import sp
from routines.es import tau
from routines.es import ts
from lib.runner import par as runpar
from lib.filesystem import orb as fsorb
from lib.filesystem import build as fbuild
from lib.filesystem import inf as finf
from lib.filesystem import read as fread
from lib.struct import tors as ptors


# Dictionary of Electronic Structure Calculations
ES_TSKS = {
    'energy': 'sp.run_energy',
    'grad': 'sp.run_gradient',
    'hess': 'sp.run_hessian',
    'vpt2': 'sp.run_vpt2'
}


def run(tsk, spc_dct, spc_name,
        thy_info, ini_thy_info,
        mr_sp_thy_info, mr_scn_thy_info,
        run_prefix, save_prefix,
        es_keyword_dct=None):
    """ run an electronic structure task
    for generating a list of conformer or tau sampling geometries
    """

    # Print the head of the task
    print('\n--------------------------------------------------------------------------------------')
    print('\nTask:', tsk, spc_name)
    print('\nOptions for electronic structure task:')
    for key, val in es_keyword_dct.items():
        print('{} = {}'.format(key, val))
    print('')

    # Set keys
    saddle = bool('ts_' in spc_name)
    # vdw = bool('vdw' in spc_name)
    spc = spc_dct[spc_name]

    # Get stuff from task
    [_, job] = tsk.split('_')

    # Run the task if an initial geom exists
    if 'init' in tsk:
        _ = run_geom_init(
            spc, thy_info, ini_thy_info,
            run_prefix, save_prefix, saddle, es_keyword_dct)
    elif 'find' in tsk:
        run_ts(spc_dct, spc_name,
               thy_info, ini_thy_info,
               mr_sp_thy_info, mref_scn_thy_info,
               run_prefix, save_prefix, es_keyword_dct)
    elif 'conf' in tsk:
        run_conformer_tsk(job, spc_dct, spc_name,
                          thy_info, ini_thy_info,
                          run_prefix, save_prefix,
                          saddle, es_keyword_dct)
    elif 'tau' in tsk:
        run_tau_tsk(job, spc, thy_info,
                    run_prefix, save_prefix,
                    es_keyword_dct)
    elif 'hr' in tsk:
        run_hr_tsk(job, spc_dct, spc_name,
                   thy_info, ini_thy_info,
                   run_prefix, save_prefix,
                   saddle, es_keyword_dct)
    elif 'irc' in tsk:
        run_irc_tsk(job, spc_dct, spc_name,
                    thy_info, ini_thy_info,
                    run_prefix, save_prefix,
                    saddle, es_keyword_dct)


# FUNCTIONS FOR SAMPLING AND SCANS #
def run_geom_init(spc, thy_info, ini_thy_info,
                  run_prefix, save_prefix, saddle, es_keyword_dct):
    """ Find the initial geometry
    """

    # Set the spc_info
    spc_info = finf.get_spc_info(spc)

    # Get es options
    [kickoff_size, kickoff_backward] = spc['kickoff']
    overwrite = es_keyword_dct['overwrite']
    dist_info = spc['dist_info'] if saddle else ()

    # Modify the theory
    mod_thy_info = fsorb.mod_orb_restrict(spc_info, thy_info)
    mod_ini_thy_info = fsorb.mod_orb_restrict(spc_info, ini_thy_info)

    # Set the filesystem objects
    thy_run_fs, thy_run_path = fbuild.spc_thy_fs_from_root(
        run_prefix, spc_info, mod_thy_info)
    thy_save_fs, thy_save_path = fbuild.spc_thy_fs_from_root(
        save_prefix, spc_info, mod_thy_info)
    ini_thy_save_fs, _ = fbuild.spc_thy_fs_from_root(
        save_prefix, spc_info, mod_ini_thy_info)
    cnf_run_fs, _ = fbuild.cnf_fs_from_thy(
        thy_run_path, saddle=saddle)
    cnf_save_fs, _ = fbuild.cnf_fs_from_thy(
        thy_save_path, saddle=saddle)

    # Set the run filesystem
    if saddle:
        _, ts_path = fbuild.ts_fs_from_thy(thy_run_path)
        run_fs = fbuild.run_fs_from_prefix(ts_path)
    else:
        run_fs = fbuild.run_fs_from_prefix(thy_run_path)

    # Get a reference geometry if one not found
    if not saddle:
        geo = geom.reference_geometry(
            spc, mod_thy_info, mod_ini_thy_info,
            thy_run_fs, thy_save_fs,
            ini_thy_save_fs,
            cnf_run_fs, cnf_save_fs, run_fs,
            kickoff_size=kickoff_size,
            kickoff_backward=kickoff_backward,
            overwrite=overwrite)
    else:
        geo = ts.sadpt_reference_geometry(
            spc, mod_thy_info, mod_ini_thy_info,
            thy_save_fs, ini_thy_save_fs,
            cnf_run_fs, cnf_save_fs, run_fs,
            dist_info=dist_info, overwrite=overwrite)

    return geo


def run_conformer_tsk(job, spc_dct, spc_name,
                      thy_info, ini_thy_info,
                      run_prefix, save_prefix,
                      saddle, es_keyword_dct):
    """ Launch tasks associated with conformers.
        Scan: Generate a set of conformer geometries and energies via
              random sampling over torsional coordinates
              following by optimization
        SP: Calculate ene, grad, ..
    """
    spc = spc_dct[spc_name]

    # Set the spc_info
    spc_info = finf.get_spc_info(spc)

    # Get es options
    overwrite = es_keyword_dct['overwrite']

    # Modify the theory
    mod_thy_info = fsorb.mod_orb_restrict(spc_info, thy_info)
    mod_ini_thy_info = fsorb.mod_orb_restrict(spc_info, ini_thy_info)

    # Set the filesystem objects
    if not saddle:
        # Build filesys for thy info
        _, thy_run_path = fbuild.spc_thy_fs_from_root(
            run_prefix, spc_info, mod_thy_info)
        thy_save_fs, thy_save_path = fbuild.spc_thy_fs_from_root(
            save_prefix, spc_info, mod_thy_info)

        # Build filesys for ini thy info
        _, ini_thy_run_path = fbuild.spc_thy_fs_from_root(
            run_prefix, spc_info, mod_ini_thy_info)
        ini_thy_save_fs, ini_thy_save_path = fbuild.spc_thy_fs_from_root(
            save_prefix, spc_info, mod_ini_thy_info)
    else:
        rxn_info = finf.rxn_info(
            spc['reacs'], spc['prods'], spc_dct)

        # Build filesys for thy info
        _, thy_run_path = fbuild.rxn_thy_fs_from_root(
            run_prefix, rxn_info, mod_thy_info)
        thy_save_fs, thy_save_path = fbuild.rxn_thy_fs_from_root(
            save_prefix, rxn_info, mod_thy_info)
        thy_save_fs, thy_save_path = fbuild.ts_fs_from_thy(thy_save_path)
        _, thy_run_path = fbuild.ts_fs_from_thy(thy_run_path)

        # Build filesys for ini thy info
        _, ini_thy_run_path = fbuild.rxn_thy_fs_from_root(
            run_prefix, rxn_info, mod_ini_thy_info)
        _, ini_thy_save_path = fbuild.rxn_thy_fs_from_root(
            save_prefix, rxn_info, mod_ini_thy_info)
        ini_thy_save_fs, ini_thy_save_path = fbuild.ts_fs_from_thy(
            ini_thy_save_path)
        _, ini_thy_run_path = fbuild.ts_fs_from_thy(
            ini_thy_run_path)

    if job == 'samp':

        # Build conformer filesys
        cnf_run_fs, _ = fbuild.cnf_fs_from_prefix(
            thy_run_path, cnf=None)
        cnf_save_fs, _ = fbuild.cnf_fs_from_prefix(
            thy_save_path, cnf=None)

        # Set up the run scripts
        _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
            *thy_info[0:2])

        # Set variables if it is a saddle
        tors_names = spc['tors_names'] if saddle else ''
        dist_info = spc['dist_info'] if saddle else ()
        two_stage = saddle
        rxn_class = spc['class'] if saddle else ''
        mc_nsamp = spc['mc_nsamp']

        # Run the sampling
        conformer.conformer_sampling(
            spc_info,
            mod_thy_info, mod_ini_thy_info,
            thy_save_fs, ini_thy_save_fs,
            cnf_run_fs, cnf_save_fs,
            opt_script_str, overwrite,
            saddle=saddle, nsamp_par=mc_nsamp,
            tors_names=tors_names, dist_info=dist_info,
            two_stage=two_stage, rxn_class=rxn_class, **opt_kwargs)

    elif job in ('energy', 'grad', 'hess', 'vpt2'):

        # Build conformer filesys
        cnf_run_fs, _ = fbuild.cnf_fs_from_prefix(
            ini_thy_run_path, cnf=None)
        cnf_save_fs, cnf_save_locs = fbuild.cnf_fs_from_prefix(
            ini_thy_save_path, cnf='min')

        # Check if locs exist, kill if it doesn't
        if not cnf_save_locs:
            print('*ERROR: No min-energy conformer found for level:',
                  ini_thy_save_path)
            sys.exit()

        # Set up the run scripts
        script_str, _, kwargs, _ = runpar.run_qchem_par(
            *thy_info[0:2])

        # Run the job over all the conformers requested by the user
        for locs in cnf_save_locs:
            geo_run_path = cnf_run_fs[-1].path([locs])
            geo_save_path = cnf_save_fs[-1].path([locs])
            cnf_run_fs[-1].create([locs])
            zma, geo = fread.get_zma_geo(cnf_save_fs, [locs])
            eval(ES_TSKS[job])(
                zma, geo, spc_info, mod_thy_info,
                cnf_save_fs, geo_run_path, geo_save_path, [locs],
                script_str, overwrite, **kwargs)


def run_tau_tsk(job, spc_dct, spc_name,
                thy_info, ini_thy_info,
                run_prefix, save_prefix,
                es_keyword_dct):
    """ Energies, gradients, and hessians,
        for set of arbitrarily sampled torsional coordinates
        with all other coordinates optimized
    """
    spc = spc_dct[spc_name]

    # Set the spc_info
    spc_info = finf.get_spc_info(spc)

    # Get es options
    overwrite = es_keyword_dct['overwrite']
    nsamp_par = spc['mc_nsamp']

    # Script
    _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
        *thy_info[0:2])

    # Modify the theory
    mod_thy_info = fsorb.mod_orb_restrict(spc_info, thy_info)
    mod_ini_thy_info = fsorb.mod_orb_restrict(spc_info, ini_thy_info)

    # Set the filesystem objects for thy info
    _, thy_run_path = fbuild.spc_thy_fs_from_root(
        run_prefix, spc_info, mod_thy_info)
    _, thy_save_path = fbuild.spc_thy_fs_from_root(
        save_prefix, spc_info, mod_thy_info)

    # Set the filesystem objects for ini thy_info
    ini_thy_save_fs, _ = fbuild.spc_thy_fs_from_root(
        save_prefix, spc_info, mod_ini_thy_info)

    # Set up tau filesystem objects
    tau_run_fs, _ = fbuild.tau_fs_from_thy(thy_run_path, tau='all')
    tau_save_fs, _ = fbuild.tau_fs_from_thy(thy_save_path, tau='all')

    if job == 'samp':
        _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
            *thy_info[0:2])
        tau.tau_sampling(
            spc_info,
            mod_thy_info, mod_ini_thy_info, ini_thy_save_fs,
            tau_run_fs, tau_save_fs,
            opt_script_str, overwrite, nsamp_par, **opt_kwargs)
    elif job in ('energy', 'grad', 'hess'):
        # Set up the run scripts
        script_str, _, kwargs, _ = runpar.run_qchem_par(
            *thy_info[0:2])
        tau_run_fs, _ = fbuild.tau_fs_from_thy(thy_run_path, tau='all')
        tau_save_fs, tau_locs = fbuild.tau_fs_from_thy(
            thy_save_path, tau='all')
        # Run the job over all the conformers requested by the user
        for locs in tau_locs:
            geo_run_path = tau_run_fs[-1].path(locs)
            geo_save_path = tau_save_fs[-1].path(locs)
            zma, geo = fread.get_zma_geo(tau_save_fs, locs)
            tau_run_fs[-1].create(locs)
            eval(ES_TSKS[job])(
                zma, geo, spc_info, mod_thy_info,
                tau_save_fs, geo_run_path, geo_save_path, locs,
                script_str, overwrite, **kwargs)


def run_hr_tsk(job, spc_dct, spc_name, thy_info, ini_thy_info,
               run_prefix, save_prefix,
               saddle, es_keyword_dct):
    """ run a scan over the specified torsional coordinates
    """

    spc = spc_dct[spc_name]

    # Set the spc_info
    spc_info = finf.get_spc_info(spc)

    # Script
    _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
        *thy_info[0:2])

    # Modify the theory
    mod_thy_info = fsorb.mod_orb_restrict(spc_info, thy_info)
    mod_ini_thy_info = fsorb.mod_orb_restrict(spc_info, ini_thy_info)

    # Set the filesystem objects
    if not saddle:
        _, thy_run_path = fbuild.spc_thy_fs_from_root(
            run_prefix, spc_info, mod_thy_info)
        _, thy_save_path = fbuild.spc_thy_fs_from_root(
            save_prefix, spc_info, mod_thy_info)
        _, ini_thy_run_path = fbuild.spc_thy_fs_from_root(
            run_prefix, spc_info, mod_ini_thy_info)
        _, ini_thy_save_path = fbuild.spc_thy_fs_from_root(
            save_prefix, spc_info, mod_ini_thy_info)
    else:
        rxn_info = finf.rxn_info(
            spc['reacs'], spc['prods'], spc_dct)

        # Build filesys for thy info
        _, thy_run_path = fbuild.rxn_thy_fs_from_root(
            run_prefix, rxn_info, mod_thy_info)
        _, thy_save_path = fbuild.rxn_thy_fs_from_root(
            save_prefix, rxn_info, mod_thy_info)
        _, thy_save_path = fbuild.ts_fs_from_thy(thy_save_path)
        _, thy_run_path = fbuild.ts_fs_from_thy(thy_run_path)

        # Build filesys for ini thy info
        _, ini_thy_run_path = fbuild.rxn_thy_fs_from_root(
            run_prefix, rxn_info, mod_ini_thy_info)
        _, ini_thy_save_path = fbuild.rxn_thy_fs_from_root(
            save_prefix, rxn_info, mod_ini_thy_info)
        _, ini_thy_save_path = fbuild.ts_fs_from_thy(
            ini_thy_save_path)
        _, ini_thy_run_path = fbuild.ts_fs_from_thy(
            ini_thy_run_path)

    # Build cnf filesys using the ini thy filesys (needed for all HR jobs)
    ini_cnf_run_fs, _ = fbuild.cnf_fs_from_prefix(
        ini_thy_run_path, cnf=None)
    ini_cnf_save_fs, ini_cnf_save_locs = fbuild.cnf_fs_from_prefix(
        ini_thy_save_path, cnf='min')
    ini_cnf_save_paths = fbuild.cnf_paths_from_locs(
        ini_cnf_save_fs, ini_cnf_save_locs)
    ini_cnf_run_paths = fbuild.cnf_paths_from_locs(
        ini_cnf_run_fs, ini_cnf_save_locs)

    # Create run fs if that directory has been deleted to run the jobs
    ini_cnf_run_fs[-1].create(ini_cnf_save_locs)

    # Get options from the dct or es options lst
    frm_bnd_key = spc['frm_bnd_key'] if saddle else ()
    brk_bnd_key = spc['brk_bnd_key'] if saddle else ()
    overwrite = es_keyword_dct['overwrite']
    ndim_tors = es_keyword_dct['ndim_tors']
    frz_all_tors = es_keyword_dct['frz_all_tors']
    scan_increment = spc['hind_inc']
    run_tors_names = spc['tors_names'] if 'tors_names' in spc else ()
    run_tors_names = [[name] for name in run_tors_names]

    # Set up the hind rot names by reading zma, geo from ini filesystem
    zma, geo = fread.get_zma_geo(ini_cnf_save_fs, ini_cnf_save_locs)
    run_tors_names, run_tors_grids = ptors.hr_prep(
        zma, geo, run_tors_names=run_tors_names,
        scan_increment=scan_increment, ndim_tors=ndim_tors,
        saddle=saddle, frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key)

    # Run the task if any torsions exist
    if run_tors_names:

        # Set constraint dct
        if not frz_all_tors:
            constraint_dct = None
        else:
            constraint_dct = ptors.build_constraint_dct(zma, run_tors_names)

        # Set up ini filesystem for scans
        ini_scn_run_fs = fbuild.scn_fs_from_cnf(
            ini_cnf_run_paths[0], constraint_dct=constraint_dct)
        ini_scn_save_fs = fbuild.scn_fs_from_cnf(
            ini_cnf_save_paths[0], constraint_dct=constraint_dct)

        if job == 'scan':

            scan.hindered_rotor_scans(
                zma, spc_info, mod_thy_info,
                ini_scn_run_fs, ini_scn_save_fs,
                run_tors_names, run_tors_grids,
                opt_script_str, overwrite,
                saddle=saddle, constraint_dct=constraint_dct, **opt_kwargs)

        elif job in ('energy', 'grad', 'hess', 'vpt2'):

            script_str, _, kwargs, _ = runpar.run_qchem_par(
                *thy_info[0:2])
            for tors_names in run_tors_names:
                scn_locs = fbuild.scn_locs_from_fs(
                    ini_scn_save_fs, tors_names, constraint_dct=constraint_dct)
                if scn_locs:
                    for locs in scn_locs:
                        geo_run_path = ini_scn_run_fs[-1].path(locs)
                        geo_save_path = ini_scn_save_fs[-1].path(locs)
                        zma, geo = fread.get_zma_geo(ini_scn_save_fs, locs)
                        ini_scn_run_fs[-1].create(locs)
                        eval(ES_TSKS[job])(
                            zma, geo, spc_info, mod_thy_info,
                            ini_scn_save_fs, geo_run_path, geo_save_path, locs,
                            script_str, overwrite, **kwargs)
                else:
                    print('*WARNING: NO SCAN INFORMATION EXISTS.',
                          'Doing scan vs cscan?')
    else:
        print('No torsional modes in the species')


def run_irc_tsk(job, spc_dct, spc_name, thy_info, ini_thy_info,
                run_prefix, save_prefix, es_keyword_dct):
    """ run a scan over the specified torsional coordinates
    """
    print('running task {}'.format('irc'))

    # Get dct for specific species task is run for
    spc = spc_dct[spc_name]

    # Set up coordinate name
    coo_name = 'RC'

    # Set the spc_info
    spc_info = finf.get_spc_info(spc)

    # Script
    _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
        *thy_info[0:2])

    # Modify the theory
    mod_thy_info = fsorb.mod_orb_restrict(spc_info, thy_info)
    mod_ini_thy_info = fsorb.mod_orb_restrict(spc_info, ini_thy_info)

    # Get options from the dct or es options lst
    overwrite = es_keyword_dct['overwrite']
    irc_idxs = spc['irc_idxs']

    # Set the filesystem objects
    rxn_info = finf.rxn_info(
        spc['reacs'], spc['prods'], spc_dct)

    _, thy_run_path = fbuild.rxn_thy_fs_from_root(
        run_prefix, rxn_info, mod_thy_info)
    _, thy_save_path = fbuild.rxn_thy_fs_from_root(
        save_prefix, rxn_info, mod_thy_info)
    _, thy_save_path = fbuild.ts_fs_from_thy(thy_save_path)
    _, thy_run_path = fbuild.ts_fs_from_thy(thy_run_path)

    _, ini_thy_run_path = fbuild.rxn_thy_fs_from_root(
        run_prefix, rxn_info, mod_ini_thy_info)
    _, ini_thy_save_path = fbuild.rxn_thy_fs_from_root(
        save_prefix, rxn_info, mod_ini_thy_info)
    _, ini_thy_save_path = fbuild.ts_fs_from_thy(ini_thy_save_path)
    _, ini_thy_run_path = fbuild.ts_fs_from_thy(ini_thy_run_path)

    ini_cnf_run_fs, _ = fbuild.cnf_fs_from_prefix(
        ini_thy_run_path, cnf=None)
    ini_cnf_save_fs, ini_cnf_save_locs = fbuild.cnf_fs_from_prefix(
        ini_thy_save_path, cnf='min')
    ini_cnf_save_paths = fbuild.cnf_paths_from_locs(
        ini_cnf_save_fs, ini_cnf_save_locs)
    ini_cnf_run_paths = fbuild.cnf_paths_from_locs(
        ini_cnf_run_fs, ini_cnf_save_locs)
    ini_cnf_run_fs[-1].create(ini_cnf_save_locs)

    ini_scn_run_fs = autofile.fs.scan(ini_cnf_run_paths[0])
    ini_scn_save_fs = autofile.fs.scan(ini_cnf_save_paths[0])

    if job == 'scan':

        # Script
        _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
            *thy_info[0:2])

        zma, geo = fread.get_zma_geo(ini_cnf_save_fs, ini_cnf_save_locs)
        ts.irc.irc_scan(
            geo, spc_info, mod_thy_info, coo_name, irc_idxs,
            ini_scn_save_fs, ini_scn_run_fs, ini_cnf_run_paths[0],
            overwrite, opt_script_str, **opt_kwargs)

    elif job in ('energy'):

        # Script
        script_str, _, kwargs, _ = runpar.run_qchem_par(
            *mod_thy_info[0:2])

        # Need to put in something with the IRC idxs
        for idx in irc_idxs:
            locs = [[coo_name], [idx]]
            geo_run_path = ini_scn_run_fs[-1].path(locs)
            geo_save_path = ini_scn_save_fs[-1].path(locs)
            zma, geo = fread.get_zma_geo(ini_scn_save_fs, locs)
            ini_scn_run_fs[-1].create(locs)
            eval(ES_TSKS[job])(
                zma, geo, spc_info, mod_thy_info,
                ini_scn_save_fs, geo_run_path, geo_save_path, locs,
                script_str, overwrite, **kwargs)


def run_ts(spc_dct, spc_name,
           thy_info, ini_thy_info,
           mr_sp_thy_info, mr_scn_thy_info,
           run_prefix, save_prefix,
           es_keyword_dct):
    """ find a transition state
    """

    # Get dct for specific species task is run for
    ts_dct = spc_dct[spc_name]

    # Build inf objects for the rxn and ts
    ts_info = ('', spc_dct[spc_name]['chg'], spc_dct[spc_name]['mul'])
    rxn_info = finf.rxn_info(
        spc_dct[spc_name]['reacs'], spc_dct[spc_name]['prods'], spc_dct)

    # Set various TS information using the dictionary
    typ = ts_dct['class']
    grid = ts_dct['grid']
    dist_name, _, update_guess, brk_name, _ = ts_dct['dist_info']

    # Get es options
    vrc_dct = {}
    overwrite = es_keyword_dct['overwrite']

    # Modify the theory
    mod_thy_info = fsorb.mod_orb_restrict(ts_info, thy_info)
    mod_ini_thy_info = fsorb.mod_orb_restrict(ts_info, ini_thy_info)
    if mr_sp_thy_info is not None:
        mod_mr_sp_thy_info = fsorb.mod_orb_restrict(ts_info, mr_sp_thy_info)
    else:
        mod_mr_sp_thy_info = None
    if mr_scn_thy_info is not None:
        mod_mr_scn_thy_info = fsorb.mod_orb_restrict(ts_info, mr_scn_thy_info)
    else:
        mod_mr_scn_thy_info = None
        
    _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
        *thy_info[0:2])

    # Build filesys for thy info for single reference
    _, thy_run_path = fbuild.rxn_thy_fs_from_root(
        run_prefix, rxn_info, mod_thy_info)
    thy_save_fs, thy_save_path = fbuild.rxn_thy_fs_from_root(
        save_prefix, rxn_info, mod_thy_info)
    
    # Build filesys for ini thy info for single reference
    _, ini_thy_run_path = fbuild.rxn_thy_fs_from_root(
        run_prefix, rxn_info, mod_ini_thy_info)
    ini_thy_save_fs, ini_thy_save_path = fbuild.rxn_thy_fs_from_root(
        save_prefix, rxn_info, mod_ini_thy_info)

    # Build the ts fs
    ts_save_fs, ts_save_path = fbuild.ts_fs_from_thy(thy_save_path)
    _, ts_run_path = fbuild.ts_fs_from_thy(thy_run_path)
    run_fs = autofile.fs.run(ts_run_path)

    # Set the cnf fs to see if TS is available or for searching
    cnf_save_fs, cnf_save_locs = fbuild.cnf_fs_from_prefix(
        thy_save_path, cnf='min')
    cnf_run_fs, _ = fbuild.cnf_fs_from_thy(
        thy_run_path, cnf=None, saddle=True)

    if cnf_save_locs and not overwrite:
        # ts_class, ts_original_zma, ts_tors_names, ts_dist_info
        # geo, zma, final_dist = check_filesys_for_ts(
        #     ts_dct, ts_zma, cnf_save_fs, overwrite,
        #     typ, dist_info, dist_name, bkp_ts_class_data)
        zma = cnf_save_fs[-1].file.zmatrix.read(cnf_save_locs)

        # Add an angle check which is added to spc dct for TS (crap code...)
        vals = automol.zmatrix.values(zma)
        final_dist = vals[dist_name]
        dist_info[1] = final_dist
        angle = ts.chk.check_angle(
            ts_dct['zma'],
            ts_dct['dist_info'],
            ts_dct['class'])
        ts_dct['dist_info'][1] = final_dist
        ts_dct['dist_info'].append(angle)

    else:
        if nobarrier(ts_dct):
            print('Running Scan for Barrierless TS:')

            # Build multireference thy info objects
            if mod_mr_scn_thy_info:
                _, thy_run_path = fbuild.rxn_thy_fs_from_root(
                    run_prefix, rxn_info, mod_mr_scn_thy_info)
                thy_save_fs, thy_save_path = fbuild.rxn_thy_fs_from_root(
                    save_prefix, rxn_info, mod_mr_scn_thy_info)
                scn_run_fs = autofile.fs.scan(thy_run_path)
                scn_save_fs = autofile.fs.scan(thy_save_path)
            else:
                print('Need mlvl specified')
            if mod_mr_sp_thy_info:
                _, thy_run_path = fbuild.rxn_thy_fs_from_root(
                    run_prefix, rxn_info, mod_mr_sp_thy_info)
                thy_save_fs, thy_save_path = fbuild.rxn_thy_fs_from_root(
                    save_prefix, rxn_info, mod_mr_sp_thy_info)
            else:
                print('Need mlvl specified')
            
            # Set up the scan filesys (need scan and cscan for rc
            scn_run_fs = autofile.fs.scan(thy_run_path)
            scn_save_fs = autofile.fs.scan(thy_save_path)

            # Run the barrierless transition state
            tsk_status, switch = ts.find.find_barrierless_transition_state(
                ts_info, ts_zma, ts_dct, spc_dct,
                grid,
                dist_name,
                rxn_run_path, rxn_save_path,
                rad_rad_ts,
                mod_ini_thy_info, mod_thy_info,
                multi_opt_info, multi_sp_info,
                run_prefix, save_prefix,
                scn_run_fs, scn_save_fs,
                overwrite, vrc_dct,
                update_guess, **opt_kwargs)

        # Run SadPt search 
        if not nobarrier(ts_dct) or (nobarrier(ts_dct) and switch):

            # Check and see if a zma is found from the filesystem
            ini_cnf_save_fs, ini_cnf_save_locs = fbuild.cnf_fs_from_prefix(
                ini_thy_save_path, cnf='min')
            guess_zma, _ = fread.get_zma_geo(ini_cnf_save_fs, ini_cnf_save_locs)

            # If no guess zma, run a TS searching algorithm
            if guess_zma is None:
                print('No Guess ZMA from inp thy. Running Scan for Fixed TS...')
                scn_run_fs = autofile.fs.scan(thy_run_path)
                scn_save_fs = autofile.fs.scan(thy_save_path)
                guess_zma = ts.find.run_sadpt_scan(
                    typ, grid, dist_name, brk_name, ts_dct['zma'], ts_info,
                    mod_thy_info,
                    scn_run_fs, scn_save_fs, opt_script_str,
                    overwrite, update_guess, **opt_kwargs)

            # Optimize the saddle point
            print('Running Scan for Fixed TS:')
            ts.find.find_sadpt_transition_state(
                opt_script_str,
                run_fs,
                guess_zma,
                ts_info,
                mod_thy_info,
                overwrite,
                ts_save_path,
                ts_save_fs,
                dist_name,
                dist_info,
                thy_save_fs,
                cnf_run_fs,
                cnf_save_fs,
                **opt_kwargs)


def nobarrier(ts_dct):
    rad_rad = bool('radical radical' in typ)
    low_spin = bool('low' in typ)
    no_elim = bool('elimination' not in ts_dct['class'])
    return rad_rad and low_sping and no_elim
