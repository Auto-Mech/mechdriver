""" eletronic structure routines modules
"""

# electronic structure routines
import autofile
from routines.es import conformer
from routines.es import geom
from routines.es import scan
from routines.es import sp
from routines.es import tau
from routines.es import ts
from routines.es import wells
from routines.es import find
from routines.es import variational
from lib.runner import par as runpar
from lib.filesystem import orb as fsorb
from lib.filesystem import build as fbuild
from lib.filesystem import inf as finf


__all__ = [
    'conformer',
    'geom',
    'scan',
    'sp',
    'tau',
    'ts',
    'wells',
    'find',
    'variational'
]


# Dictionary of Electronic Structure Calculations
ES_TSKS = {
    'energy': 'sp.run_energy',
    'grad': 'sp.run_gradient',
    'hess': 'sp.run_hessian',
    'vpt2': 'sp.run_vpt2'
}


def run_tsk(tsk, spc_dct, spc_name,
            thy_info, ini_thy_info,
            run_prefix, save_prefix,
            es_options=None):
    """ run an electronic structure task
    for generating a list of conformer or tau sampling geometries
    """

    print('Task:', tsk, spc_name)

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
            run_prefix, save_prefix, saddle, es_options)
    elif 'find' in tsk:
        run_ts(spc_dct, spc_name, thy_info, ini_thy_info,
               run_prefix, save_prefix, es_options)
    elif 'conf' in tsk:
        run_conformer_tsk(job, spc_dct, spc_name, thy_info,
                          run_prefix, save_prefix,
                          saddle, es_options)
    elif 'tau' in tsk:
        run_tau_tsk(job, spc, thy_info,
                    run_prefix, save_prefix,
                    es_options)
    elif 'hr' in tsk:
        run_hr_tsk(job, spc_dct, spc_name, thy_info,
                   run_prefix, save_prefix,
                   saddle, es_options)
    elif 'irc' in tsk:
        pass
        # run_irc_tsk(job, spc, thy_info, es_options)


# FUNCTIONS FOR SAMPLING AND SCANS #
def run_geom_init(spc, thy_info, ini_thy_info,
                  run_prefix, save_prefix, saddle, es_options):
    """ Find the initial geometry
    """
    # Set the spc_info
    spc_info = finf.get_spc_info(spc)

    # Get es options
    [kickoff_size, kickoff_backward] = spc['kickoff']
    overwrite = bool('overwrite' in es_options)
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


def run_conformer_tsk(job, spc_dct, spc_name, thy_info,
                      run_prefix, save_prefix,
                      saddle, es_options):
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
    overwrite = bool('overwrite' in es_options)

    # Modify the theory
    mod_thy_info = fsorb.mod_orb_restrict(spc_info, thy_info)

    # Set the filesystem objects
    if not saddle:
        _, thy_run_path = fbuild.spc_thy_fs_from_root(
            run_prefix, spc_info, mod_thy_info)
        thy_save_fs, thy_save_path = fbuild.spc_thy_fs_from_root(
            save_prefix, spc_info, mod_thy_info)
    else:
        rxn_info = finf.rxn_info(
            spc['reacs'], spc['prods'], spc_dct)
        _, thy_run_path = fbuild.rxn_thy_fs_from_root(
            run_prefix, rxn_info, mod_thy_info)
        thy_save_fs, thy_save_path = fbuild.rxn_thy_fs_from_root(
            save_prefix, rxn_info, mod_thy_info)
        thy_save_fs, thy_save_path = fbuild.ts_fs_from_thy(thy_save_path)
        thy_run_fs, thy_run_path = fbuild.ts_fs_from_thy(thy_run_path)
    cnf_run_fs, _ = fbuild.cnf_fs_from_prefix(
        thy_run_path, cnf=None)
    cnf_save_fs, cnf_save_locs = fbuild.cnf_fs_from_prefix(
        thy_save_path, cnf=None)

    if job == 'samp':

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
            spc_info, mod_thy_info,
            thy_save_fs, cnf_run_fs, cnf_save_fs,
            opt_script_str, overwrite,
            saddle=saddle, nsamp_par=mc_nsamp,
            tors_names=tors_names, dist_info=dist_info,
            two_stage=two_stage, rxn_class=rxn_class, **opt_kwargs)

    elif job in ('energy', 'grad', 'hess', 'vpt2'):
        cnf_run_fs, _ = fbuild.cnf_fs_from_prefix(
            thy_run_path, cnf=None)
        cnf_save_fs, cnf_save_locs = fbuild.cnf_fs_from_prefix(
            thy_save_path, cnf='min')

        # Set up the run scripts
        script_str, _, kwargs, _ = runpar.run_qchem_par(
            *thy_info[0:2])

        # Run the job over all the conformers requested by the user
        print('running task {}'.format('gradient'))
        for locs in cnf_save_locs:
            geo_run_path = cnf_run_fs[-1].path([locs])
            geo_save_path = cnf_save_fs[-1].path([locs])
            # if cnf_save_fs[-1].file.zmatrix.exists([locs]):
            if False:
                geo = cnf_save_fs[-1].file.zmatrix.read([locs])
            else:
                geo = cnf_save_fs[-1].file.geometry.read([locs])
            eval(ES_TSKS[job])(
                geo, spc_info, mod_thy_info,
                cnf_save_fs, geo_run_path, geo_save_path, [locs],
                script_str, overwrite, **kwargs)


def run_tau_tsk(job, spc, thy_info, run_prefix, save_prefix,
                es_options):
    """ Energies, gradients, and hessians,
        for set of arbitrarily sampled torsional coordinates
        with all other coordinates optimized
    """
    print('running task {}'.format('tau'))

    # Set the spc_info
    spc_info = finf.get_spc_info(spc)

    # Get es options
    overwrite = bool('overwrite' in es_options)
    nsamp_par = spc['mc_nsamp']

    # Script
    _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
        *thy_info[0:2])

    # Modify the theory
    mod_thy_info = fsorb.mod_orb_restrict(spc_info, thy_info)

    # Set the filesystem objects
    _, thy_run_path = fbuild.spc_thy_fs_from_root(
        run_prefix, spc_info, mod_thy_info)
    thy_save_fs, thy_save_path = fbuild.spc_thy_fs_from_root(
        save_prefix, spc_info, mod_thy_info)
    tau_run_fs, _ = fbuild.tau_fs_from_thy(thy_run_path, tau='all')
    tau_save_fs, _ = fbuild.tau_fs_from_thy(thy_save_path, tau='all')

    if job == 'samp':
        _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
            *thy_info[0:2])
        tau.tau_sampling(
            spc_info, mod_thy_info, thy_save_fs, tau_run_fs, tau_save_fs,
            opt_script_str, overwrite, nsamp_par, **opt_kwargs)
    elif job in ('energy', 'grad', 'hess'):
        # Set up the run scripts
        script_str, _, kwargs, _ = runpar.run_qchem_par(
            *thy_info[0:2])
        print('running task {}'.format(job))
        tau_run_fs, _ = fbuild.tau_fs_from_thy(thy_run_path, tau='all')
        tau_save_fs, tau_locs = fbuild.tau_fs_from_thy(
            thy_save_path, tau='all')
        # Run the job over all the conformers requested by the user
        for locs in tau_locs:
            geo_run_path = tau_run_fs[-1].path(locs)
            geo_save_path = tau_save_fs[-1].path(locs)
            if tau_save_fs[-1].file.zmatrix.exists(locs):
                geo = tau_save_fs[-1].file.zmatrix.read(locs)
            else:
                geo = tau_save_fs[-1].file.geometry.read(locs)
            eval(ES_TSKS[job])(
                geo, spc_info, mod_thy_info,
                tau_save_fs, geo_run_path, geo_save_path, locs,
                script_str, overwrite, **kwargs)


def run_hr_tsk(job, spc_dct, spc_name, thy_info,
               run_prefix, save_prefix,
               saddle, es_options):
    """ run a scan over the specified torsional coordinates
    """
    print('running task {}'.format('hr'))

    spc = spc_dct[spc_name]

    # Set the spc_info
    spc_info = finf.get_spc_info(spc)

    # Script
    _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
        *thy_info[0:2])

    # Modify the theory
    mod_thy_info = fsorb.mod_orb_restrict(spc_info, thy_info)

    # Set the filesystem objects
    if not saddle:
        _, thy_run_path = fbuild.spc_thy_fs_from_root(
            run_prefix, spc_info, mod_thy_info)
        _, thy_save_path = fbuild.spc_thy_fs_from_root(
            save_prefix, spc_info, mod_thy_info)
    else:
        rxn_info = finf.rxn_info(
            spc['reacs'], spc['prods'], spc_dct)
        _, thy_run_path = fbuild.rxn_thy_fs_from_root(
            run_prefix, rxn_info, mod_thy_info)
        thy_save_fs, thy_save_path = fbuild.rxn_thy_fs_from_root(
            save_prefix, rxn_info, mod_thy_info)
        thy_save_fs, thy_save_path = fbuild.ts_fs_from_thy(thy_save_path)
        _, thy_run_path = fbuild.ts_fs_from_thy(thy_run_path)
    cnf_run_fs, _ = fbuild.cnf_fs_from_prefix(
        thy_run_path, cnf=None)
    cnf_save_fs, cnf_save_locs = fbuild.cnf_fs_from_prefix(
        thy_save_path, cnf='min')
    cnf_save_paths = fbuild.cnf_paths_from_locs(
        cnf_save_fs, cnf_save_locs)
    cnf_run_paths = fbuild.cnf_paths_from_locs(
        cnf_run_fs, cnf_save_locs)

    # Set the filesystem objects
    # _, thy_run_path = fbuild.spc_thy_fs_from_root(
    #     run_prefix, spc_info, mod_thy_info)
    # _, thy_save_path = fbuild.spc_thy_fs_from_root(
    #     save_prefix, spc_info, mod_thy_info)
    # cnf_run_fs, _ = fbuild.cnf_fs_from_thy(
    #     thy_run_path, cnf=None, saddle=saddle)
    # cnf_save_fs, cnf_save_locs = fbuild.cnf_fs_from_thy(
    #     thy_save_path, cnf='min', saddle=saddle)
    # cnf_save_paths = fbuild.cnf_paths_from_locs(
    #     cnf_save_fs, cnf_save_locs)
    # cnf_run_paths = fbuild.cnf_paths_from_locs(
    #     cnf_run_fs, cnf_save_locs)

    # Get options from the dct or es options lst
    if saddle:
        frm_bnd_key = spc['frm_bnd_key']
        brk_bnd_key = spc['brk_bnd_key']
    else:
        frm_bnd_key = ()
        brk_bnd_key = ()
    overwrite = bool('overwrite' in es_options)
    ndim_tors = 'mdhr' if 'mdhr' in es_options else '1dhr'
    frz_all_tors = bool('frz_all_tors' in es_options)
    scan_increment = spc['hind_inc']
    run_tors_names = spc['tors_names'] if 'tors_names' in spc else ()
    run_tors_names = [[name] for name in run_tors_names]
    print('tors', run_tors_names)

    geo = cnf_save_fs[-1].file.geometry.read(cnf_save_locs)
    zma = cnf_save_fs[-1].file.zmatrix.read(cnf_save_locs)

    # Set up the hind rot names
    run_tors_names, run_tors_grids = scan.hr_prep(
        zma, geo, run_tors_names=run_tors_names,
        scan_increment=scan_increment, ndim_tors=ndim_tors,
        saddle=saddle, frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key)

    print('tors', run_tors_names)
    print('run paths', cnf_run_paths[0])
    print('save paths', cnf_save_paths[0])

    # Get a list of the other tors coords to freeze and set the filesystem
    if frz_all_tors:
        constraint_dct = scan.build_constraint_dct(zma, run_tors_names)
        scn_run_fs = autofile.fs.cscan(cnf_run_paths[0])
        scn_save_fs = autofile.fs.cscan(cnf_save_paths[0])
    else:
        constraint_dct = None
        scn_run_fs = autofile.fs.scan(cnf_run_paths[0])
        scn_save_fs = autofile.fs.scan(cnf_save_paths[0])

    if job == 'scan':
        scan.hindered_rotor_scans(
            zma, spc_info, thy_info, scn_run_fs, scn_save_fs,
            run_tors_names, run_tors_grids,
            opt_script_str, overwrite,
            saddle=saddle, constraint_dct=constraint_dct, **opt_kwargs)
    elif job in ('energy', 'grad', 'hess', 'vpt2'):
        if run_tors_names:
            script_str, _, kwargs, _ = runpar.run_qchem_par(
                *thy_info[0:2])
            for tors_names in run_tors_names:
                if not frz_all_tors:
                    scn_locs = fbuild.scn_locs_from_fs(
                        scn_save_fs, tors_names)
                else:
                    scn_locs = fbuild.cscn_locs_from_fs(
                        scn_save_fs, tors_names)
                for locs in scn_locs:
                    geo_run_path = scn_run_fs[-1].path(locs)
                    geo_save_path = scn_save_fs[-1].path(locs)
                    geo = scn_save_fs[-1].file.geometry.read(locs)
                    eval(ES_TSKS[job])(
                        geo, spc_info, mod_thy_info,
                        scn_save_fs, geo_run_path, geo_save_path, locs,
                        script_str, overwrite, **kwargs)


def run_irc_tsk(job, spc, thy_info,
                run_prefix, save_prefix,
                saddle, es_options):
    """ run a scan over the specified torsional coordinates
    """
    print('running task {}'.format('irc'))

    # Set up coordinate name
    coo_name = 'RC'

    # Set the spc_info
    spc_info = finf.get_spc_info(spc)

    # Script
    _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
        *thy_info[0:2])

    # Modify the theory
    mod_thy_info = fsorb.mod_orb_restrict(spc_info, thy_info)

    # Get options from the dct or es options lst
    overwrite = bool('overwrite' in es_options)
    irc_idxs = spc['irc_idxs']

    # Set the filesystem objects
    _, thy_run_path = fbuild.thy_fs_from_root(
        run_prefix, spc_info, mod_thy_info)
    _, thy_save_path = fbuild.thy_fs_from_root(
        save_prefix, spc_info, mod_thy_info)
    cnf_run_fs = fbuild.cnf_fs_from_thy(
        thy_run_path, saddle=saddle)
    cnf_save_fs = fbuild.cnf_fs_from_thy(
        thy_save_path, saddle=saddle)
    cnf_locs, cnf_save_paths = fbuild.cnf_locs_from_fs(
        cnf_save_fs, cnf='min')
    _, cnf_run_paths = fbuild.cnf_locs_from_fs(
        cnf_run_fs, cnf='min')
    run_fs = autofile.fs.run(cnf_run_paths[0])
    scn_run_fs = autofile.fs.scan(cnf_run_paths[0])
    scn_save_fs = autofile.fs.scan(cnf_save_paths[0])

    if job == 'scan':
        if cnf_save_fs[-1].file.zmatrix.exists(cnf_locs):
            geo = cnf_save_fs[-1].file.zmatrix.read(cnf_locs)
        else:
            geo = cnf_save_fs[-1].file.geometry.read(cnf_locs)
        variational.irc.irc_scan(
            geo, spc_info, mod_thy_info, coo_name, irc_idxs,
            scn_save_fs, scn_run_fs, run_fs,
            overwrite, opt_script_str, **opt_kwargs)
    elif job in ('energy'):
        script_str, _, kwargs, _ = runpar.run_qchem_par(
            *mod_thy_info[0:2])
        scn_locs = fbuild.cscn_locs_from_fs(
            scn_save_fs, [coo_name])
        for locs in scn_locs:
            geo_run_path = scn_run_fs[-1].path(locs)
            geo_save_path = scn_save_fs[-1].path(locs)
            geo = scn_save_fs[-1].file.geometry.read(locs)
            eval(ES_TSKS[job])(
                geo, spc_info, mod_thy_info,
                scn_save_fs, geo_run_path, geo_save_path, locs,
                script_str, overwrite, **kwargs)


def run_ts(spc_dct, spc_name,
           thy_info, ini_thy_info,
           run_prefix, save_prefix,
           es_options):
    """ find a transition state
    """
    #        for sadpt in ts_dct:
    #            if not ts_dct[sadpt]['class']:
    #                print('skipping reaction because type =',
    #                      ts_dct[sadpt]['class'])
    mod_multi_opt_info, mod_multi_sp_info = [], []
    vrc_dct = {}
    overwrite = bool('overwrite' in es_options)

    # Filesystem
    rxn_info = finf.rxn_info(
        spc_dct[spc_name]['reacs'], spc_dct[spc_name]['prods'], spc_dct)
    ts_info = ('', spc_dct[spc_name]['chg'], spc_dct[spc_name]['mul'])
    mod_thy_info = fsorb.mod_orb_restrict(ts_info, thy_info)
    mod_ini_thy_info = fsorb.mod_orb_restrict(ts_info, ini_thy_info)
    _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
        *thy_info[0:2])

    # set up fs
    if mod_multi_opt_info:
        _, thy_run_path = fbuild.rxn_thy_fs_from_root(
            run_prefix, rxn_info, mod_multi_opt_info)
        thy_save_fs, thy_save_path = fbuild.rxn_thy_fs_from_root(
            save_prefix, rxn_info, mod_multi_opt_info)
    else:
        _, thy_run_path = fbuild.rxn_thy_fs_from_root(
            run_prefix, rxn_info, mod_thy_info)
        thy_save_fs, thy_save_path = fbuild.rxn_thy_fs_from_root(
            save_prefix, rxn_info, mod_thy_info)
    ts_save_fs, ts_save_path = fbuild.ts_fs_from_thy(thy_save_path)
    ts_run_fs, ts_run_path = fbuild.ts_fs_from_thy(thy_run_path)
    cnf_run_fs, _ = fbuild.cnf_fs_from_thy(
        thy_run_path, cnf=None, saddle=True)
    cnf_save_fs, _ = fbuild.cnf_fs_from_thy(
        thy_save_path, cnf=None, saddle=True)
    scn_run_fs = autofile.fs.scan(thy_run_path)
    scn_save_fs = autofile.fs.scan(thy_save_path)
    run_fs = autofile.fs.run(ts_run_path)

    # Run
    find.find_ts(
        spc_dct, spc_dct[spc_name],
        spc_dct[spc_name]['zma'], ts_info,
        mod_ini_thy_info, mod_thy_info,
        mod_multi_opt_info, mod_multi_sp_info,
        thy_save_fs,
        cnf_run_fs, cnf_save_fs,
        scn_run_fs, scn_save_fs,
        run_fs,
        ts_save_fs, ts_save_path,
        run_prefix, save_prefix,
        opt_script_str,
        overwrite, vrc_dct,
        rad_rad_ts='ptst',
        **opt_kwargs)
#
#
# def run_vdw():
#     """ find a van der waals well
#     """
#     routines.es.wells.find_vdw(
#         spc, spc_dct, thy_info, ini_thy_info,
#         vdw_params,
#         thy_dct[es_run_key]['mc_nsamp'], run_prefix,
#         save_prefix, 0.1, False,
#         overwrite)
# def set_geom_analysis(tsk, ini_filesys, saddle):
#     """ set geometry analysis things
#     """
#     avail, selection, avail = False, None
#     if 'conf' in tsk and not saddle:
#         ini_save_fs = ini_filesys[3]
#         avail = fcheck.check_save(ini_save_fs, tsk, 'conf')
#         selection = 'min'
#     elif 'tau' in tsk and not saddle:
#         ini_save_fs = ini_filesys[5]
#         avail = fcheck.check_save(ini_save_fs, tsk, 'tau')
#         selection = 'all'
#     elif 'scan' in tsk and not saddle:
#         ini_save_fs = ini_filesys[7]
#         avail = fcheck.check_save(ini_save_fs, tsk, 'scan')
#         selection = 'all'
#
#     return avail, selection
