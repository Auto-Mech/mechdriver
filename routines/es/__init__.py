""" eletronic structure routines modules
"""

# electronic structure routines
from routines.es import conformer
from routines.es import geom
from routines.es import scan
from routines.es import sp
from routines.es import tau
from routines.es import ts
from routines.es import wells
from routines.es import find
from routines.es import variational
from lib.phydat import phycon
from lib.runner import par as runpar
from lib.filesystem import minc as fsmin
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
    # Conformers
    'conf_energy': 'run_ene',
    'conf_grad': 'run_grad',
    'conf_hess': 'run_hess',
    'conf_vpt2': 'run_vpt2',
    'conf_reopt': 'run_reopt',
    'conf_samp': 'run_conf_samp',
    # Tau Sampling
    'tau_energy': 'run_ene',
    'tau_grad': 'run_grad',
    'tau_hess': 'run_hess',
    'tau_vpt2': 'run_vpt2',
    'tau_reopt': 'run_reopt',
    'tau_samp': 'run_tau_samp',
    # Hindered Rotors
    'hr_energy': 'run_ene',
    'hr_grad': 'run_grad',
    'hr_hess': 'run_hess',
    'hr_reopt': 'run_reopt',
    'hr_scan': 'run_hr_scan',
    # MEP
    'mep_energy': 'run_ene',
    'mep_grad': 'run_grad',
    'mep_hess': 'run_hess',
    'mep_reopt': 'run_reopt'
}


def run_tsk(tsk, spc_name, spc
            ini_thy_level, thy_level,
            es_options=None):
    """ run an electronic structure task
    for generating a list of conformer or tau sampling geometries
    """
    
    print('Task:', tsk, spc_name)

    # Set keys
    # saddle = bool('ts_' in spc_name)
    # vdw = bool('vdw' in spc_name)

    # Get stuff from task
    [calc, job] = tsk.split('_')

    # Get an initial reference geom for the caluclations
    geom = run_geom_init(spc, es_options)

    # Run the task if an initial geom exists
    if geom:
        if 'init' in tsk:
            pass
        elif 'find' in tsk:
            run_ts()
        elif 'conf' in tsk:
            run_conformer_tsk(job, spc, thy_info, es_options)
        elif 'tau' in tsk:
            run_tau_tsk(job, spc, thy_info, es_options)
        elif 'hr' in tsk:
            run_hr_tsk(job, spc, thy_info, es_options)
        elif 'irc' in tsk:
            run_irc_tsk(job, spc, thy_info, es_options)


# FUNCTIONS FOR SAMPLING AND SCANS #
def run_geom_init(spc, es_options):
    """ Find the initial geometry
    """
    [kickoff_size, kickoff_backward] = spc['kickoff']
    # Get a reference geometry if one not found
    if not saddle:
        geo = geom.reference_geometry(
            spc, thy_level, ini_thy_level, filesys, ini_filesys,
            kickoff_size=kickoff_size,
            kickoff_backward=kickoff_backward,
            overwrite=overwrite)
    else:
        geo = ts.sadpt_reference_geometry(
            spc, thy_level, ini_thy_level, filesys, ini_filesys,
            spc['dist_info'], overwrite)


def run_conformer_tsk(job, spc, thy_info, es_options):
    """ Launch tasks associated with conformers. 
        Scan: Generate a set of conformer geometries and energies via
              random sampling over torsional coordinates following by optimization
        SP: Calculate ene, grad, ..
    """

    # Set the spc_info
    spc_info = finf.get_spc_info(spc)

    # Set up the run scripts
    _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
        *thy_level[0:2])

    # Set the filesystem objects
    thy_save_fs, thy_save_path = fbuild.thy_fs_from_root(
        save_prefix, spc_info, thy_info)
    _, thy_run_path = fbuild.thy_fs_from_root(
        run_prefix, spc_info, thy_info)
    cnf_run_fs, _ = cnf_fs_from_thy(thy_run_path, spc_info, thy_info, cnf='min')
    cnf_save_fs, _ = cnf_fs_from_thy(thy_save_path, spc_info, thy_info, cnf='min')

    if job == 'samp': 

        # Set variables if it is a saddle 
        tors_names = spc['tors_names'] if saddle else ''
        dist_info = spc['dist_info'] if saddle else ()
        two_stage = True if saddle else False
        rxn_class = spc['rxn_class'] if saddle else ''

        # Run the sampling
        conformer.conformer_sampling(
            spc_info, thy_level, thy_save_fs, cnf_run_fs, cnf_save_fs, script_str,
            overwrite, saddle=saddle, nsamp_par=mc_nsamp,
            tors_names=tors_names, dist_info=dist_info,
             two_stage=two_stage, rxn_class=rxn_class, **kwargs)

    elif job in ('energy', 'gradient', 'hessian', 'vpt2'):
        
        print('running task {}'.format('gradient'))
        # Run the job over all the conformers requested by the user 
        for locs in cnf_locs:
            eval(ES_TSKS[tsk])(
                spc_info, thy_level, cnf_run_fs, cnf_save_fs, locs,
                script_str, overwrite, **kwargs)


def run_tau_tsk(filesys, params, opt_kwargs):
    """ energies, gradients, and hessians,
    for set of arbitrarily sampled torsional coordinates
    with all other coordinates optimized
    """
    print('running task {}'.format('tau'))

    _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
        *thy_level[0:2])

    # Set the filesystem objects
    thy_save_fs, thy_save_path = fbuild.thy_fs_from_root(
        save_prefix, spc_info, thy_info)
    _, thy_run_path = fbuild.thy_fs_from_root(
        run_prefix, spc_info, thy_info)
    tau_run_fs, _ = tau_fs_from_thy(thy_run_path, spc_info, thy_info, tau='min')
    tau_save_fs, _ = tau_fs_from_thy(thy_save_path, spc_info, thy_info, tau='min')

    if job == 'samp': 

        # Set variables if it is a saddle 
        tors_names = spc['tors_names'] if saddle else ''
        dist_info = spc['dist_info'] if saddle else ()
        two_stage = True if saddle else False
        rxn_class = spc['rxn_class'] if saddle else ''

        tau.tau_sampling(
            spc_info, thy_level, thy_save_fs, tau_run_fs, tau_save_fs,
            script_str, overwrite, nsamp_par, **opt_kwargs):
    
    elif job in ('energy', 'gradient', 'hessian', 'vpt2'):

        print('running task {}'.format('gradient'))
        # Run the job over all the conformers requested by the user 
        for locs in cnf_locs:
            eval(ES_TSKS[tsk])(
                spc_info, thy_level, cnf_run_fs, cnf_save_fs, locs,
                script_str, overwrite, **kwargs)


def run_hr_tsk(filesys, params, opt_kwargs):
    """ run a scan over the specified torsional coordinates
    """
    print('running task {}'.format('hr'))


    _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
        *thy_level[0:2])

    # Set the filesystem objects
    thy_save_fs, thy_save_path = fbuild.thy_fs_from_root(
        save_prefix, spc_info, thy_info)
    _, thy_run_path = fbuild.thy_fs_from_root(
        run_prefix, spc_info, thy_info)
    tau_run_fs, _ = tau_fs_from_thy(thy_run_path, spc_info, thy_info, tau='min')
    tau_save_fs, _ = tau_fs_from_thy(thy_save_path, spc_info, thy_info, tau='min')

    # Get options from the dct or es options lst
    params['tors_model'] = tors_model
    # if 'hind_def' in spc:
    #     params['run_tors_names'] = spc['hind_def']
    # if 'hind_inc' in spc:
    #     params['scan_increment'] = spc['hind_inc'] * phycon.DEG2RAD
    # else:
    #     params['scan_increment'] = 30. * phycon.DEG2RAD
    if saddle:
        frm_bnd_key = spc['frm_bnd_key']
        brk_bnd_key = spc['brk_bnd_key']
    overwrite = True if 'overwrite' in es_options else False
    ndims = 'mdhr' if 'mdhr' in es_options else '1dhr'
    frz_tors = True if 'frz_tors' in es_options else False
    adiab_tors = False
    tors_model=(ndims, frz_tors, adiab_tors)

    if job == 'scan':
        scan.hindered_rotor_scans(
            spc_info, thy_level, cnf_run_fs, cnf_save_fs, script_str, overwrite,
            scan_increment=30.0, saddle=False,
            run_tors_names=(), frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key,
            gradient=False, hessian=False,
            tors_model=('1dhr', False, False), **opt_kwargs)


def run_irc_tsk(filesys, params, opt_kwargs):
    """ run a scan over the specified torsional coordinates
    """
    print('running task {}'.format('irc'))
    _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
        *scn_thy_info[0:2])

    if job == 'scan':
        variational.irc.irc_opt(
            ts_dct,
            thy_info,
            irc_idxs,
            overwrite)
    elif:
        routines.es.variational.irc.irc_sp(
            ts_dct[sadpt],
            thy_info,
            sp_thy_info,
            irc_idxs,
            es_options=es_options)


def run_ts():
    """ find a transition state
    """
    routines.es.find.find_ts(
        spc_dct, ts_dct[sadpt],
        ts_dct[sadpt]['zma'],
        ini_thy_info, thy_info,
        multi_opt_info, multi_sp_info,
        run_prefix, save_prefix,
        overwrite, vrc_dct,
        rad_rad_ts=rad_rad_ts,
        es_options=es_options)


def run_vdw():
    """ find a van der waals well
    """
    routines.es.wells.find_vdw(
        spc, spc_dct, thy_info, ini_thy_info,
        vdw_params,
        thy_dct[es_run_key]['mc_nsamp'], run_prefix,
        save_prefix, 0.1, False,
        overwrite)


def cnf_set():
    """ a """
    # Set the filesystem locs
    if selection == 'all':
        locs_lst = save_dir[-1].existing()
    elif selection == 'min':
        locs_lst = [fsmin.min_energy_conformer_locators(save_dir)]
    elif selection == 'subset':
        min_locs = fsmin.min_energy_conformer_locators(save_dir)
        min_ene = energy.read(min_lcs)  
        ini_locs_lst = save_dir[-1].existing()
        locs_lst = []
        for locs in ini_locs_lst:
            ene = energy.read(locs) - min_ene
            if ene <= ene_cut_off:
                locs_lst.append(locs)



def set_geom_analysis(tsk, ini_filesys, saddle):
    """ set geometry analysis things
    """
    avail, selection, avail = False, None
    if 'conf' in tsk and not saddle:
        ini_save_fs = ini_filesys[3]
        avail = fcheck.check_save(ini_save_fs, tsk, 'conf')
        selection = 'min'
    elif 'tau' in tsk and not saddle:
        ini_save_fs = ini_filesys[5]
        avail = fcheck.check_save(ini_save_fs, tsk, 'tau')
        selection = 'all'
    elif 'scan' in tsk and not saddle:
        ini_save_fs = ini_filesys[7]
        avail = fcheck.check_save(ini_save_fs, tsk, 'scan')
        selection = 'all'

    return avail, selection


if any(string in tsk for string in ('samp', 'scan', 'geom')):
    if not vdw:
        pass
    else:
        routines.es.wells.fake_geo_gen(tsk, thy_level, filesys)

