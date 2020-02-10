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
    'find'
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


def run_ene(params, kwargs):
    """ energy for geometry in fiven fs directory
    """
    print('running task {}'.format('energy'))
    sp.run_energy(**params, **kwargs)


def run_grad(params, kwargs):
    """ gradient for geometry in given fs directory
    """
    print('running task {}'.format('gradient'))
    sp.run_gradient(**params, **kwargs)


def run_hess(params, kwargs):
    """ hessian for geometry in given fs directory
    """
    print('running task {}'.format('hessian'))
    sp.run_hessian(**params, **kwargs)


def run_conf_samp(filesys, params, opt_kwargs):
    """ generate a set of conformer geometries and energies via
    random sampling over torsional coordinates following by optimization
    """
    if params['nsamp_par'][0]:
        print('running task {} with abcd of {}'.format(
            'conformer sampling', ' '.join(
                [str(x) for x in params['nsamp_par'][1:]])))
    else:
        print('running task {} for {:g} points'.format(
            'conformer sampling', params['nsamp_par'][5]))
    params['thy_save_fs'] = filesys[3]
    params['cnf_run_fs'] = filesys[4]
    params['cnf_save_fs'] = filesys[5]
    conformer.conformer_sampling(**params, **opt_kwargs)


def run_hr_scan(filesys, params, opt_kwargs):
    """ run a scan over the specified torsional coordinates
    """
    print('running task {}'.format('hr'))
    params['cnf_run_fs'] = filesys[4]
    params['cnf_save_fs'] = filesys[5]
    scan.hindered_rotor_scans(**params, **opt_kwargs)


def run_tau_samp(filesys, params, opt_kwargs):
    """ energies, gradients, and hessians,
    for set of arbitrarily sampled torsional coordinates
    with all other coordinates optimized
    """

    params['tau_run_fs'] = filesys[6]
    params['tau_save_fs'] = filesys[7]
    params['thy_save_fs'] = filesys[3]
    tau.tau_sampling(**params, **opt_kwargs)
    # del params['thy_save_fs']
    # del params['nsamp_par']
    # moldr.tau.run_tau_gradients(**params, **opt_kwargs)
    # moldr.tau.run_tau_hessians(**params, **opt_kwargs)


def geometry_generation(tsk, spc, mc_nsamp,
                        ini_thy_level, thy_level, ini_filesys, filesys,
                        overwrite, saddle=False, kickoff=(0.1, False),
                        tors_model=('1dhr', False)):
    """ run an electronic structure task
    for generating a list of conformer or tau sampling geometries
    """
    # Separate the kickoff keyword for now
    [kickoff_size, kickoff_backward] = kickoff

    spc_info = finf.get_spc_info(spc)

    # Get a reference geometry
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

    if geo:
        print('Task:', tsk)
        _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
            *thy_level[0:2])
        params = {'spc_info': spc_info,
                  'thy_level': thy_level,
                  'script_str': opt_script_str,
                  'overwrite': overwrite}
        if saddle:
            params['saddle'] = True
            params['tors_names'] = spc['tors_names']
        if tsk in ['conf_samp', 'tau_samp']:
            params['nsamp_par'] = mc_nsamp
            if saddle:
                params['dist_info'] = spc['dist_info']
                params['two_stage'] = True
                params['rxn_class'] = spc['class']
        elif tsk in ['hr_scan']:
            params['tors_model'] = tors_model
            if 'hind_def' in spc:
                params['run_tors_names'] = spc['hind_def']
            if 'hind_inc' in spc:
                params['scan_increment'] = spc['hind_inc'] * phycon.DEG2RAD
            else:
                params['scan_increment'] = 30. * phycon.DEG2RAD
            if saddle:
                params['frm_bnd_key'] = spc['frm_bnd_key']
                params['brk_bnd_key'] = spc['brk_bnd_key']
                print('key test in ts_geom gen:',
                      params['frm_bnd_key'],
                      params['brk_bnd_key'])

        if tsk in ES_TSKS:
            eval(ES_TSKS[tsk])(filesys, params, opt_kwargs)


def geometry_analysis(tsk, thy_level, ini_filesys,
                      spc, overwrite,
                      saddle=False, selection='min'):
    """ run the specified electronic structure task
    for a set of geometries
    """
    params = {}
    spc_info = finf.get_spc_info(spc)

    print('Task:', tsk)
    if 'conf' in tsk:
        run_dir = ini_filesys[2]
        save_dir = ini_filesys[3]
    elif 'tau' in tsk:
        run_dir = ini_filesys[4]
        save_dir = ini_filesys[5]
    elif 'hr' in tsk:
        run_dir = ini_filesys[6]
        save_dir = ini_filesys[7]
        if saddle:
            params['frm_bnd_key'] = spc['frm_bnd_key']
            params['brk_bnd_key'] = spc['brk_bnd_key']
            print('key test in ts_geom anal:',
                  params['frm_bnd_key'],
                  params['brk_bnd_key'])
    else:
        return
    if isinstance(selection, str):
        if selection == 'all':
            locs_lst = save_dir.leaf.existing()
        elif selection == 'min':
            locs_lst = [fsmin.min_energy_conformer_locators(save_dir)]
    else:
        locs_lst = selection

    sp_script_str, _, kwargs, _ = runpar.run_qchem_par(
        *thy_level[0:2])
    params['spc_info'] = spc_info
    params['thy_level'] = thy_level
    params['script_str'] = sp_script_str
    params['overwrite'] = overwrite

    # cycle over the locations
    if tsk in ES_TSKS:
        task_call = eval(ES_TSKS[tsk])
        for locs in locs_lst:
            if locs:
                params['geo_run_fs'] = run_dir
                params['geo_save_fs'] = save_dir
                params['locs'] = locs
                task_call(params, kwargs)
            else:
                print('No initial geometry available for {} on {}'.format(
                    spc_info[0], '/'.join(thy_level[1:3])))
