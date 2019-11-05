""" cycle over electronic structure calls
for species, TSs, and vdw species
"""
import os
import numpy
import automol
import elstruct
import autofile
import varecof_io
import moldr
from datalibs import phycon


def run_energy(params, kwargs):
    """ energy for geometry in fiven fs directory
    """
    print('running task {}'.format('energy'))
    moldr.sp.run_energy(**params, **kwargs)


def run_grad(params, kwargs):
    """ gradient for geometry in given fs directory
    """
    print('running task {}'.format('gradient'))
    moldr.sp.run_gradient(**params, **kwargs)


def run_hess(params, kwargs):
    """ hessian for geometry in given fs directory
    """
    print('running task {}'.format('hessian'))
    moldr.sp.run_hessian(**params, **kwargs)


def run_vpt2(params, kwargs):
    """ second order vibrational perturbation theory for minimum energy conformer
    """
    print('running task {}'.format('vpt2'))
    print('anharm may not be working in moldr right now')
    moldr.sp.run_vpt2(**params, **kwargs)


def run_irc(params, kwargs):
    """ irc for geometry in given fs directory
    """
    print('running task {}'.format('irc'))
    moldr.ts.run_irc(**params, **kwargs)


def run_single_conformer(spc_info, thy_level, fs, overwrite, saddle=False, dist_info=[]):
    """ generate single optimized geometry for randomly sampled initial torsional angles
    """
    mc_nsamp = [False, 0, 0, 0, 0, 1]
    sp_script_str, _, kwargs, _ = moldr.util.run_qchem_par(*thy_level[0:2])
    thy_save_fs = fs[3]
    two_stage = False
    if saddle:
        two_stage = True
    # if saddle:
        # thy_save_fs=fs[11]
    moldr.conformer.conformer_sampling(
        spc_info=spc_info,
        thy_level=thy_level,
        thy_save_fs=thy_save_fs,
        cnf_run_fs=fs[4],
        cnf_save_fs=fs[5],
        script_str=sp_script_str,
        overwrite=overwrite,
        nsamp_par=mc_nsamp,
        saddle=saddle,
        dist_info=dist_info,
        two_stage=two_stage,
        **kwargs,
    )


def run_conf_samp(fs, params, opt_kwargs):
    """ generate a set of conformer geometries and energies via
    random sampling over torsional coordinates following by optimization
    """
    if params['nsamp_par'][0]:
        print('running task {} with abcd of {}'.format('conformer sampling', ' '.join(
            [str(x) for x in params['nsamp_par'][1:]])))
    else:
        print('running task {} for {:g} points'.format(
            'conformer sampling', params['nsamp_par'][5]))
    #params['nsamp_par'] = spcdic['mc_nsamp']
    params['thy_save_fs'] = fs[3]
    params['cnf_run_fs'] = fs[4]
    params['cnf_save_fs'] = fs[5]
    moldr.conformer.conformer_sampling(**params, **opt_kwargs)


def run_hr_scan(fs, params, opt_kwargs):
    """ run a scan over the specified torsional coordinates
    """
    print('running task {}'.format('hr'))
    params['cnf_run_fs'] = fs[4]
    params['cnf_save_fs'] = fs[5]
    #print('frm_bnd_key in run_hr_scan:', params['frm_bnd_key'])
    moldr.scan.hindered_rotor_scans(**params, **opt_kwargs)


def run_tau_sampling(fs, params, opt_kwargs):
    """ energies, gradients, and hessians,
    for set of arbitrarily sampled torsional coordinates
    with all other coordinates optimized
    """

#    nsamp_par = opt_kwargs['mc_nsamp']  #turn to tau specific at some point
    params['tau_run_fs'] = fs[6]
    params['tau_save_fs'] = fs[7]
    params['thy_save_fs'] = fs[3]
#    params['nsamp_par'] = nsamp_par
    moldr.tau.tau_sampling(**params, **opt_kwargs)

    del params['thy_save_fs']
    del params['nsamp_par']

    moldr.tau.run_tau_gradients(**params, **opt_kwargs)
    moldr.tau.run_tau_hessians(**params, **opt_kwargs)


def geometry_generation(tsk, spcdic, es_dct, thy_level, fs,
    spc_info, overwrite):
    """ run an electronic structure task
    for generating a list of conformer or tau sampling geometries
    """
    print('Task in geometry_generation:', tsk)
    _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(*thy_level[0:2])
    params = {'spc_info': spc_info,
              'thy_level': thy_level,
              'script_str': opt_script_str,
              'overwrite': overwrite}
    choose_function = {'conf_samp': 'run_conf_samp',
                       'tau_samp': 'run_tau_samp',
                       'hr_scan': 'run_hr_scan'}

    if tsk in ['conf_samp', 'tau_samp']:
        params['nsamp_par'] = es_dct['mc_nsamp']
    elif tsk in ['hr_scan']:
        if 'hind_inc' in spcdic:
            params['scan_increment'] = spcdic['hind_inc']
        else:
            params['scan_increment'] = 30. * phycon.DEG2RAD

    if tsk in choose_function:
        eval(choose_function[tsk])(fs, params, opt_kwargs)


def ts_geometry_generation(tsk, spcdic, es_dct, thy_level, fs, spc_info, overwrite):
    """ run an electronic structure task
    for generating a list of conformer or tau sampling geometries
    """
    # fs[3] = fs[11]
    _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(*thy_level[0:2], saddle=True)
    print('tsk test in ts_geometry_generation:', tsk)
    params = {'spc_info': spc_info,
              'thy_level': thy_level,
              'script_str': opt_script_str,
              'saddle' :  True,
              'tors_names': spcdic['tors_names'],
              'overwrite': overwrite}
    choose_function = {'conf_samp': 'run_conf_samp',
                       'tau_samp': 'run_tau_samp',
                       'hr_scan': 'run_hr_scan'}

    # if 'insertion' and 'low' in spcdct['class']:
    #     opt_kwargs['scf_options'] = [
    #         elstruct.option.specify(
    #             elstruct.Option.Scf.Guess.MIX)
    #     ]
    if tsk in ['conf_samp', 'tau_samp']:
        params['nsamp_par'] = es_dct['mc_nsamp']
        params['dist_info'] = spcdic['dist_info']
        params['two_stage'] = True
    if tsk in ['conf_samp']:
        params['rxn_class'] = spcdic['class']
    elif tsk in ['hr_scan']:
        if 'hind_inc' in spcdic:
            params['scan_increment'] = spcdic['hind_inc']
        else:
            params['scan_increment'] = 30. * phycon.DEG2RAD
        params['frm_bnd_key'] = spcdic['frm_bnd_key']
        params['brk_bnd_key'] = spcdic['brk_bnd_key']
        print('key test in ts_geom gen:', params['frm_bnd_key'], params['brk_bnd_key'])

    if tsk in choose_function:
        eval(choose_function[tsk])(fs, params, opt_kwargs)


def geometry_analysis(tsk, thy_level, ini_fs, selection, spc_info, overwrite):
    """ run the specified electronic structure task
    for a set of geometries
    """

    print('Task in geometry_analysis:', tsk)
    # specify the fs for the runs
    if 'conf' in tsk:
        run_dir = ini_fs[2]
        save_dir = ini_fs[3]
    elif 'tau' in tsk:
        run_dir = ini_fs[4]
        save_dir = ini_fs[5]
    elif 'hr' in tsk:
        run_dir = ini_fs[6]
        save_dir = ini_fs[7]
    else:
        return
    # still need to setup mep
#    elif 'mep' in tsk:
#        run_dir = ini_fs[6]
#        save_dir = ini_fs[7]
    if isinstance(selection, str):
        if selection == 'all':
            locs_lst = save_dir.leaf.existing()
        elif selection == 'min':
            locs_lst = [moldr.util.min_energy_conformer_locators(save_dir)]
    elif isinstance(selection, int):
        locs_lst = moldr.util.geom_sort(save_dir)[0:selection]
    else:
        locs_lst = selection

    sp_script_str, _, kwargs, _ = moldr.util.run_qchem_par(*thy_level[0:2])
    params = {'spc_info': spc_info,
              'thy_level': thy_level,
              'script_str': sp_script_str,
              'overwrite': overwrite}
    choose_function = {'conf_energy': 'run_energy',
                       'tau_energy': 'run_energy',
                       'hr_energy': 'run_energy',
                       'mep_energy': 'run_energy',
                       'conf_grad': 'run_grad',
                       'tau_grad': 'run_grad',
                       'hr_grad': 'run_grad',
                       'mep_grad': 'run_grad',
                       'conf_hess': 'run_hess',
                       'tau_hess': 'run_hess',
                       'hr_hess': 'run_hess',
                       'mep_hess': 'run_hess',
                       'conf_vpt2': 'run_vpt2',
                       'tau_vpt2': 'run_vpt2',
                       'conf_reopt': 'run_reopt',
                       'tau_reopt': 'run_reopt',
                       'hr_reopt': 'run_reopt',
                       'mep_reopt': 'run_reopt'}

    # cycle over the locations
    if tsk in choose_function:
        task_call = eval(choose_function[tsk])
        for locs in locs_lst:
            if locs:
                params['geo_run_fs'] = run_dir
                params['geo_save_fs'] = save_dir
                params['locs'] = locs
                task_call(params, kwargs)
            else:
                print('No initial geometry available for {} on {}'.format(
                    spc_info[0], '/'.join(thy_level[1:3])))


def ts_geometry_analysis(tsk, thy_level, ini_fs, selection, spc_info, spc_dic, overwrite):
    """ run the specified electronic structure task
    for a set of geometries
    """

    print('Task in ts geometry_analysis:', tsk)
    params = {'spc_info': spc_info,
              'thy_level': thy_level,
              'overwrite': overwrite}
    # specify the fs for the runs
    if 'conf' in tsk:
        run_dir = ini_fs[2]
        save_dir = ini_fs[3]
    elif 'tau' in tsk:
        run_dir = ini_fs[4]
        save_dir = ini_fs[5]
    elif 'hr' in tsk:
        run_dir = ini_fs[6]
        save_dir = ini_fs[7]
        params['frm_bnd_key'] = spc_dic['frm_bnd_key']
        params['brk_bnd_key'] = spc_dic['brk_bnd_key']
        print('key test in ts_geom anal:', params['frm_bnd_key'], params['brk_bnd_key'])
    else:
        return
    # still need to setup mep
#    elif 'mep' in tsk:
#        run_dir = ini_fs[6]
#        save_dir = ini_fs[7]
    if isinstance(selection, str):
        if selection == 'all':
            locs_lst = save_dir.leaf.existing()
        elif selection == 'min':
            locs_lst = [moldr.util.min_energy_conformer_locators(save_dir)]
    elif isinstance(selection, int):
        locs_lst = moldr.util.geom_sort(save_dir)[0:selection]
    else:
        locs_lst = selection
    sp_script_str, _, kwargs, _ = moldr.util.run_qchem_par(*thy_level[0:2], saddle=True)
    params['script_str'] = sp_script_str
    choose_function = {'conf_energy': 'run_energy',
                       'tau_energy': 'run_energy',
                       'hr_energy': 'run_energy',
                       'mep_energy': 'run_energy',
                       'conf_grad': 'run_grad',
                       'tau_grad': 'run_grad',
                       'hr_grad': 'run_grad',
                       'mep_grad': 'run_grad',
                       'conf_hess': 'run_hess',
                       'tau_hess': 'run_hess',
                       'hr_hess': 'run_hess',
                       'mep_hess': 'run_hess',
                       'conf_vpt2': 'run_vpt2',
                       'tau_vpt2': 'run_vpt2',
                       'conf_reopt': 'run_reopt',
                       'tau_reopt': 'run_reopt',
                       'hr_reopt': 'run_reopt',
                       'mep_reopt': 'run_reopt'}

    # cycle over the locations

    if tsk in choose_function:
        task_call = eval(choose_function[tsk])
        for locs in locs_lst:
            if locs:
                params['geo_run_fs'] = run_dir
                params['geo_save_fs'] = save_dir
                params['locs'] = locs
                task_call(params, kwargs)
            else:
                print('No initial geometry available for {} on {}'.format(
                    spc_info[0], '/'.join(thy_level[1:3])))


def get_spc_run_path(run_prefix, spc_info):
    """ create species run path
    """
    spc_run_fs = autofile.fs.species(run_prefix)
    spc_run_fs.leaf.create(spc_info)
    spc_run_path = spc_run_fs.leaf.path(spc_info)
    return spc_run_path


def get_spc_save_path(save_prefix, spc_info):
    """ create species save path
    """
    spc_save_fs = autofile.fs.species(save_prefix)
    spc_save_fs.leaf.create(spc_info)
    spc_save_path = spc_save_fs.leaf.path(spc_info)
    return spc_save_path


def get_thy_run_path(run_prefix, spc_info, thy_info):
    """ create theory run path
    """
    orb_restr = moldr.util.orbital_restriction(
        spc_info, thy_info)
    thy_lvl = thy_info[0:3]
    thy_lvl.append(orb_restr)
    spc_run_path = get_spc_run_path(run_prefix, spc_info)
    thy_run_fs = autofile.fs.theory(spc_run_path)
    thy_run_fs.leaf.create(thy_lvl)
    thy_run_path = thy_run_fs.leaf.path(thy_lvl)
    return thy_run_path


def get_thy_save_fs(save_prefix, spc_info, thy_info):
    """ create theory save filesystem
    """
    orb_restr = moldr.util.orbital_restriction(
        spc_info, thy_info)
    thy_lvl = thy_info[0:3]
    thy_lvl.append(orb_restr)
    spc_save_path = get_spc_save_path(save_prefix, spc_info)
    thy_save_fs = autofile.fs.theory(spc_save_path)
    return thy_save_fs, thy_lvl


def get_thy_save_path(save_prefix, spc_info, thy_info):
    """ create theory save path
    """
    orb_restr = moldr.util.orbital_restriction(
        spc_info, thy_info)
    thy_lvl = thy_info[0:3]
    thy_lvl.append(orb_restr)
    spc_save_path = get_spc_save_path(save_prefix, spc_info)
    thy_save_fs = autofile.fs.theory(spc_save_path)
    thy_save_fs.leaf.create(thy_lvl)
    thy_save_path = thy_save_fs.leaf.path(thy_lvl)
    return thy_save_path


def rxn_info(run_prefix, save_prefix, ts, spc_dct, thy_info, ini_thy_info=None):
    """ prepare rxn info and reverse reactants and products if reaction is endothermic
    """
    rxn_ichs = [[], []]
    rxn_chgs = [[], []]
    rxn_muls = [[], []]
    print('\n TS for {}: {} = {}'.format(
        ts, '+'.join(spc_dct[ts]['reacs']), '+'.join(spc_dct[ts]['prods'])))
    reacs = spc_dct[ts]['reacs']
    prods = spc_dct[ts]['prods']
    for spc in reacs:
        rxn_ichs[0].append(spc_dct[spc]['ich'])
        rxn_chgs[0].append(spc_dct[spc]['chg'])
        rxn_muls[0].append(spc_dct[spc]['mul'])
    for spc in prods:
        rxn_ichs[1].append(spc_dct[spc]['ich'])
        rxn_chgs[1].append(spc_dct[spc]['chg'])
        rxn_muls[1].append(spc_dct[spc]['mul'])
    # check direction of reaction
    try:
        rxn_exo = moldr.util.reaction_energy(
            save_prefix, rxn_ichs, rxn_chgs, rxn_muls, thy_info)
    except:
        rxn_exo = moldr.util.reaction_energy(
            save_prefix, rxn_ichs, rxn_chgs, rxn_muls, ini_thy_info)
    print('reaction is {:.2f} endothermic'.format(rxn_exo*phycon.EH2KCAL))
    if rxn_exo > 0:
        rxn_ichs = rxn_ichs[::-1]
        rxn_chgs = rxn_chgs[::-1]
        rxn_muls = rxn_muls[::-1]
        spc_dct[ts]['reacs'] = prods
        spc_dct[ts]['prods'] = reacs
        print('Reaction will proceed as {}: {} = {}'.format(
            ts, '+'.join(spc_dct[ts]['reacs']), '+'.join(spc_dct[ts]['prods'])))
        # print('ts search will be performed in reverse direction')

    # set up the filesystem
    rxn_ichs, rxn_chgs, rxn_muls = autofile.system.sort_together(
        rxn_ichs, rxn_chgs, rxn_muls)
    low_mul = min(
        automol.mult.ts._low(rxn_muls[0]), automol.mult.ts._low(rxn_muls[1]))
    high_mul = max(
        automol.mult.ts._high(rxn_muls[0]), automol.mult.ts._high(rxn_muls[1]))

    return rxn_ichs, rxn_chgs, rxn_muls, low_mul, high_mul


def get_rxn_fs(run_prefix, save_prefix, ts):
    """ get filesystems for a reaction
    """
    rxn_ichs = ts['rxn_ichs']
    rxn_chgs = ts['rxn_chgs']
    rxn_muls = ts['rxn_muls']
    ts_mul = ts['mul']

    rxn_run_fs = autofile.fs.reaction(run_prefix)
    rxn_run_fs.leaf.create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
    rxn_run_path = rxn_run_fs.leaf.path(
        [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

    rxn_ichs = tuple(map(tuple, rxn_ichs))
    rxn_chgs = tuple(map(tuple, rxn_chgs))
    rxn_muls = tuple(map(tuple, rxn_muls))
    rxn_save_fs = autofile.fs.reaction(save_prefix)
    rxn_save_fs.leaf.create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
    rxn_save_path = rxn_save_fs.leaf.path(
        [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

    return rxn_run_fs, rxn_save_fs, rxn_run_path, rxn_save_path


def get_zmas(
        reacs, prods, spc_dct, ini_thy_info, save_prefix, run_prefix,
        kickoff_size, kickoff_backward, projrot_script_str):
    """get the zmats for reactants and products using the initial level of theory
    """
    if len(reacs) > 2:
        ich = spc_dct[reacs[-1]]['ich']
        ichgeo = automol.inchi.geometry(ich)
        ichzma = automol.geom.zmatrix(ichgeo)
        reacs = reacs[:-1]
    elif len(prods) > 2:
        ich = spc_dct[prods[-1]]['ich']
        ichgeo = automol.inchi.geometry(ich)
        ichzma = automol.geom.zmatrix(ichgeo)
        prods = prods[:-1]
    rct_geos, rct_cnf_save_fs_lst = get_geos(
        reacs, spc_dct, ini_thy_info, save_prefix, run_prefix, kickoff_size,
        kickoff_backward, projrot_script_str)
    prd_geos, prd_cnf_save_fs_lst = get_geos(
        prods, spc_dct, ini_thy_info, save_prefix, run_prefix, kickoff_size,
        kickoff_backward, projrot_script_str)
    rct_zmas = list(map(automol.geom.zmatrix, rct_geos))
    prd_zmas = list(map(automol.geom.zmatrix, prd_geos))
    for geo in prd_geos:
        xyzs = automol.geom.coordinates(geo)
    if len(rct_zmas) > 2:
        rct_zmas.append(ichzma)
    if len(prd_zmas) > 2:
        prd_zmas.append(ichzma)
    return rct_zmas, prd_zmas, rct_cnf_save_fs_lst, prd_cnf_save_fs_lst


def get_geos(
        spcs, spc_dct, ini_thy_info, save_prefix, run_prefix, kickoff_size,
        kickoff_backward, projrot_script_str):
    """get geos for reactants and products using the initial level of theory
    """
    spc_geos = []
    cnf_save_fs_lst = []
    for spc in spcs:
        spc_info = [spc_dct[spc]['ich'], spc_dct[spc]['chg'], spc_dct[spc]['mul']]
        orb_restr = moldr.util.orbital_restriction(spc_info, ini_thy_info)
        ini_thy_level = ini_thy_info[0:3]
        ini_thy_level.append(orb_restr)
        spc_save_fs = autofile.fs.species(save_prefix)
        spc_save_fs.leaf.create(spc_info)
        spc_save_path = spc_save_fs.leaf.path(spc_info)
        spc_run_fs = autofile.fs.species(run_prefix)
        spc_run_fs.leaf.create(spc_info)
        spc_run_path = spc_run_fs.leaf.path(spc_info)
        ini_thy_save_fs = autofile.fs.theory(spc_save_path)
        ini_thy_save_path = ini_thy_save_fs.leaf.path(ini_thy_level[1:4])
        ini_thy_run_fs = autofile.fs.theory(spc_run_path)
        ini_thy_run_path = ini_thy_run_fs.leaf.path(ini_thy_level[1:4])
        cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)
        cnf_save_fs_lst.append(cnf_save_fs)
        cnf_run_fs = autofile.fs.conformer(ini_thy_run_path)
        min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
        if min_cnf_locs:
            geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        else:
            run_fs = autofile.fs.run(ini_thy_run_path)
            run_fs.trunk.create()
            tmp_ini_fs = [None, ini_thy_save_fs]
            tmp_fs = [spc_save_fs, spc_run_fs, ini_thy_save_fs, ini_thy_run_fs,
                      cnf_save_fs, cnf_run_fs, run_fs]
            geo = moldr.geom.reference_geometry(
                spc_dct[spc], ini_thy_level, ini_thy_level, tmp_fs, tmp_ini_fs,
                kickoff_size, kickoff_backward, projrot_script_str,
                overwrite=False)
        spc_geos.append(geo)
    return spc_geos, cnf_save_fs_lst


def ts_class(rct_zmas, prd_zmas, rad_rad, ts_mul, low_mul, high_mul, rct_cnf_save_fs_lst, prd_cnf_save_fs_lst):
    """ determine type of reaction and related ts info from the reactant and product z-matrices.
    Returns the type, the transition state z-matrix, the name of the coordinate to optimize,
    the grid of values for the initial grid search, the torsion names and symmetries, and
    whether or not to update the guess on successive steps.
    These parameters are set for both the initial and a backup evaluation for if the initial ts
    search fails.
    """

    # Force a trimolecular reaction to behave like a bimolecular.
    # Termolecular species generally arise from the direct decomposition of some inital product.
    # We need to be able to find the TS for the channel preceding that direct decomposition.
    rct_tors_names = []
    if len(rct_zmas) > 2:
        ret = automol.zmatrix.ts.addition(rct_zmas[1:-1], [prd_zmas[-1]])
        new_zma, dist_name, rct_tors_names = ret
        new_zma = automol.zmatrix.standard_form(new_zma)
        babs2 = automol.zmatrix.get_babs2(new_zma, dist_name)
        new_zma = automol.zmatrix.set_values(new_zma, {dist_name: 2.2, babs2: 180.})
        rct_zmas = [rct_zmas[0], new_zma]
    elif len(prd_zmas) > 2:
        ret = automol.zmatrix.ts.addition(prd_zmas[1:-1], [prd_zmas[-1]])
        new_zma, dist_name, rct_tors_names = ret
        new_zma = automol.zmatrix.standard_form(new_zma)
        babs1 = automol.zmatrix.get_babs1(new_zma, dist_name)
        aabs1 = babs1.replace('D', 'A')
        new_zma = automol.zmatrix.set_values(
            new_zma, {dist_name: 2.2, aabs1: 170. * phycon.DEG2RAD})
        prd_zmas = [prd_zmas[0], new_zma]

    typ = None
    bkp_typ = ''
    brk_name = []
    frm_bnd_key = []
    brk_bnd_key = []
    # cycle through each of the possible reaction types checking if the reaction is in that class
    # Check both orders of reactants and products
    for direction in ('forward', 'reverse'):

        print('direction')
        print(direction)

        # Set proper cnf filesystem and flip reactants and products for second check
        if direction == 'forward':
            cnf_save_fs_lst = rct_cnf_save_fs_lst
        elif direction == 'reverse':
            zmas = [rct_zmas, prd_zmas]
            rct_zmas, prd_zmas = zmas[1], zmas[0]
            cnf_save_fs_lst = prd_cnf_save_fs_lst

        # Check for addition
        ret = automol.zmatrix.ts.addition(rct_zmas, prd_zmas, rct_tors_names)
        if ret:
            typ = 'addition'
            ts_zma, dist_name, tors_names = ret
            if ts_mul == high_mul:
                typ += ': high spin'
            elif ts_mul == low_mul:
                typ += ': low spin'
            # set up beta scission as a fall back option for failed addition TS search
            ret2 = automol.zmatrix.ts.beta_scission(rct_zmas, prd_zmas)
            if ret2:
                bkp_typ = 'beta scission'
                bkp_ts_zma, bkp_dist_name, bkp_tors_names = ret2

        # Check for beta-scission
        if typ is None:
            ret = automol.zmatrix.ts.beta_scission(rct_zmas, prd_zmas)
            if ret:
                typ = 'beta scission'
                ts_zma, dist_name, tors_names = ret
                ret2 = automol.zmatrix.ts.addition(prd_zmas, rct_zmas, rct_tors_names)
                if ret2:
                    bkp_typ = 'addition'
                    bkp_ts_zma, bkp_dist_name, bkp_tors_names = ret2

        # Check for hydrogen migration
        if typ is None:
            orig_dist = automol.zmatrix.ts.min_hyd_mig_dist(rct_zmas, prd_zmas)
            if orig_dist:
                rct_zmas = moldr.util.min_dist_conformer_zma_geo(orig_dist, cnf_save_fs_lst[0])
                ret = automol.zmatrix.ts.hydrogen_migration(rct_zmas, prd_zmas)
                if ret:
                    typ = 'hydrogen migration'
                    ts_zma, dist_name, frm_bnd_key, brk_bnd_key, tors_names = ret

        # Check for hydrogen abstraction
        if typ is None:
            ret = automol.zmatrix.ts.hydrogen_abstraction(rct_zmas, prd_zmas, sigma=False)
            #print('abstraction ret test in ts_class:', ret)
            if ret:
                typ = 'hydrogen abstraction'
                ts_zma, dist_name, frm_bnd_key, brk_bnd_key, tors_names = ret
                if ts_mul == high_mul:
                    typ += ': high spin'
                elif ts_mul == low_mul:
                    typ += ': low spin'
                brk_name = automol.zmatrix.bond_key_from_idxs(ts_zma, brk_bnd_key)
            #print('key test in ts_class:', frm_bnd_key, brk_bnd_key)
                    
        # Need special cases for (i) hydrogen abstraction where the radical is a sigma radical
        # and (ii) for abstraction of a heavy atom rather than a hydrogen atom. 
        # add these later
        # ret = automol.zmatrix.ts.hydrogen_abstraction(rct_zmas, prd_zmas, sigma=True)

        # Check for insertion
        if typ is None:
            ret = automol.zmatrix.ts.insertion(rct_zmas, prd_zmas)
            if ret:
                typ = 'insertion'
                ts_zma, dist_name, tors_names = ret
                if ts_mul == high_mul:
                    typ += ': high spin'
                elif ts_mul == low_mul:
                    typ += ': low spin'

        # Check for subsitution
        if typ is None:
            ret = automol.zmatrix.ts.substitution(rct_zmas, prd_zmas)
            if ret:
                typ = 'substitution'
                ts_zma, dist_name, tors_names = ret
                if ts_mul == high_mul:
                    typ += ': high spin'
                elif ts_mul == low_mul:
                    typ += ': low spin'

        # Check for elimination
        if typ is None:
            orig_dist = automol.zmatrix.ts.min_unimolecular_elimination_dist(rct_zmas, prd_zmas)
            if orig_dist:
                rct_zmas = moldr.util.min_dist_conformer_zma_geo(orig_dist, cnf_save_fs_lst[0])
                ret = automol.zmatrix.ts.concerted_unimolecular_elimination(rct_zmas, prd_zmas)
                if ret:
                    typ = 'elimination'
                    ts_zma, dist_name, brk_name, frm_bnd_key, tors_names = ret
                    if ts_mul == high_mul:
                        typ += ': high spin'
                    elif ts_mul == low_mul:
                        typ += ': low spin'

        # Break if reaction found
        if typ is not None:
            break

    # Nothing was found
    if typ is None:
        print("Failed to classify reaction.")
        return [], []

    # set up back up options for any radical radical case
    if rad_rad:
        typ = 'radical radical ' + typ
    print("Type: {}".format(typ))
    if bkp_typ:
        if rad_rad:
            bkp_typ = 'radical radical ' + bkp_typ

    # determine the grid for the preliminary grid search for all the different reaction types
    dist_coo, = automol.zmatrix.coordinates(ts_zma)[dist_name]
    syms = automol.zmatrix.symbols(ts_zma)
    bnd_len_key = tuple(sorted(map(syms.__getitem__, dist_coo)))

    bnd_len_dct = {
        ('C', 'C'): 1.54 * phycon.ANG2BOHR,
        ('C', 'H'): 1.09 * phycon.ANG2BOHR,
        ('H', 'H'): 0.74 * phycon.ANG2BOHR,
        ('N', 'N'): 1.45 * phycon.ANG2BOHR,
        ('O', 'O'): 1.48 * phycon.ANG2BOHR,
        ('C', 'N'): 1.47 * phycon.ANG2BOHR,
        ('C', 'O'): 1.43 * phycon.ANG2BOHR,
        ('H', 'O'): 0.95 * phycon.ANG2BOHR,
    }

    npoints = 8
    npoints1 = 4
    npoints2 = 4
    if 'beta scission' in bkp_typ:
        rmin = 1.4 * phycon.ANG2BOHR
        rmax = 2.0 * phycon.ANG2BOHR
        if bnd_len_key in bnd_len_dct:
            bnd_len = bnd_len_dct[bnd_len_key]
            npoints = 14
            rmin = bnd_len + 0.1 * phycon.ANG2BOHR
            rmax = bnd_len + 0.8 * phycon.ANG2BOHR
        bkp_grid = numpy.linspace(rmin, rmax, npoints)
        bkp_update_guess = False
    elif 'addition' in bkp_typ:
        rmin = 1.6 * phycon.ANG2BOHR
        rmax = 2.8 * phycon.ANG2BOHR
        if bnd_len_key in bnd_len_dct:
            npoints = 14
            bnd_len = bnd_len_dct[bnd_len_key]
            rmin = bnd_len + 0.1 * phycon.ANG2BOHR
            rmax = bnd_len + 1.4 * phycon.ANG2BOHR
        bkp_grid = numpy.linspace(rmin, rmax, npoints)
        bkp_update_guess = False

    if 'beta scission' in typ:
        rmin = 1.4 * phycon.ANG2BOHR
        rmax = 2.0 * phycon.ANG2BOHR
        if bnd_len_key in bnd_len_dct:
            npoints = 14
            bnd_len = bnd_len_dct[bnd_len_key]
            rmin = bnd_len + 0.1 * phycon.ANG2BOHR
            rmax = bnd_len + 0.8 * phycon.ANG2BOHR
        grid = numpy.linspace(rmin, rmax, npoints)
        update_guess = False

    elif 'radical radical addition' in typ:
        rstart = 2.4 * phycon.ANG2BOHR
        rend1 = 1.8 * phycon.ANG2BOHR
        rend2 = 3.0 * phycon.ANG2BOHR
        grid1 = numpy.linspace(rstart, rend1, npoints1)
        grid2 = numpy.linspace(rstart, rend2, npoints2)
        grid2 = numpy.delete(grid2, 0)
        grid = [grid1, grid2]
        update_guess = True

    elif 'addition' in typ:
        npoints = 14
        rmin = 1.6 * phycon.ANG2BOHR
        rmax = 2.8 * phycon.ANG2BOHR
        if bnd_len_key in bnd_len_dct:
            bnd_len = bnd_len_dct[bnd_len_key]
            rmin = bnd_len + 0.1 * phycon.ANG2BOHR
            rmax = bnd_len + 1.2 * phycon.ANG2BOHR
        #grid = numpy.linspace(rmin, rmax, npoints)
        gfac = 1.1
        grid = [rmin]
        rstp = 0.05
        rgrid = rmin
        for idx in range(npoints):
            rgrid += rstp
            if rgrid == rmax:
                break
            grid.append(rgrid)
            rstp = rstp * gfac
        grid = numpy.array(grid)
        #print('grid test:', grid)
        update_guess = False

    elif 'hydrogen migration' in typ:

        interval = 0.3*phycon.ANG2BOHR
        # get rmax from ts_zma
        rmax = automol.zmatrix.values(ts_zma)[dist_name]
        rmin1 = 2.*phycon.ANG2BOHR
        rmin2 = 1.3*phycon.ANG2BOHR
        if bnd_len_key in bnd_len_dct:
            bnd_len = bnd_len_dct[bnd_len_key]
            rmin2 = bnd_len + 0.05 * phycon.ANG2BOHR
            #rmin2 = bnd_len + 0.1 * phycon.ANG2BOHR
        if rmax > rmin1:
            npoints = (rmax-rmin1)/interval
            grid1 = numpy.linspace(rmax, rmin1, npoints)
        else:
            grid1 = []
        grid2 = numpy.linspace(rmin1, rmin2, 18)
        grid = numpy.concatenate((grid1, grid2), axis=None)
        update_guess = True

    elif 'elimination' in typ:

        brk_coo, = automol.zmatrix.coordinates(ts_zma)[brk_name]
        brk_len_key = tuple(sorted(map(syms.__getitem__, brk_coo)))

        interval = 0.2*phycon.ANG2BOHR
        rmin = 1.4 * phycon.ANG2BOHR
        rmax = 2.8 * phycon.ANG2BOHR
        if bnd_len_key in bnd_len_dct:
            bnd_len = bnd_len_dct[bnd_len_key]
            brk_len = bnd_len_dct[brk_len_key]
            r1min = bnd_len + 0.2 * phycon.ANG2BOHR
            r1max = bnd_len + 1.4 * phycon.ANG2BOHR
            r2min = brk_len + 0.2 * phycon.ANG2BOHR
            r2max = brk_len + 0.8 * phycon.ANG2BOHR
            grid1 = numpy.linspace(r1min, r1max, 8)
            grid2 = numpy.linspace(r2min, r2max, 4)
            grid = [grid1, grid2]
            update_guess = False

    elif 'radical radical hydrogen abstraction' in typ:
        rstart = 2.4 * phycon.ANG2BOHR
        rend1 = 1.4 * phycon.ANG2BOHR
        rend2 = 3.0 * phycon.ANG2BOHR
        grid1 = numpy.linspace(rstart, rend1, npoints1)
        grid2 = numpy.linspace(rstart, rend2, npoints2)
        grid2 = numpy.delete(grid2, 0)
        grid = [grid1, grid2]
        update_guess = True

    elif 'hydrogen abstraction' in typ:
        npoints = 16
        rmin = 0.7 * phycon.ANG2BOHR
        rmax = 2.2 * phycon.ANG2BOHR
        if bnd_len_key in bnd_len_dct:
            bnd_len = bnd_len_dct[bnd_len_key]
            rmin = bnd_len
            rmax = bnd_len + 1.0 * phycon.ANG2BOHR
        grid = numpy.linspace(rmin, rmax, npoints)
        update_guess = False

    elif 'substitution' in typ:
        npoints = 14
        rmin = 0.7 * phycon.ANG2BOHR
        rmax = 2.4 * phycon.ANG2BOHR
        if bnd_len_key in bnd_len_dct:
            bnd_len = bnd_len_dct[bnd_len_key]
            rmin = bnd_len
            rmax = bnd_len + 1.4 * phycon.ANG2BOHR
        grid = numpy.linspace(rmin, rmax, npoints)
        update_guess = False

    elif 'insertion' in typ:
        npoints = 16
        rmin = 1.4 * phycon.ANG2BOHR
        rmax = 2.4 * phycon.ANG2BOHR
        if bnd_len_key in bnd_len_dct:
            bnd_len = bnd_len_dct[bnd_len_key]
            rmin = bnd_len
            rmax = bnd_len + 1.4 * phycon.ANG2BOHR
        grid = numpy.linspace(rmin, rmax, npoints)
        update_guess = False

    elif 'radical radical' in typ:
        grid = None
        update_guess = True

    if typ:
        ts_class_data = [
            typ, ts_zma, dist_name, brk_name, grid, frm_bnd_key, brk_bnd_key,
            tors_names, update_guess]
    else:
        ts_class_data = []
    if bkp_typ:
        bkp_ts_class_data = [
            bkp_typ, bkp_ts_zma, bkp_dist_name, bkp_grid, bkp_tors_names, bkp_update_guess]
    else:
        bkp_ts_class_data = []

    return ts_class_data, bkp_ts_class_data


def find_ts(
        spc_dct, ts_dct, ts_info, ts_zma, typ, dist_info, grid,
        bkp_ts_class_data, ini_thy_info, thy_info, run_prefix, save_prefix,
        rxn_run_path, rxn_save_path, overwrite, attempt=1):
    """ find the ts geometry
    """
    print('prepping ts scan for:', typ)

    _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(*thy_info[0:2], saddle=True)

    orb_restr = moldr.util.orbital_restriction(ts_info, thy_info)
    ref_level = thy_info[0:3]
    ref_level.append(orb_restr)

    thy_run_fs = autofile.fs.theory(rxn_run_path)
    thy_run_fs.leaf.create(ref_level[1:4])
    thy_run_path = thy_run_fs.leaf.path(ref_level[1:4])

    thy_save_fs = autofile.fs.theory(rxn_save_path)
    thy_save_fs.leaf.create(ref_level[1:4])
    thy_save_path = thy_save_fs.leaf.path(ref_level[1:4])

    scn_run_fs = autofile.fs.scan(thy_run_path)
    scn_save_fs = autofile.fs.scan(thy_save_path)

    ts_run_fs = autofile.fs.ts(thy_run_path)
    ts_run_fs.trunk.create()
    ts_run_path = ts_run_fs.trunk.path()
    run_fs = autofile.fs.run(ts_run_path)

    ts_save_fs = autofile.fs.ts(thy_save_path)
    ts_save_fs.trunk.create()
    ts_save_path = ts_save_fs.trunk.path()

    cnf_run_fs = autofile.fs.conformer(ts_run_path)
    cnf_save_fs = autofile.fs.conformer(ts_save_path)
    cnf_save_fs.trunk.create()

    dist_name = dist_info[0]
    update_guess = dist_info[2]
    brk_name = dist_info[3]

    # Check if TS already is found, and determine if it fits original guess
    min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
    # added a check for the presence ts in run directory but not in save directory,
    # in case the save had to be removed for some reason
    # once the code is working cleanly this should not be needed
    #if not min_cnf_locs:
        #opt_ret = moldr.driver.read_job(
        #    job='optimization',
        #    run_fs=run_fs,
        #)
        #if opt_ret is not None:
        #    inf_obj, _, out_str = opt_ret
        #    prog = inf_obj.prog
        #    method = inf_obj.method
        #    ene = elstruct.reader.energy(prog, method, out_str)
        #    geo = elstruct.reader.opt_geometry(prog, out_str)
        #    zma = elstruct.reader.opt_zmatrix(prog, out_str)
#
#            print(" - Saving...")
#            print(" - Save path: {}".format(ts_save_path))
#
#            ts_save_fs.trunk.file.energy.write(ene)
#            ts_save_fs.trunk.file.geometry.write(geo)
#            ts_save_fs.trunk.file.zmatrix.write(zma)
#
#            vals = automol.zmatrix.values(zma)
#            final_dist = vals[dist_name]
#            dist_info[1] = final_dist
#            # run_single_conformer(ts_info, ref_level, fs, overwrite, saddle=True, dist_info=dist_info)
#            moldr.conformer.save_conformers(
#                cnf_run_fs, cnf_save_fs, saddle=True, dist_info=dist_info)
#            min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
#            print('min_cnf_locs test in run_save:', min_cnf_locs)
    # end of run path checkingA
    #
    # check to see if rxn class for already found ts is of expected class
    # do this by comparing names
    if min_cnf_locs and not overwrite:
        cnf_path = cnf_save_fs.trunk.path()
        print('Found TS at {}'.format(cnf_path))
        geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        zma = cnf_save_fs.leaf.file.zmatrix.read(min_cnf_locs)
        chk_bkp = False
        #print('z-matrix test:')
        #print(automol.zmatrix.string(zma))
        #print(automol.zmatrix.string(ts_zma))
        if automol.zmatrix.names(zma) == automol.zmatrix.names(ts_zma):
            # print('zmatrix and ts-zmatrix are not identical')
            if not automol.zmatrix.almost_equal(zma, ts_zma, 4e-1, True):
                # print('zmatrix and ts-zmatrix are not nearly identical')
                if 'babs1' in automol.zmatrix.names(ts_zma):
                    babs1 = 170. * phycon.DEG2RAD
                    if automol.zmatrix.values(ts_zma)['babs1'] == babs1:
                        babs1 = 85. * phycon.DEG2RAD
                    ts_zma = automol.zmatrix.set_value(
                        ts_zma, {'babs1': babs1})
                    ts_dct['original_zma'] = ts_zma
                    if not automol.zmatrix.almost_equal(zma, ts_zma, 4e-1):
                        #print('check true in 3:')
                        chk_bkp = True
                else:
                    #print('check true in 2:')
                    chk_bkp = True
        else:
            #print('check true in 1:')
            chk_bkp = True

        is_bkp = False
        #print('bkp test:', chk_bkp, bkp_ts_class_data, is_bkp)
        if chk_bkp and bkp_ts_class_data:
            [bkp_typ, bkp_ts_zma, bkp_dist_name, bkp_grid, bkp_tors_names,
             bkp_update_guess] = bkp_ts_class_data
            if automol.zmatrix.names(zma) == automol.zmatrix.names(bkp_ts_zma):
                if automol.zmatrix.almost_equal(zma, bkp_ts_zma, 4e-1, True):
                    is_bkp = True
                elif 'babs1' in automol.zmatrix.names(bkp_ts_zma):
                    babs1 = 170. * phycon.DEG2RAD
                    if automol.zmatrix.values(bkp_ts_zma)['babs1'] == babs1:
                        babs1 = 85. * phycon.DEG2RAD
                    bkp_ts_zma = automol.zmatrix.set_value(bkp_ts_zma, {'babs1': babs1})
                    if not automol.zmatrix.almost_equal(zma, bkp_ts_zma, 4e-1):
                        is_bkp = True
        if not chk_bkp:
            print("TS is type {}".format(typ))
        elif is_bkp:
            print('updating reaction class to {}'.format(bkp_typ))
            ts_dct['class'] = bkp_typ
            ts_dct['original_zma'] = bkp_ts_zma
            bkp_dist_info = [bkp_dist_name, 0., bkp_update_guess]
            ts_dct['dist_info'] = bkp_dist_info
            ts_dct['tors_names'] = bkp_tors_names
            print("TS is backup type {}".format(bkp_typ))
        else:
            print("TS may not be original type or backup type")
            print("Some part of the z-matrices have changed")
        print('class test:', ts_dct['class'])
        vals = automol.zmatrix.values(zma)
        final_dist = vals[dist_name]
        dist_info[1] = final_dist

    # Find TS
    else:
        fs = [None, None, ts_run_fs, ts_save_fs,
              cnf_run_fs, cnf_save_fs, None, None,
              scn_run_fs, scn_save_fs, run_fs]

        print('running ts scan:')
        if 'radical radical addition' in typ or 'radical radical hydrogen abstraction' in typ:
            # run mep scan
            multi_info = ['molpro2015', 'caspt2', 'cc-pvdz', 'RR']

            orb_restr = moldr.util.orbital_restriction(ts_info, multi_info)
            multi_level = multi_info[0:3]
            multi_level.append(orb_restr)

            thy_run_fs = autofile.fs.theory(rxn_run_path)
            thy_run_fs.leaf.create(multi_level[1:4])
            thy_run_path = thy_run_fs.leaf.path(multi_level[1:4])

            thy_save_fs = autofile.fs.theory(rxn_save_path)
            thy_save_fs.leaf.create(multi_level[1:4])
            thy_save_path = thy_save_fs.leaf.path(multi_level[1:4])

            scn_run_fs = autofile.fs.scan(thy_run_path)
            scn_save_fs = autofile.fs.scan(thy_save_path)

            ts_formula = automol.geom.formula(automol.zmatrix.geometry(ts_zma))
            grid1 = grid[0]
            grid2 = grid[1]
            grid = numpy.append(grid[0], grid[1])
            high_mul = ts_dct['high_mul']
            print('starting multiref scan:', scn_run_fs.trunk.path())
            vtst = True
            if vtst:
                gradient = False
                hessian = True
            moldr.scan.run_multiref_rscan(
                formula=ts_formula,
                high_mul=high_mul,
                zma=ts_zma,
                spc_info=ts_info,
                multi_level=multi_level,
                dist_name=dist_name,
                grid1=grid1,
                grid2=grid2,
                scn_run_fs=scn_run_fs,
                scn_save_fs=scn_save_fs,
                script_str=opt_script_str,
                overwrite=overwrite,
                update_guess=update_guess,
                gradient=gradient,
                hessian=hessian,
                **opt_kwargs
            )

            moldr.scan.save_scan(
                scn_run_fs=scn_run_fs,
                scn_save_fs=scn_save_fs,
                coo_names=[dist_name],
                gradient=gradient,
                hessian=hessian,
            )

            locs = [[dist_name], [grid1[0]]]
            # calculate and save the infinite seperation energy
            ts_zma = scn_save_fs.leaf.file.zmatrix.read(locs)
            rcts = ts_dct['reacs']
            spc_1_info = [spc_dct[rcts[0]]['ich'], spc_dct[rcts[0]]['chg'], spc_dct[rcts[0]]['mul']]
            spc_2_info = [spc_dct[rcts[1]]['ich'], spc_dct[rcts[1]]['chg'], spc_dct[rcts[1]]['mul']]

            inf_sep_ene = moldr.scan.infinite_separation_energy(
                spc_1_info, spc_2_info, ts_info, high_mul, ts_zma, ini_thy_info, thy_info,
                multi_info, run_prefix, save_prefix, scn_run_fs, scn_save_fs, locs)

            inf_locs = [[dist_name], [1000.]]
            scn_save_fs.leaf.create(inf_locs)
            scn_save_fs.leaf.file.energy.write(inf_sep_ene, inf_locs)

            geo = automol.zmatrix.geometry(ts_zma)
            zma = ts_zma
            final_dist = grid1[0]

            vrctst = False
            if vrctst:

                # continue on to finish setting up the correction potential


                # now call vrctst which sets up all the vrctst files
                input_strs = moldr.vrctst.input_prep(ts_zma, dist_name)
                [struct_inp_str, lr_divsur_inp_str, sr_divsur_inp_str, tst_inp_str,
                 els_inp_str, mc_flux_inp_str, conv_inp_str] = input_strs

                # generate the molpro template file
                      # Write the *.tml input string
                memory = 4.0
                basis = 'cc-pvdz'
                num_act_elc = high_mul
                num_act_orb = num_act_elc

                ts_formula = automol.geom.formula(automol.zmatrix.geometry(ts_zma))

                _, wfn_str = moldr.ts.cas_options(
                    ts_info, ts_formula, num_act_elc, num_act_orb, high_mul)
                method = '{rs2c, shift=0.25}'
                # inf_sep_energy = -78.137635
                tml_inp_str = varecof_io.writer.input_file.tml(
                    memory, basis, wfn_str, method, inf_sep_ene)
                print('\n\nmol.tml:')
                print(tml_inp_str)

                vrc_path = os.path.join(os.getcwd(), 'vrc')
                scr_path = os.path.join(vrc_path, 'scratch')
                os.makedirs(vrc_path, exist_ok=True)
                os.makedirs(scr_path, exist_ok=True)
                machines = ['b450:8', 'b451:8', 'b452:8', 'b453:8']

                print('vrc_path test:', vrc_path)
                with open(os.path.join(vrc_path, 'machines'), 'w') as machine_file:
                    for machine in machines:
                        machine_file.writelines(machine + '\n')
                with open(os.path.join(vrc_path, 'structure.inp'), 'w') as inp_file:
                    inp_file.write(struct_inp_str)
                with open(os.path.join(vrc_path, 'sr_divsur.inp'), 'w') as inp_file:
                    inp_file.write(sr_divsur_inp_str)
                with open(os.path.join(vrc_path, 'lr_divsur.inp'), 'w') as inp_file:
                    inp_file.write(lr_divsur_inp_str)
                with open(os.path.join(vrc_path, 'tst.inp'), 'w') as inp_file:
                    inp_file.write(tst_inp_str)
                with open(os.path.join(vrc_path, 'molpro.inp'), 'w') as inp_file:
                    inp_file.write(els_inp_str)
                with open(os.path.join(vrc_path, 'mc_flux.inp'), 'w') as inp_file:
                    inp_file.write(mc_flux_inp_str)
                with open(os.path.join(vrc_path, 'convert.inp'), 'w') as inp_file:
                    inp_file.write(conv_inp_str)
                with open(os.path.join(vrc_path, 'mol.tml'), 'w') as tml_file:
                    tml_file.write(tml_inp_str)

                geo = automol.zmatrix.geometry(ts_zma)
                zma = ts_zma
                final_dist = grid1[0]

        else:
            if 'elimination' in typ:
                grid1, grid2 = grid
                grid_dct = {dist_name: grid1, brk_name: grid2}
            else:
                grid_dct = {dist_name: grid}
            moldr.scan.run_scan(
                zma=ts_zma,
                spc_info=ts_info,
                thy_level=ref_level,
                grid_dct=grid_dct,
                scn_run_fs=scn_run_fs,
                scn_save_fs=scn_save_fs,
                script_str=opt_script_str,
                saddle=False,
                overwrite=overwrite,
                update_guess=update_guess,
                reverse_sweep=False,
                fix_failures=False,
                **opt_kwargs,
                )
            if 'elimination' in typ:
                moldr.scan.save_scan(
                    scn_run_fs=scn_run_fs,
                    scn_save_fs=scn_save_fs,
                    coo_names=[dist_name, brk_name],
                    )
            else:    
                moldr.scan.save_scan(
                    scn_run_fs=scn_run_fs,
                    scn_save_fs=scn_save_fs,
                    coo_names=[dist_name],
                    )

            print('typ test:', typ)
            if 'elimination' in typ:
                enes_lst = []
                locs_lst_lst = []
                for grid_val_j in grid2:
                    locs_list = []
                    for grid_val_i in grid1:
                        locs_list.append([[dist_name, brk_name], [grid_val_i, grid_val_j]])
                    print('locs_lst', locs_list)    
                    enes = []
                    locs_lst = []
                    for locs in locs_list:
                        if scn_save_fs.leaf.exists(locs):
                            print(scn_save_fs.leaf.path(locs))
                            enes.append(scn_save_fs.leaf.file.energy.read(locs))
                            locs_lst.append(locs)
                    locs_lst_lst.append(locs_lst)
                    enes_lst.append(enes)
                    print('enes_lst', enes_lst)    
                max_enes = []  
                max_locs = []
                for idx_j, enes in enumerate(enes_lst):
                    max_ene = -10000.
                    max_loc = ''
                    for idx_i, ene in enumerate(enes):
                        if ene > max_ene:
                            max_ene = ene
                            max_loc = locs_lst_lst[idx_j][idx_i]
                            print('new max', max_ene, max_loc)    
                    max_enes.append(max_ene)
                    max_locs.append(max_loc)
                min_ene = 10000.    
                locs = []
                for idx_j, ene in enumerate(max_enes):
                    if ene < min_ene:
                        min_ene = ene
                        locs = max_locs[idx_j]
                max_locs = locs
                max_ene = min_ene
                print('min max loc', max_ene, max_locs)
                print('min max loc', scn_save_fs.leaf.path(max_locs))
                max_zma = scn_save_fs.leaf.file.zmatrix.read(max_locs)
            else:
                locs_list = []
                locs_lst = []
                enes = []
                for grid_val_i in grid:
                    locs_list.append([[dist_name], [grid_val_i]])
                for locs in locs_list:
                    if scn_save_fs.leaf.exists(locs):
                        enes.append(scn_save_fs.leaf.file.energy.read(locs))
                        locs_lst.append(locs)
                max_ene = max(enes)
                max_idx = enes.index(max(enes))
                if 'migration' in typ:
                    max_grid_val = grid[max_idx]
                    max_zma = automol.zmatrix.set_values(
                        ts_zma, {dist_name: max_grid_val})
                else:
                    max_locs = locs_lst[max_idx]
                    max_zma = scn_save_fs.leaf.file.zmatrix.read(max_locs)
                print('enes test in vtst:', enes)

            print('geometry for maximum along scan:', max_zma)
            print('energy for maximum along scan:', max_ene)
            print('optimizing ts')
            # find saddlepoint from maximum on the grid opt scan

            moldr.driver.run_job(
                job='optimization',
                script_str=opt_script_str,
                run_fs=run_fs,
                geom=max_zma,
                spc_info=ts_info,
                thy_level=ref_level,
                saddle=True,
                overwrite=overwrite,
                **opt_kwargs,
                )

            opt_ret = moldr.driver.read_job(
                job='optimization',
                run_fs=run_fs,
            )
            if opt_ret is not None:
                inf_obj, _, out_str = opt_ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)
                geo = elstruct.reader.opt_geometry(prog, out_str)
                zma = elstruct.reader.opt_zmatrix(prog, out_str)

                print(" - Saving...")
                print(" - Save path: {}".format(ts_save_path))

                ts_save_fs.trunk.file.energy.write(ene)
                ts_save_fs.trunk.file.geometry.write(geo)
                ts_save_fs.trunk.file.zmatrix.write(zma)

                vals = automol.zmatrix.values(zma)
                final_dist = vals[dist_name]
                dist_info[1] = final_dist
                print('Test final distance for reactant coordinate', final_dist)
                run_single_conformer(ts_info, ref_level, fs, overwrite, saddle=True, dist_info=dist_info)

            elif 'addition' in typ and bkp_ts_class_data and attempt > 2:
                bkp_typ, bkp_ts_zma, bkp_dist_name, bkp_grid, bkp_tors_names, bkp_update_guess = bkp_ts_class_data
                print('TS find failed. Attempting to find with new reaction class: {}'.format(bkp_typ))
                bkp_dist_info = [bkp_dist_name, 0., bkp_update_guess]
                ts_dct['class'] = bkp_typ
                ts_dct['original_zma'] = bkp_ts_zma
                ts_dct['dist_info'] = bkp_dist_info
                ts_dct['tors_names'] = bkp_tors_names
                attempt += 1
                geo, zma, final_dist = find_ts(
                    spc_dct, ts_dct, ts_info, bkp_ts_zma, bkp_typ, bkp_dist_info,
                    bkp_grid, None, ini_thy_info, thy_info, run_prefix,
                    save_prefix, rxn_run_path, rxn_save_path, overwrite=True,
                    attempt=attempt)
            elif ('addition ' in typ or 'abstraction' in typ) and attempt < 3:
                babs1 = 170. * phycon.DEG2RAD
                if automol.zmatrix.values(ts_zma)['babs1'] == babs1:
                    babs1 = 85. * phycon.DEG2RAD
                print('TS find failed. Attempting to find with new angle of attack: {:.1f}'.format(babs1))
                ts_zma = automol.zmatrix.set_value(ts_zma, {'babs1': babs1})
                ts_dct['original_zma'] = ts_zma
                attempt += 1
                geo, zma, final_dist = find_ts(
                    spc_dct, ts_dct, ts_info, ts_zma, typ, dist_info, grid,
                    bkp_ts_class_data, ini_thy_info, thy_info, run_prefix,
                    save_prefix, rxn_run_path, rxn_save_path, overwrite=True,
                    attempt=attempt)
            elif 'beta scission' in typ and bkp_ts_class_data and attempt < 2:
                [bkp_typ, bkp_ts_zma, bkp_dist_name, bkp_grid, bkp_tors_names, bkp_update_guess] = bkp_ts_class_data
                print('TS find failed. Attempting to find with new reaction class: {}'.format(bkp_typ))
                bkp_dist_info = [bkp_dist_name, 0., bkp_update_guess]
                ts_dct['class'] = bkp_typ
                ts_dct['original_zma'] = bkp_ts_zma
                ts_dct['dist_info'] = bkp_dist_info
                ts_dct['tors_names'] = bkp_tors_names
                attempt += 1
                geo, zma, final_dist = find_ts(
                    spc_dct, ts_dct, ts_info, bkp_ts_zma, bkp_typ, bkp_dist_info,
                    bkp_grid, None, ini_thy_info, thy_info, run_prefix,
                    save_prefix, rxn_run_path, rxn_save_path, overwrite=True,
                    attempt=attempt)
            else:
                geo = 'failed'
                zma = 'failed'
                final_dist = 0.
    return geo, zma, final_dist


#def variational_data():
#    """ Perform the calculations to do variational calculations: VTST and VRCTST
#    """
#
#    if rxn_has_sadpt:
#        if method == 'vtst':
#            run_irc()
#    else:
#        if method == 'vtst':
#            #do vtst stuff
#        elif method == 'vrctst':
#            run_vrctst()


def find_vdw(ts_name, spc_dct, thy_info, ini_thy_info, vdw_params,
             nsamp_par, run_prefix, save_prefix, kickoff_size, kickoff_backward,
             projrot_script_str, overwrite):
    """ find van der Waals structures for all the pairs of species in a reaction list
    """
    new_vdws = []
    _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(*thy_info[:2])
    mul = spc_dct[ts_name]['low_mul']
    vdw_names_lst = []
    if vdw_params[0]:
        vdw_names_lst.append([sorted(spc_dct[ts_name]['reacs']), mul, 'r'])
    if vdw_params[1]:
        vdw_names_lst.append([sorted(spc_dct[ts_name]['prods']), mul, 'p'])

    for names, ts_mul, label in vdw_names_lst:
        if len(names) < 2:
            print("Cannot find van der Waals well for unimolecular reactant or product")
        ichs = list(map(lambda name: spc_dct[name]['ich'], names))
        chgs = list(map(lambda name: spc_dct[name]['chg'], names))
        muls = list(map(lambda name: spc_dct[name]['mul'], names))

        # theory
        prog = thy_info[0]
        method = thy_info[1]
        _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(prog, method)

        geos = []
        ntaudof = 0.
        for name, ich, chg, mul in zip(names, ichs, chgs, muls):
            spc_info = [ich, chg, mul]
            orb_restr = moldr.util.orbital_restriction(spc_info, ini_thy_info)
            ini_thy_level = ini_thy_info[0:3]
            ini_thy_level.append(orb_restr)
            orb_restr = moldr.util.orbital_restriction(spc_info, thy_info)
            thy_level = thy_info[0:3]
            thy_level.append(orb_restr)
            spc_run_fs = autofile.fs.species(run_prefix)
            spc_run_fs.leaf.create(spc_info)
            spc_run_path = spc_run_fs.leaf.path(spc_info)
            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs.leaf.create(spc_info)
            spc_save_path = spc_save_fs.leaf.path(spc_info)

            thy_run_fs = autofile.fs.theory(spc_run_path)
            thy_run_fs.leaf.create(thy_level[1:4])
            thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])
            thy_save_fs = autofile.fs.theory(spc_save_path)
            thy_save_fs.leaf.create(thy_level[1:4])
            thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])
            run_fs = autofile.fs.run(thy_run_path)

            ini_thy_save_fs = autofile.fs.theory(spc_save_path)
            ini_thy_save_fs.leaf.create(ini_thy_level[1:4])

            cnf_run_fs = autofile.fs.conformer(thy_run_path)
            cnf_save_fs = autofile.fs.conformer(thy_save_path)

            ini_fs = [None, ini_thy_save_fs]
            fs = [spc_run_fs, spc_save_fs, thy_run_fs, thy_save_fs,
                  cnf_run_fs, cnf_save_fs, None, None,
                  None, None, run_fs]
    # fs = [None, None, thy_run_fs, thy_save_fs,
          # cnf_run_fs, cnf_save_fs, None, None,
            geo = moldr.geom.reference_geometry(
                spc_dct[name], thy_level, ini_thy_level, fs, ini_fs,
                kickoff_size=kickoff_size,
                kickoff_backward=kickoff_backward,
                projrot_script_str=projrot_script_str,
                overwrite=overwrite)
            geos.append(geo)
            gra = automol.geom.graph(geo)
            ntaudof += len(automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
        nsamp = moldr.util.nsamp_init(nsamp_par, ntaudof)
        geo1, geo2 = geos
        geo1 = automol.geom.mass_centered(geo1)
        geo2 = automol.geom.mass_centered(geo2)
        min_ene = 0.
        for idx in range(int(nsamp)):
            print('Optimizing vdw geometry {}/{}'.format(idx+1, nsamp))
            angs1 = numpy.multiply(
                numpy.random.rand(3), [1*numpy.pi, 2*numpy.pi, 2*numpy.pi])
            angs2 = numpy.multiply(
                numpy.random.rand(3), [1*numpy.pi, 2*numpy.pi, 2*numpy.pi])
            angs12 = numpy.multiply(
                numpy.random.rand(2), [1*numpy.pi, 2*numpy.pi])
            geo1 = automol.geom.euler_rotated(geo1, *angs1)
            geo2 = automol.geom.euler_rotated(geo2, *angs2)
            dist_cutoff = 3.0 * phycon.ANG2BOHR

            geo = automol.geom.join(geo1, geo2, dist_cutoff, *angs12)
            print("Species: {}".format('+'.join(names)))
            print('vdw starting geometry')
            print(automol.geom.xyz_string(geo))

   #  set up the filesystem
            ich = automol.inchi.recalculate(automol.inchi.join(ichs))
            chg = sum(chgs)
            mul = ts_mul
            spc_info = (ich, chg, mul)
            #orb_restr = moldr.util.orbital_restriction(mul, thy_info[0:3] restrict_open_shell)
            #orb_restr = restrict_open_shell
            spc_run_fs = autofile.fs.species(run_prefix)
            spc_run_fs.leaf.create(spc_info)
            spc_run_path = spc_run_fs.leaf.path(spc_info)
            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs.leaf.create(spc_info)
            spc_save_path = spc_save_fs.leaf.path(spc_info)
            orb_restr = moldr.util.orbital_restriction(spc_info, thy_info)
            thy_level = thy_info[0:3]
            thy_level.append(orb_restr)
            thy_run_fs = autofile.fs.theory(spc_run_path)
            thy_run_fs.leaf.create(thy_level[1:4])
            thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])
            thy_save_fs = autofile.fs.theory(spc_save_path)
            thy_save_fs.leaf.create(thy_level[1:4])
            thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])
            run_fs = autofile.fs.run(thy_run_path)
   #  generate reference geometry
   #  generate the z-matrix and sampling ranges

            moldr.driver.run_job(
                job=elstruct.Job.OPTIMIZATION,
                geom=geo,
                spc_info=spc_info,
                thy_level=thy_level,
                run_fs=run_fs,
                script_str=opt_script_str,
                overwrite=overwrite,
                **opt_kwargs,
            )

   #  save info for the initial geometry (from inchi or from save directory)
            ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if ret:
                print('Saving reference geometry')
                print(" - Save path: {}".format(thy_save_path))

                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                geo = elstruct.reader.opt_geometry(prog, out_str)
                print('vdw ending geometry')
                print(automol.geom.xyz_string(geo))
                thy_save_fs.leaf.file.geometry.write(geo, thy_level[1:4])
                ene = elstruct.reader.energy(prog, method, out_str)
                if ene < min_ene:
                    min_ene = ene
                    print('ene test in vdw')
                    print(ene)
                    thy_save_fs.leaf.file.energy.write(ene, thy_level[1:4])
                    print('Saving reference geometry')
                    print(" - Save path: {}".format(thy_save_path))
                    vdw_name = label + ts_name.replace('ts', 'vdw')
                    spc_dct[vdw_name] = spc_dct[ts_name].copy()
                    spc_dct[vdw_name]['ich'] = ich
                    spc_dct[vdw_name]['mul'] = mul
                    spc_dct[vdw_name]['chg'] = chg
                    spc_dct[vdw_name]['dist_info'][1] = dist_cutoff
                    fs = [spc_run_fs, spc_save_fs, thy_run_fs, thy_save_fs,
                          cnf_run_fs, cnf_save_fs, None, None,
                          None, None, run_fs]
                    #Make a fake conformer
                    cnf_save_fs = autofile.fs.conformer(thy_save_path)
                    cnf_run_fs = autofile.fs.conformer(thy_run_path)
                    cnf_save_fs.trunk.create()
                    cnf_run_fs.trunk.create()
                    tors_range_dct = {}
                    cinf_obj = autofile.system.info.conformer_trunk(0, tors_range_dct)
                    cinf_obj.nsamp = 1
                    cnf_save_fs.trunk.file.info.write(cinf_obj)
                    locs_lst = cnf_save_fs.leaf.existing()
                    if not locs_lst:
                        cid = autofile.system.generate_new_conformer_id()
                        locs = [cid]
                    else:
                        locs = locs_lst[0]
                    cnf_save_fs.leaf.create(locs)
                    cnf_run_fs.leaf.create(locs)
                    cnf_save_fs.leaf.file.geometry_info.write(
                        inf_obj, locs)
                    cnf_save_fs.leaf.file.geometry_input.write(
                        inp_str, locs)
                    cnf_save_fs.leaf.file.energy.write(ene, locs)
                    cnf_save_fs.leaf.file.geometry.write(geo, locs)
        if min_ene:
            new_vdws.append(vdw_name)

    return new_vdws


def fake_conf(thy_level, fs, inf=[]):
    cnf_save_fs = fs[5]
    cnf_run_fs = fs[4]
    thy_save_fs = fs[3]
    run_fs = fs[-1]
    thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])
    geo = thy_save_fs.leaf.file.geometry.read(thy_level[1:4])
    if inf:
        inf_obj, ene = inf
    else:
        ene = thy_save_fs.leaf.file.energy.read(thy_level[1:4])
        inf_obj = run_fs.trunk.file.info.read()
    tors_range_dct = {}
    cinf_obj = autofile.system.info.conformer_trunk(0, tors_range_dct)
    cinf_obj.nsamp = 1
    cnf_save_fs = autofile.fs.conformer(thy_save_path)
    cnf_save_fs.trunk.create()
    cnf_run_fs.trunk.create()
    cnf_save_fs.trunk.file.info.write(cinf_obj)
    cnf_run_fs.trunk.file.info.write(cinf_obj)
    locs_lst = cnf_save_fs.leaf.existing()
    if not locs_lst:
        cid = autofile.system.generate_new_conformer_id()
        locs = [cid]
    else:
        locs = locs_lst[0]
    cnf_save_fs.leaf.create(locs)
    cnf_run_fs.leaf.create(locs)
    cnf_save_fs.leaf.file.geometry_info.write(
        inf_obj, locs)
    cnf_run_fs.leaf.file.geometry_info.write(
        inf_obj, locs)
    method = inf_obj.method
    cnf_save_fs.leaf.file.energy.write(ene, locs)
    cnf_run_fs.leaf.file.energy.write(ene, locs)
    cnf_save_fs.leaf.file.geometry.write(geo, locs)
    cnf_run_fs.leaf.file.geometry.write(geo, locs)


def fake_geo_gen(tsk, spcdic, es_dct, thy_level, fs,
       spc_info, overwrite):
    if 'conf' in tsk:
        fake_conf(thy_level, fs)

    if 'scan' in tsk:
        pass
    if 'tau' in tsk:
        pass


def get_thy_info(lvldic):
    """ convert theory level dictionary to theory information array
    """
    err_msg = ''
    info = ['program', 'method', 'basis', 'orb_res']
    for i, inf in enumerate(info):
        if inf in lvldic:
            info[i] = lvldic[inf]
        else:
            err_msg = inf
    if err_msg:
        print('ERROR: No {} found'.format(err_msg))
    return info


def get_spc_info(spc_dct_i):
    """ convert species dictionary to species_info array
    """
    err_msg = ''
    props = ['ich', 'chg', 'mul']
    for i, prop in enumerate(props):
        if prop in spc_dct_i:
            props[i] = spc_dct_i[prop]
        else:
            err_msg = prop
    if err_msg:
        print('ERROR: No {} found'.format(err_msg))
    return props
