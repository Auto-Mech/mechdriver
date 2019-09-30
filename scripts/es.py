""" cycle over electronic structure calls
for species, TSs, and vdw species
"""
import numpy
from qcelemental import constants as qcc
import automol
import elstruct
import thermo
import autofile
import moldr

ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')
WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

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


def run_single_conformer(spc_info, thy_level, fs, overwrite, saddle=False, dist_info=[]):
    """ generate single optimized geometry for randomly sampled initial torsional angles
    """
    mc_nsamp = [False, 0, 0, 0, 0, 1]
    sp_script_str, _, kwargs, _ = moldr.util.run_qchem_par(*thy_level[0:2])
    thy_save_fs = fs[3]
    two_stage = False
    if saddle:
        two_stage=True
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
            params['scan_increment'] = 30. * qcc.conversion_factor('degree', 'radian')

    if tsk in choose_function:
        eval(choose_function[tsk])(fs, params, opt_kwargs)


def ts_geometry_generation(tsk, spcdic, es_dct, thy_level, fs, spc_info, overwrite):
    """ run an electronic structure task
    for generating a list of conformer or tau sampling geometries
    """
    # fs[3] = fs[11]
    _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(*thy_level[0:2])
    params = {'spc_info': spc_info,
              'thy_level': thy_level,
              'script_str': opt_script_str,
              'saddle' :  True,
              'tors_names': spcdic['tors_names'],
              'overwrite': overwrite}
    choose_function = {'conf_samp': 'run_conf_samp',
                       'tau_samp': 'run_tau_samp',
                       'hr_scan': 'run_hr_scan'}

    if tsk in ['conf_samp', 'tau_samp']:
        params['nsamp_par'] = es_dct['mc_nsamp']
        params['dist_info'] = spcdic['dist_info']
        params['two_stage'] = True
    elif tsk in ['hr_scan']:
        if 'hind_inc' in spcdic:
            params['scan_increment'] = spcdic['hind_inc']
        else:
            params['scan_increment'] = 30. * qcc.conversion_factor('degree', 'radian')

    if tsk in choose_function:
        eval(choose_function[tsk])(fs, params, opt_kwargs)


def geometry_analysis(tsk, thy_level, ini_fs, selection, spc_info,
        overwrite):
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


def ts_geometry_analysis(tsk, thy_level, ini_fs, selection, spc_info, overwrite):
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


def get_spc_run_path(run_prefix, spc_info):
    spc_run_fs = autofile.fs.species(run_prefix)
    spc_run_fs.leaf.create(spc_info)
    spc_run_path = spc_run_fs.leaf.path(spc_info)
    return spc_run_path


def get_spc_save_path(save_prefix, spc_info):
    spc_save_fs = autofile.fs.species(save_prefix)
    spc_save_fs.leaf.create(spc_info)
    spc_save_path = spc_save_fs.leaf.path(spc_info)
    return spc_save_path


def get_thy_save_fs(save_prefix, spc_info, thy_info):
    orb_restr = moldr.util.orbital_restriction(
        spc_info, thy_info)
    thy_lvl = thy_info[0:3]
    thy_lvl.append(orb_restr)
    spc_save_path = get_spc_save_path(save_prefix, spc_info)
    thy_save_fs = autofile.fs.theory(spc_save_path)
    return thy_save_fs, thy_lvl


def get_thy_run_path(run_prefix, spc_info, thy_info):
    orb_restr = moldr.util.orbital_restriction(
        spc_info, thy_info)
    thy_lvl = thy_info[0:3]
    thy_lvl.append(orb_restr)
    spc_run_path = get_spc_run_path(run_prefix, spc_info)
    thy_run_fs = autofile.fs.theory(spc_run_path)
    thy_run_fs.leaf.create(thy_lvl)
    thy_run_path = thy_run_fs.leaf.path(thy_lvl)
    return thy_run_path


def get_thy_save_path(save_prefix, spc_info, thy_info):
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
    rxn_ichs = [[], []]
    rxn_chgs = [[], []]
    rxn_muls = [[], []]
    print('ts test:', ts)
    print('spc_dct reacs test:',spc_dct[ts]['reacs'])
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
    print('checking exothermicity of reaction')
    try:
        rxn_exo = moldr.util.reaction_energy(
            save_prefix, rxn_ichs, rxn_chgs, rxn_muls, thy_info)
    except:
        rxn_exo = moldr.util.reaction_energy(
            save_prefix, rxn_ichs, rxn_chgs, rxn_muls, ini_thy_info)
    print('reaction is {:.2f}'.format(rxn_exo))
    if rxn_exo > 0:
        rxn_ichs = rxn_ichs[::-1]
        rxn_chgs = rxn_chgs[::-1]
        rxn_muls = rxn_muls[::-1]
        print('ts search will be performed in reverse direction')

    # set up the filesystem
    is_rev = autofile.system.reaction_is_reversed(
        rxn_ichs, rxn_chgs, rxn_muls)
    rxn_ichs, rxn_chgs, rxn_muls = autofile.system.sort_together(
        rxn_ichs, rxn_chgs, rxn_muls)
    print("The reaction direction is {}"
          .format('backward' if is_rev else 'forward'))

    low_mul = automol.mult.ts._low(rxn_muls[0])
    high_mul = automol.mult.ts._high(rxn_muls[0])

    return rxn_ichs, rxn_chgs, rxn_muls, low_mul, high_mul


def get_rxn_fs(run_prefix, save_prefix, ts):
    """get filesystems for a reaction
    """
    rxn_ichs = ts['rxn_ichs']
    rxn_chgs = ts['rxn_chgs']
    rxn_muls = ts['rxn_muls']
    # print('ts_mul test:', ts_mul)
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
    rct_geos, cnf_save_fs_lst = get_geos(
        reacs, spc_dct, ini_thy_info, save_prefix, run_prefix, kickoff_size,
        kickoff_backward, projrot_script_str)
    prd_geos, _ = get_geos(
        prods, spc_dct, ini_thy_info, save_prefix, run_prefix, kickoff_size,
        kickoff_backward, projrot_script_str)
    rct_zmas = list(map(automol.geom.zmatrix, rct_geos))
    prd_zmas = list(map(automol.geom.zmatrix, prd_geos))
    return rct_zmas, prd_zmas, cnf_save_fs_lst


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


def ts_class(rct_zmas, prd_zmas, rad_rad, rct_cnf_save_fs):
    typ = None
    ret = automol.zmatrix.ts.beta_scission(rct_zmas, prd_zmas)
    if ret and typ is None:
        typ = 'beta scission'
        ts_zma, dist_name, tors_names = ret
        print('beta scission')
        print('ts zma:', ts_zma)
        print('dist name:', dist_name)
        print('tors names:', tors_names)

    ret = automol.zmatrix.ts.addition(rct_zmas, prd_zmas)
    if ret and typ is None:
        typ = 'addition'
        ts_zma, dist_name, tors_names = ret
        print('addn')
        print('ts zma:')
        print(ts_zma)
        print('dist name:')
        print(dist_name)
        print('tors names:')
        print(tors_names)

    ret = automol.zmatrix.ts.hydrogen_migration(rct_zmas, prd_zmas)
    if ret and typ is None:
        typ = 'hydrogen migration'
        ts_zma, dist_name, tors_names = ret
        # rct_zmas = list(moldr.util.min_dist_conformer_zma(dist_name, rct_cnf_save_fs[0]))
        rct_zmas = moldr.util.min_dist_conformer_zma(dist_name, rct_cnf_save_fs[0])
        ret = automol.zmatrix.ts.hydrogen_migration(rct_zmas, prd_zmas)
        ts_zma, dist_name, tors_names = ret
        print('H migration')
        print('ts zma:', ts_zma)
        print('dist name:', dist_name)
        print('tors names:', tors_names)

    # fix this later
    # ret = automol.zmatrix.ts.hydrogen_abstraction(rct_zmas, prd_zmas,
    #                                               sigma=True)
    ret = automol.zmatrix.ts.hydrogen_abstraction(rct_zmas, prd_zmas,
                                                  sigma=False)
    if ret and typ is None:
        typ = 'hydrogen abstraction'
        ts_zma, dist_name, tors_names = ret
        print('H abs')
        print('ts zma:', ts_zma)
        print('dist name:', dist_name)
        print('tors names:', tors_names)

    if typ is None:
        print("Failed to classify reaction.")
    else:
        if rad_rad:
            typ = 'radical radical ' + typ
        print("Type: {}".format(typ))

        # determine the grid
        dist_coo, = automol.zmatrix.coordinates(ts_zma)[dist_name]
        syms = automol.zmatrix.symbols(ts_zma)
        bnd_len_key = tuple(sorted(map(syms.__getitem__, dist_coo)))

        bnd_len_dct = {
            ('C', 'C'): 1.54 * ANG2BOHR,
            ('C', 'H'): 1.09 * ANG2BOHR,
            ('H', 'H'): 0.74 * ANG2BOHR,
            ('N', 'N'): 1.45 * ANG2BOHR,
            ('O', 'O'): 1.48 * ANG2BOHR,
            ('C', 'N'): 1.47 * ANG2BOHR,
            ('C', 'O'): 1.43 * ANG2BOHR,
            ('H', 'O'): 1.20 * ANG2BOHR,
        }

        npoints = 8
        npoints1 = 4
        npoints2 = 4
        if typ in ('beta scission'):
            rmin = 1.4 * ANG2BOHR
            rmax = 2.0 * ANG2BOHR
            if bnd_len_key in bnd_len_dct:
                bnd_len = bnd_len_dct[bnd_len_key]
                rmin = bnd_len + 0.1 * ANG2BOHR
                rmax = bnd_len + 0.5 * ANG2BOHR
            grid = numpy.linspace(rmin, rmax, npoints)
            update_guess = False

        elif typ in ('addition'):
            rmin = 1.6 * ANG2BOHR
            rmax = 2.8 * ANG2BOHR
            if bnd_len_key in bnd_len_dct:
                bnd_len = bnd_len_dct[bnd_len_key]
                rmin = bnd_len + 0.2 * ANG2BOHR
                rmax = bnd_len + 1.4 * ANG2BOHR
            grid = numpy.linspace(rmin, rmax, npoints)
            update_guess = False

        elif typ == 'hydrogen migration':

            interval = 0.2*ANG2BOHR
            rmin = 1.4 * ANG2BOHR
            rmax = 2.8 * ANG2BOHR
            if bnd_len_key in bnd_len_dct:
                rmax = bnd_len_dct[bnd_len_key]
                rmin1 = 2.*ANG2BOHR
                rmin2 = 1.1*ANG2BOHR
                if rmax > rmin:
                    npoints = (rmax-rmin)/interval
                    grid1 = numpy.linspace(rmax, rmin, npoints)
                else:
                    grid1 = []
                grid2 = numpy.linspace(rmin1, rmin2, 9)
                grid = numpy.concatenate((grid1, grid2), axis=None)
                update_guess = True

        elif typ == 'hydrogen abstraction':
            rmin = 0.7 * ANG2BOHR
            rmax = 2.2 * ANG2BOHR
            if bnd_len_key in bnd_len_dct:
                bnd_len = bnd_len_dct[bnd_len_key]
                rmin = bnd_len
                rmax = bnd_len + 1.0 * ANG2BOHR
            grid = numpy.linspace(rmin, rmax, npoints)
            update_guess = False

        elif typ == 'radical radical addition':
            rstart = 2.4 * ANG2BOHR
            rend1 = 3.0 * ANG2BOHR
            rend2 = 1.8 * ANG2BOHR
            grid1 = numpy.linspace(rstart, rend1, npoints1)
            grid2 = numpy.linspace(rstart, rend2, npoints2)
            grid2 = numpy.delete(grid2, 0)
            grid = [grid1, grid2]
            update_guess = False
        elif 'radical radical' in typ:
            grid = None
            update_guess = False

        return typ, ts_zma, dist_name, grid, tors_names, update_guess


def find_ts(ts_dct, ts_info, ts_zma, typ, dist_info, grid, thy_info, rxn_run_path, rxn_save_path, overwrite):
    """ find the ts geometry
    """
    print('prepping ts scan for:', typ)
    if 'radical radical' in typ or not typ:
        print('skipping reaction because type =:', typ)
        return 'Failure', None, None

    _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(*thy_info[0:2])

    print('rxn_run_path test in find_ts:', rxn_run_path)
    print('ts_info test in find_ts:', ts_info)
    print('thy_info test in find_ts:', thy_info)
    orb_restr = moldr.util.orbital_restriction(ts_info, thy_info)
    ref_level = thy_info[0:3]
    ref_level.append(orb_restr)
    print('ref_level test in find_ts:', ref_level)

    thy_run_fs = autofile.fs.theory(rxn_run_path)
    thy_run_fs.leaf.create(ref_level[1:4])
    thy_run_path = thy_run_fs.leaf.path(ref_level[1:4])

    thy_save_fs = autofile.fs.theory(rxn_save_path)
    thy_save_fs.leaf.create(ref_level[1:4])
    thy_save_path = thy_save_fs.leaf.path(ref_level[1:4])

    scn_run_fs = autofile.fs.scan(thy_run_path)
    scn_save_fs = autofile.fs.scan(thy_save_path)

    print('thy_run_path in ts_opt:', thy_run_path)

    ts_run_fs = autofile.fs.ts(thy_run_path)
    ts_run_fs.trunk.create()
    ts_run_path = ts_run_fs.trunk.path()
    run_fs = autofile.fs.run(ts_run_path)

    print('ts_run_path:', ts_run_path)

    ts_save_fs = autofile.fs.ts(thy_save_path)
    ts_save_fs.trunk.create()
    ts_save_path = ts_save_fs.trunk.path()
    print('ts_save_path:', ts_save_path)

    cnf_run_fs = autofile.fs.conformer(ts_run_path)
    cnf_save_fs = autofile.fs.conformer(ts_save_path)

    # fs = [None, None, thy_run_fs, thy_save_fs,
          # cnf_run_fs, cnf_save_fs, None, None,
          # scn_run_fs, scn_save_fs, ts_run_fs, ts_save_fs, run_fs]

    fs = [None, None, ts_run_fs, ts_save_fs,
          cnf_run_fs, cnf_save_fs, None, None,
          scn_run_fs, scn_save_fs, run_fs]

    dist_name = dist_info[0]
    update_guess = dist_info[2]

    print('running ts scan:')
    if typ == 'radical radical addition':
        ts_formula = automol.geom.formula(automol.zmatrix.geometry(ts_zma))
        grid = numpy.append(grid[0], grid[1])
        high_mul = ts_dct['high_mul']
        moldr.scan.run_multiref_rscan(
            formula=ts_formula,
            high_mul=high_mul,
            zma=ts_zma,
            spc_info=ts_info,
            thy_level=ref_level,
            dist_name=dist_name,
            grid1=grid1,
            grid2=grid2,
            scn_run_fs=scn_run_fs,
            scn_save_fs=scn_save_fs,
            script_str=opt_script_str,
            overwrite=overwrite,
            update_guess=update_guess,
            **opt_kwargs
        )

        moldr.scan.save_scan(
            scn_run_fs=scn_run_fs,
            scn_save_fs=scn_save_fs,
            coo_names=[dist_name],
        )

        nsamp_max = 2000
        nsamp_min = 500
        flux_err = 5
        pes_size = 1
        tst_inp_str = varecof_io.writer.write_tst_input(
            nsamp_max, nsamp_min, flux_err, pes_size)

        print('\ntst.inp:')
        print(tst_inp_str)
        # Write the divsur input file string; distances in Angstrom
        #distances = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
        #divsur_inp_str = varecof_io.writer.write_divsur_input(
        #    distances)
        #print('\ndivsur.inp:')
        #print(divsur_inp_str)
    else:
        print('grid dct test:', dist_name, grid)
        print('scan fs test:', scn_run_fs.trunk.path())
        moldr.scan.run_scan(
            zma=ts_zma,
            spc_info=ts_info,
            thy_level=ref_level,
            grid_dct={dist_name: grid},
            scn_run_fs=scn_run_fs,
            scn_save_fs=scn_save_fs,
            script_str=opt_script_str,
            overwrite=overwrite,
            update_guess=update_guess,
            reverse_sweep=False,
            **opt_kwargs
        )

    moldr.scan.save_scan(
        scn_run_fs=scn_run_fs,
        scn_save_fs=scn_save_fs,
        coo_names=[dist_name],
    )

    locs_lst = [
        locs for locs in scn_save_fs.leaf.existing([[dist_name]])
        if scn_save_fs.leaf.file.energy.exists(locs)]
    enes = [scn_save_fs.leaf.file.energy.read(locs)
            for locs in locs_lst]
    max_locs = locs_lst[enes.index(max(enes))]
    max_ene = max(enes)
    max_zma = scn_save_fs.leaf.file.zmatrix.read(max_locs)
    print('geometry for maximum along scan:', max_zma)
    print('energy for maximum along scan:', max_ene)

    print('optimizing ts')
    # find saddlepoint from maximum on the grid opt scan

    print('starting ts optimization')
    print('theory_level=:', thy_info)
    print('ts_run_path=:', ts_run_path)
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

    else:
        geo = 'failed'
        zma = 'failed'

    return geo, zma, final_dist


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
            dist_cutoff = 3.*qcc.conversion_factor('angstrom', 'bohr')

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
                    print(vdw_name)
                    print(ich)
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

#############
def species_qchem(
        spc_names, spc_info, run_opt_levels, ref_high_level,
        run_high_levels, geom_dct, run_prefix, save_prefix, qchem_flags,
        nsamp_pars, scan_increment, kickoff_pars, overwrite,
        ):
    """ run specified electronic structure calls for the given set of species
    """
    idx = 0
    run_ini_geom = qchem_flags[idx]
    idx += 1
    run_remove_imag = qchem_flags[idx]
    idx += 1
    run_conf_samp = qchem_flags[idx]
    idx += 1
    run_min_grad = qchem_flags[idx]
    idx += 1
    run_min_hess = qchem_flags[idx]
    idx += 1
    run_min_vpt2 = qchem_flags[idx]
    idx += 1
    run_conf_scan = qchem_flags[idx]
    idx += 1
    run_conf_grad = qchem_flags[idx]
    idx += 1
    run_conf_hess = qchem_flags[idx]
    idx += 1
    run_tau_samp = qchem_flags[idx]
    idx += 1
    run_tau_grad = qchem_flags[idx]
    idx += 1
    run_tau_hess = qchem_flags[idx]
    idx += 1
    run_hl_min_ene = qchem_flags[idx]

    if run_ini_geom:
        print('The initial geometries will be checked for imaginary frequencies')

    if run_remove_imag:
        print('If there is an imaginary frequency for any species it will be removed by kicking off from the saddle point')

    if run_conf_samp:
        nsamp_conf_par = nsamp_pars[0]
        print('The optimal conformers will be found through Monte Carlo sampling of torsions')

        if run_min_grad:
            print('The gradient will be determined for the minimum energy conformers')

        if run_min_hess:
            print('The hessian will be determined for the minimum energy conformers')

        if run_min_vpt2:
            print('Second order vibrational perturbation theory will be performed for the minimum energy conformers')

        if run_conf_scan:
            print('One-dimensional torsional scans will be performed starting from the minimum energy conformers')

        if run_conf_grad:
            print('The gradient will be determined for each point along the one-dimensional torsional scans')

        if run_conf_hess:
            print('The hessian will be determined for each point along the one-dimensional torsional scans')

    if run_tau_samp:
        print('Random sampling of the torsional geometries will be performed')
        nsamp_tau_par = nsamp_pars[1]

        if run_tau_grad:
            print('The gradient will be determined for each randomly sampled torsional geometry')

        if run_tau_hess:
            print('The hessian will be determined for each randomly sampled torsional geometry')

    if run_hl_min_ene:
        print('Higher level energies will be evaluated at the minimum conformer geometry')


    kickoff_backward, kickoff_size = kickoff_pars

    for name in spc_names:
        # species
        print("Species: {}".format(name))
        ich = spc_info[name][0]
        smi = automol.inchi.smiles(ich)
        print("smiles: {}".format(smi), "inchi: {}".format(ich))

        for opt_level_idx, _ in enumerate(run_opt_levels):
            # theory
            prog = run_opt_levels[opt_level_idx][0]
            method = run_opt_levels[opt_level_idx][1]
            SP_SCRIPT_STR, OPT_SCRIPT_STR, KWARGS, OPT_KWARGS = moldr.util.run_qchem_par(prog, method)

            # set up the file systems
            orb_restr = moldr.util.orbital_restriction(
                spc_info[name], run_opt_levels[opt_level_idx])
            thy_level = run_opt_levels[opt_level_idx][0:3]
            thy_level.append(orb_restr)

            spc_run_fs = autofile.fs.species(run_prefix)
            spc_run_fs.leaf.create(spc_info[name])
            spc_run_path = spc_run_fs.leaf.path(spc_info[name])

            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs.leaf.create(spc_info[name])
            spc_save_path = spc_save_fs.leaf.path(spc_info[name])

            thy_run_fs = autofile.fs.theory(spc_run_path)
            thy_run_fs.leaf.create(thy_level[1:4])
            thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])

            thy_save_fs = autofile.fs.theory(spc_save_path)
            thy_save_fs.leaf.create(thy_level[1:4])
            thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])

            run_fs = autofile.fs.run(thy_run_path)

            cnf_run_fs = autofile.fs.conformer(thy_run_path)
            cnf_save_fs = autofile.fs.conformer(thy_save_path)

            tau_run_fs = autofile.fs.tau(thy_run_path)
            tau_save_fs = autofile.fs.tau(thy_save_path)

            # this uses theory run path - should start with a check in save path to see if initial geometry has already been saved
            # eventually theory data will be removed
            # also may need to remove hessian etc from saved geometry ...
            if run_ini_geom:                                                        

                geo_init = moldr.util.reference_geometry(
                    spc_info=spc_info[name],
                    thy_level=thy_level,
                    thy_fs=thy_save_fs,
                    geom_dct=geom_dct)

#                    geo_init = moldr.util.reference_geometry(
#                    spc_info=spc_info[name],
#                    theory_level=run_opt_levels[opt_level_idx],
#                    prefix=save_prefix,
#                    geom_dct=geom_dct)

                geo = moldr.geom.run_initial_geometry_opt(
                    spc_info=spc_info[name],
                    thy_level=thy_level,
                    run_fs=run_fs,
                    thy_run_fs=thy_run_fs,
                    thy_save_fs=thy_save_fs,
                    script_str=opt_script_str,
                    overwrite=overwrite,
                    geo_init=geo_init,
                    **OPT_KWARGS,
                )

#                geo = moldr.driver.run_initial_geometry_opt(
#                    spc_info=spc_info[name],
#                    theory_level=run_opt_levels[opt_level_idx],
#                    run_prefix=spc_run_path,
#                    save_prefix=spc_save_path,
#                    script_str=OPT_SCRIPT_STR,
#                    overwrite=overwrite,
#                    geo_init=geo_init,
#                    **OPT_KWARGS,
#                )

            if run_remove_imag:
                imag, geo, disp_xyzs = moldr.geom.run_check_imaginary(
                    spc_info=spc_info[name],
                    thy_level=thy_level,
                    thy_run_fs=thy_run_fs,
                    thy_save_fs=thy_save_fs,
                    script_str=opt_script_str,
                    overwrite=overwrite,
                    **KWARGS,
                )
#                imag, geo, disp_xyzs = moldr.driver.run_check_imaginary(
#                    spc_info=spc_info[name],
#                    theory_level=run_opt_levels[opt_level_idx],
#                    run_prefix=spc_run_path,
#                    save_prefix=spc_save_path,
#                    script_str=OPT_SCRIPT_STR,
#                    overwrite=overwrite,
#                    **KWARGS,
#                )
                if imag:
                    moldr.geom.run_kickoff_saddle(
                        geo, disp_xyzs,
                        spc_info=spc_info[name],
                        thy_level=thy_level,
                        run_fs=run_fs,
                        thy_run_fs=thy_run_fs,
                        script_str=OPT_SCRIPT_STR,
                        kickoff_size=kickoff_size,
                        kickoff_backward=kickoff_backward,
                        opt_cart=False,
                        **OPT_KWARGS)
#                    moldr.driver.run_kickoff_saddle(
#                        geo, disp_xyzs,
#                        spc_info=spc_info[name],
#                        theory_level=run_opt_levels[opt_level_idx],
#                        run_path=thy_run_path,
#                        script_str=OPT_SCRIPT_STR,
#                        kickoff_backward=kickoff_backward,
#                        kickoff_size=kickoff_size,
#                        opt_cart=False,
#                        **OPT_KWARGS)
                    print('removing saddlepoint hessian')

                    run_fs.leaf.remove([elstruct.Job.HESSIAN])
#                    save_fs.leaf.remove('hess')

                    moldr.geom.save_initial_geometry(
                        spc_info=spc_info[name],
                        thy_level=thy_level,
                        run_fs=run_fs,
                        thy_run_fs=thy_run_fs,
                        thy_save_fs=thy_save_fs,
                    )

#                    run_fs = autofile.fs.run(thy_run_path)
#                    run_fs.leaf.remove([elstruct.Job.HESSIAN])
#                    save_fs = autofile.fs.save(thy_save_path)
#                    save_fs.leaf.remove('hess')
#
#                    moldr.driver.save_initial_geometry(
#                        spc_info=spc_info[name],
#                        theory_level=run_opt_levels[opt_level_idx],
#                        run_prefix=spc_run_path,
#                        save_prefix=spc_save_path,
#                    )

            if run_conf_samp:
                moldr.conformer.conformer_sampling(
                    spc_info=spc_info[name],
                    thy_level=thy_level,
                    thy_save_fs=thy_save_fs,
                    cnf_run_fs=cnf_run_fs,
                    cnf_save_fs=cnf_save_fs,
                    script_str=OPT_SCRIPT_STR,
                    overwrite=overwrite,
                    nsamp_par=nsamp_conf_par,
                    **OPT_KWARGS,
                )

            if run_min_grad:
                moldr.sp.run_minimum_energy_gradient(
                    spc_info=spc_info[name],
                    thy_level=thy_level,
                    cnf_run_fs=cnf_run_fs,
                    cnf_save_fs=cnf_save_fs,
                    script_str=OPT_SCRIPT_STR,
                    overwrite=overwrite,
                    **KWARGS,
                )

            if run_min_hess:
                moldr.sp.run_minimum_energy_hessian(
                    spc_info=spc_info[name],
                    thy_level=thy_level,
                    cnf_run_fs=cnf_run_fs,
                    cnf_save_fs=cnf_save_fs,
                    script_str=OPT_SCRIPT_STR,
                    overwrite=overwrite,
                    **KWARGS,
                )

            if run_min_vpt2:
                moldr.sp.run_minimum_energy_vpt2(
                    spc_info=spc_info[name],
                    thy_level=thy_level,
                    cnf_run_fs=cnf_run_fs,
                    cnf_save_fs=cnf_save_fs,
                    script_str=OPT_SCRIPT_STR,
                    overwrite=overwrite,
                    **KWARGS,
                )

            if run_conf_scan:
                moldr.scan.hindered_rotor_scans(
                    spc_info=spc_info[name],
                    thy_level=thy_level,
                    cnf_run_fs=cnf_run_fs,
                    cnf_save_fs=cnf_save_fs,
                    script_str=OPT_SCRIPT_STR,
                    overwrite=overwrite,
                    scan_increment=scan_increment,
                    **OPT_KWARGS,
                )

            if run_conf_grad:
                moldr.conformer.run_conformer_gradients(
                    spc_info=spc_info[name],
                    thy_level=thy_level,
                    cnf_run_fs=cnf_run_fs,
                    cnf_save_fs=cnf_save_fs,
                    script_str=OPT_SCRIPT_STR,
                    overwrite=overwrite,
                    **KWARGS,
                )

            if run_conf_hess:
                moldr.conformer.run_conformer_hessians(
                    spc_info=spc_info[name],
                    thy_level=thy_level,
                    cnf_run_fs=cnf_run_fs,
                    cnf_save_fs=cnf_save_fs,
                    script_str=OPT_SCRIPT_STR,
                    overwrite=overwrite,
                    **KWARGS,
                )

            if run_tau_samp:
                moldr.tau.tau_sampling(
                    spc_info=spc_info[name],
                    thy_level=thy_level,
                    thy_save_fs=thy_save_fs,
                    tau_run_fs=tau_run_fs,
                    tau_save_fs=tau_save_fs,
                    script_str=OPT_SCRIPT_STR,
                    overwrite=overwrite,
                    nsamp_par=nsamp_tau_par,
                    **OPT_KWARGS,
                )

                if run_tau_grad:
                    moldr.tau.run_tau_gradients(
                        spc_info=spc_info[name],
                        thy_level=thy_level,
                        tau_run_fs=tau_run_fs,
                        tau_save_fs=tau_save_fs,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=overwrite,
                        **KWARGS,
                    )

                if run_tau_hess:
                    moldr.tau.run_tau_hessians(
                        spc_info=spc_info[name],
                        thy_level=thy_level,
                        tau_run_fs=tau_run_fs,
                        tau_save_fs=tau_save_fs,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=overwrite,
                        **KWARGS,
                    )

        for high_level_idx, _ in enumerate(run_high_levels):

            orb_restr = moldr.util.orbital_restriction(
                spc_info[name], ref_high_level)
            ref_level = ref_high_level[0:3]
            ref_level.append(orb_restr)

            spc_run_fs = autofile.fs.species(run_prefix)
            spc_run_fs.leaf.create(spc_info[name])
            spc_run_path = spc_run_fs.leaf.path(spc_info[name])
            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs.leaf.create(spc_info[name])
            spc_save_path = spc_save_fs.leaf.path(spc_info[name])

            ref_run_fs = autofile.fs.theory(spc_run_path)
            ref_run_fs.leaf.create(ref_level[1:4])
            ref_run_path = ref_run_fs.leaf.path(ref_level[1:4])
            ref_save_fs = autofile.fs.theory(spc_save_path)
            ref_save_fs.leaf.create(ref_level[1:4])
            ref_save_path = ref_save_fs.leaf.path(ref_level[1:4])

            cnf_run_fs = autofile.fs.conformer(ref_run_path)
            cnf_save_fs = autofile.fs.conformer(ref_save_path)
            min_cnf_locs = moldr.util.min_energy_conformer_locators(
                cnf_save_fs)
            cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
            cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
            min_cnf_geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)

            orb_restr = moldr.util.orbital_restriction(
                spc_info[name], run_high_levels[high_level_idx])
            thy_level = run_high_levels[high_level_idx][0:3]
            thy_level.append(orb_restr)
#            print('thy_level test:', thy_level)
#            print('conf_run_path test:', cnf_run_path)

            sp_run_fs = autofile.fs.single_point(cnf_run_path)
            sp_save_fs = autofile.fs.single_point(cnf_save_path)

            # evaluate the high level energy and save it

            prog = run_high_levels[high_level_idx][0]
            method = run_high_levels[opt_level_idx][1]
            SP_SCRIPT_STR, OPT_SCRIPT_STR, KWARGS, OPT_KWARGS = (
                moldr.util.run_qchem_par(prog, method))
            if run_hl_min_ene:
                moldr.sp.run_single_point_energy(
                    geo=min_cnf_geo,
                    spc_info=spc_info[name],
                    thy_level=thy_level,
                    sp_run_fs=sp_run_fs,
                    sp_save_fs=sp_save_fs,
                    script_str=SP_SCRIPT_STR,
                    overwrite=overwrite,
                    **KWARGS,
                )

                # add cycle over conformer geometries
                cnf_run_fs = autofile.fs.conformer(ref_run_path)
                cnf_save_fs = autofile.fs.conformer(ref_save_path)

                cnf_locs_lst = cnf_save_fs.leaf.existing()
                for locs in cnf_locs_lst:
                    cnf_run_path = cnf_run_fs.leaf.path(locs)
                    cnf_save_path = cnf_save_fs.leaf.path(locs)
                    cnf_geo = cnf_save_fs.leaf.file.geometry.read(locs)

                    moldr.sp.run_single_point_energy(
                        geo=min_cnf_geo,
                        spc_info=spc_info[name],
                        thy_level=thy_level,
                        sp_run_fs=sp_run_fs,
                        sp_save_fs=sp_save_fs,
                        script_str=SP_SCRIPT_STR,
                        overwrite=overwrite,
                        **KWARGS,
                    )

            run_hl_conf_opt = False
            if run_hl_conf_opt:
                # cycle over conformer geometries
                cnf_run_fs = autofile.fs.conformer(run_prefix)
                cnf_save_fs = autofile.fs.conformer(save_prefix)

                cnf_locs_lst = cnf_save_fs.leaf.existing()
                for locs in cnf_locs_lst:
                    cnf_run_path = cnf_run_fs.leaf.path(locs)
                    cnf_save_path = cnf_save_fs.leaf.path(locs)
                    cnf_geo = cnf_save_fs.leaf.file.geometry.read(locs)
                    cnf_geo_zma = cnf_save_fs.leaf.file.zmatrix.read(locs)

                    moldr.driver.run_job(
                        job='optimization',
                        script_str=SCRIPT_STR,
                        run_fs=cnf_run_fs,
                        geo=cnf_geo_zma,
                        spc_info=spc_info[name],
                        thy_level=run_high_levels[high_level_idx],
                        overwrite=overwrite,
                        **KWARGS,
                    )
                    # need to add a save


def ts_qchem(
        rxn_info_lst, smi_dct, chg_dct, mul_dct,
        run_opt_levels, ref_high_level, run_high_levels, geom_dct, run_prefix,
        save_prefix, qchem_flags, nsamp_pars, scan_increment, kickoff_pars,
        overwrite,
        ):
    """ run specified electronic structure calls for the transition states of
        the specified set of reactions
    """
    idx = 0
    run_ts_conf_samp = qchem_flags[idx]
    idx += 1
    run_ts_min_grad = qchem_flags[idx]
    idx += 1
    run_ts_min_hess = qchem_flags[idx]
    idx += 1
    run_ts_min_vpt2 = qchem_flags[idx]
    idx += 1
    run_ts_conf_scan = qchem_flags[idx]
    idx += 1
    run_ts_conf_grad = qchem_flags[idx]
    idx += 1
    run_ts_conf_hess = qchem_flags[idx]
    idx += 1
    run_ts_tau_samp = qchem_flags[idx]
    idx += 1
    run_ts_tau_grad = qchem_flags[idx]
    idx += 1
    run_ts_tau_hess = qchem_flags[idx]
    idx += 1
    run_ts_hl_min_ene = qchem_flags[idx]
    idx += 1
    run_ts_kicks_qchem = qchem_flags[idx]

    if run_ts_conf_samp:
        nsamp_ts_conf_par = nsamp_pars[0]
        print('The optimal conformers of the ts will be found through Monte Carlo sampling of torsions')

    if run_ts_min_grad:
        print('The TS gradient will be determined for the minimum energy conformers')

    if run_ts_min_hess:
        print('The TS hessian will be determined for the minimum energy conformers')

    if run_ts_min_vpt2:
        print('Second order vibrational perturbation theory will be performed for the TS at the minimum energy conformers')

    if run_ts_conf_scan:
        print('One-dimensional torsional scans will be performed for the TS starting from the minimum energy conformers')

    if run_ts_conf_grad:
        print('The gradient will be determined for each point along the one-dimensional torsional scans at the TS')

    if run_ts_conf_hess:
        print('The hessian will be determined for each point along the one-dimensional torsional scans at the TS')

    if run_ts_tau_samp:
        print('Random sampling of the torsional geometries in the TS will be performed') 
        nsamp_ts_tau_par = nsamp_pars[1]

    if run_ts_tau_grad:
        print('The gradient will be determined for each randomly sampled torsional geometry at the TS')

    if run_ts_tau_hess:
        print('The TS hessian will be determined for each randomly sampled torsional geometry at the TS')

    if run_ts_hl_min_ene:
        print('Higher level energies will be evaluated for the TS')

    if run_ts_kicks_qchem:
        print('Reactants and products will be determined by kicking off from the TS')

    nsamp_ts_conf_par = nsamp_pars[0]
    nsamp_ts_tau_par = nsamp_pars[1]
    for _, rct_names, prd_names, rxn_name in rxn_info_lst:
#    for rct_names, prd_names in zip(rct_names_lst, prd_names_lst):
    # print the CHEMKIN reaction name for reference
#        rxn_name = '='.join(['+'.join(rct_names), '+'.join(prd_names)])
        print()
        print("Reaction: {}".format(rxn_name))

        # determine inchis, charges, and multiplicities

        rct_smis = list(map(smi_dct.__getitem__, rct_names))
        prd_smis = list(map(smi_dct.__getitem__, prd_names))
        rct_ichs = list(map(automol.smiles.inchi, rct_smis))
        prd_ichs = list(map(automol.smiles.inchi, prd_smis))
        rct_chgs = list(map(chg_dct.__getitem__, rct_names))
        prd_chgs = list(map(chg_dct.__getitem__, prd_names))
        rct_muls = list(map(mul_dct.__getitem__, rct_names))
        prd_muls = list(map(mul_dct.__getitem__, prd_names))

        # determine the transition state multiplicity
        ts_low_mul = automol.mult.ts.low(rct_muls, prd_muls)
        ts_high_mul = automol.mult.ts._high(rct_muls)
        ts_mul = ts_low_mul
#        print('high_mul in es:', ts_high_mul)
        ts_chg = sum(rct_chgs)
#        print('ts_chg test:',ts_chg)
        ts_info = ('', ts_chg, ts_mul)

        # theory
        for opt_level_idx, _ in enumerate(run_opt_levels):

            # check direction of reaction
            rxn_ichs = [rct_ichs, prd_ichs]
            rxn_chgs = [rct_chgs, prd_chgs]
            rxn_muls = [rct_muls, prd_muls]
            rxn_exo = moldr.util.reaction_energy(
                save_prefix, rxn_ichs, rxn_chgs, rxn_muls, run_opt_levels[opt_level_idx])
            print(rxn_exo)
            if rxn_exo > 0:
                rct_ichs, prd_ichs = prd_ichs, rct_ichs
                rct_chgs, prd_chgs = prd_chgs, rct_chgs
                rct_muls, prd_muls = prd_muls, rct_muls
                print('ts search will be performed in reverse direction')

            # obtain geometries from a hierachy of (i) data directory and (ii)
            # previous species calculation
            rct_geos = []
            for ich, chg, mul in zip(rct_ichs, rct_chgs, rct_muls):
                rct_info = [ich, chg, mul]
                orb_restr = moldr.util.orbital_restriction(
                    rct_info, run_opt_levels[opt_level_idx])
                thy_level = run_opt_levels[opt_level_idx][0:3]
                thy_level.append(orb_restr)
                spc_save_fs = autofile.fs.species(save_prefix)
                spc_save_fs.leaf.create(rct_info)
                spc_save_path = spc_save_fs.leaf.path(rct_info)
                thy_save_fs = autofile.fs.theory(spc_save_path)
                thy_save_fs.leaf.create(thy_level[1:4])
                thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])

                geo = moldr.util.reference_geometry(
                    rct_info, thy_level=thy_level, thy_fs=thy_save_fs,
                    geom_dct=geom_dct)
                # if conformer search already done then replace with conformational minimum geometry
                cnf_save_fs = autofile.fs.conformer(thy_save_path)
                min_cnf_locs = moldr.util.min_energy_conformer_locators(
                   cnf_save_fs)
                if min_cnf_locs:
                    cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
                    geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
                rct_geos.append(geo)

            prd_geos = []
            for ich, chg, mul in zip(prd_ichs, prd_chgs, prd_muls):
                prd_info = [ich, chg, mul]
                orb_restr = moldr.util.orbital_restriction(
                    prd_info, run_opt_levels[opt_level_idx])
                thy_level = run_opt_levels[opt_level_idx][0:3]
                thy_level.append(orb_restr)
                spc_save_fs = autofile.fs.species(save_prefix)
                spc_save_fs.leaf.create(prd_info)
                spc_save_path = spc_save_fs.leaf.path(prd_info)
                thy_save_fs = autofile.fs.theory(spc_save_path)
                thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])
                geo = moldr.util.reference_geometry(
                    prd_info, thy_level=thy_level, thy_fs=thy_save_fs,
                    geom_dct=geom_dct)
                cnf_save_fs = autofile.fs.conformer(thy_save_path)
                min_cnf_locs = moldr.util.min_energy_conformer_locators(
                   cnf_save_fs)
                if min_cnf_locs:
                    cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
                    geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
                prd_geos.append(geo)

            # determine the transition state z-matrix
            # replace this with save values if they are available
            rct_zmas = list(map(automol.geom.zmatrix, rct_geos))
            prd_zmas = list(map(automol.geom.zmatrix, prd_geos))

            typ = None

            ret = automol.zmatrix.ts.beta_scission(rct_zmas, prd_zmas)
            if ret and typ is None:
                typ = 'beta scission'
                ts_zma, dist_name, tors_names = ret
                print('beta scission')
                print('ts zma:', ts_zma)
                print('dist name:', dist_name)
                print('tors names:', tors_names)

            ret = automol.zmatrix.ts.addition(rct_zmas, prd_zmas)
            if ret and typ is None:
                typ = 'addition'
                ts_zma, dist_name, tors_names = ret
                print('addn')
                print('ts zma:', ts_zma)
                print('dist name:', dist_name)
                print('tors names:', tors_names)
                if rct_muls[0] != 1 and rct_muls[1] != 1:
                   typ = 'radical radical addition' 

            # fix this later
            # ret = automol.zmatrix.ts.hydrogen_abstraction(rct_zmas, prd_zmas,
            #                                               sigma=True)
            ret = automol.zmatrix.ts.hydrogen_abstraction(rct_zmas, prd_zmas,
                                                          sigma=False)
            if ret and typ is None:
                typ = 'hydrogen abstraction'
                ts_zma, dist_name, tors_names = ret
                print('H abs')
                print('ts zma:', ts_zma)
                print('dist name:', dist_name)
                print('tors names:', tors_names)

            if typ is None:
                print("Failed to classify reaction.")
            else:
                print("Type: {}".format(typ))

                # determine the grid
                dist_coo, = automol.zmatrix.coordinates(ts_zma)[dist_name]
                syms = automol.zmatrix.symbols(ts_zma)
                bnd_len_key = tuple(sorted(map(syms.__getitem__, dist_coo)))

                bnd_len_dct = {
                    ('C', 'C'): 1.54 * ANG2BOHR,
                    ('C', 'H'): 1.09 * ANG2BOHR,
                    ('H', 'H'): 0.74 * ANG2BOHR,
                    ('N', 'N'): 1.45 * ANG2BOHR,
                    ('O', 'O'): 1.48 * ANG2BOHR,
                    ('C', 'N'): 1.47 * ANG2BOHR,
                    ('C', 'O'): 1.43 * ANG2BOHR,
                    ('H', 'O'): 1.20 * ANG2BOHR,
                    ('H', 'N'): 0.99 * ANG2BOHR,
                }

                npoints = 8
                npoints1 = 4
                npoints2 = 4
                if typ in ('beta scission', 'addition'):
                    rmin = 1.4 * ANG2BOHR
                    rmax = 2.8 * ANG2BOHR
                    if bnd_len_key in bnd_len_dct:
                        bnd_len = bnd_len_dct[bnd_len_key]
                        rmin = bnd_len + 0.2 * ANG2BOHR
                        rmax = bnd_len + 1.6 * ANG2BOHR
                    grid = numpy.linspace(rmin, rmax, npoints)
                elif typ == 'hydrogen abstraction':
                    rmin = 0.7 * ANG2BOHR
                    rmax = 2.2 * ANG2BOHR
                    if bnd_len_key in bnd_len_dct:
                        bnd_len = bnd_len_dct[bnd_len_key]
                        rmin = bnd_len
                        rmax = bnd_len + 1.0 * ANG2BOHR
                    grid = numpy.linspace(rmin, rmax, npoints)
                elif typ == 'radical radical addition':
                    rstart = 2.4 * ANG2BOHR
                    rend1 = 3.0 * ANG2BOHR
                    rend2 = 1.8 * ANG2BOHR
                    grid1 = numpy.linspace(rstart, rend1, npoints1)
                    grid2 = numpy.linspace(rstart, rend2, npoints2)
                    grid2 = numpy.delete(grid2, 0)

                # set up the file systems
                prog = run_opt_levels[opt_level_idx][0]
                method = run_opt_levels[opt_level_idx][1]
                SCRIPT_STR, OPT_SCRIPT_STR, KWARGS, OPT_KWARGS = (
                    moldr.util.run_qchem_par(prog, method))

    #            ts_orb_restr = moldr.util.orbital_restriction(
    #                ts_info, run_opt_levels[opt_level_idx])
    #            thy_level = run_opt_levels[opt_level_idx][0:3]
    #            thy_level.append(ts_orb_restr)

    #            thy_run_fs = autofile.fs.theory(spc_run_path)
    #            thy_run_fs.leaf.create(thy_level[1:4])
    #            thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])

    #            thy_save_fs = autofile.fs.theory(spc_save_path)
    #            thy_save_fs.leaf.create(thy_level[1:4])
    #            thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])

    #            run_fs = autofile.fs.run(thy_run_path)

    #            cnf_run_fs = autofile.fs.conformer(thy_run_path)
    #            cnf_save_fs = autofile.fs.conformer(thy_save_path)

    #            tau_run_fs = autofile.fs.tau(thy_run_path)
    #            tau_save_fs = autofile.fs.tau(thy_save_path)
                # construct the filesystem
                rxn_ichs = [rct_ichs, prd_ichs]
                rxn_chgs = [rct_chgs, prd_chgs]
                rxn_muls = [rct_muls, prd_muls]

                # set up the filesystem
                is_rev = autofile.system.reaction_is_reversed(
                    rxn_ichs, rxn_chgs, rxn_muls)
                rxn_ichs, rxn_chgs, rxn_muls = autofile.system.sort_together(
                    rxn_ichs, rxn_chgs, rxn_muls)
                print(" - The reaction direction is {}"
                      .format('backward' if is_rev else 'forward'))


                rxn_run_fs = autofile.fs.reaction(run_prefix)
                rxn_run_fs.leaf.create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
                rxn_run_path = rxn_run_fs.leaf.path(
                    [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

                rxn_ichs = tuple(map(tuple, rxn_ichs))
                rxn_chgs = tuple(map(tuple, rxn_chgs))
                rxn_muls = tuple(map(tuple, rxn_muls))
#                print('rxn_save test0', rxn_ichs, rxn_chgs, rxn_muls, ts_mul)
#                print(save_prefix)
                rxn_save_fs = autofile.fs.reaction(save_prefix)
                rxn_save_fs.leaf.create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
                rxn_save_path = rxn_save_fs.leaf.path(
                    [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
                print(rxn_save_path)

                ts_orb_restr = moldr.util.orbital_restriction(
                    ts_info, run_opt_levels[opt_level_idx])
                ref_level = run_opt_levels[opt_level_idx][0:3]
                ref_level.append(ts_orb_restr)
#                print('ref level test:', ref_level[1:4])
#                print('rxn run path:', rxn_run_path)

                thy_run_fs = autofile.fs.theory(rxn_run_path)
                thy_run_fs.leaf.create(ref_level[1:4])
                thy_run_path = thy_run_fs.leaf.path(
                        ref_level[1:4])

                thy_save_fs = autofile.fs.theory(rxn_save_path)
                thy_save_fs.leaf.create(ref_level[1:4])
                thy_save_path = thy_save_fs.leaf.path(ref_level[1:4])

                scn_run_fs = autofile.fs.scan(thy_run_path)
                scn_save_fs = autofile.fs.scan(thy_save_path)

                print('entering run_scan:')

                if typ == 'radical radical addition':
                    ts_formula = ''
                    for ich in rct_ichs:
                        formula_i = thermo.util.inchi_formula(ich)
                        formula_i_dict = thermo.util.get_atom_counts_dict(formula_i)
                        ts_formula = automol.formula._formula.join(ts_formula, formula_i_dict)

                    grid = numpy.append(grid1, grid2)
                    high_mul = automol.mult.ts._high(rct_muls)
                    moldr.scan.run_multiref_rscan(
                        formula=ts_formula,
                        high_mul=high_mul,
                        zma=ts_zma,
                        spc_info=ts_info,
                        thy_level=ref_level,
                        dist_name=dist_name,
                        grid1=grid1,
                        grid2=grid2,
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=SCRIPT_STR,
                        overwrite=overwrite,
                        update_guess=True,
                        **OPT_KWARGS
                    )

                    moldr.scan.save_scan(
                        scn_run_fs=scn_run_fs,
                        scn_save_fs=scn_save_fs,
                        coo_names=[dist_name],
                    )

                    ref_ene = -40.

                    nsamp_max = 2000
                    nsamp_min = 500
                    flux_err = 5
                    pes_size = 1
                    tst_inp_str = varecof_io.writer.write_tst_input(
                        nsamp_max, nsamp_min, flux_err, pes_size)

                    print('\ntst.inp:')
                    print(tst_inp_str)

                    # Write the divsur input file string; distances in Angstrom
                    distances = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
                    divsur_inp_str = varecof_io.writer.write_divsur_input(
                        distances)
                    print('\ndivsur.inp:')
                    print(divsur_inp_str)

                    with open(os.path.join(varecof_path, 'divsur.inp'), 'w') as divsur_file:
                        divsur_file.write(divsur_inp_str)

                    # Write the els input string
                    exe_path = molpro_path_str
                    base_name = formula
                    els_inp_str = varecof_io.writer.write_els_input(
                        exe_path, base_name)
                    print('\nels.inp:')
                    print(els_inp_str)

                    with open(os.path.join(varecof_path, 'els.inp'), 'w') as els_file:
                        els_file.write(els_inp_str)

                    # Write the structure input string
                    struct_inp_str = varecof_io.writer.write_structure_input(
                        rct_geos[0], rct_geos[1])
                    print('\nstructure.inp:')
                    print(struct_inp_str)

                    with open(os.path.join(varecof_path, 'structure.inp'), 'w') as struct_file:
                        struct_file.write(tml_inp_str)


                    # Write the *.tml input string
                    electron_count = automol.formula._formula.electron_count(formula)
                    # this is only for 2e,2o case
                    closed_orb = electron_count//2 - 1
                    occ_orb = electron_count//2 + 1
                    # end of 2e,2o case
                    two_spin = spc_info[2]-1
                    chg = spc_info[1]
                    wfn_info = [electron_count, closed_orb, occ_orb, chg, high_spin, low_spin]
                    method = method
                    shift_key = '0.2'
                    ipea_shift_key = '0.25'
                    inf_sep_energy = -654.3210123456
                    tml_inp_str = varecof_io.writer.write_tml_input(
                        prog = prog,
                        method = method,
                        basis = basis,
                        memory = kwargs['memory'],
                        wfn_info = wfn_info,
                        inf_sep_energy = inf_sep_energy
                        )
                    print('\nmol.tml:')
                    print(tml_inp_str)

                    with open(os.path.join(varecof_path, basename+'.tml'), 'w') as tml_file:
                        tml_file.write(tml_inp_str)

                else:

                    moldr.scan.run_scan(
                        zma=ts_zma,
                        spc_info=ts_info,
                        thy_level=ref_level,
                        grid_dct={dist_name: grid},
                        scn_run_fs=scn_run_fs,
                        scn_save_fs=scn_save_fs,
                        script_str=SCRIPT_STR,
                        overwrite=overwrite,
                        update_guess=False,
                        reverse_sweep=False,
                        **OPT_KWARGS
                    )

                    moldr.scan.save_scan(
                        scn_run_fs=scn_run_fs,
                        scn_save_fs=scn_save_fs,
                        coo_names=[dist_name],
                    )

                    scn_save_fs = autofile.fs.scan(thy_save_path)
                    locs_lst = [
                        locs for locs in scn_save_fs.leaf.existing([[dist_name]])
                        if scn_save_fs.leaf.file.energy.exists(locs)]
    #                print(locs_lst)
                    enes = [scn_save_fs.leaf.file.energy.read(locs)
                            for locs in locs_lst]
                    max_locs = locs_lst[enes.index(max(enes))]
                    max_ene = max(enes)
                    max_zma = scn_save_fs.leaf.file.zmatrix.read(max_locs)
                    print('geometry for maximum along scan:', max_zma)
                    print('energy for maximum along scan:', max_ene)

                    print('optimizing ts')
                    # find saddlepoint from maximum on the grid opt scan
                    print('thy_run_path in ts_opt:', thy_run_path)
                    ts_run_fs = autofile.fs.ts(thy_run_path)
                    ts_run_fs.trunk.create()
                    ts_run_path = ts_run_fs.trunk.path()
                    run_fs = autofile.fs.run(ts_run_path)

                    print('ts_run_path:', ts_run_path)

                    ts_save_fs = autofile.fs.ts(thy_save_path)
                    ts_save_fs.trunk.create()
                    ts_save_path = ts_save_fs.trunk.path()
                    print('ts_save_path:', ts_save_path)
                    scn_run_fs = autofile.fs.scan(ts_run_path)
                    scn_save_fs = autofile.fs.scan(ts_save_path)

                    print('starting ts optimization')
                    print('theory_level=:', run_opt_levels[opt_level_idx])
                    print('ts_run_path=:', ts_run_path)
                    moldr.driver.run_job(
                        job='optimization',
                        script_str=SCRIPT_STR,
                        run_fs=run_fs,
                        geom=max_zma,
                        spc_info=ts_info,
                        thy_level=ref_level,
                        saddle=True,
                        overwrite=overwrite,
                        **OPT_KWARGS,
                    )
                    opt_ret = moldr.driver.read_job(
                        job='optimization',
                        run_fs=run_fs,
                    )
                    if opt_ret is not None:
                        inf_obj, inp_str, out_str = opt_ret
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

                    if run_ts_conf_scan:
                        zma = ts_save_fs.trunk.file.zmatrix.read()
                        val_dct = automol.zmatrix.values(zma)
                        tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
                            zma, tors_names, scan_increment)
                        tors_grids = [
                            numpy.linspace(*linspace) + val_dct[name]
                            for name, linspace in zip(tors_names, tors_linspaces)]
                        for tors_name, tors_grid in zip(tors_names, tors_grids):
                            moldr.scan.run_scan(
                                zma=zma,
                                spc_info=ts_info,
                                thy_level=ref_level,
                                grid_dct={tors_name: tors_grid},
                                scn_run_fs=scn_run_fs,
                                scn_save_fs=scn_save_fs,
                                script_str=SCRIPT_STR,
                                overwrite=overwrite,
                                saddle=True,
                                **OPT_KWARGS,
                            )

                            moldr.scan.save_scan(
                                scn_run_fs=scn_run_fs,
                                scn_save_fs=scn_save_fs,
                                coo_names=[tors_name],
                            )

                        hind_rot_dct = {}

                        min_ene = ts_save_fs.trunk.file.energy.read()
                        for tors_name in tors_names:
                            enes = [scn_save_fs.leaf.file.energy.read(locs)
                                    for locs in scn_save_fs.leaf.existing([tors_name])]
                            enes = numpy.subtract(enes, min_ene)
                            hind_rot_dct[tors_name] = enes*EH2KCAL

                        print('ts hindered rotor potential')
                        print(hind_rot_dct)

                    if run_ts_conf_samp:
                        # set up cnf_run_fs and cnf_save_fs
#                        print('ts_conf_test:')
#                        print(ts_info)
#                        print(run_opt_levels[opt_level_idx])
#                        print(ts_run_path)
#                        print(ts_save_path)
#                        print(rxn_run_path)
#                        print(rxn_save_path)
#                        print(tors_names)
                        cnf_run_fs = autofile.fs.conformer(thy_run_path)
                        cnf_save_fs = autofile.fs.conformer(thy_save_path)
                        moldr.ts.ts_conformer_sampling(
                            spc_info=ts_info,
                            geo=geo,
                            zma=zma,
                            tors_names=tors_names,
                            thy_level=ref_level,
                            thy_save_fs=thy_save_fs,
                            cnf_run_fs=cnf_run_fs,
                            cnf_save_fs=cnf_save_fs,
                            script_str=SCRIPT_STR,
                            overwrite=overwrite,
                            saddle=True,
                            nsamp_par=nsamp_ts_conf_par,
                            **OPT_KWARGS,
                        )

                        if run_ts_min_grad:
                            moldr.sp.run_minimum_energy_gradient(
                                spc_info=ts_info,
                                thy_level=ref_level,
                                cnf_run_fs=cnf_run_fs,
                                cnf_save_fs=cnf_save_fs,
                                script_str=SCRIPT_STR,
                                overwrite=overwrite,
                                **KWARGS,
                            )

                        if run_ts_min_hess:
                            moldr.sp.run_minimum_energy_hessian(
                                spc_info=ts_info,
                                thy_level=ref_level,
                                cnf_run_fs=cnf_run_fs,
                                cnf_save_fs=cnf_save_fs,
                                script_str=SCRIPT_STR,
                                overwrite=overwrite,
                                **KWARGS,
                            )

                        if run_ts_min_vpt2:
                            moldr.sp.run_minimum_energy_vpt2(
                                spc_info=ts_info,
                                thy_level=ref_level,
                                cnf_run_fs=cnf_run_fs,
                                cnf_save_fs=cnf_save_fs,
                                script_str=SCRIPT_STR,
                                overwrite=overwrite,
                                saddle=True,
                                **KWARGS,
                            )

                        if run_ts_tau_samp:

                            moldr.tau.save_tau(
                                tau_run_fs=ts_run_fs,
                                tau_save_fs=ts_save_fs,
                            )

                            zma = ts_save_fs.trunk.file.zmatrix.read()
                            tors_ranges = automol.zmatrix.torsional_sampling_ranges(
                                zma, tors_names)
                            tors_range_dct = dict(zip(tors_names, tors_ranges))

                            moldr.tau.run_tau(
                                zma=zma,
                                spc_info=ts_info,
                                thy_level=ref_level,
                                nsamp=nsamp_ts_tau_par,
                                tors_range_dct=tors_range_dct,
                                run_fs=ts_run_fs,
                                save_fs=ts_save_fs,
                                script_str=SCRIPT_STR,
                                overwrite=overwrite,
                                saddle=True,
                                **OPT_KWARGS,
                            )
        #                        saddle=True,
        # used to have saddle=True in call, but this is not used. Probably a bug.
        # Need to pass saddle information somewhere and use it

                            moldr.tau.save_tau(
                                run_prefix=thy_run_path,
                                save_prefix=thy_save_path,
                            )

                        if run_ts_kicks_qchem:
                            ret = moldr.ts.read_job(
                                job=elstruct.Job.HESSIAN, prefix=ts_run_path)
                            if ret:
                                inf_obj, _, out_str = ret
                                prog = inf_obj.prog
                                hess = elstruct.reader.hessian(prog, out_str)
                                freqs = elstruct.util.harmonic_frequencies(geo, hess, project=True)
                                norm_coos = elstruct.util.normal_coordinates(geo, hess, project=True)
                                assert freqs[0] < -100

                                print('Kicking off from saddle in forward direction')
                                im_norm_coo = numpy.array(norm_coos)[:, 0]
                                disp_xyzs = numpy.reshape(im_norm_coo, (-1, 3))
                                dir_afs = autofile.fs.direction()
                                fwd_run_path = dir_afs.direction.dir.path(thy_run_path, [True])
                                dir_afs.direction.dir.create(thy_run_path, [True])
                                fwd_save_path = dir_afs.direction.dir.path(thy_save_path, [True])
                                dir_afs.direction.dir.create(thy_save_path, [True])
                                print(automol.geom.string(geo))
                                print(disp_xyzs)
                                moldr.geom.run_kickoff_saddle(
                                    geo, disp_xyzs, ts_info, method, basis, orb_restr,
                                    fwd_run_path, SCRIPT_STR, prog, overwrite,
                                    kickoff_size=kickoff_size, kickoff_backward=False,
                                    opt_cart=True, **OPT_KWARGS)
                                print('Saving product of kick off from saddle in forward direction')
                                ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, prefix=fwd_run_path)
                                if ret:
                                    inf_obj, inp_str, out_str = ret
                                    prog = inf_obj.prog
                                    method = inf_obj.method
                                    ene = elstruct.reader.energy(prog, method, out_str)
                                    geo = elstruct.reader.opt_geometry(prog, out_str)
                                    fwd_save_path = dir_afs.direction.dir.path(thy_save_path, [True])
                                    print('save path', fwd_save_path)
                                    dir_afs.direction.file.geometry_info.write(inf_obj, thy_save_path, [True])
                                    dir_afs.direction.file.geometry_input.write(inp_str, thy_save_path, [True])
                                    dir_afs.direction.file.geometry.write(geo, thy_save_path, [True])
                                    dir_afs.direction.file.energy.write(ene, thy_save_path, [True])

                                print('Kicking off from saddle in backward direction')
                                bwd_run_path = dir_afs.direction.dir.path(thy_run_path, [False])
                                dir_afs.direction.dir.create(thy_run_path, [False])
                                bwd_save_path = dir_afs.direction.dir.path(thy_save_path, [False])
                                dir_afs.direction.dir.create(thy_save_path, [False])
                                moldr.geom.run_kickoff_saddle(
                                    geo, disp_xyzs, chg, mul, method, basis,
                                    orb_restr, bwd_run_path, SCRIPT_STR, prog,
                                    overwrite, kickoff_size=kickoff_size,
                                    kickoff_backward=True, **OPT_KWARGS)
                                print('Saving product of kick off from saddle in backward direction')
                                ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, prefix=bwd_run_path)
                                if ret:
                                    inf_obj, inp_str, out_str = ret
                                    prog = inf_obj.prog
                                    method = inf_obj.method
                                    ene = elstruct.reader.energy(prog, method, out_str)
                                    geo = elstruct.reader.opt_geometry(prog, out_str)
                                    bwd_save_path = dir_afs.direction.dir.path(thy_save_path, [False])
                                    print('save path', bwd_save_path)
                                    dir_afs.direction.file.geometry_info.write(inf_obj, thy_save_path, [False])
                                    dir_afs.direction.file.geometry_input.write(inp_str, thy_save_path, [False])
                                    dir_afs.direction.file.geometry.write(geo, thy_save_path, [False])
                                    dir_afs.direction.file.energy.write(ene, thy_save_path, [False])

                    for high_level_idx, _ in enumerate(run_high_levels):

                        orb_restr = moldr.util.orbital_restriction(
                            ts_info, ref_high_level)
                        ref_level = ref_high_level[0:3]
                        ref_level.append(orb_restr)
#                        print('ref level test:', ref_level)

                        ref_run_fs = autofile.fs.theory(rxn_run_path)
                        ref_run_fs.leaf.create(ref_level[1:4])
                        ref_run_path = ref_run_fs.leaf.path(ref_level[1:4])
                        ref_save_fs = autofile.fs.theory(rxn_save_path)
                        ref_save_fs.leaf.create(ref_level[1:4])
                        ref_save_path = ref_save_fs.leaf.path(ref_level[1:4])

                        thy_run_fs = autofile.fs.theory(ref_run_path)
                        thy_run_fs.leaf.create(ref_level[1:4])
                        thy_run_path = thy_run_fs.leaf.path(
                                ref_level[1:4])

                        thy_save_fs = autofile.fs.theory(ref_save_path)
                        thy_save_fs.leaf.create(ref_level[1:4])
                        thy_save_path = thy_save_fs.leaf.path(ref_level[1:4])

                        print('thy_run_path in ts_opt:', thy_run_path)
                        ts_run_fs = autofile.fs.ts(thy_run_path)
                        ts_run_fs.trunk.create()
                        ts_run_path = ts_run_fs.trunk.path()
                        print('ts_run_path:', ts_run_path)

                        ts_save_fs = autofile.fs.ts(thy_save_path)
                        ts_save_fs.trunk.create()
                        ts_save_path = ts_save_fs.trunk.path()
                        print('ts_save_path:', ts_save_path)
                        # evaluate the high level energy and save it

                        cnf_run_fs = autofile.fs.conformer(ts_run_path)
                        cnf_save_fs = autofile.fs.conformer(ts_save_path)
                        min_cnf_locs = moldr.util.min_energy_conformer_locators(
                            cnf_save_fs)
                        cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
                        min_cnf_geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)

                        prog = run_high_levels[high_level_idx][0]
                        method = run_high_levels[opt_level_idx][1]
                        SP_SCRIPT_STR, OPT_SCRIPT_STR, KWARGS, OPT_KWARGS = (
                            moldr.util.run_qchem_par(prog, method))
                        if run_ts_hl_min_ene:
                            moldr.sp.run_single_point_energy(
                                geo=min_cnf_geo,
                                spc_info=ts_info,
                                thy_level=ref_level,
                                run_fs=cnf_run_fs,
                                save_fs=cnf_save_fs,
                                script_str=SP_SCRIPT_STR,
                                overwrite=overwrite,
                                **KWARGS,
                            )

def vdw_qchem(
        rct_names_lst, prd_names_lst, smi_dct, chg_dct, mul_dct, run_opt_levels, ref_high_level,
        run_high_levels, geom_dct, run_prefix, save_prefix, qchem_flags,
        nsamp_pars, scan_increment, kickoff_pars, overwrite,
        ):

    ntaudof = 5
    nsamp_par = nsamp_pars[2]
    nsamp = nsamp_init(nsamp_par, ntaudof)

    VDW_NAMES_LST = []
    for rct_names, prd_names in zip(rct_names_lst, prd_names_lst):
        rct_muls = list(map(mul_dct.__getitem__, rct_names))
        prd_muls = list(map(mul_dct.__getitem__, prd_names))
        ts_mul = automol.mult.ts.low(rct_muls, prd_muls)
        if len(rct_names) == 2:
            if sorted(rct_names) not in VDW_NAMES_LST:
                VDW_NAMES_LST.append([sorted(rct_names), ts_mul])
        if len(prd_names) == 2:
            if sorted(prd_names) not in VDW_NAMES_LST:
                VDW_NAMES_LST.append([sorted(prd_names), ts_mul])

    for names, ts_mul in VDW_NAMES_LST:
        smis = list(map(smi_dct.__getitem__, names))
        ichs = list(map(automol.smiles.inchi, smis))
        chgs = list(map(chg_dct.__getitem__, names))
        muls = list(map(mul_dct.__getitem__, names))

        for opt_level_idx, _ in enumerate(run_opt_levels):
            # theory
            prog = run_opt_levels[opt_level_idx][0]
            method = run_opt_levels[opt_level_idx][1]
            SP_SCRIPT_STR, OPT_SCRIPT_STR, KWARGS, OPT_KWARGS = moldr.util.run_qchem_par(prog, method)

            geos = []
            for ich, chg, mul in zip(ichs, chgs, muls):
                spc_info = [ich, chg, mul]
                orb_restr = moldr.util.orbital_restriction(spc_info, run_opt_levels[opt_level_idx][0])
                geo = moldr.util.reference_geometry(
                    spc_info=spc_info, thy_level=run_opt_levels[opt_level_idx],
                    prefix=save_prefix, geom_dct=GEOM_DCT)
                geos.append(geo)
               
            geo1, geo2 = geos
            geo1 = automol.geom.mass_centered(geo1)
            geo2 = automol.geom.mass_centered(geo2)
            for idx in range(nsamp):
                print('Optimizing vdw geometry {}/{}'.format(idx+1, nsamp))
                angs1 = numpy.multiply(
                    numpy.random.rand(3), [1*numpy.pi, 2*numpy.pi, 2*numpy.pi])
                angs2 = numpy.multiply(
                    numpy.random.rand(3), [1*numpy.pi, 2*numpy.pi, 2*numpy.pi])
                angs12 = numpy.multiply(
                    numpy.random.rand(2), [1*numpy.pi, 2*numpy.pi])
                geo1 = automol.geom.euler_rotated(geo1, *angs1)
                geo2 = automol.geom.euler_rotated(geo2, *angs2)
                dist_cutoff = 3.*qcc.conversion_factor('angstrom', 'bohr')

                geo = automol.geom.join(geo1, geo2, dist_cutoff, *angs12)
                print("Species: {}".format('+'.join(names)))
                print('vdw starting geometry')
                print(automol.geom.xyz_string(geo))

        # set up the filesystem
                ich = automol.inchi.recalculate(automol.inchi.join(ichs))
                chg = sum(chgs)
                mul = ts_mul
                orb_restr = moldr.util.orbital_restriction(mul, RESTRICT_OPEN_SHELL)
                spc_run_path = moldr.util.species_path(ich, chg, mul, run_prefix)
                spc_save_path = moldr.util.species_path(ich, chg, mul, save_prefix)
                thy_run_path = moldr.util.theory_path(method, basis, orb_restr, spc_run_path)
                thy_save_path = moldr.util.theory_path(method, basis, orb_restr, spc_save_path)

        # generate reference geometry
        # generate the z-matrix and sampling ranges

                moldr.driver.run_job(
                    job=elstruct.Job.OPTIMIZATION,
                    geom=geo,
                    spc_info=ts_info,
                    thy_level=run_opt_levels[opt_level_idx],
                    prefix=thy_run_path,
                    script_str=SCRIPT_STR,
                    overwrite=overwrite,
                    **OPT_KWARGS,
                )

        # save info for the initial geometry (from inchi or from save directory)
                ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, prefix=thy_run_path)
                if ret:
                    print('Saving reference geometry')
                    print(" - Save path: {}".format(thy_save_path))

                    inf_obj, inp_str, out_str = ret
                    prog = inf_obj.prog
                    method = inf_obj.method
                    geo = elstruct.reader.opt_geometry(prog, out_str)
                    thy_afs = autofile.fs.theory()
                    thy_afs.theory.file.geometry.write(geo, spc_save_path, [method, basis, orb_restr])
                    ene = elstruct.reader.energy(prog, method, out_str)
                    print('ene test in vdw')
                    print(ene)
                    thy_afs.theory.file.energy.write(ene, spc_save_path, [method, basis, orb_restr])


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