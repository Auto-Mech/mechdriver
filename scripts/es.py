""" reaction list test
"""
import numpy
from qcelemental import constants as qcc
import automol
import elstruct
import autofile
import moldr
import scripts

ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')
WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

import logging
log   = logging.getLogger(__name__)



def run_freq(save_prefix, spcdic, params, OPT_KWARGS):
    log.info('    | running task {}'.format('freq'))
    moldr.driver.run_minimum_energy_hessian(**params, **OPT_KWARGS)
    moldr.driver.run_minimum_energy_gradient(**params, **OPT_KWARGS)
    return

def run_anharm(prefixes, spcdic, params, OPT_KWARGS):
    log.info('    | running task {}'.format('anharm'))
    print('anharm is not working in moldr right now')
    #moldr.driver.run_minimum_energy_vpt2(**params, **OPT_KWARGS)
    return

def run_mc(prefixes, spcdic, params, OPT_KWARGS):
    if spcdic['mc_nsamp'][0]:
       log.info('    | running task {} with abcd of {}'.format('mc', ' '.join([str(x) for x in spcdic['mc_nsamp'][1:]])))
    else:
        log.info('    | running task {} for {:g} points'.format('mc', spcdic['mc_nsamp'][5]))
    spc_info      = get_spc_info(spcdic)
    spc_run_path  = get_spc_run_path( prefixes[0], spc_info) 
    spc_save_path = get_spc_save_path(prefixes[1], spc_info) 
    params[  'run_prefix'] = spc_run_path
    params[ 'save_prefix'] = spc_save_path
    params[   'nsamp_par'] = spcdic['mc_nsamp']
    moldr.driver.conformer_sampling(**params, **OPT_KWARGS)
    return

def run_hr(prefixes, spcdic, params, OPT_KWARGS):
    log.info('    | running task {}'.format('hr'))
    params['scan_increment'] = spcdic['hind_inc']
    moldr.driver.hindered_rotor_scans(**params, **OPT_KWARGS)
    return

def run_task(tsk, spcdic, thydic, initial_thy_info, run_prefix, save_prefix, overwrite):

    spc_info      = get_spc_info(spcdic)
    thy_info      = get_thy_info(thydic)
    spc_run_path  = get_spc_run_path(run_prefix, spc_info) 
    thy_run_path  = get_thy_run_path(run_prefix, spc_info, thy_info)
    spc_save_path = get_spc_save_path(save_prefix, spc_info) 
    thy_save_path = get_thy_save_path(save_prefix, spc_info, thy_info) 
    sp_script_str, opt_script_str, KWARGS, OPT_KWARGS = moldr.util.run_qchem_par(thy_info[0])
    prefixes = [run_prefix, save_prefix]
 
    params =  {    'spc_info': spc_info, 
               'theory_level': thy_info, 
                 'run_prefix': thy_run_path, 
                'save_prefix': thy_save_path, 
                 'script_str': opt_script_str,
                  'overwrite': overwrite} 

    choose_function = {'freq':   'run_freq',
                       'anharm': 'run_anharm', 
                       'mc':     'run_mc',
                       'hr':     'run_hr'}

    if tsk in choose_function:
        eval( choose_function[tsk] )(prefixes, spcdic, params, OPT_KWARGS)

    elif tsk == 'sp':  #I still have to make a function for this one
        log.info('    | initializing...')
        ref_run_path  = get_thy_run_path(run_prefix, spc_info, initial_thy_info) 
        ref_save_path = get_thy_save_path(save_prefix, spc_info, initial_thy_info)
        #Make sure there is at least one CONF for this level of theory
        params[  'run_prefix'] = spc_run_path
        params[ 'save_prefix'] = spc_save_path
        params['theory_level'] = initial_thy_info
        params[   'nsamp_par'] = [False, 0, 0, 0, 0, 1]
        moldr.driver.conformer_sampling(**params, **OPT_KWARGS)
        min_cnf_locs  = moldr.util.min_energy_conformer_locators(
            ref_save_path)
        cnf_run_fs   = autofile.fs.conformer(ref_run_path)
        cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
        cnf_save_fs  = autofile.fs.conformer(ref_save_path)
        cnf_save_path= cnf_save_fs.leaf.path(min_cnf_locs)
        min_cnf_geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        
        log.info('    | running sp on geo')
        log.info(min_cnf_geo)
        del params['nsamp_par']
        params[        'geo'] = min_cnf_geo
        params[ 'run_prefix'] = cnf_run_path
        params['save_prefix'] = cnf_save_path
        params[ 'script_str'] = sp_script_str

        moldr.driver.run_single_point_energy(**params, **KWARGS)
    return


def geo_init(params):
    geo, msg =  moldr.util.reference_geometry(**params)
    return geo, msg

def geo_ref(spcdic, initial_thy_info, running_thy_info, run_prefix, save_prefix, overwrite):


    log.info('    | initializing...')
    spc_info       = get_spc_info(spcdic)

    spc_run_path  = get_spc_run_path(run_prefix, spc_info) 
    spc_save_path = get_spc_save_path(save_prefix, spc_info)
 
    #Check to see if geo already exists at running_theroy
    params =  {    'spc_info': spc_info ,
               'theory_level': running_thy_info, 
                     'prefix': save_prefix, 
                 'return_msg': True} 
    geom_obj = spcdic['geoobj']
    params['geom_obj'] = geom_obj
    geo_init_, msg = geo_init(params)

    #If not, Compute geo at running_theory, using geo from iniitial_level as the starting point
    if 'getting reference geometry from directory'  in msg:
        geo = geo_init_

    else:
        #Obtain initial geometry
        if 'placeholder' in initial_thy_info: #placeholders are put if geo is to be read in from spc input file
            geo_init_ = geom_obj
            msg      = 'found initial geometry in species input file'
        else:
            params['theory_level'] = initial_thy_info
            params[    'geom_obj'] = geom_obj
            geo_init_, msg = geo_init(params)
        log.info('    | ' + msg)

        #Optimize from initial geometry to get reference geometry
        _, script_str, _, OPT_KWARGS = moldr.util.run_qchem_par(running_thy_info[0])
        params['theory_level'] = running_thy_info
        params[  'script_str'] = script_str
        params[   'overwrite'] = overwrite
        params[  'run_prefix'] = spc_run_path
        params[ 'save_prefix'] = spc_save_path
        params[    'geo_init'] = geo_init_
        del params['geom_obj']
        del params[  'prefix']
        geo, msg = moldr.driver.run_initial_geometry_opt(**params, **OPT_KWARGS)

    log.info('    | ' + msg)
    return geo 

def remove_imag(spcdic, thydic, run_prefix, save_prefix, overwrite):
    
    log.info('    | the initial geometries will be checked for imaginary frequencies')

    KICKOFF_SIZE = 0.1
    KICKOFF_BACKWARD = False

    spc_info      = get_spc_info(spcdic)
    thy_info      = get_thy_info(thydic)
    spc_run_path  = get_spc_run_path(run_prefix, spc_info) 
    thy_run_path  = get_thy_run_path(run_prefix, spc_info, thy_info) 
    spc_save_path = get_spc_save_path(save_prefix, spc_info) 
    thy_save_path = get_thy_save_path(save_prefix, spc_info, thy_info) 
    _, script_str, KWARGS, OPT_KWARGS = moldr.util.run_qchem_par(thy_info[0])
    
    params =  {    'spc_info': spc_info, 
               'theory_level': thy_info, 
                 'run_prefix': spc_run_path, 
                'save_prefix': spc_save_path, 
                 'script_str': script_str,
                  'overwrite': overwrite} 

    imag, geo, disp_xyzs = moldr.driver.run_check_imaginary(**params, **KWARGS)
    if imag:
        log.info('  | imaginary frequency detected, attempting to kick off')
        del params['run_prefix']
        del params['save_prefix']
        del params['overwrite']
        params[        'run_path'] = thy_run_path
        params['kickoff_backward'] = KICKOFF_BACKWARD
        params[    'kickoff_size'] = KICKOFF_SIZE
        params[        'opt_cart'] = False

        moldr.driver.run_kickoff_saddle(geo, disp_xyzs, **params, **OPT_KWARGS)
        log.info('  | removing saddlepoint hessian')

        run_fs = autofile.fs.run(thy_run_path)
        run_fs.leaf.remove([elstruct.Job.HESSIAN])
        moldr.driver.save_initial_geometry(
            spc_info=spc_info,
            theory_level=thy_info,
            run_prefix=spc_run_path,
            save_prefix=spc_save_path,
            )

    return

def ts_geo_ref(spc, spcs, reacs, prods, initial_thy_info, running_thy_info, run_prefix, save_prefix, overwrite):

    #I will do a better way of getting ts_info
    reaczmats = [automol.geom.zmatrix(spcs[x]['geoobj']) for x in reacs]
    prodzmats = [automol.geom.zmatrix(spcs[x]['geoobj']) for x in prods]
    ts_zma, dist_name, grid = ts_params(reaczmats, prodzmats)

    #Find TS filesystems
    ini_rxn_run_path, ini_rxn_save_path, ini_ts_info = rxn_file_paths(run_prefix, save_prefix, reacs, prods, spcs, initial_thy_info)
    run_rxn_run_path, run_rxn_save_path, run_ts_info = rxn_file_paths(run_prefix, save_prefix, reacs, prods, spcs, running_thy_info, initial_thy_info)
    orb_restr = moldr.util.orbital_restriction(run_ts_info, running_thy_info)
    ref_level = running_thy_info[1:3]
    ref_level.append(orb_restr)
    
    run_thy_run_fs = autofile.fs.theory(run_rxn_run_path)
    run_thy_run_fs.leaf.create(ref_level)
    run_thy_run_path = run_thy_run_fs.leaf.path(ref_level)
    
    run_thy_save_fs = autofile.fs.theory(run_rxn_save_path)
    run_thy_save_fs.leaf.create(ref_level)
    run_thy_save_path = run_thy_save_fs.leaf.path(ref_level)
    
    orb_restr = moldr.util.orbital_restriction(ini_ts_info, initial_thy_info)
    ref_level = initial_thy_info[1:3]
    ref_level.append(orb_restr)

    ini_thy_run_fs = autofile.fs.theory(ini_rxn_run_path)
    ini_thy_run_fs.leaf.create(ref_level)
    ini_thy_run_path = ini_thy_run_fs.leaf.path(ref_level)
    
    ini_thy_save_fs = autofile.fs.theory(ini_rxn_save_path)
    ini_thy_save_fs.leaf.create(ref_level)
    ini_thy_save_path = ini_thy_save_fs.leaf.path(ref_level)
    
    run_ts_run_path, run_ts_save_path, run_ts_save_fs = ts_paths(run_thy_run_path, run_thy_save_path)
    ini_ts_run_path, ini_ts_save_path, ini_ts_save_fs = ts_paths(ini_thy_run_path, ini_thy_save_path)
    
    log.info('   | running ts scan:')
 
    spc_save_fs = autofile.fs.species(save_prefix)

    try:
         zmat_init_ = run_ts_save_fs.trunk.file.zmatrix.read()
         geo_init_  = run_ts_save_fs.trunk.file.geometry.read()
         msg = 'getting reference geometry from directory {}'.format(run_ts_save_path)
         return geo

    except:
         zmat_init_ = ini_ts_save_fs.trunk.file.zmatrix.read()
         geo_init_  = ini_ts_save_fs.trunk.file.geometry.read()
         msg = 'initial zmat found in {}'.format(ini_ts_save_path)

         _, script_str, _, OPT_KWARGS = moldr.util.run_qchem_par(running_thy_info[0])
         params =  {    'spc_info': run_ts_info ,
                  'theory_level': running_thy_info, 
                     'prefix': run_ts_run_path,
                 'script_str': script_str,
                 'overwrite': overwrite,
                       'geom':  geo_init_,
                        'job':  'optimization',
                     'saddle': True}
         moldr.driver.run_job(**params)
         opt_ret = moldr.driver.read_job(job='optimization', prefix=run_ts_run_path)
         if opt_ret is not None:
             inf_obj, inp_str, out_str = opt_ret
             prog = inf_obj.prog
             method = inf_obj.method
             ene = elstruct.reader.energy(prog, method, out_str)
             geo = elstruct.reader.opt_geometry(prog, out_str)
             zma = elstruct.reader.opt_zmatrix(prog, out_str)

             print(" - Saving...")
             print(" - Save path: {}".format(run_ts_save_path))
             msg = 'ts optimization successful at running level' 

             run_ts_save_fs.trunk.file.energy.write(ene)
             run_ts_save_fs.trunk.file.geometry.write(geo)
             run_ts_save_fs.trunk.file.zmatrix.write(zma)
    return geo

def ts_run_task(tsk, spc, spcs, reacs, prods, thydic, initial_thy_info, run_prefix, save_prefix, overwrite):
    return

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

def get_thy_run_path(run_prefix, spc_info, thy_info):
    orb_restr = moldr.util.orbital_restriction(
        spc_info, thy_info)
    thy_lvl = thy_info[1:3]
    thy_lvl.append(orb_restr)
    spc_run_path = get_spc_run_path(run_prefix, spc_info)
    thy_run_fs   = autofile.fs.theory(spc_run_path)
    thy_run_fs.leaf.create(thy_lvl)
    thy_run_path = thy_run_fs.leaf.path(thy_lvl)
    return thy_run_path

def get_thy_save_path(save_prefix, spc_info, thy_info):
    orb_restr = moldr.util.orbital_restriction(
        spc_info, thy_info)
    thy_lvl = thy_info[1:3]
    thy_lvl.append(orb_restr)
    spc_save_path = get_spc_save_path(save_prefix, spc_info)
    thy_save_fs = autofile.fs.theory(spc_save_path)
    thy_save_fs.leaf.create(thy_lvl)
    thy_save_path = thy_save_fs.leaf.path(thy_lvl)
    return thy_save_path

def rxn_file_paths(run_prefix, save_prefix, reacs, prods, spcs, thy_info, initial_thy_info=None):
    ts_mul = automol.mult.ts.low([spcs[spc]['mult'] for spc in reacs], [spcs[spc]['mult'] for spc in reacs])
    ts_chg = sum([spcs[spc]['charge'] for spc in reacs])
    ts_info = ('', ts_chg, ts_mul)

    rxn_ichs = [[],[]] 
    rxn_chgs = [[],[]]
    rxn_muls = [[],[]]
    for spc in reacs:
         rxn_ichs[0].append(spcs[spc][ 'inchi'])
         rxn_chgs[0].append(spcs[spc]['charge'])
         rxn_muls[0].append(spcs[spc][  'mult'])
    for spc in prods:
         rxn_ichs[1].append(spcs[spc][ 'inchi'])
         rxn_chgs[1].append(spcs[spc]['charge'])
         rxn_muls[1].append(spcs[spc][  'mult'])

    # check direction of reaction
    log.info('    | checking exothermicity of reaction')
    try:
        rxn_exo = moldr.util.reaction_energy(
            save_prefix, rxn_ichs, rxn_chgs, rxn_muls, thy_info)
    except:
        rxn_exo = moldr.util.reaction_energy(
            save_prefix, rxn_ichs, rxn_chgs, rxn_muls, initial_thy_info)
    log.info('    | reaction is {:.2f}'.format(rxn_exo))
    if rxn_exo > 0:
        rxn_ichs =  rxn_ichs[::-1]
        rxn_chgs =  rxn_chgs[::-1]
        rxn_muls =  rxn_muls[::-1]
        log.info('    | ts search will be performed in reverse direction')
    
    # set up the filesystem
    is_rev = autofile.system.reaction_is_reversed(
        rxn_ichs, rxn_chgs, rxn_muls)
    rxn_ichs, rxn_chgs, rxn_muls = autofile.system.sort_together(
        rxn_ichs, rxn_chgs, rxn_muls)
    log.info("    | The reaction direction is {}"
          .format('backward' if is_rev else 'forward'))

    ts_info = ['', ts_chg, ts_mul]

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
    return rxn_run_path, rxn_save_path, ts_info

def ts_paths(thy_run_path, thy_save_path):

    print('thy_run_path in ts_opt:', thy_run_path)
    ts_run_fs = autofile.fs.ts(thy_run_path)
    ts_run_fs.trunk.create()
    ts_run_path = ts_run_fs.trunk.path()
    print('ts_run_path:', ts_run_path)
    
    ts_save_fs = autofile.fs.ts(thy_save_path)
    ts_save_fs.trunk.create()
    ts_save_path = ts_save_fs.trunk.path()
    print('ts_save_path:', ts_save_path)
    
    return ts_run_path, ts_save_path, ts_save_fs


def ts_params(rct_zmas, prd_zmas):
    typ = None
    ret = automol.zmatrix.ts.beta_scission(rct_zmas, prd_zmas)
    if ret and typ is None:
        typ = 'beta scission'
        ts_zma, dist_name, tors_names = ret
        log.info('    | beta scission')
        log.info('    | ts zma:', ts_zma)
        log.info('    | dist name:', dist_name)
        log.info('    | tors names:', tors_names)

    ret = automol.zmatrix.ts.addition(rct_zmas, prd_zmas)
    if ret and typ is None:
        typ = 'addition'
        ts_zma, dist_name, tors_names = ret
        log.info('    | addn')
        log.info('    | ts zma:')
        log.info(ts_zma)
        log.info('    | dist name:')
        log.info(dist_name)
        log.info('    | tors names:')
        log.info(tors_names)

    # fix this later
    # ret = automol.zmatrix.ts.hydrogen_abstraction(rct_zmas, prd_zmas,
    #                                               sigma=True)
    ret = automol.zmatrix.ts.hydrogen_abstraction(rct_zmas, prd_zmas,
                                                  sigma=False)
    if ret and typ is None:
        typ = 'hydrogen abstraction'
        ts_zma, dist_name, tors_names = ret
        log.info('    | H abs')
        log.info('    | ts zma:', ts_zma)
        log.info('    | dist name:', dist_name)
        log.info('    | tors names:', tors_names)

    if typ is None:
        log.info("    | Failed to classify reaction.")
    else:
        log.info("    | Type: {}".format(typ))

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

        if typ in ('beta scission', 'addition'):
            rmin = 1.4 * ANG2BOHR
            rmin = 2.8 * ANG2BOHR
            if bnd_len_key in bnd_len_dct:
                bnd_len = bnd_len_dct[bnd_len_key]
                rmin = bnd_len + 0.2 * ANG2BOHR
                rmax = bnd_len + 1.6 * ANG2BOHR
        elif typ == 'hydrogen abstraction':
            rmin = 0.7 * ANG2BOHR
            rmax = 2.2 * ANG2BOHR
            if bnd_len_key in bnd_len_dct:
                bnd_len = bnd_len_dct[bnd_len_key]
                rmin = bnd_len
                rmax = bnd_len + 1.0 * ANG2BOHR

        npoints = 8
        grid = numpy.linspace(rmin, rmax, npoints)
        return ts_zma, dist_name, grid

def find_ts(run_prefix, save_prefix, reacs, prods, spcs, thy_info, overwrite):

    log.info('   | prepping ts scan:')
    script_str, opt_script_str, KWARGS, OPT_KWARGS = moldr.util.run_qchem_par(thy_info[0])
    
    
    rxn_run_path, rxn_save_path, ts_info = rxn_file_paths(run_prefix, save_prefix, reacs, prods, spcs, thy_info)
    reaczmats = [automol.geom.zmatrix(spcs[x]['geoobj']) for x in reacs]
    prodzmats = [automol.geom.zmatrix(spcs[x]['geoobj']) for x in prods]
    ts_zma, dist_name, grid = ts_params(reaczmats, prodzmats)
    orb_restr = moldr.util.orbital_restriction(ts_info, thy_info)
    ref_level = thy_info[1:3]
    ref_level.append(orb_restr)
    
    thy_run_fs = autofile.fs.theory(rxn_run_path)
    thy_run_fs.leaf.create(ref_level)
    thy_run_path = thy_run_fs.leaf.path(ref_level)
    
    thy_save_fs = autofile.fs.theory(rxn_save_path)
    thy_save_fs.leaf.create(ref_level)
    thy_save_path = thy_save_fs.leaf.path(ref_level)
    
    log.info('   | running ts scan:')
 
    moldr.driver.run_scan(
        zma=ts_zma,
        spc_info=ts_info,
        theory_level=thy_info,
        grid_dct={dist_name: grid},
        run_prefix=thy_run_path,
        save_prefix=thy_save_path,
        script_str=script_str,
        overwrite=overwrite,
        update_guess=False,
        reverse_sweep=False,
        **OPT_KWARGS
    )
    
    moldr.driver.save_scan(
        run_prefix=thy_run_path,
        save_prefix=thy_save_path,
        coo_names=[dist_name],
    )
    
    scn_save_fs = autofile.fs.scan(thy_save_path)
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

    ts_run_path, ts_save_path, ts_save_fs = ts_paths(thy_run_path, thy_save_path)
    print('starting ts optimization')
    print('theory_level=:', thy_info)
    print('ts_run_path=:', ts_run_path)

    moldr.driver.run_job(
        job='optimization',
        script_str=script_str,
        prefix=ts_run_path,
        geom=max_zma,
        spc_info=ts_info,
        theory_level=thy_info,
        saddle=True,
        overwrite=overwrite,
        **OPT_KWARGS,
    )
    opt_ret = moldr.driver.read_job(
        job='optimization',
        prefix=ts_run_path,
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

    else:                                                      
        geo = 'failed'                                         
        zma = 'failed'
    return geo, zma
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
            SP_SCRIPT_STR, OPT_SCRIPT_STR, KWARGS, OPT_KWARGS = moldr.util.run_qchem_par(prog)

            orb_restr = moldr.util.orbital_restriction(
                spc_info[name], run_opt_levels[opt_level_idx])
            thy_level = run_opt_levels[opt_level_idx][1:3]
            thy_level.append(orb_restr)

            # a. conformer sampling
            spc_run_fs = autofile.fs.species(run_prefix)
            spc_run_fs.leaf.create(spc_info[name])
            spc_run_path = spc_run_fs.leaf.path(spc_info[name])

            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs.leaf.create(spc_info[name])
            spc_save_path = spc_save_fs.leaf.path(spc_info[name])

            thy_run_fs = autofile.fs.theory(spc_run_path)
            thy_run_fs.leaf.create(thy_level)
            thy_run_path = thy_run_fs.leaf.path(thy_level)

            thy_save_fs = autofile.fs.theory(spc_save_path)
            thy_save_fs.leaf.create(thy_level)
            thy_save_path = thy_save_fs.leaf.path(thy_level)

            # this uses theory run path - should start with a check in save path to see if initial geometry has already been saved
            # eventually theory data will be removed
            # also may need to remove hessian etc from saved geometry ...
            if run_ini_geom:                                                        

                geo_init = moldr.util.reference_geometry(
                    spc_info=spc_info[name],
                    theory_level=run_opt_levels[opt_level_idx],
                    prefix=save_prefix,
                    geom_dct=geom_dct)

                geo = moldr.driver.run_initial_geometry_opt(
                    spc_info=spc_info[name],
                    theory_level=run_opt_levels[opt_level_idx],
                    run_prefix=spc_run_path,
                    save_prefix=spc_save_path,
                    script_str=OPT_SCRIPT_STR,
                    overwrite=overwrite,
                    geo_init=geo_init,
                    **OPT_KWARGS,
                )

            if run_remove_imag:
                imag, geo, disp_xyzs = moldr.driver.run_check_imaginary(
                    spc_info=spc_info[name],
                    theory_level=run_opt_levels[opt_level_idx],
                    run_prefix=spc_run_path,
                    save_prefix=spc_save_path,
                    script_str=OPT_SCRIPT_STR,
                    overwrite=overwrite,
                    **KWARGS,
                )
                if imag:
                    moldr.driver.run_kickoff_saddle(
                        geo, disp_xyzs,
                        spc_info=spc_info[name],
                        theory_level=run_opt_levels[opt_level_idx],
                        run_path=thy_run_path,
                        script_str=OPT_SCRIPT_STR,
                        kickoff_backward=kickoff_backward,
                        kickoff_size=kickoff_size,
                        opt_cart=False,
                        **OPT_KWARGS)
                    print('removing saddlepoint hessian')

                    run_fs = autofile.fs.run(thy_run_path)
                    run_fs.leaf.remove([elstruct.Job.HESSIAN])
                    save_fs = autofile.fs.save(thy_save_path)
                    save_fs.leaf.remove('hess')

                    moldr.driver.save_initial_geometry(
                        spc_info=spc_info[name],
                        theory_level=run_opt_levels[opt_level_idx],
                        run_prefix=spc_run_path,
                        save_prefix=spc_save_path,
                    )

            if run_conf_samp:
                moldr.driver.conformer_sampling(
                    spc_info=spc_info[name],
                    theory_level=run_opt_levels[opt_level_idx],
                    run_prefix=spc_run_path,
                    save_prefix=spc_save_path,
                    script_str=OPT_SCRIPT_STR,
                    overwrite=overwrite,
                    nsamp_par=nsamp_conf_par,
                    **OPT_KWARGS,
                )

                if run_min_grad:
                    moldr.driver.run_minimum_energy_gradient(
                        spc_info=spc_info[name],
                        theory_level=run_opt_levels[opt_level_idx],
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=overwrite,
                        **KWARGS,
                    )

                if run_min_hess:
                    moldr.driver.run_minimum_energy_hessian(
                        spc_info=spc_info[name],
                        theory_level=run_opt_levels[opt_level_idx],
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=overwrite,
                        **KWARGS,
                    )

                if run_min_vpt2:
                    moldr.driver.run_minimum_energy_vpt2(
                        spc_info=spc_info[name],
                        theory_level=run_opt_levels[opt_level_idx],
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=overwrite,
                        **KWARGS,
                    )

                if run_conf_scan:
                    moldr.driver.hindered_rotor_scans(
                        spc_info=spc_info[name],
                        theory_level=run_opt_levels[opt_level_idx],
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=overwrite,
                        scan_increment=scan_increment,
                        **OPT_KWARGS,
                    )

                if run_conf_grad:
                    moldr.driver.run_conformer_gradients(
                        spc_info=spc_info[name],
                        theory_level=run_opt_levels[opt_level_idx],
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=overwrite,
                        **KWARGS,
                    )

                if run_conf_hess:
                    moldr.driver.run_conformer_hessians(
                        spc_info=spc_info[name],
                        theory_level=run_opt_levels[opt_level_idx],
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=overwrite,
                        **KWARGS,
                    )

            if run_tau_samp:
                moldr.driver.tau_sampling(
                    spc_info=spc_info[name],
                    theory_level=run_opt_levels[opt_level_idx],
                    run_prefix=spc_run_path,
                    save_prefix=spc_save_path,
                    script_str=OPT_SCRIPT_STR,
                    overwrite=overwrite,
                    nsamp_par=nsamp_tau_par,
                    **OPT_KWARGS,
                )

                if run_tau_grad:
                    moldr.driver.run_tau_gradients(
                        spc_info=spc_info[name],
                        theory_level=run_opt_levels[opt_level_idx],
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=overwrite,
                        **KWARGS,
                    )

                if run_tau_hess:
                    moldr.driver.run_tau_hessians(
                        spc_info=spc_info[name],
                        theory_level=run_opt_levels[opt_level_idx],
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=overwrite,
                        **KWARGS,
                    )

        for high_level_idx, _ in enumerate(run_high_levels):

            orb_restr = moldr.util.orbital_restriction(
                spc_info[name], ref_high_level)
            ref_level = ref_high_level[1:3]
            ref_level.append(orb_restr)

            spc_run_fs = autofile.fs.species(run_prefix)
            spc_run_fs.leaf.create(spc_info[name])
            spc_run_path = spc_run_fs.leaf.path(spc_info[name])
            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs.leaf.create(spc_info[name])
            spc_save_path = spc_save_fs.leaf.path(spc_info[name])

            ref_run_fs = autofile.fs.theory(spc_run_path)
            ref_run_fs.leaf.create(ref_level)
            ref_run_path = ref_run_fs.leaf.path(ref_level)
            ref_save_fs = autofile.fs.theory(spc_save_path)
            ref_save_fs.leaf.create(ref_level)
            ref_save_path = ref_save_fs.leaf.path(ref_level)

            min_cnf_locs = moldr.util.min_energy_conformer_locators(
                ref_save_path)
            cnf_run_fs = autofile.fs.conformer(ref_run_path)
            cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
            cnf_save_fs = autofile.fs.conformer(ref_save_path)
            cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
            min_cnf_geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)

            # evaluate the high level energy and save it

            prog = run_high_levels[high_level_idx][0]
            SP_SCRIPT_STR, OPT_SCRIPT_STR, KWARGS, OPT_KWARGS = (
                moldr.util.run_qchem_par(prog))
            if run_hl_min_ene:
                moldr.driver.run_single_point_energy(
                    geo=min_cnf_geo,
                    spc_info=spc_info[name],
                    theory_level=run_high_levels[high_level_idx],
                    run_prefix=cnf_run_path,
                    save_prefix=cnf_save_path,
                    script_str=SP_SCRIPT_STR,
                    overwrite=overwrite,
                    **KWARGS,
                )

            if run_hl_conf_ene:
                # add cycle over conformer geometries
                cnf_run_fs = autofile.fs.conformer(run_prefix)
                cnf_save_fs = autofile.fs.conformer(save_prefix)

                cnf_locs_lst = cnf_save_fs.leaf.existing()
                for locs in cnf_locs_lst:
                    cnf_run_path = cnf_run_fs.leaf.path(locs)
                    cnf_save_path = cnf_save_fs.leaf.path(locs)
                    cnf_geo = cnf_save_fs.leaf.file.geometry.read(locs)

                    moldr.driver.run_single_point_energy(
                        geo=min_cnf_geo,
                        spc_info=spc_info[name],
                        theory_level=run_high_levels[high_level_idx],
                        run_prefix=cnf_run_path,
                        save_prefix=cnf_save_path,
                        script_str=SP_SCRIPT_STR,
                        overwrite=overwrite,
                        **KWARGS,
                    )

            if run_hl_conf_opt:
                # add cycle over conformer geometries
                cnf_run_fs = autofile.fs.conformer(run_prefix)
                cnf_save_fs = autofile.fs.conformer(save_prefix)

                cnf_locs_lst = cnf_save_fs.leaf.existing()
                for locs in cnf_locs_lst:
                    cnf_run_path = cnf_run_fs.leaf.path(locs)
                    cnf_save_path = cnf_save_fs.leaf.path(locs)
                    cnf_geo = cnf_save_fs.leaf.file.geometry.read(locs)

                    moldr.driver.run_optimization(
                        geo=min_cnf_geo,
                        spc_info=spc_info[name],
                        theory_level=run_high_levels[high_level_idx],
                        run_prefix=cnf_run_path,
                        save_prefix=cnf_save_path,
                        script_str=SP_SCRIPT_STR,
                        overwrite=overwrite,
                        **KWARGS,
                    )



def ts_qchem(
        rct_names_lst, prd_names_lst, smi_dct, chg_dct, mul_dct,
        run_opt_levels, ref_high_level,
        run_high_levels, geom_dct, run_prefix, save_prefix, qchem_flags,
        nsamp_pars, scan_increment, kickoff_pars, overwrite,
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
    for rct_names, prd_names in zip(rct_names_lst, prd_names_lst):
        # print the CHEMKIN reaction name for reference
        rxn_name = '='.join(['+'.join(rct_names), '+'.join(prd_names)])
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
        ts_mul = automol.mult.ts.low(rct_muls, prd_muls)
        ts_chg = sum(rct_chgs)
        print('ts_chg test:', ts_chg)
        ts_info = ('', ts_chg, ts_mul)

        # theory
        for opt_level_idx, _ in enumerate(run_opt_levels):
            prog = run_opt_levels[opt_level_idx][0]
            SCRIPT_STR, OPT_SCRIPT_STR, KWARGS, OPT_KWARGS = (
                moldr.util.run_qchem_par(prog))

            ts_orb_restr = moldr.util.orbital_restriction(
                ts_info, run_opt_levels[opt_level_idx])
            thy_level = run_opt_levels[opt_level_idx][1:3]
            thy_level.append(ts_orb_restr)

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
                geo = moldr.util.reference_geometry(
                    rct_info, run_opt_levels[opt_level_idx], save_prefix,
                    geom_dct)
                rct_geos.append(geo)

            prd_geos = []
            for ich, chg, mul in zip(prd_ichs, prd_chgs, prd_muls):
                prd_info = [ich, chg, mul]
                geo = moldr.util.reference_geometry(
                    prd_info, run_opt_levels[opt_level_idx], save_prefix,
                    geom_dct)
                prd_geos.append(geo)

            # determine the transition state z-matrix
            # replace this with save values if they are available
            rct_zmas = list(map(automol.geom.zmatrix, rct_geos))
            prd_zmas = list(map(automol.geom.zmatrix, prd_geos))

            typ = None

            # # (migrations are not yet implemented)
            # ret = automol.zmatrix.ts.hydrogen_migration(rct_zmas, prd_zmas)
            # if ret and typ is None:
            #     typ = 'hydrogen migration'
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

                if typ in ('beta scission', 'addition'):
                    rmin = 1.4 * ANG2BOHR
                    rmin = 2.8 * ANG2BOHR
                    if bnd_len_key in bnd_len_dct:
                        bnd_len = bnd_len_dct[bnd_len_key]
                        rmin = bnd_len + 0.2 * ANG2BOHR
                        rmax = bnd_len + 1.6 * ANG2BOHR
                elif typ == 'hydrogen abstraction':
                    rmin = 0.7 * ANG2BOHR
                    rmax = 2.2 * ANG2BOHR
                    if bnd_len_key in bnd_len_dct:
                        bnd_len = bnd_len_dct[bnd_len_key]
                        rmin = bnd_len
                        rmax = bnd_len + 1.0 * ANG2BOHR

                npoints = 8
                grid = numpy.linspace(rmin, rmax, npoints)

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

                ts_chg = 0
                for rct_chg in rct_chgs:
                    ts_chg += rct_chg
                ts_info = ['', ts_chg, ts_mul]

                rxn_run_fs = autofile.fs.reaction(run_prefix)
                rxn_run_fs.leaf.create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
                rxn_run_path = rxn_run_fs.leaf.path(
                    [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

                rxn_ichs = tuple(map(tuple, rxn_ichs))
                rxn_chgs = tuple(map(tuple, rxn_chgs))
                rxn_muls = tuple(map(tuple, rxn_muls))
                print('rxn_save test0', rxn_ichs, rxn_chgs, rxn_muls, ts_mul)
                print(save_prefix)
                rxn_save_fs = autofile.fs.reaction(save_prefix)
                rxn_save_fs.leaf.create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
                rxn_save_path = rxn_save_fs.leaf.path(
                    [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
                print(rxn_save_path)

                orb_restr = moldr.util.orbital_restriction(
                    ts_info, run_opt_levels[opt_level_idx])
                ref_level = run_opt_levels[opt_level_idx][1:3]
                ref_level.append(orb_restr)
                print('ref level test:', ref_level)

                thy_run_fs = autofile.fs.theory(rxn_run_path)
                thy_run_fs.leaf.create(ref_level)
                thy_run_path = thy_run_fs.leaf.path(
                    ref_level)

                thy_save_fs = autofile.fs.theory(rxn_save_path)
                thy_save_fs.leaf.create(ref_level)
                thy_save_path = thy_save_fs.leaf.path(ref_level)
                print(thy_save_path)

                print('entering run_scan:')

                moldr.driver.run_scan(
                    zma=ts_zma,
                    spc_info=ts_info,
                    theory_level=run_opt_levels[opt_level_idx],
                    grid_dct={dist_name: grid},
                    run_prefix=thy_run_path,
                    save_prefix=thy_save_path,
                    script_str=SCRIPT_STR,
                    overwrite=overwrite,
                    update_guess=False,
                    reverse_sweep=False,
                    **OPT_KWARGS
                )

                moldr.driver.save_scan(
                    run_prefix=thy_run_path,
                    save_prefix=thy_save_path,
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
                print('ts_run_path:', ts_run_path)

                ts_save_fs = autofile.fs.ts(thy_save_path)
                ts_save_fs.trunk.create()
                ts_save_path = ts_save_fs.trunk.path()
                print('ts_save_path:', ts_save_path)

                print('starting ts optimization')
                print('theory_level=:', run_opt_levels[opt_level_idx])
                print('ts_run_path=:', ts_run_path)
                moldr.driver.run_job(
                    job='optimization',
                    script_str=SCRIPT_STR,
                    prefix=ts_run_path,
                    geom=max_zma,
                    spc_info=ts_info,
                    theory_level=run_opt_levels[opt_level_idx],
                    saddle=True,
                    overwrite=overwrite,
                    **OPT_KWARGS,
                )
                opt_ret = moldr.driver.read_job(
                    job='optimization',
                    prefix=ts_run_path,
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
                        moldr.driver.run_scan(
                            zma=zma,
                            spc_info=ts_info,
                            theory_level=run_opt_levels[opt_level_idx],
                            grid_dct={tors_name: tors_grid},
                            run_prefix=ts_run_path,
                            save_prefix=ts_save_path,
                            script_str=SCRIPT_STR,
                            overwrite=overwrite,
                            saddle=True,
                            **OPT_KWARGS,
                        )

                        moldr.driver.save_scan(
                            run_prefix=ts_run_path,
                            save_prefix=ts_save_path,
                            coo_names=[tors_name],
                        )

                    hind_rot_dct = {}
                    scn_run_fs = autofile.fs.scan(ts_run_path)
                    scn_save_fs = autofile.fs.scan(ts_save_path)

                    min_ene = ts_save_fs.trunk.file.energy.read()
                    for tors_name in tors_names:
                        enes = [scn_save_fs.leaf.file.energy.read(locs)
                                for locs in scn_save_fs.leaf.existing([tors_name])]
                        enes = numpy.subtract(enes, min_ene)
                        hind_rot_dct[tors_name] = enes*EH2KCAL

                    print('ts hindered rotor potential')
                    print(hind_rot_dct)

                if run_ts_conf_samp:
                    print('ts_conf_test:')
                    print(ts_info)
                    print(run_opt_levels[opt_level_idx])
                    print(ts_run_path)
                    print(ts_save_path)
                    print(rxn_run_path)
                    print(rxn_save_path)
                    print(tors_names)
                    moldr.driver.ts_conformer_sampling(
                        spc_info=ts_info,
                        geo=geo,
                        zma=zma,
                        tors_names=tors_names,
                        theory_level=run_opt_levels[opt_level_idx],
                        run_prefix=rxn_run_path,
                        save_prefix=rxn_save_path,
                        script_str=SCRIPT_STR,
                        overwrite=overwrite,
                        saddle=True,
                        nsamp_par=nsamp_ts_conf_par,
                        **OPT_KWARGS,
                    )

                    if run_ts_min_grad:
                        moldr.driver.run_minimum_energy_gradient(
                            spc_info=ts_info,
                            theory_level=run_opt_levels[opt_level_idx],
                            run_prefix=ts_run_path,
                            save_prefix=ts_save_path,
                            script_str=SCRIPT_STR,
                            overwrite=overwrite,
                            **KWARGS,
                        )

                    if run_ts_min_hess:
                        moldr.driver.run_minimum_energy_hessian(
                            spc_info=ts_info,
                            theory_level=run_opt_levels[opt_level_idx],
                            run_prefix=ts_run_path,
                            save_prefix=ts_save_path,
                            script_str=SCRIPT_STR,
                            overwrite=overwrite,
                            **KWARGS,
                        )

                    if run_ts_min_vpt2:
                        moldr.driver.run_minimum_energy_vpt2(
                            spc_info=ts_info,
                            theory_level=run_opt_levels[opt_level_idx],
                            run_prefix=ts_run_path,
                            save_prefix=ts_save_path,
                            script_str=SCRIPT_STR,
                            overwrite=overwrite,
                            saddle=True,
                            **KWARGS,
                        )

                    if run_ts_tau_samp:

                        moldr.driver.save_tau(
                            run_prefix=thy_run_path,
                            save_prefix=thy_save_path,
                        )

                        zma = ts_save_fs.trunk.file.zmatrix.read()
                        tors_ranges = automol.zmatrix.torsional_sampling_ranges(
                            zma, tors_names)
                        tors_range_dct = dict(zip(tors_names, tors_ranges))

                        moldr.driver.run_tau(
                            zma=zma,
                            spc_info=ts_info,
                            theory_level=run_opt_levels[opt_level_idx],
                            nsamp=nsamp_ts_tau_par,
                            tors_range_dct=tors_range_dct,
                            run_prefix=thy_run_path,
                            save_prefix=thy_save_path,
                            script_str=SCRIPT_STR,
                            overwrite=overwrite,
                            saddle=True,
                            **OPT_KWARGS,
                        )
    #                        saddle=True,
    # used to have saddle=True in call, but this is not used. Probably a bug.
    # Need to pass saddle information somewhere and use it

                        moldr.driver.save_tau(
                            run_prefix=thy_run_path,
                            save_prefix=thy_save_path,
                        )

                    if run_ts_kicks_qchem:
                        ret = moldr.driver.read_job(
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
                            moldr.driver.run_kickoff_saddle(
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
                            moldr.driver.run_kickoff_saddle(
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
                    ref_level = ref_high_level[1:3]
                    ref_level.append(orb_restr)
                    print('ref level test:', ref_level)

                    ref_run_fs = autofile.fs.theory(rxn_run_path)
                    ref_run_fs.leaf.create(ref_level)
                    ref_run_path = ref_run_fs.leaf.path(ref_level)
                    ref_save_fs = autofile.fs.theory(rxn_save_path)
                    ref_save_fs.leaf.create(ref_level)
                    ref_save_path = ref_save_fs.leaf.path(ref_level)

                    thy_run_fs = autofile.fs.theory(rxn_run_path)
                    thy_run_fs.leaf.create(ref_level)
                    thy_run_path = thy_run_fs.leaf.path(
                        ref_level)

                    thy_save_fs = autofile.fs.theory(rxn_save_path)
                    thy_save_fs.leaf.create(ref_level)
                    thy_save_path = thy_save_fs.leaf.path(ref_level)

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

                    min_cnf_locs = moldr.util.min_energy_conformer_locators(
                        ts_save_path)
                    cnf_run_fs = autofile.fs.conformer(ts_run_path)
                    cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
                    cnf_save_fs = autofile.fs.conformer(ts_save_path)
                    cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
                    min_cnf_geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)

                    prog = run_high_levels[high_level_idx][0]
                    SP_SCRIPT_STR, OPT_SCRIPT_STR, KWARGS, OPT_KWARGS = (
                        moldr.util.run_qchem_par(prog))
                    if run_ts_hl_min_ene:
                        moldr.driver.run_single_point_energy(
                            geo=min_cnf_geo,
                            spc_info=ts_info,
                            theory_level=run_high_levels[high_level_idx],
                            run_prefix=cnf_run_path,
                            save_prefix=cnf_save_path,
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
            SP_SCRIPT_STR, OPT_SCRIPT_STR, KWARGS, OPT_KWARGS = moldr.util.run_qchem_par(prog)

            geos = []
            for ich, chg, mul in zip(ichs, chgs, muls):
                spc_info = [ich, chg, mul]
                orb_restr = moldr.util.orbital_restriction(spc_info, run_opt_levels[opt_level_idx][0])
                geo = moldr.util.reference_geometry(
                    spc_info=spc_info, theory_level=run_opt_levels[opt_level_idx],
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
                    theory_level=run_opt_levels[opt_level_idx],
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


def get_spc_info(spcdic):
    err_msg = ''
    props   = ['inchi', 'charge', 'mult']
    for i, prop in enumerate(props):
        if prop in spcdic:
            props[i] = spcdic[prop]
        else:
            err_msg = prop
    if err_msg:
         log.error('No {} found'.format(err_msg))
    return props

def get_thy_info(lvldic):
    err_msg = ''
    info    = ['program', 'method', 'basis', 'orb_res']
    for i, inf in enumerate(info):
        if inf in lvldic:
            info[i] = lvldic[inf]
        else:
            err_msg = inf
    if err_msg:
         log.error('No {} found'.format(err_msg))
    return info
 

