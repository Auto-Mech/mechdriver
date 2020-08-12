""" eletronic structure routines modules
"""

import sys
import itertools
import importlib
import autofile
import automol
from routines.es import runner as es_runner
from routines.es.findts import run as runts
# from routines.es.newts import run as runts
from routines.es._routines import conformer
from routines.es._routines import geom
from routines.es._routines import hr
from routines.es._routines import tau
from routines.es._routines import irc
from lib import filesys
from lib import structure
from lib.structure import instab


# Dictionary of Electronic Structure Calculators
SP_MODULE = importlib.import_module('routines.es._routines.sp')
ES_TSKS = {
    'energy': SP_MODULE.run_energy,
    'grad': SP_MODULE.run_gradient,
    'hess': SP_MODULE.run_hessian,
    'vpt2': SP_MODULE.run_vpt2
}


def run_tsk(tsk, spc_dct, spc_name,
            thy_dct, es_keyword_dct,
            run_prefix, save_prefix):
    """ run an electronic structure task
    for generating a list of conformer or tau sampling geometries
    """

    # Print the head of the task
    print(('\n------------------------------------------------' +
           '--------------------------------------'))
    print('\nTask:', tsk, spc_name)
    print('\nOptions for electronic structure task:')
    for key, val in es_keyword_dct.items():
        print('{} = {}'.format(key, val))
    print('')

    # If species is unstable, set task to 'none'
    ini_thy_info = filesys.inf.get_es_info(
        es_keyword_dct['inplvl'], thy_dct)
    stable = instab.check_unstable_species(
        tsk, spc_dct, spc_name, ini_thy_info, save_prefix)
    print()

    if stable:

        # Set keys
        saddle = bool('ts_' in spc_name)
        # vdw = bool('vdw' in spc_name)
        spc = spc_dct[spc_name]

        # Get stuff from task
        job = tsk.split('_', 1)[1]

        # Run the task if an initial geom exists
        if 'init' in tsk:
            _ = geom_init(
                spc, thy_dct, es_keyword_dct,
                run_prefix, save_prefix, saddle)
        elif 'conf' in tsk:
            conformer_tsk(
                job, spc_dct, spc_name, thy_dct, es_keyword_dct,
                run_prefix, save_prefix, saddle)
        elif 'tau' in tsk:
            tau_tsk(
                job, spc_dct, spc_name, thy_dct, es_keyword_dct,
                run_prefix, save_prefix, saddle)
        elif 'hr' in tsk:
            hr_tsk(
                job, spc_dct, spc_name, thy_dct, es_keyword_dct,
                run_prefix, save_prefix, saddle)
        elif 'irc' in tsk:
            irc_tsk(
                job, spc_dct, spc_name, thy_dct, es_keyword_dct,
                run_prefix, save_prefix, saddle)
        elif 'find' in tsk:
            runts(
                job, spc_dct, spc_name, thy_dct, es_keyword_dct,
                run_prefix, save_prefix)

    else:
        print('\nSkipping task for unstable species...')


# FUNCTIONS FOR SAMPLING AND SCANS #
def geom_init(spc, thy_dct, es_keyword_dct,
                  run_prefix, save_prefix, saddle):
    """ Find the initial geometry
    """

    # Set the spc_info
    spc_info = filesys.inf.get_spc_info(spc)

    # Get es options
    [kickoff_size, kickoff_backward] = spc['kickoff']
    overwrite = es_keyword_dct['overwrite']
    # retryfail = es_keyword_dct['retryfail']

    # Get the thery info
    ini_thy_info = filesys.inf.get_es_info(
        es_keyword_dct['inplvl'], thy_dct)
    thy_info = filesys.inf.get_es_info(
        es_keyword_dct['runlvl'], thy_dct)
    mod_thy_info = filesys.inf.modify_orb_restrict(spc_info, thy_info)
    mod_ini_thy_info = filesys.inf.modify_orb_restrict(spc_info, ini_thy_info)

    # Set the filesystem objects
    thy_run_fs, thy_run_path = filesys.build.spc_thy_fs_from_root(
        run_prefix, spc_info, mod_thy_info)
    thy_save_fs, thy_save_path = filesys.build.spc_thy_fs_from_root(
        save_prefix, spc_info, mod_thy_info)
    _, ini_thy_save_path = filesys.build.spc_thy_fs_from_root(
        save_prefix, spc_info, mod_ini_thy_info)
    cnf_run_fs, _ = filesys.build.cnf_fs_from_thy(
        thy_run_path, mod_thy_info, saddle=saddle)
    cnf_save_fs, _ = filesys.build.cnf_fs_from_thy(
        thy_save_path, mod_thy_info, saddle=saddle)
    ini_cnf_save_fs, ini_min_cnf_locs = filesys.build.cnf_fs_from_thy(
        ini_thy_save_path, mod_ini_thy_info, cnf='min')

    # Set the run filesystem
    if saddle:
        _, ts_path = filesys.build.ts_fs_from_thy(thy_run_path)
        run_fs = filesys.build.run_fs_from_prefix(ts_path)
    else:
        run_fs = filesys.build.run_fs_from_prefix(thy_run_path)

    # Set up the script
    _, opt_script_str, _, opt_kwargs = es_runner.qchem_params(
        *thy_info[0:2])

    # Get a reference geometry if one not found
    geo = geom.reference_geometry(
        spc, spc_info, mod_thy_info,
        thy_run_fs, thy_save_fs,
        cnf_run_fs, cnf_save_fs,
        ini_cnf_save_fs, ini_min_cnf_locs,
        run_fs,
        opt_script_str, overwrite,
        kickoff_size=kickoff_size,
        kickoff_backward=kickoff_backward,
        **opt_kwargs)

    return geo


def conformer_tsk(job, spc_dct, spc_name,
                      thy_dct, es_keyword_dct,
                      run_prefix, save_prefix, saddle):
    """ Launch tasks associated with conformers.
        Scan: Generate a set of conformer geometries and energies via
              random sampling over torsional coordinates
              following by optimization
        SP: Calculate ene, grad, ..
    """
    spc = spc_dct[spc_name]

    # Set the spc_info
    spc_info = filesys.inf.get_spc_info(spc)

    # Get es options
    overwrite = es_keyword_dct['overwrite']
    retryfail = es_keyword_dct['retryfail']

    # Modify the theory
    ini_thy_info = filesys.inf.get_es_info(
        es_keyword_dct['inplvl'], thy_dct)
    thy_info = filesys.inf.get_es_info(
        es_keyword_dct['runlvl'], thy_dct)
    mod_thy_info = filesys.inf.modify_orb_restrict(spc_info, thy_info)
    mod_ini_thy_info = filesys.inf.modify_orb_restrict(spc_info, ini_thy_info)

    # Set the filesystem objects
    if not saddle:
        _, ini_thy_run_path = filesys.build.spc_thy_fs_from_root(
            run_prefix, spc_info, mod_ini_thy_info)
        ini_thy_fs = filesys.build.spc_thy_fs_from_root(
            save_prefix, spc_info, mod_ini_thy_info)
        [_, ini_thy_save_path] = ini_thy_fs

    else:
        rxn_info = filesys.inf.rxn_info(
            spc['reacs'], spc['prods'], spc_dct)

        _, ini_thy_run_path = filesys.build.rxn_thy_fs_from_root(
            run_prefix, rxn_info, mod_ini_thy_info)
        _, ini_thy_save_path = filesys.build.rxn_thy_fs_from_root(
            save_prefix, rxn_info, mod_ini_thy_info)
        _, ini_thy_save_path = filesys.build.ts_fs_from_thy(
            ini_thy_save_path)
        _, ini_thy_run_path = filesys.build.ts_fs_from_thy(
            ini_thy_run_path)

    if job == 'samp':

        # Build the thy filesyste
        if not saddle:
            _, thy_run_path = filesys.build.spc_thy_fs_from_root(
                run_prefix, spc_info, mod_thy_info)
            thy_save_fs, thy_save_path = filesys.build.spc_thy_fs_from_root(
                save_prefix, spc_info, mod_thy_info)
        else:
            _, thy_run_path = filesys.build.rxn_thy_fs_from_root(
                run_prefix, rxn_info, mod_thy_info)
            thy_save_fs, thy_save_path = filesys.build.rxn_thy_fs_from_root(
                save_prefix, rxn_info, mod_thy_info)
            thy_save_fs, thy_save_path = filesys.build.ts_fs_from_thy(
                thy_save_path)
            _, thy_run_path = filesys.build.ts_fs_from_thy(thy_run_path)

        # Build conformer filesys
        cnf_run_fs, _ = filesys.build.cnf_fs_from_prefix(
            thy_run_path, mod_thy_info, cnf=None)
        cnf_save_fs, _ = filesys.build.cnf_fs_from_prefix(
            thy_save_path, mod_thy_info, cnf=None)

        # Build the ini zma filesys
        ini_cnf_save_fs, ini_cnf_save_locs = filesys.build.cnf_fs_from_prefix(
            ini_thy_save_path, mod_ini_thy_info, cnf='min')
        ini_cnf_save_paths = filesys.build.cnf_paths_from_locs(
            ini_cnf_save_fs, ini_cnf_save_locs)
        ini_zma_save_fs, _ = filesys.build.zma_fs_from_prefix(
            ini_cnf_save_paths[0], zma_idxs=[0])

        # Set up the run scripts
        _, opt_script_str, _, opt_kwargs = es_runner.qchem_params(
            *thy_info[0:2])

        # Set variables if it is a saddle
        two_stage = saddle
        rxn_class = spc['class'] if saddle else ''
        mc_nsamp = spc['mc_nsamp']

        # Read the geometry and zma from the ini file system
        if not saddle:
            geo = ini_cnf_save_fs[-1].file.geometry.read(ini_cnf_save_locs)
            zma = ini_zma_save_fs[-1].file.zmatrix.read([0])
            tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
            geo_path = ini_cnf_save_fs[-1].path(ini_cnf_save_locs)
        else:
            print('ini path', ini_thy_save_path)
            geo = ini_cnf_save_fs[-1].file.geometry.read(ini_cnf_save_locs)
            # geo = ini_thy_save_fs[0].file.geometry.read()
            zma = ini_zma_save_fs[-1].file.zmatrix.read([0])
            tors_names = spc['amech_ts_tors_names']
            # geo_path = thy_save_fs[0].path()
            geo_path = ini_cnf_save_fs[-1].path(ini_cnf_save_locs)

        print('Sampling done using geom from {}'.format(geo_path))

        # Run the sampling
        _ = conformer.conformer_sampling(
            zma, spc_info,
            mod_thy_info, thy_save_fs,
            cnf_run_fs, cnf_save_fs,
            opt_script_str, overwrite,
            saddle=saddle, nsamp_par=mc_nsamp,
            tors_names=tors_names,
            two_stage=two_stage, retryfail=retryfail,
            rxn_class=rxn_class, **opt_kwargs)

    elif job in ('energy', 'grad', 'hess', 'vpt2'):

        # Build conformer filesys
        cnf_run_fs, _ = filesys.build.cnf_fs_from_prefix(
            ini_thy_run_path, mod_ini_thy_info, cnf=None)
        cnf_save_fs, cnf_save_locs = filesys.build.cnf_fs_from_prefix(
            ini_thy_save_path, mod_ini_thy_info, cnf='min')

        # Check if locs exist, kill if it doesn't
        if not cnf_save_locs:
            print('*ERROR: No min-energy conformer found for level:',
                  ini_thy_save_path)
            sys.exit()

        # Set up the run scripts
        script_str, _, kwargs, _ = es_runner.qchem_params(
            *thy_info[0:2])

        # Run the job over all the conformers requested by the user
        for locs in cnf_save_locs:
            geo_run_path = cnf_run_fs[-1].path([locs])
            geo_save_path = cnf_save_fs[-1].path([locs])
            cnf_run_fs[-1].create([locs])
            zma, geo = filesys.inf.cnf_fs_zma_geo(cnf_save_fs, [locs])
            ES_TSKS[job](
                zma, geo, spc_info, mod_thy_info,
                cnf_save_fs, geo_run_path, geo_save_path, [locs],
                script_str, overwrite,
                retryfail=retryfail, **kwargs)


def tau_tsk(job, spc_dct, spc_name,
                thy_dct, es_keyword_dct,
                run_prefix, save_prefix, saddle):
    """ Energies, gradients, and hessians,
        for set of arbitrarily sampled torsional coordinates
        with all other coordinates optimized
    """
    spc = spc_dct[spc_name]

    # Set the spc_info
    spc_info = filesys.inf.get_spc_info(spc)

    # Get es options
    overwrite = es_keyword_dct['overwrite']
    retryfail = es_keyword_dct['retryfail']
    scan_increment = spc['hind_inc']
    nsamp_par = spc['tau_nsamp']

    # Script
    _, opt_script_str, _, opt_kwargs = es_runner.qchem_params(
        *thy_info[0:2])

    # Modify the theory
    ini_thy_info = filesys.inf.get_es_info(
        es_keyword_dct['inplvl'], thy_dct)
    thy_info = filesys.inf.get_es_info(
        es_keyword_dct['runlvl'], thy_dct)
    mod_thy_info = filesys.inf.modify_orb_restrict(spc_info, thy_info)
    mod_ini_thy_info = filesys.inf.modify_orb_restrict(spc_info, ini_thy_info)

    # Set the filesystem objects for thy info
    _, thy_run_path = filesys.build.spc_thy_fs_from_root(
        run_prefix, spc_info, mod_thy_info)
    _, thy_save_path = filesys.build.spc_thy_fs_from_root(
        save_prefix, spc_info, mod_thy_info)

    # Set the filesystem objects for ini thy_info
    _, ini_thy_save_path = filesys.build.spc_thy_fs_from_root(
        save_prefix, spc_info, mod_ini_thy_info)

    # Get the geom and energy of reference species
    ini_cnf_save_fs, ini_cnf_locs = filesys.build.cnf_fs_from_prefix(
        ini_thy_save_path, mod_ini_thy_info, cnf='min')
    zma, geo = filesys.inf.cnf_fs_zma_geo(
        ini_cnf_save_fs, ini_cnf_locs)
    ref_ene = ini_cnf_save_fs[-1].file.energy.read(ini_cnf_locs)

    # Bond key stuff
    if saddle:
        frm_bnd_keys, brk_bnd_keys = structure.ts.rxn_bnd_keys(
            ini_cnf_save_fs, ini_cnf_locs, zma_locs=[0])
    else:
        frm_bnd_keys, brk_bnd_keys = (), ()

    # Set up the torsion info
    dct_tors_names, _ = structure.tors.names_from_dct(
        spc, '1dhr')
    amech_spc_tors_names = structure.tors.names_from_geo(
        geo, '1dhr', saddle=saddle)
    if dct_tors_names:
        run_tors_names = dct_tors_names
    else:
        run_tors_names = amech_spc_tors_names
        print('Using tors names generated by AutoMech...')

    run_tors_names, _, _ = structure.tors.hr_prep(
        zma, tors_name_grps=run_tors_names,
        scan_increment=scan_increment, tors_model='1dhr',
        frm_bnd_keys=frm_bnd_keys, brk_bnd_keys=brk_bnd_keys)

    # Run the task if any torsions exist
    if run_tors_names:

        # Set up tau filesystem objects
        tau_run_fs, _ = filesys.build.tau_fs_from_thy(
            thy_run_path, tau='all')
        tau_save_fs, tau_save_locs = filesys.build.tau_fs_from_thy(
            thy_save_path, tau='all')

        if job == 'samp':

            # Set up the script
            _, opt_script_str, _, opt_kwargs = es_runner.qchem_params(
                *thy_info[0:2])

            # Run sampling
            tau.tau_sampling(
                zma, ref_ene,
                spc_info, run_tors_names, nsamp_par,
                mod_ini_thy_info,
                tau_run_fs, tau_save_fs,
                opt_script_str, overwrite,
                saddle=saddle, **opt_kwargs)

        elif job in ('energy', 'grad'):

            # Set up the run scripts
            script_str, _, kwargs, _ = es_runner.qchem_params(
                *thy_info[0:2])
            # Run the job over all the conformers requested by the user
            for locs in tau_save_locs:
                geo_run_path = tau_run_fs[-1].path(locs)
                geo_save_path = tau_save_fs[-1].path(locs)
                geo = tau_save_fs[-1].file.geometry.read(locs)
                zma = None
                tau_run_fs[-1].create(locs)
                ES_TSKS[job](
                    zma, geo, spc_info, mod_thy_info,
                    tau_save_fs, geo_run_path, geo_save_path, locs,
                    script_str, overwrite,
                    retryfail=retryfail, **kwargs)
                print('\n')

        elif job == 'hess':

            # Add the hessian max
            hessmax = es_keyword_dct['hessmax']

            # Set up the run scripts
            script_str, _, kwargs, _ = es_runner.qchem_params(
                *thy_info[0:2])
            # Run the job over all the conformers requested by the user
            hess_cnt = 0
            for locs in tau_save_locs:
                print('\nHESS Number {}'.format(hess_cnt+1))
                geo_run_path = tau_run_fs[-1].path(locs)
                geo_save_path = tau_save_fs[-1].path(locs)
                if not tau_save_fs[-1].file.hessian.exists(locs):
                    geo = tau_save_fs[-1].file.geometry.read(locs)
                    zma = None
                    tau_run_fs[-1].create(locs)
                    ES_TSKS[job](
                        zma, geo, spc_info, mod_thy_info,
                        tau_save_fs, geo_run_path, geo_save_path, locs,
                        script_str, overwrite,
                        retryfail=retryfail, **kwargs)
                    hess_cnt += 1
                else:
                    print('Hessian found and saved previously at {}'.format(
                        geo_save_path))

                    hess_cnt += 1
                if hess_cnt == hessmax:
                    break

    else:
        print('No torsional modes in the species')


def hr_tsk(job, spc_dct, spc_name,
               thy_dct, es_keyword_dct,
               run_prefix, save_prefix, saddle):
    """ run a scan over the specified torsional coordinates
    """

    spc = spc_dct[spc_name]

    # Set the spc_info
    spc_info = filesys.inf.get_spc_info(spc)

    # Modify the theory
    ini_thy_info = filesys.inf.get_es_info(
        es_keyword_dct['inplvl'], thy_dct)
    thy_info = filesys.inf.get_es_info(
        es_keyword_dct['runlvl'], thy_dct)
    mod_thy_info = filesys.inf.modify_orb_restrict(spc_info, thy_info)
    mod_ini_thy_info = filesys.inf.modify_orb_restrict(spc_info, ini_thy_info)
    
    # Script
    _, opt_script_str, _, opt_kwargs = es_runner.qchem_params(
        *mod_thy_info[0:2])

    # Set the filesystem objects
    if not saddle:
        _, ini_thy_run_path = filesys.build.spc_thy_fs_from_root(
            run_prefix, spc_info, mod_ini_thy_info)
        inifs = filesys.build.spc_thy_fs_from_root(
            save_prefix, spc_info, mod_ini_thy_info)
        [ini_thy_save_fs, ini_thy_save_path] = inifs
    else:
        rxn_info = filesys.inf.rxn_info(
            spc['reacs'], spc['prods'], spc_dct)

        _, ini_thy_run_path = filesys.build.rxn_thy_fs_from_root(
            run_prefix, rxn_info, mod_ini_thy_info)
        _, ini_thy_save_path = filesys.build.rxn_thy_fs_from_root(
            save_prefix, rxn_info, mod_ini_thy_info)
        inifs = filesys.build.ts_fs_from_thy(
            ini_thy_save_path)
        [ini_thy_save_fs, ini_thy_save_path] = inifs
        _, ini_thy_run_path = filesys.build.ts_fs_from_thy(
            ini_thy_run_path)

    # Build cnf filesys using the ini thy filesys (needed for all HR jobs)
    ini_cnf_run_fs, _ = filesys.build.cnf_fs_from_prefix(
        ini_thy_run_path, mod_ini_thy_info, cnf=None)
    ini_cnf_save_fs, ini_cnf_save_locs = filesys.build.cnf_fs_from_prefix(
        ini_thy_save_path, mod_ini_thy_info, cnf='min')
    ini_cnf_save_paths = filesys.build.cnf_paths_from_locs(
        ini_cnf_save_fs, ini_cnf_save_locs)
    ini_cnf_run_paths = filesys.build.cnf_paths_from_locs(
        ini_cnf_run_fs, ini_cnf_save_locs)

    # Create run fs if that directory has been deleted to run the jobs
    ini_cnf_run_fs[-1].create(ini_cnf_save_locs)

    # Get options from the dct or es options lst
    overwrite = es_keyword_dct['overwrite']
    retryfail = es_keyword_dct['retryfail']
    tors_model = es_keyword_dct['tors_model']
    resamp_min = es_keyword_dct['resamp_min']
    scan_increment = spc['hind_inc']

    # Bond key stuff
    if saddle:
        frm_bnd_keys, brk_bnd_keys = structure.ts.rxn_bnd_keys(
            ini_cnf_save_fs, ini_cnf_save_locs, zma_locs=[0])
    else:
        frm_bnd_keys, brk_bnd_keys = (), ()

    # Read fs for zma and geo
    zma, geo = filesys.inf.cnf_fs_zma_geo(ini_cnf_save_fs, ini_cnf_save_locs)

    # Set up the torsion info
    dct_tors_names, amech_sadpt_tors_names = structure.tors.names_from_dct(
        spc, tors_model)
    amech_spc_tors_names = structure.tors.names_from_geo(
        geo, tors_model, saddle=saddle)
    if dct_tors_names:
        run_tors_names = dct_tors_names
    else:
        run_tors_names = amech_spc_tors_names
        print('Using tors names generated by AutoMech...')

    run_tors_names, run_tors_grids, _ = structure.tors.hr_prep(
        zma, tors_name_grps=run_tors_names,
        scan_increment=scan_increment, tors_model=tors_model,
        frm_bnd_keys=frm_bnd_keys, brk_bnd_keys=brk_bnd_keys)

    # Run the task if any torsions exist
    if run_tors_names:

        # Set constraints
        if tors_model in ('1dhr', 'mdhr', 'mdhrv'):
            const_names = ()
        else:
            if tors_model == '1dhrf':
                if saddle:
                    const_names = tuple(
                        itertools.chain(*amech_sadpt_tors_names))
                else:
                    const_names = tuple(
                        itertools.chain(*amech_spc_tors_names))
            elif tors_model == '1dhrfa':
                coords = list(automol.zmatrix.coordinates(zma))
                const_names = tuple(coord for coord in coords)

        # Set if scan is rigid or relaxed
        scn_typ = 'relaxed' if tors_model != '1dhrfa' else 'rigid'

        # Set up ini filesystem for scans
        _, ini_zma_run_path = filesys.build.zma_fs_from_prefix(
            ini_cnf_run_paths[0], zma_idxs=[0])
        _, ini_zma_save_path = filesys.build.zma_fs_from_prefix(
            ini_cnf_save_paths[0], zma_idxs=[0])

        if job == 'scan':

            # Run the scan
            hr.hindered_rotor_scans(
                zma, spc_info, mod_thy_info, ini_thy_save_fs,
                ini_zma_run_path, ini_zma_save_path,
                run_tors_names, run_tors_grids,
                opt_script_str, overwrite,
                scn_typ=scn_typ,
                saddle=saddle, const_names=const_names,
                retryfail=retryfail, **opt_kwargs)

            # Read and print the potential
            ref_ene = ini_cnf_save_fs[-1].file.energy.read(ini_cnf_save_locs)
            tors_pots, tors_zmas = {}, {}
            for tors_names, tors_grids in zip(run_tors_names, run_tors_grids):
                constraint_dct = structure.tors.build_constraint_dct(
                    zma, const_names, tors_names)
                pot, _, _, _, zmas = structure.tors.read_hr_pot(
                    tors_names, tors_grids,
                    ini_cnf_save_paths[0],
                    mod_ini_thy_info, ref_ene,
                    constraint_dct,
                    read_zma=True)
                tors_pots[tors_names] = pot
                tors_zmas[tors_names] = zmas

            # Print potential
            structure.tors.print_hr_pot(tors_pots)

            # Launch new conformer sampling from negative if requested
            if resamp_min:

                # Check for new minimum conformer
                new_min_zma = structure.tors.check_hr_pot(tors_pots, tors_zmas)
            
                if new_min_zma is not None:
                    print('Finding new low energy conformer...')
                    #     _ = conformer.conformer_sampling(
                    #         zma, spc_info,
                    #         mod_thy_info, thy_save_fs,
                    #         cnf_run_fs, cnf_save_fs,
                    #         opt_script_str, overwrite,
                    #         saddle=saddle, nsamp_par=mc_nsamp,
                    #         tors_names=run_tors_names,
                    #         two_stage=two_stage, retryfail=retryfail,
                    #         rxn_class=rxn_class, **opt_kwargs)

        elif job in ('energy', 'grad', 'hess', 'vpt2'):

            script_str, _, kwargs, _ = es_runner.qchem_params(
                *thy_info[0:2])
            for tors_names in run_tors_names:

                # Set the constraint dct and filesys for the scan
                constraint_dct = structure.tors.build_constraint_dct(
                    zma, const_names, tors_names)
                ini_scn_run_fs = filesys.build.scn_fs_from_cnf(
                    ini_zma_run_path, constraint_dct=constraint_dct)
                ini_scn_save_fs = filesys.build.scn_fs_from_cnf(
                    ini_zma_save_path, constraint_dct=constraint_dct)

                scn_locs = filesys.build.scn_locs_from_fs(
                    ini_scn_save_fs, tors_names, constraint_dct=constraint_dct)
                if scn_locs:
                    for locs in scn_locs:
                        geo_run_path = ini_scn_run_fs[-1].path(locs)
                        geo_save_path = ini_scn_save_fs[-1].path(locs)
                        geo = ini_scn_save_fs[-1].file.geometry.read(locs)
                        zma = ini_scn_save_fs[-1].file.zmatrix.read(locs)
                        ini_scn_run_fs[-1].create(locs)
                        ES_TSKS[job](
                            zma, geo, spc_info, mod_thy_info,
                            ini_scn_save_fs, geo_run_path, geo_save_path, locs,
                            script_str, overwrite,
                            retryfail=retryfail, **kwargs)
                        print('\n')
                else:
                    print('*WARNING: NO SCAN INFORMATION EXISTS.',
                          'Doing scan vs cscan?')
    else:
        print('No torsional modes in the species')


def irc_tsk(job, spc_dct, spc_name,
                thy_dct, es_keyword_dct,
                run_prefix, save_prefix, saddle):
    """ run a scan over the specified torsional coordinates
    """

    # Get dct for specific species task is run for
    spc = spc_dct[spc_name]

    # Set up coordinate name
    coo_name = 'RC'

    # Set the spc_info
    spc_info = filesys.inf.get_spc_info(spc)

    # Script
    _, opt_script_str, _, opt_kwargs = es_runner.qchem_params(
        *thy_info[0:2])

    # Modify the theory
    mod_thy_info = filesys.inf.modify_orb_restrict(spc_info, thy_info)
    mod_ini_thy_info = filesys.inf.modify_orb_restrict(spc_info, ini_thy_info)

    # Get options from the dct or es options lst
    overwrite = es_keyword_dct['overwrite']
    # retryfail = es_keyword_dct['retryfail']
    irc_idxs = spc['irc_idxs']

    # Set the filesystem objects
    rxn_info = filesys.inf.rxn_info(
        spc['reacs'], spc['prods'], spc_dct)

    _, thy_run_path = filesys.build.rxn_thy_fs_from_root(
        run_prefix, rxn_info, mod_thy_info)
    _, thy_save_path = filesys.build.rxn_thy_fs_from_root(
        save_prefix, rxn_info, mod_thy_info)
    _, thy_save_path = filesys.build.ts_fs_from_thy(thy_save_path)
    _, thy_run_path = filesys.build.ts_fs_from_thy(thy_run_path)

    _, ini_thy_run_path = filesys.build.rxn_thy_fs_from_root(
        run_prefix, rxn_info, mod_ini_thy_info)
    _, ini_thy_save_path = filesys.build.rxn_thy_fs_from_root(
        save_prefix, rxn_info, mod_ini_thy_info)
    _, ini_thy_save_path = filesys.build.ts_fs_from_thy(ini_thy_save_path)
    _, ini_thy_run_path = filesys.build.ts_fs_from_thy(ini_thy_run_path)

    ini_cnf_run_fs, _ = filesys.build.cnf_fs_from_prefix(
        ini_thy_run_path, mod_ini_thy_info, cnf=None)
    ini_cnf_save_fs, ini_cnf_save_locs = filesys.build.cnf_fs_from_prefix(
        ini_thy_save_path, mod_ini_thy_info, cnf='min')
    ini_cnf_save_paths = filesys.build.cnf_paths_from_locs(
        ini_cnf_save_fs, ini_cnf_save_locs)
    ini_cnf_run_paths = filesys.build.cnf_paths_from_locs(
        ini_cnf_run_fs, ini_cnf_save_locs)
    ini_cnf_run_fs[-1].create(ini_cnf_save_locs)

    ini_scn_run_fs = autofile.fs.scan(ini_cnf_run_paths[0])
    ini_scn_save_fs = autofile.fs.scan(ini_cnf_save_paths[0])

    if job == 'scan':

        # Script
        _, opt_script_str, _, opt_kwargs = es_runner.qchem_params(
            *thy_info[0:2])

        zma, geo = filesys.inf.get_zma_geo(ini_cnf_save_fs, ini_cnf_save_locs)
        irc.scan(
            geo, spc_info, mod_thy_info, coo_name, irc_idxs,
            ini_scn_save_fs, ini_scn_run_fs, ini_cnf_run_paths[0],
            overwrite, opt_script_str, **opt_kwargs)

    elif job in ('energy', 'grad', 'hess'):

        # Script
        script_str, _, kwargs, _ = es_runner.qchem_params(
            *mod_thy_info[0:2])

        # Need to put in something with the IRC idxs
        for idx in irc_idxs:
            locs = [[coo_name], [idx]]
            geo_run_path = ini_scn_run_fs[-1].path(locs)
            geo_save_path = ini_scn_save_fs[-1].path(locs)
            zma, geo = filesys.inf.get_zma_geo(ini_scn_save_fs, locs)
            ini_scn_run_fs[-1].create(locs)
            ES_TSKS[job](
                zma, geo, spc_info, mod_thy_info,
                ini_scn_save_fs, geo_run_path, geo_save_path, locs,
                script_str, overwrite, **kwargs)
            print('\n')
