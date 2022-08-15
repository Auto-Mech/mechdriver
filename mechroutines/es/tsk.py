""" eletronic structure routines modules
"""

import importlib
import autofile
import automol
import elstruct
from phydat import phycon
from mechanalyzer.inf import thy as tinfo
from mechanalyzer.inf import rxn as rinfo
from mechanalyzer.inf import spc as sinfo
from mechlib import filesys
from mechlib.filesys import build_fs
from mechlib.filesys import root_locs
from mechlib.amech_io import printer as ioprinter
from mechroutines.es._routines import conformer
from mechroutines.es._routines import hr
from mechroutines.es._routines import tau
from mechroutines.es.ts import findts
from mechroutines.es.ts import ts_zma_locs

from mechroutines.es.runner import scan
from mechroutines.es.runner import qchem_params


# Dictionary of Electronic Structure Calculators
SP_MODULE = importlib.import_module('mechroutines.es._routines.sp')
ES_TSKS = {
    'energy': SP_MODULE.run_energy,
    'grad': SP_MODULE.run_gradient,
    'hess': SP_MODULE.run_hessian,
    'vpt2': SP_MODULE.run_vpt2,
    'prop': SP_MODULE.run_prop
}


def run_tsk(tsk, spc_dct, spc_name,
            thy_dct, es_keyword_dct,
            run_prefix, save_prefix,
            print_debug=False):
    """ Execute the specified electronic structure task.

        :param tsk: name of task
        :type tsk: str
        :param spc_dct:
        :type spc_dct:
        :param spc_name: name of species
        :type spc_name: str
        :param thy_dct:
        :type thy_dct:
        :param es_keyword_dct: keyword-value pairs for task
        :type es_keyword_dct: dict[str:str]
        :param run_prefix: root-path to the run-filesystem
        :type run_prefix: str
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
        :param print_debug: option to print extra debug information
        :type print_debug: bool
    """

    ioprinter.task_header(tsk, spc_name)
    ioprinter.keyword_list(es_keyword_dct, thy_dct)

    skip = skip_task(tsk, spc_dct, spc_name,
                     thy_dct, es_keyword_dct, save_prefix)
    if not skip:
        # Get stuff from task
        job = tsk.split('_', 1)[1]

        # Run the task if an initial geom exists
        if 'init' in tsk:
            _ = geom_init(
                spc_dct, spc_name, thy_dct, es_keyword_dct,
                run_prefix, save_prefix)
        elif 'conf' in tsk:
            conformer_tsk(
                job, spc_dct, spc_name, thy_dct, es_keyword_dct,
                run_prefix, save_prefix, print_debug=print_debug)
        elif 'tau' in tsk:
            tau_tsk(
                job, spc_dct, spc_name, thy_dct, es_keyword_dct,
                run_prefix, save_prefix)
        elif 'hr' in tsk:
            hr_tsk(
                job, spc_dct, spc_name, thy_dct, es_keyword_dct,
                run_prefix, save_prefix)
        elif 'find' in tsk or 'rpath' in tsk:
            findts(
                tsk, spc_dct, spc_name, thy_dct, es_keyword_dct,
                run_prefix, save_prefix)

    ioprinter.task_footer()


# FUNCTIONS FOR SAMPLING AND SCANS #
def geom_init(spc_dct, spc_name, thy_dct, es_keyword_dct,
              run_prefix, save_prefix):
    """ Execute the task for a species used to seed the
        filesystem with a reliable initial conformer.

        :param spc_dct:
        :type spc_dct:
        :param spc_name: name of species
        :type spc_name: str
        :param thy_dct:
        :type thy_dct:
        :param es_keyword_dct: keyword-val pairs for electronic structure task
        :type es_keyword_dct: dict[str:str]
        :param run_prefix: root-path to the run-filesystem
        :type run_prefix: str
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
    """

    spc_dct_i = spc_dct[spc_name]
    spc_info = sinfo.from_dct(spc_dct_i, canonical=True)
    
    # Get the theory info
    method_dct = thy_dct.get(es_keyword_dct['runlvl'])
    ini_method_dct = thy_dct.get(es_keyword_dct['inplvl'])
    thy_info = tinfo.from_dct(method_dct)
    ini_thy_info = tinfo.from_dct(ini_method_dct)
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
    mod_ini_thy_info = tinfo.modify_orb_label(
        ini_thy_info, spc_info)

    # Set the filesystem objects
    _, ini_cnf_save_fs = build_fs(
        run_prefix, save_prefix, 'CONFORMER',
        spc_locs=spc_info, thy_locs=mod_ini_thy_info[1:])
    cnf_run_fs, cnf_save_fs = build_fs(
        run_prefix, save_prefix, 'CONFORMER',
        spc_locs=spc_info, thy_locs=mod_thy_info[1:])

    # Get a reference geometry if one not found
    success = conformer.initial_conformer(
        spc_dct_i, spc_info, ini_method_dct, method_dct,
        ini_cnf_save_fs, cnf_run_fs, cnf_save_fs,
        es_keyword_dct)

    return success


def conformer_tsk(job, spc_dct, spc_name,
                  thy_dct, es_keyword_dct,
                  run_prefix, save_prefix, print_debug=False):
    """ Prepares and executes all electronic structure tasks that
        generate information for species and transition state conformers.
        This includes sampling and optimization procedures to generate
        conformer structures, as well as __ calculations using some
        saved conformer as input.

        :param job(subtask): calculatiion(s) to perform for conformer
        :type job: str
        :param spc_dct:
        :type spc_dct:
        :param spc_name: name of species
        :type spc_name: str
        :param thy_dct:
        :type thy_dct:
        :param es_keyword_dct: keyword-values for electronic structure task
        :type es_keyword_dct: dict[str:str]
        :param run_prefix: root-path to the run-filesystem
        :type run_prefix: str
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
    """

    saddle = bool('ts_' in spc_name)

    spc_dct_i = spc_dct[spc_name]

    # Set the spc_info
    if not saddle:
        spc_info = sinfo.from_dct(spc_dct_i, canonical=True)
    else:
        spc_info = rinfo.ts_info(spc_dct_i['canon_rxn_info'])
    zrxn = spc_dct_i.get('zrxn', None)

    overwrite = es_keyword_dct['overwrite']
    retryfail = es_keyword_dct['retryfail']
    nprocs = 1
    # Modify the theory
    method_dct = thy_dct.get(es_keyword_dct['runlvl'])
    ini_method_dct = thy_dct.get(es_keyword_dct['inplvl'])
    thy_info = tinfo.from_dct(method_dct)
    ini_thy_info = tinfo.from_dct(ini_method_dct)
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
    mod_ini_thy_info = tinfo.modify_orb_label(
        ini_thy_info, spc_info)

    # New filesystem objects
    _root = root_locs(spc_dct_i, saddle=saddle, name=spc_name)
    ini_cnf_run_fs, ini_cnf_save_fs = build_fs(
        run_prefix, save_prefix, 'CONFORMER',
        thy_locs=mod_ini_thy_info[1:],
        **_root)
    cnf_run_fs, cnf_save_fs = build_fs(
        run_prefix, save_prefix, 'CONFORMER',
        thy_locs=mod_thy_info[1:],
        **_root)

    if job == 'samp':

        # Build the ini zma filesys
        user_conf_ids = spc_dct_i.get('conf_id')
        if user_conf_ids is None:
            ini_loc_info = filesys.mincnf.min_energy_conformer_locators(
                ini_cnf_save_fs, mod_ini_thy_info)
            ini_locs, ini_min_cnf_path = ini_loc_info
        else:
            print(f'Using user specified conformer IDs: {user_conf_ids}')
            ini_locs = user_conf_ids

        if any(ini_locs):
            ini_zma_save_fs = autofile.fs.zmatrix(ini_min_cnf_path)

            # Set up the run scripts
            script_str, kwargs = qchem_params(
                method_dct, elstruct.Job.OPTIMIZATION)

            # Set variables if it is a saddle
            two_stage = saddle
            mc_nsamp = spc_dct_i['mc_nsamp']
            resave = es_keyword_dct['resave']

            # Read the geometry and zma from the ini file system
            geo = ini_cnf_save_fs[-1].file.geometry.read(ini_locs)
            zma_locs = (0,)
            if saddle:
                zma_locs = ts_zma_locs(spc_dct, spc_name, ini_zma_save_fs)
            zma = ini_zma_save_fs[-1].file.zmatrix.read(zma_locs)

            # Read the torsions from the ini file sys
            if ini_zma_save_fs[-1].file.torsions.exists(zma_locs):
                tors_dct = ini_zma_save_fs[-1].file.torsions.read(zma_locs)
                rotors = automol.rotor.from_data(zma, tors_dct)
                tors_names = automol.rotor.names(rotors, flat=True)
            else:
                tors_names = ()

            geo_path = ini_cnf_save_fs[-1].path(ini_locs)
            ioprinter.initial_geom_path('Sampling started', geo_path)

            # Check runsystem for equal ring CONF make conf_fs
            # Else make new ring conf directory
            rid = conformer.rng_loc_for_geo(geo, cnf_save_fs)

            if rid is None:
                conformer.single_conformer(
                    zma, spc_info, mod_thy_info,
                    cnf_run_fs, cnf_save_fs,
                    script_str, overwrite,
                    retryfail=retryfail, zrxn=zrxn,
                    **kwargs)

                rid = conformer.rng_loc_for_geo(
                    geo, cnf_save_fs)

            # Run the sampling
            conformer.conformer_sampling(
                zma, spc_info, mod_thy_info,
                cnf_run_fs, cnf_save_fs, rid,
                script_str, overwrite,
                nsamp_par=mc_nsamp,
                tors_names=tors_names,
                zrxn=zrxn, two_stage=two_stage,
                retryfail=retryfail, resave=resave,
                repulsion_thresh=40.0, print_debug=print_debug,
                **kwargs)
        else:
            ioprinter.info_message(
                'Missing conformers. Skipping task...')

    elif job == 'pucker':

        # Build the ini zma filesys
        ini_loc_info = filesys.mincnf.min_energy_conformer_locators(
            ini_cnf_save_fs, mod_ini_thy_info)
        ini_min_locs, ini_min_cnf_path = ini_loc_info
        ini_zma_save_fs = autofile.fs.zmatrix(ini_min_cnf_path)

        # Set up the run scripts
        script_str, kwargs = qchem_params(
            method_dct, elstruct.Job.OPTIMIZATION)

        # Set variables if it is a saddle
        two_stage = saddle
        mc_nsamp = spc_dct_i['mc_nsamp']

        # Read the geometry and zma from the ini file system
        geo = ini_cnf_save_fs[-1].file.geometry.read(ini_min_locs)
        zma_locs = (0,)
        if saddle:
            zma_locs = ts_zma_locs(spc_dct, spc_name, ini_zma_save_fs)
        zma = ini_zma_save_fs[-1].file.zmatrix.read(zma_locs)

        # Read the torsions from the ini file sys
        if ini_zma_save_fs[-1].file.ring_torsions.exists(zma_locs):
            ring_tors_dct = ini_zma_save_fs[-1].file.ring_torsions.read(zma_locs)
        else:
            ring_tors_dct = {}

        geo_path = ini_cnf_save_fs[-1].path(ini_min_locs)
        ioprinter.initial_geom_path('Sampling started', geo_path)

        # Run the sampling
        conformer.ring_conformer_sampling(
            zma, spc_info, mod_thy_info,
            cnf_run_fs, cnf_save_fs,
            script_str, overwrite,
            nsamp_par=mc_nsamp,
            ring_tors_dct=ring_tors_dct, zrxn=zrxn,
            two_stage=two_stage, retryfail=retryfail,
            **kwargs)

    elif job == 'opt':

        cnf_range = es_keyword_dct['cnf_range']
        hbond_cutoffs = spc_dct_i['hbond_cutoffs']
        cnf_sort_info_lst = _sort_info_lst(
            es_keyword_dct['sort'], thy_dct, spc_info)
        resave = es_keyword_dct['resave']

        # Set up the run scripts
        script_str, kwargs = qchem_params(
            method_dct, elstruct.Job.OPTIMIZATION)

        rng_cnf_locs_lst, _ = filesys.mincnf.conformer_locators(
            cnf_save_fs, mod_thy_info,
            cnf_range='all', nprocs=nprocs)

        ini_rng_cnf_locs_lst, _ = filesys.mincnf.conformer_locators(
            ini_cnf_save_fs, mod_ini_thy_info,
            cnf_range=cnf_range, sort_info_lst=cnf_sort_info_lst,
            hbond_cutoffs=hbond_cutoffs,
            print_enes=True, nprocs=nprocs)

        # Truncate the list of the ini confs
        uni_rng_locs_lst, uni_cnf_locs_lst = conformer.unique_fs_ring_confs(
            cnf_save_fs, rng_cnf_locs_lst,
            ini_cnf_save_fs, ini_rng_cnf_locs_lst)
        # ioprinter.debug_message(
        #    'uni lst that has no similar ring', uni_rng_locs_lst)
        # ioprinter.debug_message(
        #    'uni lst that has similar ring', uni_cnf_locs_lst)

        for locs in uni_rng_locs_lst:
            # Obtain the zma from ini loc
            ini_cnf_save_path = ini_cnf_save_fs[-1].path(locs)
            ini_zma_save_fs = autofile.fs.zmatrix(ini_cnf_save_path)
            zma_locs = (0,)
            if saddle:
                zma_locs = ts_zma_locs(spc_dct, spc_name, ini_zma_save_fs)
            zma = ini_zma_save_fs[-1].file.zmatrix.read(zma_locs)

            # Make the ring filesystem
            conformer.single_conformer(
                zma, spc_info, mod_thy_info,
                cnf_run_fs, cnf_save_fs,
                script_str, overwrite,
                retryfail=retryfail, zrxn=zrxn,
                use_locs=locs, resave=resave,
                **kwargs)

        for locs in uni_cnf_locs_lst:
            ini_locs, rid = locs
            # Obtain the zma from ini loc
            ini_cnf_save_path = ini_cnf_save_fs[-1].path(ini_locs)
            ini_zma_save_fs = autofile.fs.zmatrix(ini_cnf_save_path)
            zma_locs = (0,)
            if saddle:
                zma_locs = ts_zma_locs(spc_dct, spc_name, ini_zma_save_fs)
            zma = ini_zma_save_fs[-1].file.zmatrix.read(zma_locs)
            # obtain conformer filesys associated with ring at the runlevel
            cid = autofile.schema.generate_new_conformer_id()
            conformer.single_conformer(
                zma, spc_info, mod_thy_info,
                cnf_run_fs, cnf_save_fs,
                script_str, overwrite,
                retryfail=retryfail, zrxn=zrxn,
                use_locs=(rid, cid), resave=resave,
                **kwargs)

        # print all geometres within cnfrange
        rng_cnf_locs_lst, _ = filesys.mincnf.conformer_locators(
            cnf_save_fs, mod_thy_info,
            cnf_range=cnf_range, sort_info_lst=cnf_sort_info_lst,
            hbond_cutoffs=hbond_cutoffs, nprocs=nprocs)
        for locs in rng_cnf_locs_lst:
            geo = cnf_save_fs[-1].file.geometry.read(locs)
            ioprinter.geometry(geo)

    elif job in ('energy', 'grad', 'hess', 'vpt2', 'prop'):

        cnf_range = es_keyword_dct['cnf_range']
        hbond_cutoffs = spc_dct_i['hbond_cutoffs']
        cnf_sort_info_lst = _sort_info_lst(
            es_keyword_dct['sort'], thy_dct, spc_info)

        user_conf_ids = spc_dct_i.get('conf_id')
        if user_conf_ids is None:
            ini_rng_cnf_locs_lst, _ = filesys.mincnf.conformer_locators(
                ini_cnf_save_fs, mod_ini_thy_info,
                cnf_range=cnf_range, sort_info_lst=cnf_sort_info_lst,
                hbond_cutoffs=hbond_cutoffs,
                print_enes=True, nprocs=nprocs)
        else:
            print(f'Using user specified conformer IDs: {user_conf_ids}')
            ini_rng_cnf_locs_lst = (user_conf_ids,)

        # Check if locs exist, kill if it doesn't
        if not ini_rng_cnf_locs_lst:
            ioprinter.error_message(
                'No min-energy conformer found for level:')

        else:

            # Grab frequencies for the reference, print ref freqs
            if job == 'hess':
                if ini_cnf_save_fs[-1].file.harmonic_frequencies.exists(
                    ini_rng_cnf_locs_lst[0]):
                    frq = ini_cnf_save_fs[-1].file.harmonic_frequencies.read(
                        ini_rng_cnf_locs_lst[0])
                    ref_val = frq
                else:
                    ref_val = None
                if ref_val is not None and zrxn is not None:
                    ref_path = cnf_save_fs[-1].path(ini_rng_cnf_locs_lst[0])
                    print('Found reference frequencies for saddle-point '
                          f'checks for conformer at\n {ref_path}')
                    ioprinter.frequencies(ref_val)
            else:
                ref_val = None

            # Run the job over all the conformers requested by the user
            print('Going over all requested conformers for task...\n')
            for ini_locs in ini_rng_cnf_locs_lst:
                ini_cnf_run_fs[-1].create(ini_locs)
                geo_save_path = ini_cnf_save_fs[-1].path(ini_locs)
                ini_zma_save_fs = autofile.fs.zmatrix(geo_save_path)
                print('Running task for geometry at ', geo_save_path)
                geo = ini_cnf_save_fs[-1].file.geometry.read(ini_locs)
                zma_locs = (0,)
                if saddle:
                    zma_locs = ts_zma_locs(spc_dct, spc_name, ini_zma_save_fs)
                zma = ini_zma_save_fs[-1].file.zmatrix.read(zma_locs)

                script_str, kwargs = qchem_params(
                    method_dct, geo=geo, spc_info=spc_info)

                ES_TSKS[job](
                    zma, geo, spc_info, mod_thy_info,
                    ini_cnf_run_fs, ini_cnf_save_fs, ini_locs, run_prefix,
                    script_str, overwrite, zrxn=zrxn,
                    retryfail=retryfail, method_dct=method_dct,
                    ref_val=ref_val,
                    **kwargs)
                print('\n === FINISHED CONF ===\n')


def tau_tsk(job, spc_dct, spc_name,
            thy_dct, es_keyword_dct,
            run_prefix, save_prefix):
    """ Energies, gradients, and hessians,
        for set of arbitrarily sampled torsional coordinates
        with all other coordinates optimized

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
    spc_dct_i = spc_dct[spc_name]

    # Set the spc_info
    spc_info = sinfo.from_dct(spc_dct_i, canonical=True)

    # Get es options
    overwrite = es_keyword_dct['overwrite']
    retryfail = es_keyword_dct['retryfail']
    # scan_increment = spc_dct_i['hind_inc']
    nsamp_par = spc_dct_i['tau_nsamp']
    zrxn = spc_dct_i.get('zrxn', None)

    # Modify the theory
    method_dct = thy_dct.get(es_keyword_dct['runlvl'])
    ini_method_dct = thy_dct.get(es_keyword_dct['inplvl'])
    thy_info = tinfo.from_dct(method_dct)
    ini_thy_info = tinfo.from_dct(ini_method_dct)
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
    mod_ini_thy_info = tinfo.modify_orb_label(ini_thy_info, spc_info)

    # Script
    script_str, kwargs = qchem_params(
        method_dct, elstruct.Job.OPTIMIZATION)

    # Set the filesystem objects for thy info
    _, ini_cnf_save_fs = build_fs(
        run_prefix, save_prefix, 'CONFORMER',
        spc_locs=spc_info, thy_locs=mod_ini_thy_info[1:])
    ini_loc_info = filesys.mincnf.min_energy_conformer_locators(
        ini_cnf_save_fs, mod_ini_thy_info)
    ini_min_cnf_locs, ini_min_cnf_path = ini_loc_info

    ini_zma_save_fs = autofile.fs.zmatrix(ini_min_cnf_path)
    geo = ini_cnf_save_fs[-1].file.geometry.read(ini_min_cnf_locs)
    zma_locs = (0,)
    zma = ini_zma_save_fs[-1].file.zmatrix.read(zma_locs)
    ini_sp_save_fs = autofile.fs.single_point(ini_min_cnf_path)
    if ini_sp_save_fs[-1].file.energy.exists(mod_ini_thy_info[1:4]):
        ref_ene = ini_sp_save_fs[-1].file.energy.read(mod_ini_thy_info[1:4])
    else:
        ref_ene = ini_cnf_save_fs[-1].file.energy.read(ini_min_cnf_locs)

    # Get the tors names
    if ini_zma_save_fs[-1].file.torsions.exists(zma_locs):
        tors_dct = ini_zma_save_fs[-1].file.torsions.read(zma_locs)
        torsions = automol.rotor.from_data(zma, tors_dct)
    else:
        torsions = ()

    saddle = bool('ts_' in spc_name)

    # Set the database style
    db_style = 'jsondb'
    # db_style = 'directory'

    # Run the task if any torsions exist
    if torsions and not saddle:

        # Set up tau filesystem objects
        tau_run_fs, tau_save_fs = build_fs(
            run_prefix, save_prefix, 'TAU',
            spc_locs=spc_info, thy_locs=mod_thy_info[1:])

        if db_style == 'jsondb':
            tau_save_fs[-1].root.create()
            tau_save_fs[-1].json_create()
            for locs in tau_save_fs[-1].existing():
                if tau_save_fs[-1].file.geometry.exists(locs):
                    geol = tau_save_fs[-1].file.geometry.read(locs)
                    tau_save_fs[-1].json.geometry.write(geol, locs)
                if tau_save_fs[-1].file.energy.exists(locs):
                    enel = tau_save_fs[-1].file.energy.read(locs)
                    tau_save_fs[-1].json.energy.write(enel, locs)
                if tau_save_fs[-1].file.geometry_info.exists(locs):
                    geo_infl = tau_save_fs[-1].file.geometry_info.read(locs)
                    tau_save_fs[-1].json.geometry_info.write(geo_infl, locs)
                if tau_save_fs[-1].file.geometry_input.exists(locs):
                    inp_strl = tau_save_fs[-1].file.geometry_input.read(locs)
                    tau_save_fs[-1].json.geometry_input.write(inp_strl, locs)
                if tau_save_fs[-1].file.gradient_input.exists(locs):
                    inp_strl = tau_save_fs[-1].file.gradient_input.read(locs)
                    tau_save_fs[-1].json.gradient_input.write(inp_strl, locs)
                if tau_save_fs[-1].file.hessian_input.exists(locs):
                    inp_strl = tau_save_fs[-1].file.hessian_input.read(locs)
                    tau_save_fs[-1].json.hessian_input.write(inp_strl, locs)
                if tau_save_fs[-1].file.gradient_info.exists(locs):
                    inf_objl = tau_save_fs[-1].file.gradient_info.read(locs)
                    tau_save_fs[-1].json.gradient_info.write(inf_objl, locs)
                if tau_save_fs[-1].file.hessian_info.exists(locs):
                    inf_objl = tau_save_fs[-1].file.hessian_info.read(locs)
                    tau_save_fs[-1].json.hessian_info.write(inf_objl, locs)
                if tau_save_fs[-1].file.gradient.exists(locs):
                    gradl = tau_save_fs[-1].file.gradient.read(locs)
                    tau_save_fs[-1].json.gradient.write(gradl, locs)
                if tau_save_fs[-1].file.hessian.exists(locs):
                    hessl = tau_save_fs[-1].file.hessian.read(locs)
                    tau_save_fs[-1].json.energy.hessian(hessl, locs)
                if tau_save_fs[-1].file.zmatrix.exists(locs):
                    zmatl = tau_save_fs[-1].file.zmatrix.read(locs)
                    tau_save_fs[-1].json.zmatrix.write(zmatl, locs)
                if tau_save_fs[-1].file.harmonic_frequencies.exists(locs):
                    hfreql = tau_save_fs[-1].file.harmonic_frequencies.read(
                        locs)
                    tau_save_fs[-1].json.harmonic_frequencies.write(
                        hfreql, locs)
                save_path = tau_save_fs[-1].path(locs)
                sp_save_fs = autofile.fs.single_point(save_path)
                sp_save_locs = sp_save_fs[-1].existing()
                save_path = tau_save_fs[-1].root.path()
                jsp_save_fs = autofile.fs.single_point(
                    save_path, json_layer=locs)
                for sp_locs in sp_save_locs:
                    if sp_save_fs[-1].file.energy.exists(sp_locs):
                        enel = sp_save_fs[-1].file.energy.read(sp_locs)
                        jsp_save_fs[-1].json.energy.write(enel, sp_locs)
                    if sp_save_fs[-1].file.input.exists(sp_locs):
                        inp_strl = sp_save_fs[-1].file.input.read(sp_locs)
                        jsp_save_fs[-1].json.input.write(inp_strl, sp_locs)
                    if sp_save_fs[-1].file.info.exists(sp_locs):
                        inf_objl = sp_save_fs[-1].file.info.read(sp_locs)
                        jsp_save_fs[-1].json.info.write(inf_objl, sp_locs)

        if job == 'samp':

            # Set up the script
            script_str, kwargs = qchem_params(
                method_dct, elstruct.Job.OPTIMIZATION)

            tors_names = automol.rotor.names(torsions, flat=True)
            resave = es_keyword_dct['resave']

            tau.tau_sampling(
                zma, ref_ene, spc_info,
                mod_ini_thy_info,
                tau_run_fs, tau_save_fs,
                script_str, overwrite,
                db_style=db_style,
                nsamp_par=nsamp_par,
                tors_names=tors_names,
                repulsion_thresh=40.0,
                zrxn=zrxn, resave=resave,
                **kwargs)

        elif job in ('energy', 'grad'):

            # Set up the run scripts
            script_str, kwargs = qchem_params(
                method_dct)

            # Run the job over all the conformers requested by the user
            if db_style == 'directory':
                tau_locs = tau_save_fs[-1].existing()
            elif db_style == 'jsondb':
                tau_locs = tau_save_fs[-1].json_existing()
            for locs in tau_locs:
                if db_style == 'jsondb':
                    geo_save_path = tau_save_fs[-1].root.path()
                    geo = tau_save_fs[-1].json.geometry.read(locs)
                elif db_style == 'directory':
                    geo_save_path = tau_save_fs[-1].path(locs)
                    geo = tau_save_fs[-1].file.geometry.read(locs)
                tau_run_fs[-1].create(locs)
                zma = None
                ES_TSKS[job](
                    zma, geo, spc_info, mod_thy_info,
                    tau_save_fs, locs, run_prefix,
                    script_str, overwrite,
                    retryfail=retryfail,
                    **kwargs)
                ioprinter.obj('vspace')

        elif job == 'hess':

            # Add the hessian max
            hessmax = es_keyword_dct['hessmax']

            # Set up the run scripts
            script_str, kwargs = qchem_params(
                method_dct)

            # Run the job over all the conformers requested by the user
            hess_cnt = 0
            if db_style == 'directory':
                tau_locs = tau_save_fs[-1].existing()
            elif db_style == 'jsondb':
                tau_locs = tau_save_fs[-1].json_existing()
            for locs in tau_locs:
                ioprinter.info_message(
                    f'HESS Number {hess_cnt+1}', newline=1)
                if db_style == 'directory':
                    geo_save_path = tau_save_fs[-1].path(locs)
                    if tau_save_fs[-1].file.hessian.exists(locs):
                        ioprinter.existing_path('Hessian', geo_save_path)
                        hess_cnt += 1
                        continue
                    geo = tau_save_fs[-1].file.geometry.read(locs)
                elif db_style == 'jsondb':
                    geo_save_path = tau_save_fs[-1].root.path()
                    if tau_save_fs[-1].json.hessian.exists(locs):
                        ioprinter.existing_path('Hessian', geo_save_path)
                        hess_cnt += 1
                        continue
                    geo = tau_save_fs[-1].json.geometry.read(locs)
                tau_run_fs[-1].create(locs)
                ES_TSKS[job](
                    None, geo, spc_info, mod_thy_info,
                    tau_run_fs, tau_save_fs, locs, run_prefix,
                    script_str, overwrite,
                    retryfail=retryfail, method_dct=method_dct,
                    correct_vals=False,
                    **kwargs)
                hess_cnt += 1
                if hess_cnt == hessmax:
                    break

    else:
        ioprinter.info_message('No torsional modes in the species')


def hr_tsk(job, spc_dct, spc_name,
           thy_dct, es_keyword_dct,
           run_prefix, save_prefix):
    """ Prepares and executes all electronic structure tasks that
        generate information for points along hindered-rotor coordinate
        scans which are launched from some conformer in the save filesystem.

        For species and transition state conformers.

        This includes scanning procedures to generate geometries
        (relaxed) or energies (rigid) points along
        conformer structures, as well as __ calculations using some
        saved conformer as input.

        :param job:
        :type job:
        :param spc_dct:
        :type spc_dct:
        :param spc_name:
        :type spc_name:
        :param thy_dct:
        :type thy_dct:
        :param es_keyword_dct: keyword-val pairs for electronic structure task
        :type es_keyword_dct: dict[str:str]
        :param run_prefix: root-path to the run-filesystem
        :type run_prefix: str
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
    """

    # 1) Set up input info
    # 1a) Set the spc_info
    spc_dct_i = spc_dct[spc_name]
    saddle = bool('ts_' in spc_name)
    if not saddle:
        spc_info = sinfo.from_dct(spc_dct_i, canonical=True)
    else:
        spc_info = rinfo.ts_info(spc_dct_i['canon_rxn_info'])

    # 1b) Get options from the dct or es options lst
    overwrite = es_keyword_dct['overwrite']
    retryfail = es_keyword_dct['retryfail']
    tors_model = es_keyword_dct['tors_model']
    # 1c) how parallelized do we wanna be
    # nprocs = es_keyword_dct['nprocs']
    nprocs = 1

    # 1d) Modify the theory  info
    method_dct = thy_dct.get(es_keyword_dct['runlvl'])
    ini_method_dct = thy_dct.get(es_keyword_dct['inplvl'])
    thy_info = tinfo.from_dct(method_dct)
    ini_thy_info = tinfo.from_dct(ini_method_dct)
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
    mod_ini_thy_info = tinfo.modify_orb_label(
        ini_thy_info, spc_info)

    # 2) Set the filesystem objects
    _root = root_locs(spc_dct_i, saddle=saddle, name=spc_name)
    ini_cnf_run_fs, ini_cnf_save_fs = build_fs(
        run_prefix, save_prefix, 'CONFORMER',
        thy_locs=mod_ini_thy_info[1:],
        **_root)
    # 3) Determine job locations
    # 3a) set the sorting parameters
    cnf_range = es_keyword_dct['cnf_range']
    hbond_cutoffs = spc_dct_i['hbond_cutoffs']
    user_conf_ids = spc_dct_i.get('conf_id')
    # 3b) find input geometry locations
    if user_conf_ids is None:
        cnf_sort_info_lst = _sort_info_lst(
            es_keyword_dct['sort'], thy_dct, spc_info)
        ini_min_locs_lst, ini_path_lst = filesys.mincnf.conformer_locators(
            ini_cnf_save_fs, mod_ini_thy_info,
            cnf_range=cnf_range, sort_info_lst=cnf_sort_info_lst,
            hbond_cutoffs=hbond_cutoffs,
            print_enes=True, nprocs=nprocs)
    else:
        ini_min_locs_lst = (user_conf_ids,)
        ini_path_lst = (ini_cnf_save_fs[-1].path(user_conf_ids),)
    # 3c) find run geometry locations, for jobs that require it
    if job == 'scan':
        cnf_run_fs, cnf_save_fs = build_fs(
            run_prefix, save_prefix, 'CONFORMER',
            thy_locs=mod_thy_info[1:],
            **_root)
        all_run_cnf_locs_lst, _ = filesys.mincnf.conformer_locators(
            cnf_save_fs, mod_thy_info,
            cnf_range='all', nprocs=nprocs)
        ini_to_run_locs_dct = filesys.mincnf.fs_confs_dict(
            cnf_save_fs, all_run_cnf_locs_lst,
            ini_cnf_save_fs, ini_min_locs_lst)

    # 4) Loop over input locations that jobs will be run on
    for ini_min_locs, ini_cnf_save_path in zip(ini_min_locs_lst, ini_path_lst):

        # 5) Read properties for input geometry
        ini_zma_save_fs = autofile.fs.zmatrix(ini_cnf_save_path)
        geo = ini_cnf_save_fs[-1].file.geometry.read(ini_min_locs)
        zma_locs = (0,)
        if saddle:
            zma_locs = ts_zma_locs(spc_dct, spc_name, ini_zma_save_fs)
        zma = ini_zma_save_fs[-1].file.zmatrix.read(zma_locs)
        if ini_zma_save_fs[-1].file.torsions.exists(zma_locs):
            tors_dct = ini_zma_save_fs[-1].file.torsions.read(zma_locs)
            rotors = automol.rotor.from_data(
                zma, tors_dct, multi='md' in tors_model)
        else:
            rotors = ()
        zrxn = spc_dct_i.get('zrxn', None)
        # If there aren't any rotors, stop hr tasks
        if not any(rotors):
            ioprinter.info_message('No torsional modes in the species')
            break

        if 'fa' in tors_model:
            scn = 'CSCAN'
        elif 'f' in tors_model:
            if len(rotors) > 1:
                scn = 'CSCAN'
            else:
                scn = 'SCAN'
        else:
            scn = 'SCAN'

        if job == 'scan':

            # scan) Find equivalent conformer in the run filesys
            # scan) if it doesn't exist
            # scan) run a single conformer to generate it
            min_locs = ini_to_run_locs_dct[tuple(ini_min_locs)]
            if min_locs is None:
                script_str, kwargs = qchem_params(
                    method_dct, elstruct.Job.OPTIMIZATION)
                rid = conformer.rng_loc_for_geo(geo, cnf_save_fs)
                if rid is None:
                    new_rid = autofile.schema.generate_new_ring_id()
                    new_cid = autofile.schema.generate_new_conformer_id()
                    conformer.single_conformer(
                        zma, spc_info, mod_thy_info,
                        cnf_run_fs, cnf_save_fs,
                        script_str, overwrite,
                        retryfail=retryfail, zrxn=zrxn,
                        use_locs=(new_rid, new_cid),
                        **kwargs)
                    min_locs = (new_rid, new_cid)
                else:
                    new_cid = autofile.schema.generate_new_conformer_id()
                    conformer.single_conformer(
                        zma, spc_info, mod_thy_info,
                        cnf_run_fs, cnf_save_fs,
                        script_str, overwrite,
                        retryfail=retryfail, zrxn=zrxn,
                        use_locs=(rid, new_cid),
                        **kwargs)
                    min_locs = (rid, new_cid)
                    save_locs = cnf_save_fs[-1].existing()
                    if min_locs not in save_locs:
                        locinf = filesys.mincnf.this_conformer_was_run_in_run(
                            zma, cnf_run_fs)
                        _, sym_locs_lst = locinf
                        for sym_locs in sym_locs_lst:
                            if sym_locs in save_locs:
                                min_locs = sym_locs
            cnf_save_path = cnf_save_fs[-1].path(min_locs)
            ioprinter.info_message(
                f'Same conformer saved at {ini_cnf_save_path} '
                f'and {cnf_save_path}')

            # scan) re-create run fs if its been deleted
            cnf_run_fs[-1].create(min_locs)
            cnf_run_path = cnf_run_fs[-1].path(min_locs)

            # scan) Get the runlvl zma and torsion info
            zma_save_fs = autofile.fs.zmatrix(cnf_save_path)
            geo = cnf_save_fs[-1].file.geometry.read(min_locs)
            zma_locs = (0,)
            if saddle:
                zma_locs = ts_zma_locs(spc_dct, spc_name, zma_save_fs)
            zma = zma_save_fs[-1].file.zmatrix.read(zma_locs)
            if zma_save_fs[-1].file.torsions.exists(zma_locs):
                tors_dct = zma_save_fs[-1].file.torsions.read(zma_locs)
                rotors = automol.rotor.from_data(
                    zma, tors_dct, multi='md' in tors_model)
            else:
                rotors = ()
            if 'fa' in tors_model:
                scn = 'CSCAN'
            elif 'f' in tors_model:
                if len(rotors) > 1:
                    scn = 'CSCAN'
                else:
                    scn = 'SCAN'
            else:
                scn = 'SCAN'
            scn_run_fs, scn_save_fs = build_fs(
                cnf_run_path, cnf_save_path, scn,
                zma_locs=zma_locs)

            increment = spc_dct_i.get('hind_inc', 30.0*phycon.DEG2RAD)
            hr.hindered_rotor_scans(
                zma, spc_info, mod_thy_info,
                scn_run_fs, scn_save_fs,
                rotors, tors_model, method_dct,
                overwrite,
                zrxn=zrxn,
                saddle=saddle,
                increment=increment,
                retryfail=retryfail)

        elif job == 'reopt':

            script_str, kwargs = qchem_params(
                method_dct, elstruct.Job.OPTIMIZATION)

            # pull stuff from dcts
            ethresh = es_keyword_dct['hrthresh']
            increment = spc_dct_i.get('hind_inc', 30.0*phycon.DEG2RAD)

            zrxn = spc_dct_i.get('zrxn', None)

            run_tors_names = automol.rotor.names(rotors)
            run_tors_grids = automol.rotor.grids(
                rotors, increment=increment)

            # Set constraints
            const_names = automol.zmat.set_constraint_names(
                zma, run_tors_names, tors_model)

            # Read and print the potential
            ini_cnf_save_path = ini_cnf_save_fs[-1].path(ini_min_locs)
            sp_fs = autofile.fs.single_point(ini_cnf_save_path)
            ref_ene = sp_fs[-1].file.energy.read(mod_thy_info[1:4])
            tors_pots, tors_zmas, tors_paths = {}, {}, {}
            for tors_names, tors_grids in zip(
                    run_tors_names, run_tors_grids):
                constraint_dct = automol.zmat.constraint_dct(
                    zma, const_names, tors_names)
                pot, _, _, _, zmas, paths = filesys.read.potential(
                    tors_names, tors_grids,
                    ini_cnf_save_path,
                    mod_thy_info, ref_ene,
                    constraint_dct,
                    read_zma=True,
                    read_energy_backstep=True,
                    remove_bad_points=False)
                tors_pots[tors_names] = pot
                tors_zmas[tors_names] = zmas
                tors_paths[tors_names] = paths

            # Check for new minimum conformer
            new_min_zma = hr.check_hr_pot(
                tors_pots, tors_zmas, tors_paths, emax=ethresh)

            if new_min_zma is not None:
                ioprinter.info_message(
                    'Finding new low energy conformer...', newline=1)
                new_min_geo = automol.zmat.geometry(new_min_zma)
                rid = conformer.rng_loc_for_geo(
                    new_min_geo, cnf_save_fs)
                if rid is None:
                    new_locs = None
                else:
                    cid = autofile.schema.generate_new_conformer_id()
                    new_locs = (rid, cid)
                conformer.single_conformer(
                    new_min_zma, spc_info, mod_thy_info,
                    cnf_run_fs, cnf_save_fs,
                    script_str, overwrite,
                    retryfail=retryfail, zrxn=zrxn,
                    use_locs=new_locs, **kwargs)

        elif job in ('energy', 'grad', 'hess', 'vpt2'):

            # Script (add energy script call)
            script_str, kwargs = qchem_params(
                method_dct)

            ini_cnf_save_path = ini_cnf_save_fs[-1].path(ini_min_locs)
            ini_cnf_run_path = ini_cnf_run_fs[-1].path(ini_min_locs)
            ini_scn_run_fs, ini_scn_save_fs = build_fs(
                ini_cnf_run_path, ini_cnf_save_path, scn,
                zma_locs=zma_locs)
            run_tors_names = automol.rotor.names(rotors, flat=True)
            for tors_names in run_tors_names:

                # Set the constraint dct and filesys for the scan
                const_names = automol.zmat.set_constraint_names(
                    zma, [run_tors_names], tors_model)
                constraint_dct = automol.zmat.constraint_dct(
                    zma, const_names, tors_names)

                # get the scn_locs, maybe get a function?
                _, scn_locs = scan.scan_locs(
                    ini_scn_save_fs, (tors_names,),
                    constraint_dct=constraint_dct)
                for locs in scn_locs:
                    geo = ini_scn_save_fs[-1].file.geometry.read(locs)
                    zma = ini_scn_save_fs[-1].file.zmatrix.read(locs)
                    ini_scn_run_fs[-1].create(locs)
                    ES_TSKS[job](
                        zma, geo, spc_info, mod_thy_info,
                        ini_scn_run_fs, ini_scn_save_fs, locs, run_prefix,
                        script_str, overwrite,
                        zrxn=zrxn,
                        retryfail=retryfail, **kwargs)
                    ioprinter.obj('vspace')


def skip_task(tsk, spc_dct, spc_name, thy_dct, es_keyword_dct, save_prefix):
    """ Determine if an electronic structure task should be skipped based on
        various parameters.

        :param spc_dct: species dictionary
        :type spc_dct: dictionary
        :param spc_name: name of species
        :type spc_name: string
        :rtype: bool
    """

    # Initialize skip to be false
    skip = False

    # Set theory info needed to find information
    ini_method_dct = thy_dct.get(es_keyword_dct['inplvl'])
    ini_thy_info = tinfo.from_dct(ini_method_dct)

    # Perform checks
    if 'ts' in spc_name:
        # Skip all tasks except find_ts
        # if rad-rad TS
        if tsk not in ('find_ts', 'rpath_scan'):  # generalize to other rpath
            rxn_info = spc_dct[spc_name]['canon_rxn_info']
            ts_mul = rinfo.value(rxn_info, 'tsmult')
            high_ts_mul = rinfo.ts_mult(rxn_info, rxn_mul='high')
            if rinfo.radrad(rxn_info) and ts_mul != high_ts_mul:
                skip = True
                ioprinter.info_message(
                    f'Skipping task because {spc_name}',
                    'is a low-spin radical radical reaction')
    else:
        spc_natoms = len(automol.chi.geometry(spc_dct[spc_name]['inchi']))
        if spc_natoms == 1:
            # Skip all tasks except init_geom and conf_energy
            # if species is an atom
            if tsk not in ('init_geom', 'conf_energy', 'conf_prop'):
                skip = True
                ioprinter.info_message(
                    'Skipping task for an atom...', newline=1)
        else:
            # Skip all tasks except ini_geom
            # if (non-TS) species is unstable (zrxn found (i.e. is not None))
            if tsk != 'init_geom':
                instab, path = filesys.read.instability_transformation(
                    spc_dct, spc_name, ini_thy_info, save_prefix)

                skip = (instab is not None)
                if skip:
                    ioprinter.info_message(
                        f'Found instability file at path {path}', newline=1)
                    ioprinter.info_message(
                        'Skipping task for unstable species...', newline=1)

    return skip


def _sort_info_lst(sort_str, thy_dct, spc_info):
    """ Return the levels to sort conformers by if zpve or sp
        levels were assigned in input

        if we ask for zpe(lvl_wbs),sp(lvl_b2t),gibbs(700)
        out sort_info_lst will be [('gaussian', 'wb97xd', '6-31*', 'RU'),
        ('gaussian', 'b2plypd3', 'cc-pvtz', 'RU'), None, None, 700.]
    """
    sort_lvls = [None, None, None, None, None]
    sort_typ_lst = ['freqs', 'sp', 'enthalpy', 'entropy', 'gibbs']
    if sort_str is not None:
        for sort_param in sort_str.split(','):
            idx = None
            for typ_idx, typ_str in enumerate(sort_typ_lst):
                if typ_str in sort_param:
                    lvl_key = sort_str.split(typ_str + '(')[1].split(')')[0]
                    idx = typ_idx
            if idx is not None:
                if idx < 2:
                    method_dct = thy_dct.get(lvl_key)
                    if method_dct is None:
                        ioprinter.warning_message(
                            f'no {lvl_key} in theory.dat, '
                            f'not using {sort_typ_lst[idx]} in sorting')
                        continue
                    thy_info = tinfo.from_dct(method_dct)
                    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
                    sort_lvls[idx] = mod_thy_info
                else:
                    sort_lvls[idx] = float(lvl_key)
    return sort_lvls
