""" eletronic structure routines modules
"""

import sys
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
from mechroutines.es._ts import findts
from mechroutines.es._routines import conformer
from mechroutines.es._routines import hr
from mechroutines.es._routines import tau
# from mechroutines.es._routines import rpath
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
                run_prefix, save_prefix)
        elif 'tau' in tsk:
            tau_tsk(
                job, spc_dct, spc_name, thy_dct, es_keyword_dct,
                run_prefix, save_prefix)
        elif 'hr' in tsk:
            hr_tsk(
                job, spc_dct, spc_name, thy_dct, es_keyword_dct,
                run_prefix, save_prefix)
        elif 'rpath' in tsk:
            rpath_tsk(
                job, spc_dct, spc_name, thy_dct, es_keyword_dct,
                run_prefix, save_prefix)
        elif 'find' in tsk:
            findts(
                spc_dct, spc_name, thy_dct, es_keyword_dct,
                run_prefix, save_prefix)

        # return success # complicated, since success means different thing


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
        :param es_keyword_dct: keyword-value pairs for electronic structure task
        :type es_keyword_dct: dict[str:str]
        :param run_prefix: root-path to the run-filesystem
        :type run_prefix: str
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
    """

    spc_dct_i = spc_dct[spc_name]
    spc_info = sinfo.from_dct(spc_dct_i)

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
                  run_prefix, save_prefix):
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
        spc_info = sinfo.from_dct(spc_dct_i)
    else:
        spc_info = rinfo.ts_info(spc_dct_i['rxn_info'])
    zrxn = spc_dct_i.get('zrxn', None)

    overwrite = es_keyword_dct['overwrite']
    retryfail = es_keyword_dct['retryfail']

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
        ini_loc_info = filesys.mincnf.min_energy_conformer_locators(
            ini_cnf_save_fs, mod_ini_thy_info)
        ini_locs, ini_min_cnf_path = ini_loc_info

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
        zma = ini_zma_save_fs[-1].file.zmatrix.read([0])

        # Read the torsions from the ini file sys
        if ini_zma_save_fs[-1].file.torsions.exists([0]):
            tors_dct = ini_zma_save_fs[-1].file.torsions.read([0])
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
            tors_names=tors_names, zrxn=zrxn,
            two_stage=two_stage, retryfail=retryfail, resave=resave,
            **kwargs)

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
        zma = ini_zma_save_fs[-1].file.zmatrix.read([0])

        # Read the torsions from the ini file sys
        if ini_zma_save_fs[-1].file.ring_torsions.exists([0]):
            ring_tors_dct = ini_zma_save_fs[-1].file.ring_torsions.read([0])
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
        cnf_sort_info_lst = _sort_info_lst(
            es_keyword_dct['sort'], thy_dct, spc_info)

        # Set up the run scripts
        script_str, kwargs = qchem_params(
            method_dct, elstruct.Job.OPTIMIZATION)

        rng_cnf_locs_lst, _ = filesys.mincnf.conformer_locators(
            cnf_save_fs, mod_thy_info,
            cnf_range='all')

        ini_rng_cnf_locs_lst, _ = filesys.mincnf.conformer_locators(
            ini_cnf_save_fs, mod_ini_thy_info,
            cnf_range=cnf_range, sort_info_lst=cnf_sort_info_lst,
            print_enes=True)

        # Truncate the list of the ini confs
        uni_rng_locs_lst, uni_cnf_locs_lst = conformer.unique_fs_ring_confs(
            cnf_save_fs, rng_cnf_locs_lst,
            ini_cnf_save_fs, ini_rng_cnf_locs_lst)
        ioprinter.debug_message(
            'uni lst that has no similar ring', uni_rng_locs_lst)
        ioprinter.debug_message(
            'uni lst that has similar ring', uni_cnf_locs_lst)

        for locs in uni_rng_locs_lst:
            rid, cid = locs
            # Obtain the zma from ini loc
            ini_cnf_save_path = ini_cnf_save_fs[-1].path(locs)
            ini_zma_save_fs = autofile.fs.zmatrix(ini_cnf_save_path)
            zma = ini_zma_save_fs[-1].file.zmatrix.read((0,))

            # Make the ring filesystem
            conformer.single_conformer(
                zma, spc_info, mod_thy_info,
                cnf_run_fs, cnf_save_fs,
                script_str, overwrite,
                retryfail=retryfail, zrxn=zrxn,
                use_locs=locs,
                **kwargs)

        for locs in uni_cnf_locs_lst:
            ini_locs, rid = locs
            # Obtain the zma from ini loc
            ini_cnf_save_path = ini_cnf_save_fs[-1].path(ini_locs)
            ini_zma_save_fs = autofile.fs.zmatrix(ini_cnf_save_path)
            zma = ini_zma_save_fs[-1].file.zmatrix.read((0,))
            # obtain conformer filesys associated with ring at the runlevel
            cid = autofile.schema.generate_new_conformer_id()
            conformer.single_conformer(
                zma, spc_info, mod_thy_info,
                cnf_run_fs, cnf_save_fs,
                script_str, overwrite,
                retryfail=retryfail, zrxn=zrxn,
                use_locs=(rid, cid),
                **kwargs)

        # print all geometres within cnfrange
        rng_cnf_locs_lst, _ = filesys.mincnf.conformer_locators(
            cnf_save_fs, mod_thy_info,
            cnf_range=cnf_range, sort_info_lst=cnf_sort_info_lst)
        for locs in rng_cnf_locs_lst:
            geo = cnf_save_fs[-1].file.geometry.read(locs)
            ioprinter.geometry(geo)

    elif job in ('energy', 'grad', 'hess', 'vpt2', 'prop'):

        cnf_range = es_keyword_dct['cnf_range']
        cnf_sort_info_lst = _sort_info_lst(
            es_keyword_dct['sort'], thy_dct, spc_info)

        ini_rng_cnf_locs_lst, _ = filesys.mincnf.conformer_locators(
            ini_cnf_save_fs, mod_ini_thy_info,
            cnf_range=cnf_range, sort_info_lst=cnf_sort_info_lst,
            print_enes=True)

        # Check if locs exist, kill if it doesn't
        if not ini_rng_cnf_locs_lst:
            ioprinter.error_message(
                'No min-energy conformer found for level:')
            sys.exit()

        # Set up the run scripts
        script_str, kwargs = qchem_params(
            method_dct)

        # Run the job over all the conformers requested by the user
        for ini_locs in ini_rng_cnf_locs_lst:
            ini_cnf_run_fs[-1].create(ini_locs)
            geo_save_path = ini_cnf_save_fs[-1].path(ini_locs)
            ini_zma_save_fs = autofile.fs.zmatrix(geo_save_path)
            print('Running task for geometry at ', geo_save_path)
            geo = ini_cnf_save_fs[-1].file.geometry.read(ini_locs)
            zma = ini_zma_save_fs[-1].file.zmatrix.read((0,))
            ES_TSKS[job](
                zma, geo, spc_info, mod_thy_info,
                ini_cnf_run_fs, ini_cnf_save_fs, ini_locs,
                script_str, overwrite,
                retryfail=retryfail, **kwargs)


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
    spc_info = sinfo.from_dct(spc_dct_i)

    # Get es options
    overwrite = es_keyword_dct['overwrite']
    retryfail = es_keyword_dct['retryfail']
    # scan_increment = spc_dct_i['hind_inc']
    nsamp_par = spc_dct_i['tau_nsamp']

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
    zma = ini_zma_save_fs[-1].file.zmatrix.read((0,))
    ini_sp_save_fs = autofile.fs.single_point(ini_min_cnf_path)
    if ini_sp_save_fs[-1].file.energy.exists(mod_ini_thy_info[1:4]):
        ref_ene = ini_sp_save_fs[-1].file.energy.read(mod_ini_thy_info[1:4])
    else:
        ref_ene = ini_cnf_save_fs[-1].file.energy.read(ini_min_cnf_locs)

    # Get the tors names
    ini_zma_save_fs = autofile.fs.zmatrix(ini_min_cnf_path)
    if ini_zma_save_fs[-1].file.torsions.exists([0]):
        tors_dct = ini_zma_save_fs[-1].file.torsions.read([0])
        torsions = automol.rotor.from_data(zma, tors_dct)
    else:
        torsions = ()

    saddle = bool('ts_' in spc_name)

    # Run the task if any torsions exist
    if torsions and not saddle:

        # Set up tau filesystem objects
        tau_run_fs, tau_save_fs = build_fs(
            run_prefix, save_prefix, 'TAU',
            spc_locs=spc_info, thy_locs=mod_thy_info[1:])
        # db_style = 'jsondb'
        db_style = 'directory'
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
            # Run sampling
            tau.tau_sampling(
                zma, ref_ene,
                spc_info, tors_names, nsamp_par,
                mod_ini_thy_info,
                tau_run_fs, tau_save_fs,
                script_str, overwrite,
                saddle=saddle, **kwargs)

        elif job in ('energy', 'grad'):

            # Set up the run scripts
            script_str, kwargs = qchem_params(
                method_dct)

            # Run the job over all the conformers requested by the user
            for locs in tau_save_fs[-1].existing():
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
                    tau_save_fs, locs,
                    script_str, overwrite,
                    retryfail=retryfail, **kwargs)
                ioprinter.obj('vspace')

        elif job == 'hess':

            # Add the hessian max
            hessmax = es_keyword_dct['hessmax']

            # Set up the run scripts
            script_str, kwargs = qchem_params(
                method_dct)

            # Run the job over all the conformers requested by the user
            hess_cnt = 0
            for locs in tau_save_fs.existing():
                ioprinter.info_message(
                    'HESS Number {}'.format(hess_cnt+1), newline=1)
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
                    tau_run_fs, tau_save_fs, locs,
                    script_str, overwrite,
                    retryfail=retryfail, **kwargs)
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
        :param es_keyword_dct: keyword-value pairs for electronic structure task
        :type es_keyword_dct: dict[str:str]
        :param run_prefix: root-path to the run-filesystem
        :type run_prefix: str
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
    """

    spc_dct_i = spc_dct[spc_name]
    saddle = bool('ts_' in spc_name)
    # Set the spc_info
    if not saddle:
        spc_info = sinfo.from_dct(spc_dct_i)
    else:
        spc_info = rinfo.ts_info(spc_dct_i['rxn_info'])

    # Get options from the dct or es options lst
    overwrite = es_keyword_dct['overwrite']
    retryfail = es_keyword_dct['retryfail']
    tors_model = es_keyword_dct['tors_model']

    # Modify the theory
    method_dct = thy_dct.get(es_keyword_dct['runlvl'])
    ini_method_dct = thy_dct.get(es_keyword_dct['inplvl'])
    thy_info = tinfo.from_dct(method_dct)
    ini_thy_info = tinfo.from_dct(ini_method_dct)
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
    mod_ini_thy_info = tinfo.modify_orb_label(
        ini_thy_info, spc_info)

    # Set the filesystem objects
    _root = root_locs(spc_dct_i, saddle=saddle, name=spc_name)
    ini_cnf_run_fs, ini_cnf_save_fs = build_fs(
        run_prefix, save_prefix, 'CONFORMER',
        thy_locs=mod_ini_thy_info[1:],
        **_root)
    cnf_run_fs, cnf_save_fs = build_fs(
        run_prefix, save_prefix, 'CONFORMER',
        thy_locs=mod_thy_info[1:],
        **_root)

    cnf_range = es_keyword_dct['cnf_range']
    cnf_sort_info_lst = _sort_info_lst(
        es_keyword_dct['sort'], thy_dct, spc_info)
    ini_min_locs_lst, ini_path_lst = filesys.mincnf.conformer_locators(
        ini_cnf_save_fs, mod_ini_thy_info,
        cnf_range=cnf_range, sort_info_lst=cnf_sort_info_lst,
        print_enes=True)
    all_run_cnf_locs_lst, _ = filesys.mincnf.conformer_locators(
        cnf_save_fs, mod_thy_info,
        cnf_range='all')
    ini_to_run_locs_dct = conformer.fs_confs_dict(
        cnf_save_fs, all_run_cnf_locs_lst,
        ini_cnf_save_fs, ini_min_locs_lst)

    # ini_loc_info = filesys.mincnf.min_energy_conformer_locators(
    #     ini_cnf_save_fs, mod_ini_thy_info)
    # ini_min_locs, ini_cnf_save_path = ini_loc_info

    for ini_min_locs, ini_cnf_save_path in zip(ini_min_locs_lst, ini_path_lst):


        # Read zma, geo, and torsions
        ini_zma_save_fs = autofile.fs.zmatrix(ini_cnf_save_path)
        geo = ini_cnf_save_fs[-1].file.geometry.read(ini_min_locs)
        zma = ini_zma_save_fs[-1].file.zmatrix.read((0,))
        if ini_zma_save_fs[-1].file.torsions.exists([0]):
            tors_dct = ini_zma_save_fs[-1].file.torsions.read([0])
            torsions = automol.rotor.from_data(zma, tors_dct,)
        else:
            torsions = ()

        # Find equivalent conformer in the run filesys, if it doesn't exist
        # run a single conformer to generate it
        min_locs = ini_to_run_locs_dct[tuple(ini_min_locs)]
        if min_locs is None:
            rid = conformer.rng_loc_for_geo(geo, cnf_save_fs)
            if rid is None:
                new_rid = autofile.schema.generate_new_ring_id()
                new_cid = autofile.schema.generate_new_conformer_id()
                conformer.single_conformer(
                    zma, spc_info, mod_thy_info,
                    cnf_run_fs, cnf_save_fs,
                    script_str, overwrite,
                    retryfail=retryfail, zrxn=zrxn,
                    use_locs = (new_rid, new_cid),
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
        cnf_save_path = cnf_save_fs[-1].path(min_locs)
        ioprinter.info_message('Same conformer saved at {} and {}'.format(
            ini_cnf_save_path, cnf_save_path))

        # Create run fs if that directory has been deleted to run the jobs
        # ini_cnf_run_fs[-1].create(ini_min_locs)
        # ini_cnf_run_path = ini_cnf_run_fs[-1].path(ini_min_locs)
        cnf_run_fs[-1].create(min_locs)
        cnf_run_path = cnf_run_fs[-1].path(min_locs)

        # Run the task if any torsions exist
        if any(torsions):
            if 'fa' in tors_model:
                scn = 'CSCAN'
            elif 'f' in tors_model:
                if len(torsions) > 1:
                    scn = 'CSCAN'
                else:
                    scn = 'SCAN'
            else:
                scn = 'SCAN'
            # ini_scn_run_fs, ini_scn_save_fs = build_fs(
            #     ini_cnf_run_path, ini_cnf_save_path, scn,
            #     zma_locs=(0,))
            scn_run_fs, scn_save_fs = build_fs(
                cnf_run_path, cnf_save_path, scn,
                zma_locs=(0,))

            if job == 'scan':

                increment = spc_dct_i.get('hind_inc', 30.0*phycon.DEG2RAD)
                hr.hindered_rotor_scans(
                    zma, spc_info, mod_thy_info,
                    scn_run_fs, scn_save_fs,
                    torsions, tors_model, method_dct,
                    overwrite,
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

                run_tors_names = automol.rotor.names(torsions)
                run_tors_grids = automol.rotor.grids(
                    torsions, increment=increment)

                # Set constraints
                const_names = automol.zmat.set_constraint_names(
                    zma, run_tors_names, tors_model)

                # Read and print the potential
                sp_fs = autofile.fs.single_point(cnf_save_path)
                ref_ene = sp_fs[-1].file.energy.read(mod_thy_info[1:4])
                tors_pots, tors_zmas, tors_paths = {}, {}, {}
                for tors_names, tors_grids in zip(
                        run_tors_names, run_tors_grids):
                    constraint_dct = automol.zmat.constraint_dct(
                        zma, const_names, tors_names)
                    pot, _, _, _, zmas, paths = filesys.read.potential(
                        tors_names, tors_grids,
                        cnf_save_path,
                        mod_thy_info, ref_ene,
                        constraint_dct,
                        read_zma=True)
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

                run_tors_names = automol.rotor.names(torsions, flat=True)
                for tors_names in run_tors_names:

                    # Set the constraint dct and filesys for the scan
                    const_names = automol.zmat.set_constraint_names(
                        zma, run_tors_names, tors_model)
                    constraint_dct = automol.zmat.constraint_dct(
                        zma, const_names, tors_names)

                    # get the scn_locs, maybe get a function?
                    _, scn_locs = scan.scan_locs(
                        scn_save_fs, tors_names,
                        constraint_dct=constraint_dct)
                    for locs in scn_locs:
                        geo = scn_save_fs[-1].file.geometry.read(locs)
                        zma = scn_save_fs[-1].file.zmatrix.read(locs)
                        scn_run_fs[-1].create(locs)
                        ES_TSKS[job](
                            zma, geo, spc_info, mod_thy_info,
                            scn_run_fs, scn_save_fs, locs,
                            script_str, overwrite,
                            retryfail=retryfail, **kwargs)
                        ioprinter.obj('vspace')
        else:
            ioprinter.info_message('No torsional modes in the species')


def rpath_tsk(job, spc_dct, spc_name,
              thy_dct, es_keyword_dct,
              run_prefix, save_prefix):
    """ run a scan over the specified torsional coordinates

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
    spc_dct_i = spc_dct[spc_name]

    # Set up coordinate name
    rxn_coord = es_keyword_dct.get('rxn_coord')
    if rxn_coord == 'auto':
        coord_name = ['Rn']  # grab from zrxn object
    else:
        coord_name = ['IRC']

    # Set the spc_info
    spc_info = sinfo.from_dct(spc_dct_i)

    # Modify the theory
    method_dct = thy_dct.get(es_keyword_dct['runlvl'])
    ini_method_dct = thy_dct.get(es_keyword_dct['inplvl'])
    thy_info = tinfo.from_dct(method_dct)
    ini_thy_info = tinfo.from_dct(ini_method_dct)
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
    mod_ini_thy_info = tinfo.modify_orb_label(
        ini_thy_info, spc_info)

    # Get options from the dct or es options lst
    overwrite = es_keyword_dct['overwrite']

    # Set up the script
    script_str, kwargs = qchem_params(
        method_dct, elstruct.Job.OPTIMIZATION)

    # Set the filesystem objects
    # rxn_info = spc_dct_i['rxn_info']
    # fs_rxn_info = rinfo.sort(rxn_info)

    # New filesystem objects
    if coord_name == 'irc':
        _root = root_locs(spc_dct_i, saddle=True)
        ini_cnf_run_fs, ini_cnf_save_fs = build_fs(
            run_prefix, save_prefix, 'CONFORMER',
            thy_locs=mod_ini_thy_info[1:],
            **_root)
        ini_loc_info = filesys.mincnf.min_energy_conformer_locators(
            ini_cnf_save_fs, mod_ini_thy_info)
        ini_min_locs, ini_pfx_save_path = ini_loc_info
        ini_cnf_run_fs[-1].create(ini_min_locs)
        ini_pfx_run_path = ini_cnf_run_fs[-1].path(ini_min_locs)

    else:
        # ts_info = (ts_num,)
        ts_info = (0,)
        ini_ts_run_fs, ini_ts_save_fs = build_fs(
            run_prefix, save_prefix, 'TS',
            thy_locs=mod_ini_thy_info[1:],
            **_root)
        ini_pfx_run_path = ini_ts_run_fs.path(ts_info)
        ini_pfx_save_path = ini_ts_save_fs.path(ts_info)

    # Get options from the dct or es options lst
    overwrite = es_keyword_dct['overwrite']

    ini_scn_run_fs, ini_scn_save_fs = build_fs(
        ini_pfx_run_path, ini_pfx_save_path, 'SCAN',
        zma_locs=(0,))

    # ini_zma_save_fs = autofile.fs.zmatrix(ini_cnf_save_path)
    geo = ini_cnf_save_fs[-1].file.geometry.read(ini_min_locs)
    # zma = ini_zma_save_fs[-1].file.zmatrix.read((0,))

    # Run job
    if job == 'scan':
        pass
        # if rxn_coord == 'auto':
        #     pass
        # elif rxn_coord == 'irc':
        #     rpath.irc_scan(
        #         geo, spc_info, coord_name,
        #         mod_ini_thy_info, ini_method_dct,
        #         ini_scn_save_fs, ini_cnf_run_path,
        #         overwrite)

    elif job in ('energy', 'grad', 'hess'):

        # Script
        script_str, kwargs = qchem_params(
            method_dct)

        # Need to put in something with the IRC idxs
        for locs in ini_scn_save_fs[-1].existing():
            geo = ini_scn_save_fs[-1].file.geometry.read(locs)
            ini_scn_run_fs[-1].create(locs)
            ES_TSKS[job](
                None, geo, spc_info, mod_thy_info,
                ini_scn_run_fs, ini_scn_save_fs, locs,
                script_str, overwrite, **kwargs)
            ioprinter.obj('vspace')

    elif Gjob == 'infene':
        Gpass
        # inf_sep_ene()


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
        if tsk != 'find_ts':
            rxn_info = spc_dct[spc_name]['rxn_info']
            ts_mul = rinfo.value(rxn_info, 'tsmult')
            high_ts_mul = rinfo.ts_mult(rxn_info, rxn_mul='high')
            if rinfo.radrad(rxn_info) and ts_mul != high_ts_mul:
                skip = True
                ioprinter.info_message(
                    'Skipping task because {}'.format(spc_name),
                    'is a low-spin radical radical reaction')
    else:
        spc_natoms = len(automol.inchi.geometry(spc_dct[spc_name]['inchi']))
        if spc_natoms == 1:
            # Skip all tasks except init_geom and conf_energy
            # if species is an atom
            if tsk not in ('init_geom', 'conf_energy'):
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
                        'Found instability file at path {}'.format(path),
                        newline=1)
                    ioprinter.info_message(
                        'Skipping task for unstable species...', newline=1)

    return skip
        

def _sort_info_lst(sort_str, thy_dct, spc_info):
    """ Return the levels to sort conformers by if zpve or sp
        levels were assigned in input
    """
    sort_lvls = [None, None]
    sort_typ_lst = ['zpe', 'sp']
    if sort_str is not None:
        for sort_param in sort_str.split(','):
            idx = None
            for typ_idx, typ_str in enumerate(sort_typ_lst):
                if typ_str in sort_param:
                    lvl_key = sort_str.split(typ_str + '(')[1].split(')')[0]
                    idx = typ_idx
            if idx is not None:
                method_dct = thy_dct.get(lvl_key)
                if method_dct is None:
                    ioprinter.warning_message(
                        'no {} in theory.dat, not using {} in sorting'.format(
                            lvl_key, sort_typ_lst[idx]))
                    continue
                thy_info = tinfo.from_dct(method_dct)
                mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
                sort_lvls[idx] = mod_thy_info
    return sort_lvls
