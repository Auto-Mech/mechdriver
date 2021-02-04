""" eletronic structure routines modules
"""

import sys
import importlib
import autofile
import automol
import elstruct
from mechanalyzer.inf import thy as tinfo
from mechanalyzer.inf import rxn as rinfo
from mechanalyzer.inf import spc as sinfo
from mechroutines.es._ts import findts
from mechroutines.es._routines import conformer
from mechroutines.es._routines import hr
from mechroutines.es._routines import tau
from mechroutines.es._routines import irc
from mechlib import filesys
from mechlib import structure
from mechlib.structure import instab
from mechlib.submission import qchem_params
from mechlib.amech_io import printer as ioprinter


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
            run_prefix, save_prefix):
    """ run an electronic structure task
    for generating a list of conformer or tau sampling geometries
    """

    # Print the head of the task
    ioprinter.task_header(tsk, spc_name)
    ioprinter.keyword_list(es_keyword_dct, thy_dct)

    # If species is unstable, set task to 'none'
    ini_method_dct = thy_dct.get(es_keyword_dct['inplvl'])
    ini_thy_info = tinfo.from_dct(ini_method_dct)
    stable = instab.check_unstable_species(
        tsk, spc_dct, spc_name, ini_thy_info, save_prefix)
    stable = True
    if stable:
        ioprinter.debug_message('- Proceeding with requested task...')

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
        elif 'irc' in tsk:
            irc_tsk(
                job, spc_dct, spc_name, thy_dct, es_keyword_dct,
                run_prefix, save_prefix)
        elif 'find' in tsk:
            findts(
                job, spc_dct, spc_name, thy_dct, es_keyword_dct,
                run_prefix, save_prefix)

    else:
        ioprinter.info_message(
            'Skipping task for unstable species...', newline=1)


# FUNCTIONS FOR SAMPLING AND SCANS #
def geom_init(spc_dct, spc_name, thy_dct, es_keyword_dct,
              run_prefix, save_prefix):
    """ Find the initial geometry
    """

    spc_dct_i = spc_dct[spc_name]

    # Set the spc_info
    spc_info = sinfo.from_dct(spc_dct_i)

    # Get es options
    [kickoff_size, kickoff_backward] = spc_dct_i['kickoff']
    overwrite = es_keyword_dct['overwrite']
    # retryfail = es_keyword_dct['retryfail']

    # Get the theory info
    method_dct = thy_dct.get(es_keyword_dct['runlvl'])
    ini_method_dct = thy_dct.get(es_keyword_dct['inplvl'])
    print('dct\n', method_dct)
    thy_info = tinfo.from_dct(method_dct)
    ini_thy_info = tinfo.from_dct(ini_method_dct)
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
    mod_ini_thy_info = tinfo.modify_orb_label(
        ini_thy_info, spc_info)

    # Set the filesystem objects
    thy_run_fs, thy_run_path = filesys.build.spc_thy_fs_from_root(
        run_prefix, spc_info, mod_thy_info)
    thy_save_fs, thy_save_path = filesys.build.spc_thy_fs_from_root(
        save_prefix, spc_info, mod_thy_info)
    _, ini_thy_save_path = filesys.build.spc_thy_fs_from_root(
        save_prefix, spc_info, mod_ini_thy_info)
    cnf_save_fs = autofile.fs.conformer(thy_save_path)

    # Set the run filesystem
    saddle = bool('ts_' in spc_name)
    if saddle:
        _, ts_path = filesys.build.ts_fs_from_thy(thy_run_path)
        run_fs = filesys.build.run_fs_from_prefix(ts_path)
    else:
        run_fs = filesys.build.run_fs_from_prefix(thy_run_path)

    # Get a reference geometry if one not found
    geo = conformer.initial_conformer(
        spc_dct_i, spc_info, method_dct,
        thy_run_fs, thy_save_fs,
        cnf_save_fs,
        ini_thy_save_path, mod_ini_thy_info,
        run_fs,
        overwrite,
        kickoff_size=kickoff_size,
        kickoff_backward=kickoff_backward)

    return geo


def conformer_tsk(job, spc_dct, spc_name,
                  thy_dct, es_keyword_dct,
                  run_prefix, save_prefix):
    """ Launch tasks associated with conformers.
        Scan: Generate a set of conformer geometries and energies via
              random sampling over torsional coordinates
              following by optimization
        SP: Calculate ene, grad, ..
    """

    saddle = bool('ts_' in spc_name)

    spc_dct_i = spc_dct[spc_name]

    # Set the spc_info
    if not saddle:
        spc_info = sinfo.from_dct(spc_dct_i)
        zrxn = None
    else:
        spc_info = rinfo.ts_info(spc_dct_i['rxn_info'])
        zrxn = spc_dct_i['zrxn']

    # Get es options
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

    # Set the filesystem objects
    if not saddle:
        _, ini_thy_run_path = filesys.build.spc_thy_fs_from_root(
            run_prefix, spc_info, mod_ini_thy_info)
        ini_thy_fs = filesys.build.spc_thy_fs_from_root(
            save_prefix, spc_info, mod_ini_thy_info)
        [_, ini_thy_save_path] = ini_thy_fs

    else:
        rxn_info = spc_dct_i['rxn_info']
        fs_rxn_info = rinfo.sort(rxn_info)

        _, ini_thy_run_path = filesys.build.rxn_thy_fs_from_root(
            run_prefix, fs_rxn_info, mod_ini_thy_info)
        _, ini_thy_save_path = filesys.build.rxn_thy_fs_from_root(
            save_prefix, fs_rxn_info, mod_ini_thy_info)
        _, ini_thy_save_path = filesys.build.ts_fs_from_thy(
            ini_thy_save_path)
        _, ini_thy_run_path = filesys.build.ts_fs_from_thy(
            ini_thy_run_path)

    # Build the thy filesyste
    if not saddle:
        _, thy_run_path = filesys.build.spc_thy_fs_from_root(
            run_prefix, spc_info, mod_thy_info)
        thy_save_fs, thy_save_path = filesys.build.spc_thy_fs_from_root(
            save_prefix, spc_info, mod_thy_info)
    else:
        _, thy_run_path = filesys.build.rxn_thy_fs_from_root(
            run_prefix, fs_rxn_info, mod_thy_info)
        thy_save_fs, thy_save_path = filesys.build.rxn_thy_fs_from_root(
            save_prefix, fs_rxn_info, mod_thy_info)
        thy_save_fs, thy_save_path = filesys.build.ts_fs_from_thy(
            thy_save_path)
        _, thy_run_path = filesys.build.ts_fs_from_thy(thy_run_path)

    if job == 'samp':

        # Build conformer filesys
        cnf_run_fs = autofile.fs.conformer(thy_run_path)
        cnf_save_fs = autofile.fs.conformer(thy_save_path)

        # Build the ini zma filesys
        # ini_cnf_run_fs = autofile.fs.conformer(ini_thy_run_path)
        ini_cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)
        ini_loc_info = filesys.mincnf.min_energy_conformer_locators(
            ini_cnf_save_fs, mod_ini_thy_info)
        ini_min_cnf_locs, ini_min_cnf_path = ini_loc_info
        ini_zma_save_fs, _ = filesys.build.zma_fs_from_prefix(
            ini_min_cnf_path, zma_idxs=[0])

        # Set up the run scripts
        script_str, kwargs = qchem_params(
            method_dct, elstruct.Job.OPTIMIZATION)

        # Set variables if it is a saddle
        two_stage = saddle
        mc_nsamp = spc_dct_i['mc_nsamp']

        # Read the geometry and zma from the ini file system
        geo = ini_cnf_save_fs[-1].file.geometry.read(ini_min_cnf_locs)
        zma = ini_zma_save_fs[-1].file.zmatrix.read([0])
        ioprinter.debug_message('ini path', ini_thy_save_path)

        # Read the torsions from the ini file sys
        if ini_zma_save_fs[-1].file.torsions.exists([0]):
            tors_dct = ini_zma_save_fs[-1].file.torsions.read([0])
            torsions = automol.rotor.from_data(zma, tors_dct)
            print('tors_path', ini_zma_save_fs[-1].file.torsions.path([0]))
            print('tors_lst', tors_dct)
            tors_names = automol.rotor.names(torsions, flat=True)
        else:
            tors_names = ()

        geo_path = ini_cnf_save_fs[-1].path(ini_min_cnf_locs)
        ioprinter.initial_geom_path('Sampling started', geo_path)

        # Run the sampling
        _ = conformer.conformer_sampling(
            zma, spc_info,
            mod_thy_info, thy_save_fs,
            cnf_run_fs, cnf_save_fs,
            script_str, overwrite,
            zrxn=zrxn, nsamp_par=mc_nsamp,
            tors_names=tors_names,
            two_stage=two_stage, retryfail=retryfail,
            **kwargs)

    elif job == 'opt':

        cnf_range = es_keyword_dct['cnf_range']
        ioprinter.debug_message('range', cnf_range)

        # Set up the run scripts
        script_str, kwargs = qchem_params(
            method_dct, elstruct.Job.OPTIMIZATION)

        # Build conformer filesys
        cnf_run_fs = autofile.fs.conformer(thy_run_path)

        cnf_save_fs = autofile.fs.conformer(thy_save_path)
        cnf_locs_lst, _ = filesys.mincnf.conformer_locators(
            cnf_save_fs, mod_thy_info, cnf_range='all')

        ini_cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)
        ini_cnf_locs_lst, _ = filesys.mincnf.conformer_locators(
            ini_cnf_save_fs, mod_ini_thy_info, cnf_range=cnf_range)

        # Truncate the list of the ini confs
        uni_locs_lst = conformer.unique_fs_confs(
            cnf_save_fs, cnf_locs_lst, ini_cnf_save_fs, ini_cnf_locs_lst)
        ioprinter.debug_message('uni lst', uni_locs_lst)
        for locs in uni_locs_lst:
            ini_zma_save_fs = autofile.fs.zmatrix(ini_min_cnf_path)
            geo = ini_cnf_save_fs[-1].file.geometry.read(ini_min_cnf_locs)
            zma = ini_zma_save_fs[-1].file.zmatrix.read((0,))
            cnf_run_fs[-1].create(locs)
            conformer.single_conformer(
                zma, spc_info, mod_thy_info,
                cnf_run_fs, cnf_save_fs,
                script_str, overwrite,
                retryfail=retryfail, saddle=saddle, **kwargs)

        print_cnf_locs_lst, _ = filesys.mincnf.conformer_locators(
            cnf_save_fs, mod_thy_info, cnf_range=cnf_range)
        for locs in print_cnf_locs_lst:
            geo = cnf_save_fs[-1].file.geometry.read(locs)
            ioprinter.geometry(geo)

    elif job in ('energy', 'grad', 'hess', 'vpt2', 'prop'):

        cnf_range = es_keyword_dct['cnf_range']
        # print('range', cnf_range)

        ini_cnf_run_fs = autofile.fs.conformer(ini_thy_run_path)

        ini_cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)
        ini_cnf_save_locs_lst, _ = filesys.mincnf.conformer_locators(
            ini_cnf_save_fs, mod_ini_thy_info, cnf_range=cnf_range)

        # Check if locs exist, kill if it doesn't
        if not ini_cnf_save_locs_lst:
            ioprinter.error_message(
                'No min-energy conformer found for level:', ini_thy_save_path)
            sys.exit()

        # Set up the run scripts
        script_str, kwargs = qchem_params(
            method_dct)

        # Run the job over all the conformers requested by the user
        for locs in ini_cnf_save_locs_lst:
            ioprinter.running('task for conformer: ', locs, newline=2)
            geo_run_path = ini_cnf_run_fs[-1].path(locs)
            geo_save_path = ini_cnf_save_fs[-1].path(locs)
            ini_cnf_run_fs[-1].create(locs)
            ini_zma_save_fs = autofile.fs.zmatrix(geo_save_path)
            geo = ini_cnf_save_fs[-1].file.geometry.read(locs)
            zma = ini_zma_save_fs[-1].file.zmatrix.read((0,))
            ES_TSKS[job](
                zma, geo, spc_info, mod_thy_info,
                ini_cnf_save_fs, geo_run_path, geo_save_path, locs,
                script_str, overwrite,
                retryfail=retryfail, **kwargs)


def tau_tsk(job, spc_dct, spc_name,
            thy_dct, es_keyword_dct,
            run_prefix, save_prefix):
    """ Energies, gradients, and hessians,
        for set of arbitrarily sampled torsional coordinates
        with all other coordinates optimized
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
    _, thy_run_path = filesys.build.spc_thy_fs_from_root(
        run_prefix, spc_info, mod_thy_info)
    _, thy_save_path = filesys.build.spc_thy_fs_from_root(
        save_prefix, spc_info, mod_thy_info)

    # Set the filesystem objects for ini thy_info
    _, ini_thy_save_path = filesys.build.spc_thy_fs_from_root(
        save_prefix, spc_info, mod_ini_thy_info)

    # Get the geom and energy of reference species
    ini_cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)
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
        tau_run_fs, _ = filesys.build.tau_fs_from_thy(
            thy_run_path, tau='all')
        tau_save_fs, tau_save_locs = filesys.build.tau_fs_from_thy(
            thy_save_path, tau='all')
        # db_style = 'jsondb'
        db_style = 'directory'
        if db_style == 'jsondb':
            tau_save_fs[-1].root.create()
            tau_save_fs[-1].json_create()
            for locs in tau_save_locs:
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
            tau_save_locs = tau_save_fs[-1].json_existing()

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
            for locs in tau_save_locs:
                geo_run_path = tau_run_fs[-1].path(locs)
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
                    tau_save_fs, geo_run_path, geo_save_path, locs,
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
            for locs in tau_save_locs:
                ioprinter.info_message(
                    'HESS Number {}'.format(hess_cnt+1), newline=1)
                geo_run_path = tau_run_fs[-1].path(locs)
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
                zma = None
                tau_run_fs[-1].create(locs)
                ES_TSKS[job](
                    zma, geo, spc_info, mod_thy_info,
                    tau_save_fs, geo_run_path, geo_save_path, locs,
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
    """ run a scan over the specified torsional coordinates
    """

    spc_dct_i = spc_dct[spc_name]
    saddle = bool('ts_' in spc_name)
    # Set the spc_info
    if not saddle:
        spc_info = sinfo.from_dct(spc_dct_i)
    else:
        spc_info = rinfo.ts_info(spc_dct_i['rxn_info'])

    # Modify the theory
    method_dct = thy_dct.get(es_keyword_dct['runlvl'])
    ini_method_dct = thy_dct.get(es_keyword_dct['inplvl'])
    thy_info = tinfo.from_dct(method_dct)
    ini_thy_info = tinfo.from_dct(ini_method_dct)
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
    mod_ini_thy_info = tinfo.modify_orb_label(ini_thy_info, spc_info)

    # Set the filesystem objects
    if not saddle:
        _, ini_thy_run_path = filesys.build.spc_thy_fs_from_root(
            run_prefix, spc_info, mod_ini_thy_info)
        inifs = filesys.build.spc_thy_fs_from_root(
            save_prefix, spc_info, mod_ini_thy_info)
        [ini_thy_save_fs, ini_thy_save_path] = inifs
    else:
        rxn_info = spc_dct_i['rxn_info']
        fs_rxn_info = rinfo.sort(rxn_info)

        _, ini_thy_run_path = filesys.build.rxn_thy_fs_from_root(
            run_prefix, fs_rxn_info, mod_ini_thy_info)
        _, ini_thy_save_path = filesys.build.rxn_thy_fs_from_root(
            save_prefix, fs_rxn_info, mod_ini_thy_info)
        inifs = filesys.build.ts_fs_from_thy(
            ini_thy_save_path)
        [ini_thy_save_fs, ini_thy_save_path] = inifs
        _, ini_thy_run_path = filesys.build.ts_fs_from_thy(
            ini_thy_run_path)

    # Build cnf filesys using the ini thy filesys (needed for all HR jobs)
    ini_cnf_run_fs = autofile.fs.conformer(ini_thy_run_path)
    ini_cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)

    ini_loc_info = filesys.mincnf.min_energy_conformer_locators(
        ini_cnf_save_fs, mod_ini_thy_info)
    ini_min_cnf_locs, ini_cnf_save_path = ini_loc_info

    ini_cnf_run_path = filesys.build.cnf_paths_from_locs(
        ini_cnf_run_fs, [ini_min_cnf_locs])[0]

    # Create run fs if that directory has been deleted to run the jobs
    ini_cnf_run_fs[-1].create(ini_min_cnf_locs)

    # Get options from the dct or es options lst
    overwrite = es_keyword_dct['overwrite']
    retryfail = es_keyword_dct['retryfail']
    tors_model = es_keyword_dct['tors_model']

    # Read zma, geo, and torsions
    ini_zma_save_fs = autofile.fs.zmatrix(ini_cnf_save_path)
    geo = ini_cnf_save_fs[-1].file.geometry.read(ini_min_cnf_locs)
    zma = ini_zma_save_fs[-1].file.zmatrix.read((0,))
    if ini_zma_save_fs[-1].file.torsions.exists([0]):
        tors_dct = ini_zma_save_fs[-1].file.torsions.read([0])
        torsions = automol.rotor.from_data(zma, tors_dct)
    else:
        torsions = ()

    # Run the task if any torsions exist
    if any(torsions):

        # Set up ini filesystem for scans
        _, ini_zma_run_path = filesys.build.zma_fs_from_prefix(
            ini_cnf_run_path, zma_idxs=[0])
        _, ini_zma_save_path = filesys.build.zma_fs_from_prefix(
            ini_cnf_save_path, zma_idxs=[0])

        if job == 'scan':

            # Run the scan
            hr.hindered_rotor_scans(
                zma, spc_info, mod_thy_info, ini_thy_save_fs,
                ini_zma_run_path, ini_zma_save_path,
                torsions, tors_model, method_dct,
                overwrite,
                saddle=saddle,
                retryfail=retryfail)

            # Read and print the potential
            # sp_fs = autofile.fs.single_point(ini_cnf_save_path)
            # ref_ene = sp_fs[-1].file.energy.read(mod_ini_thy_info[1:4])
            # ref_ene = ini_cnf_save_fs[-1].file.energy.read(ini_min_cnf_locs)
            # tors_pots, tors_zmas = {}, {}
            # for tors_names, tors_grids in zip(run_tors_names, run_tors_grids):
            #     constraint_dct = structure.tors.build_constraint_dct(
            #         zma, const_names, tors_names)
            #     pot, _, _, _, zmas, _ = structure.tors.read_hr_pot(
            #         tors_names, tors_grids,
            #         ini_cnf_save_path,
            #         mod_ini_thy_info, ref_ene,
            #         constraint_dct,
            #         read_zma=True)
            #     tors_pots[tors_names] = pot
            #     tors_zmas[tors_names] = zmas

            # # Print potential
            # structure.tors.print_hr_pot(tors_pots)

        # elif job == 'reopt':

        #     # pull stuff from dcts
        #     two_stage = saddle
        #     rxn_class = spc_dct_i['class'] if saddle else ''
        #     mc_nsamp = spc_dct_i['mc_nsamp']
        #     ethresh = es_keyword_dct['hrthresh']

        #     # Read and print the potential
        #     sp_fs = autofile.fs.single_point(ini_cnf_save_path)
        #     ref_ene = sp_fs[-1].file.energy.read(mod_ini_thy_info[1:4])
        #     # ref_ene = ini_cnf_save_fs[-1].file.energy.read(ini_min_cnf_locs)
        #     tors_pots, tors_zmas, tors_paths = {}, {}, {}
        #     for tors_names, tors_grids in zip(run_tors_names, run_tors_grids):
        #         constraint_dct = structure.tors.build_constraint_dct(
        #             zma, const_names, tors_names)
        #         pot, _, _, _, zmas, paths = structure.tors.read_hr_pot(
        #             tors_names, tors_grids,
        #             ini_cnf_save_path,
        #             mod_ini_thy_info, ref_ene,
        #             constraint_dct,
        #             read_zma=True)
        #         tors_pots[tors_names] = pot
        #         tors_zmas[tors_names] = zmas
        #         tors_paths[tors_names] = paths

        #     # Check for new minimum conformer
        #     new_min_zma = structure.tors.check_hr_pot(
        #         tors_pots, tors_zmas, tors_paths, emax=ethresh)

        #     if new_min_zma is not None:
        #         ioprinter.info_message(
        #             'Finding new low energy conformer...', newline=1)
        #         conformer.single_conformer(
        #             zma, spc_info, mod_thy_info,
        #             ini_cnf_run_fs, ini_cnf_save_fs,
        #             script_str, overwrite,
        #             retryfail=retryfail, rxn=rxn, **kwargs)

        elif job in ('energy', 'grad', 'hess', 'vpt2'):

            # Script (add energy script call)
            script_str, kwargs = qchem_params(
                method_dct)

            run_tors_names = automol.rotor.names(torsions, flat=True)
            for tors_names in run_tors_names:

                # Set the constraint dct and filesys for the scan
                const_names = structure.tors.set_constraint_names(
                    zma, run_tors_names, tors_model)
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
                        ioprinter.obj('vspace')
                else:
                    ioprinter.warning_message(
                        'NO SCAN INFORMATION EXISTS.',
                        'Doing scan vs cscan?')
    else:
        ioprinter.info_message('No torsional modes in the species')


def irc_tsk(job, spc_dct, spc_name,
            thy_dct, es_keyword_dct,
            run_prefix, save_prefix):
    """ run a scan over the specified torsional coordinates
    """

    # Get dct for specific species task is run for
    spc_dct_i = spc_dct[spc_name]

    # Set up coordinate name
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
    # retryfail = es_keyword_dct['retryfail']

    # Set up the script
    script_str, kwargs = qchem_params(
        method_dct, elstruct.Job.OPTIMIZATION)

    # Set the filesystem objects
    rxn_info = spc_dct_i['rxn_info']
    fs_rxn_info = rinfo.sort(rxn_info)

    _, thy_run_path = filesys.build.rxn_thy_fs_from_root(
        run_prefix, fs_rxn_info, mod_thy_info)
    _, thy_save_path = filesys.build.rxn_thy_fs_from_root(
        save_prefix, fs_rxn_info, mod_thy_info)
    _, thy_save_path = filesys.build.ts_fs_from_thy(thy_save_path)
    _, thy_run_path = filesys.build.ts_fs_from_thy(thy_run_path)

    _, ini_thy_run_path = filesys.build.rxn_thy_fs_from_root(
        run_prefix, fs_rxn_info, mod_ini_thy_info)
    _, ini_thy_save_path = filesys.build.rxn_thy_fs_from_root(
        save_prefix, fs_rxn_info, mod_ini_thy_info)
    # _, ini_ts_save_path = filesys.build.ts_fs_from_thy(ini_thy_save_path)
    # _, ini_ts_run_path = filesys.build.ts_fs_from_thy(ini_thy_run_path)

    # Build cnf filesys using the ini thy filesys (needed for all HR jobs)
    ini_cnf_run_fs = autofile.fs.conformer(ini_thy_run_path)
    ini_cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)

    ini_min_cnf_locs = filesys.mincnf.min_energy_conformer_locators(
        ini_cnf_save_fs, mod_ini_thy_info)

    ini_cnf_save_paths = filesys.build.cnf_paths_from_locs(
        ini_cnf_save_fs, ini_min_cnf_locs)
    ini_cnf_run_paths = filesys.build.cnf_paths_from_locs(
        ini_cnf_run_fs, ini_min_cnf_locs)
    ini_cnf_run_fs[-1].create(ini_min_cnf_locs)

    # Set up ini filesystem for scans
    _, ini_zma_save_path = filesys.build.zma_fs_from_prefix(
        ini_cnf_save_paths[0], zma_idxs=[0])
    ini_scn_save_fs = filesys.build.scn_fs_from_cnf(
        ini_zma_save_path, constraint_dct=None)
    _, ini_zma_run_path = filesys.build.zma_fs_from_prefix(
        ini_cnf_run_paths[0], zma_idxs=[0])
    ini_scn_run_fs = filesys.build.scn_fs_from_cnf(
        ini_zma_run_path, constraint_dct=None)

    # Run job
    if job == 'scan':
        pass
        # irc.scan(geo, spc_info, coord_name,
        #          mod_ini_thy_info, ini_method_dct,
        #          ini_scn_save_fs, ini_cnf_run_path,
        #          overwrite)

    elif job in ('energy', 'grad', 'hess'):

        # Script
        script_str, kwargs = qchem_params(
            method_dct)

        # Need to put in something with the IRC idxs
        scn_locs = filesys.build.scn_locs_from_fs(
            ini_scn_save_fs, coord_name, constraint_dct=None)
        for locs in scn_locs:
            geo_run_path = ini_scn_run_fs[-1].path(locs)
            geo_save_path = ini_scn_save_fs[-1].path(locs)
            geo = ini_scn_save_fs[-1].file.geometry.read(locs)
            # zma = ini_scn_save_fs[-1].file.zmatrix.read(locs)
            zma = None
            # FIX
            # zma, geo = mechanalyzer.inf.cnf_fs_zma_geo(
            #     ini_scn_save_fs, locs)
            ini_scn_run_fs[-1].create(locs)
            ES_TSKS[job](
                zma, geo, spc_info, mod_thy_info,
                ini_scn_save_fs, geo_run_path, geo_save_path, locs,
                script_str, overwrite, **kwargs)
            ioprinter.obj('vspace')
