""" Functions for sadpt
"""

import automol
import autofile
import elstruct
import autorun
from mechanalyzer.inf import rxn as rinfo
from mechanalyzer.inf import thy as tinfo
from mechlib.reaction import grid as rxngrid
from mechlib.amech_io import printer as ioprinter
from mechlib import filesys
from mechroutines.es import runner as es_runner
from mechroutines.es.runner import qchem_params


# SADPT FINDER FUNCTIONS
def generate_guess_structure(ts_dct, method_dct, es_keyword_dct,
                             runfs_dct, savefs_dct, zma_locs=(0,)):
    """ Checks the filesystem for z-matrices at some input
        level of theory and if nothing existsm it will
        launch a scan to find the transition state.

        :param ts_dct: dict of information for the TS
        :param method_dct:
    """

    guess_zmas = _check_filesys_for_guess(savefs_dct, zma_locs, es_keyword_dct)
    if not guess_zmas:
        guess_zmas = scan_for_guess(
            ts_dct, method_dct, runfs_dct, savefs_dct, es_keyword_dct)

    return guess_zmas


def obtain_saddle_point(guess_zmas, ts_dct, method_dct,
                        runfs_dct, savefs_dct,
                        es_keyword_dct):
    """ Given the saddle point guess structure, obtain a
        proper saddle point
    """

    # Get info (move later)
    ts_info = rinfo.ts_info(ts_dct['rxn_info'])
    mod_thy_info = tinfo.modify_orb_label(tinfo.from_dct(method_dct), ts_info)
    script_str, kwargs = qchem_params(
        method_dct, job=elstruct.Job.OPTIMIZATION)

    overwrite = es_keyword_dct['overwrite']
    ts_info = rinfo.ts_info(ts_dct['rxn_info'])
    zrxn = ts_dct['zrxn']

    # runlvl_cnf_run_fs = runfs_dct['runlvl_cnf_fs']
    # cid = [autofile.schema.generate_new_conformer_id()]
    # run_fs = autofile.fs.run(runlvl_cnf_run_fs[-1].path(cid))

    runlvl_cnf_run_fs = runfs_dct['runlvl_cnf_fs']
    rid = autofile.schema.generate_new_ring_id()
    cid = autofile.schema.generate_new_conformer_id()
    locs = (rid, cid)
    run_fs = autofile.fs.run(runlvl_cnf_run_fs[-1].path(locs))

    # Optimize the saddle point
    script_str, kwargs = qchem_params(
        method_dct)
    opt_ret = optimize_saddle_point(
        guess_zmas, ts_info, mod_thy_info,
        run_fs, script_str, overwrite, **kwargs)

    # Calculate the Hessian for the optimized structure
    if opt_ret is not None:
        # Get the Hessian and check the saddle point
        # (maybe just have remove imag do this?)
        hess_ret, _, imags = saddle_point_hessian(
            opt_ret, ts_info, method_dct,
            run_fs, overwrite)

        sadpt_status = saddle_point_checker(imags)

        # Assess saddle point, save it if viable
        # if sadpt_status == 'kickoff':
        #     opt_inf_obj, _, opt_out_str = opt_ret
        #     opt_prog = opt_inf_obj.prog
        #     geo = elstruct.reader.opt_geometry(opt_prog, opt_out_str)

        #     # Need to return an opt_ret
        #     geo, _ = geom.remove_imag(
        #         geo, ts_info, mod_thy_info, ts_run_fs, run_fs,
        #         kickoff_size=0.1, kickoff_backward=False, kickoff_mode=1,
        #         overwrite=False)

        #     sadpt_status = saddle_point_checker(imags)

        if sadpt_status == 'success':
            runlvl_cnf_save_fs, _ = savefs_dct['runlvl_cnf_fs']
            filesys.save.conformer(
                opt_ret, hess_ret, runlvl_cnf_save_fs, mod_thy_info[1:],
                zrxn=zrxn, rng_locs=(rid,), tors_locs=(cid,), zma_locs=None)
    else:
        ioprinter.warning_message(
            '\n TS optimization failed. No geom to check and save.')


def _check_filesys_for_guess(savefs_dct, zma_locs, es_keyword_dct):
    """ Check if the filesystem for any TS structures at the input
        level of theory
    """

    ioprinter.info_message(
        '\nSearching save filesys for guess Z-Matrix calculated',
        'at {} level...'.format(es_keyword_dct['inplvl']))

    ini_zma_fs = savefs_dct['inilvl_zma_fs']

    guess_zmas = []
    if ini_zma_fs is not None:
        if ini_zma_fs[-1].file.zmatrix.exists(zma_locs):
            geo_path = ini_zma_fs[-1].file.zmatrix.exists(zma_locs)
            ioprinter.info_message(' - Z-Matrix found.')
            ioprinter.info_message(
                ' - Reading Z-Matrix from path {}'.format(geo_path))
            guess_zmas.append(
                ini_zma_fs[-1].file.zmatrix.read(zma_locs))

    return guess_zmas


def scan_for_guess(ts_dct, method_dct, runfs_dct, savefs_dct,
                   es_keyword_dct):
    """ saddle point scan code
    """

    print(' - No Z-Matrix is found in save filesys.')
    print('\nRunning scan to generate guess Z-Matrix for opt...')

    # Get es keyword info
    overwrite = es_keyword_dct['overwrite']

    # Get info from the dct
    ts_zma = ts_dct['zma']
    zrxn = ts_dct['zrxn']   # convert to zrxn
    ts_info = rinfo.ts_info(ts_dct['rxn_info'])

    scn_save_fs = savefs_dct['runlvl_scn_fs']

    # Build grid and names appropriate for reaction type
    scan_inf = automol.reac.build_scan_info(zrxn, ts_zma)
    coord_names, constraint_dct, coord_grids, update_guess = scan_inf

    # Get filesystem information
    if constraint_dct is None:
        scn_run_fs = runfs_dct['runlvl_scn_fs']
        scn_save_fs = savefs_dct['runlvl_scn_fs']
    else:
        scn_run_fs = runfs_dct['runlvl_cscn_fs']
        scn_save_fs = savefs_dct['runlvl_cscn_fs']

    # Set up script string and kwargs
    mod_thy_info = tinfo.modify_orb_label(tinfo.from_dct(method_dct), ts_info)
    script_str, kwargs = qchem_params(
        method_dct, job=elstruct.Job.OPTIMIZATION)

    es_runner.scan.execute_scan(
        zma=ts_zma,
        spc_info=ts_info,
        mod_thy_info=mod_thy_info,
        coord_names=coord_names,
        coord_grids=coord_grids,
        scn_run_fs=scn_run_fs,
        scn_save_fs=scn_save_fs,
        scn_typ='relaxed',
        script_str=script_str,
        overwrite=overwrite,
        update_guess=update_guess,
        reverse_sweep=False,
        saddle=False,
        constraint_dct=constraint_dct,
        retryfail=False,
        **kwargs,
        )

    guess_zmas = rxngrid.grid_maximum_zmatrices(
        zrxn.class_, ts_zma, coord_grids, coord_names, scn_save_fs,
        mod_thy_info, constraint_dct)

    return guess_zmas


def optimize_saddle_point(guess_zmas, ts_info, mod_thy_info,
                          run_fs, opt_script_str, overwrite,
                          **opt_kwargs):
    """ Optimize the transition state structure obtained from the grid search
    """

    ioprinter.info_message(
        '\nOptimizing guess Z-Matrix obtained from scan or filesys...')

    if len(guess_zmas) == 1:
        ioprinter.info_message(
            'There is 1 guess Z-Matrix',
            'to attempt to find saddle point.',
            newline=1)
    else:
        ioprinter.info_message(
            'There are {} guess Z-Matrices'.format(len(guess_zmas)),
            'to attempt to find saddle point.', newline=1)

    # Loop over all the guess zmas to find a TS
    opt_ret = None
    for idx, zma in enumerate(guess_zmas):
        ioprinter.info_message(
            '\nOptimizing guess Z-Matrix {}...'.format(idx+1))

        # Run the transition state optimization
        opt_success, opt_ret = es_runner.execute_job(
            job='optimization',
            script_str=opt_script_str,
            run_fs=run_fs,
            geo=zma,
            spc_info=ts_info,
            thy_info=mod_thy_info,
            saddle=True,
            overwrite=overwrite,
            **opt_kwargs,
            )

        # If successful, break loop. If not, try constrained opt+full opt seq
        if opt_success:
            break

        # frozen_coords_lst = ((), tors_names)
        # success, opt_ret = es_runner.multi_stage_optimization(
        #     script_str=opt_script_str,
        #     run_fs=run_fs,
        #     geo=inp_geom,
        #     spc_info=spc_info,
        #     thy_info=thy_info,
        #     frozen_coords_lst=frozen_coords_lst,
        #     overwrite=overwrite,
        #     saddle=saddle,
        #     retryfail=retryfail,
        #     **kwargs
        # )

    return opt_ret


def saddle_point_hessian(opt_ret, ts_info, method_dct,
                         run_fs, overwrite):
    """ run things for checking Hessian
    """

    mod_thy_info = tinfo.modify_orb_label(tinfo.from_dct(method_dct), ts_info)
    script_str, kwargs = qchem_params(
        method_dct)

    # Obtain geometry from optimization
    opt_inf_obj, _, opt_out_str = opt_ret
    opt_prog = opt_inf_obj.prog
    geo = elstruct.reader.opt_geometry(opt_prog, opt_out_str)

    # Run a Hessian
    hess_success, hess_ret = es_runner.execute_job(
        job='hessian',
        script_str=script_str,
        run_fs=run_fs,
        geo=geo,
        spc_info=ts_info,
        thy_info=mod_thy_info,
        overwrite=overwrite,
        **kwargs,
        )

    # If successful, Read the geom and energy from the optimization
    if hess_success:
        hess_inf_obj, _, hess_out_str = hess_ret
        hess = elstruct.reader.hessian(hess_inf_obj.prog, hess_out_str)
        freq_run_path = run_fs[-1].path(['hessian'])
        run_fs[-1].create(['hessian'])
        script_str = autorun.SCRIPT_DCT['projrot']
        freqs, _, imags, _ = autorun.projrot.frequencies(
            script_str, freq_run_path, [geo], [[]], [hess])
    else:
        freqs, imags = [], []

    return hess_ret, freqs, imags


def saddle_point_checker(imags):
    """ run things for checking Hessian
    """

    big_imag, kick_imag = 0, 0

    ioprinter.checking('the imaginary frequencies of the saddle point...')
    if len(imags) < 1:
        ioprinter.warning_message('No imaginary modes for geometry')
        status = 'fail'
    else:
        if len(imags) > 1:
            ioprinter.warning_message(
                'More than one imaginary mode for geometry')
        for idx, imag in enumerate(imags):
            if imag <= 50.0:
                ioprinter.warning_message(
                    'Mode {} {} cm-1 is low,'.format(str(idx+1), imag))
            elif 50.0 < imag <= 200.0:
                lowstr = 'Mode {} {} cm-1 is low,'.format(str(idx+1), imag)
                ioprinter.debug_message(
                    lowstr + ' check mode and see if it should be corrected')
                big_imag += 1
                # Adding to the kick counter kills code for good TSs
                # Some addditions of big species have low mode of this
                # ioprinter.warning_message(
                #     lowstr + 'need a kickoff procedure to remove')
                # kick_imag += 1
            else:
                ioprinter.debug_message(
                    'Mode {} {} cm-1 likely fine,'.format(str(idx+1), imag))
                big_imag += 1

        if big_imag > 1:
            ioprinter.warning_message(
                'More than one imaginary mode for geometry')
            if kick_imag >= 1:
                ioprinter.debug_message('Will kickoff to get saddle point')
                status = 'kickoff'
        elif big_imag == 1:
            status = 'success'

    return status
