""" Functions for sadpt
"""

import automol
import autofile
import elstruct
from mechanalyzer.inf import rxn as rinfo
from mechanalyzer.inf import thy as tinfo
# from mechroutines.es._routines import _geom as geom
from mechroutines.es import runner as es_runner
from mechlib import structure
from mechlib.reaction import grid as rxngrid
from mechlib.submission import qchem_params


# SADPT FINDER FUNCTIONS
def generate_guess_structure(ts_dct, method_dct, es_keyword_dct,
                             runfs_dct, savefs_dct):
    """ Checks the filesystem for z-matrices at some input
        level of theory and if nothing existsm it will
        launch a scan to find the transition state.

        :param ts_dct: dict of information for the TS
        :param method_dct:
    """

    guess_zmas = _check_filesys_for_guess(savefs_dct, (0,), es_keyword_dct)
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

    # Optimize the saddle point
    script_str, kwargs = qchem_params(
        method_dct)
    opt_ret = optimize_saddle_point(
        guess_zmas, ts_info, mod_thy_info,
        runfs_dct, script_str, overwrite, **kwargs)

    # Calculate the Hessian for the optimized structure
    if opt_ret is not None:
        # Get the Hessian and check the saddle point
        # (maybe just have remove imag do this?)
        hess_ret, freqs, imags = saddle_point_hessian(
            opt_ret, ts_info, method_dct,
            runfs_dct, overwrite)

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
            save_saddle_point(
                zrxn, opt_ret, hess_ret, freqs, imags,
                mod_thy_info, savefs_dct,
                zma_locs=[0])

    else:
        print('\n TS optimization failed. No geom to check and save.')


def _check_filesys_for_guess(savefs_dct, zma_locs, es_keyword_dct):
    """ Check if the filesystem for any TS structures at the input
        level of theory
    """

    print('\nSearching save filesys for guess Z-Matrix calculated',
          'at {} level...'.format(es_keyword_dct['inplvl']))

    ini_zma_fs = savefs_dct['inilvl_zma_fs']

    guess_zmas = []
    if ini_zma_fs is not None:
        if ini_zma_fs[-1].file.zmatrix.exists(zma_locs):
            geo_path = ini_zma_fs[-1].file.zmatrix.exists(zma_locs)
            print(' - Z-Matrix found.')
            print(' - Reading Z-Matrix from path {}'.format(geo_path))
            guess_zmas.append(
                ini_zma_fs[-1].file.zmatrix.read(zma_locs))

    return guess_zmas


def scan_for_guess(ts_dct, method_dct, runfs_dct, savefs_dct,
                   es_keyword_dct):
    """ saddle point scan code
    """

    print(' - No Z-Matrix is found in save filesys.')
    print('\nRunning scan to generate guess Z-Matrix for opt...')

    # Get filesystem information
    thy_save_fs = ()
    scn_run_fs = runfs_dct['runlvl_scn_fs']
    scn_save_fs = savefs_dct['runlvl_scn_fs']

    # Get es keyword info
    overwrite = es_keyword_dct['overwrite']

    # Get info from the dct
    ts_zma = ts_dct['zma']
    zrxn = ts_dct['zrxn']   # convert to zrxn
    ts_info = rinfo.ts_info(ts_dct['rxn_info'])

    # Build grid and names appropriate for reaction type
    scan_inf = automol.reac.build_scan_info(zrxn, ts_zma)
    coord_names, constraint_dct, coord_grids, update_guess = scan_inf

    # Set up script string and kwargs
    mod_thy_info = tinfo.modify_orb_label(tinfo.from_dct(method_dct), ts_info)
    script_str, kwargs = qchem_params(
        method_dct, job=elstruct.Job.OPTIMIZATION)

    es_runner.scan.execute_scan(
        zma=ts_zma,
        spc_info=ts_info,
        mod_thy_info=mod_thy_info,
        thy_save_fs=thy_save_fs,
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
        chkstab=False,
        **kwargs,
        )

    # Find the structure at the maximum on the grid opt scan
    # if 'elimination' in rxn_typ:
    #     [grid1, grid2] = coord_grids
    #     max_zma = rxngrid.find_max_2d(
    #         grid1, grid2, frm_name, brk_name, scn_save_fs,
    #         mod_thy_info, constraint_dct)
    #     guess_zmas = [max_zma]
    # else:
    guess_zmas = rxngrid.find_max_1d(
        zrxn.class_, coord_grids[0], ts_zma, coord_names[0], scn_save_fs,
        mod_thy_info, constraint_dct)

    return guess_zmas


def optimize_saddle_point(guess_zmas, ts_info, mod_thy_info,
                          runfs_dct, opt_script_str, overwrite,
                          **opt_kwargs):
    """ Optimize the transition state structure obtained from the grid search
    """

    _, ts_run_path = runfs_dct['runlvl_ts_fs']
    run_fs = autofile.fs.run(ts_run_path)

    print('\nOptimizing guess Z-Matrix obtained from scan or filesys...')

    if len(guess_zmas) == 1:
        print('\nThere is 1 guess Z-Matrix',
              'to attempt to find saddle point.')
    else:
        print('\nThere are {} guess Z-Matrices'.format(len(guess_zmas)),
              'to attempt to find saddle point.')

    # Loop over all the guess zmas to find a TS
    opt_ret = None
    for idx, zma in enumerate(guess_zmas):
        print('\nOptimizing guess Z-Matrix {}...'.format(idx+1))

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
                         runfs_dct, overwrite):
    """ run things for checking Hessian
    """

    _, ts_run_path = runfs_dct['runlvl_ts_fs']
    run_fs = autofile.fs.run(ts_run_path)

    mod_thy_info = tinfo.modify_orb_label(tinfo.from_dct(method_dct), ts_info)
    script_str, kwargs = qchem_params(
        method_dct)

    print('\nCalculating Hessian for the optimized geometry...')

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
        freqs, _, imags, _ = structure.vib.projrot_freqs(
            [geo], [hess], freq_run_path)
    else:
        freqs, imags = [], []

    return hess_ret, freqs, imags


def saddle_point_checker(imags):
    """ run things for checking Hessian
    """

    print('Assessing the saddle point...')

    big_imag, kick_imag = 0, 0

    print('\nChecking the imaginary frequencies of the saddle point...')
    if len(imags) < 1:
        print('WARNING: No imaginary modes for geometry')
        status = 'fail'
    else:
        if len(imags) > 1:
            print('WARNING: More than one imaginary mode for geometry')
        for idx, imag in enumerate(imags):
            if imag <= 50.0:
                print('Mode {} {} cm-1 is low,'.format(str(idx+1), imag))
            elif 50.0 < imag <= 200.0:
                lowstr = 'Mode {} {} cm-1 is low,'.format(str(idx+1), imag)
                print(lowstr + 'need a kickoff procedure to remove')
                kick_imag += 1
            else:
                print('Mode {} {} cm-1 likely fine,'.format(str(idx+1), imag))
                big_imag += 1

        if big_imag > 1:
            print('WARNING: More than one imaginary mode for geometry')
            if kick_imag >= 1:
                print('Will kickoff to get saddle point')
                status = 'kickoff'
        elif big_imag == 1:
            status = 'success'

    return status


def save_saddle_point(zrxn, opt_ret, hess_ret, freqs, imags,
                      mod_thy_info, savefs_dct,
                      zma_locs=(0,)):
    """ Optimize the transition state structure obtained from the grid search
    """
    print('savefs keys', list(savefs_dct.keys()))
    cnf_save_fs, _ = savefs_dct['runlvl_cnf_fs']
    ts_save_fs, ts_save_path = savefs_dct['runlvl_ts_fs']

    print('Saving saddle point...')

    # Read the geom, energy, and Hessian from output
    opt_inf_obj, opt_inp_str, opt_out_str = opt_ret
    opt_prog = opt_inf_obj.prog
    opt_method = opt_inf_obj.method
    ene = elstruct.reader.energy(opt_prog, opt_method, opt_out_str)
    geo = elstruct.reader.opt_geometry(opt_prog, opt_out_str)
    zma = elstruct.reader.opt_zmatrix(opt_prog, opt_out_str)
    print('TS Geometry:')
    print(automol.geom.string(geo))
    print()
    # Build new zma using x2z and new torsion coordinates
    print(" - Reading hessian from output...")
    # zma, tors_names, frm_bnd_keys, brk_bnd_keys = _new_zma(
    #     geo, frm_bnd_keys, brk_bnd_keys)

    print(" - Reading hessian from output...")
    hess_inf_obj, hess_inp_str, hess_out_str = hess_ret
    hess_prog = hess_inf_obj.prog
    hess = elstruct.reader.hessian(hess_prog, hess_out_str)
    freqs = sorted([-1.0*val for val in imags] + freqs)
    print('TS freqs: {}'.format(' '.join(str(freq) for freq in freqs)))

    # Save the information into the filesystem
    print(" - Saving...")
    print(" - Save path: {}".format(ts_save_path))

    # Save geom in the upper theory/TS layer
    ts_save_fs[0].file.geometry.write(geo)

    # Save this structure as first conformer
    locs = [autofile.schema.generate_new_conformer_id()]
    cnf_save_fs[-1].create(locs)
    cnf_save_fs[-1].file.geometry_info.write(opt_inf_obj, locs)
    cnf_save_fs[-1].file.geometry_input.write(opt_inp_str, locs)
    cnf_save_fs[-1].file.hessian_info.write(hess_inf_obj, locs)
    cnf_save_fs[-1].file.hessian_input.write(hess_inp_str, locs)
    cnf_save_fs[-1].file.energy.write(ene, locs)
    cnf_save_fs[-1].file.geometry.write(geo, locs)
    cnf_save_fs[-1].file.hessian.write(hess, locs)
    cnf_save_fs[-1].file.harmonic_frequencies.write(freqs, locs)
    cnf_save_path = cnf_save_fs[-1].path(locs)

    # Save the zmatrix information in a zma filesystem
    cnf_save_path = cnf_save_fs[-1].path(locs)
    zma_save_fs = autofile.fs.zmatrix(cnf_save_path)
    zma_save_fs[-1].create(zma_locs)
    zma_save_fs[-1].file.geometry_info.write(opt_inf_obj, zma_locs)
    zma_save_fs[-1].file.geometry_input.write(opt_inp_str, zma_locs)
    zma_save_fs[-1].file.zmatrix.write(zma, zma_locs)
    zma_save_fs[-1].file.reaction.write(zrxn, zma_locs)

    # Save the torsions
    # geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    # zrxn = automol.reac.insert_dummy_atoms(rxn, gdummy_key_dct)
    rotors = automol.rotor.from_zma(zma)  # Add a graph for the TS zma
    zma_save_fs[-1].file.torsions.write(rotors, zma_locs)

    # Save the energy in a single-point filesystem
    print(" - Saving energy...")
    sp_save_fs = autofile.fs.single_point(cnf_save_path)
    sp_save_fs[-1].create(mod_thy_info[1:4])
    sp_save_fs[-1].file.input.write(opt_inp_str, mod_thy_info[1:4])
    sp_save_fs[-1].file.info.write(opt_inf_obj, mod_thy_info[1:4])
    sp_save_fs[-1].file.energy.write(ene, mod_thy_info[1:4])
