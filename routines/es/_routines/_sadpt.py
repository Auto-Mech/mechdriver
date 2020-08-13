""" Functions for sadpt
"""

import automol
import autofile
from autofile import fs
import elstruct
from routines.es._routines import _scan as scan
from routines.es import runner as es_runner
from lib import structure
from lib import filesys
from lib.reaction import grid as rxngrid


def check_filesys_for_guess(ini_ts_save_path, mod_ini_thy_info, zma_locs=(0,)):
    """ Check if the filesystem for any TS structures at the input
        level of theory
    """
    guess_zmas = []

    # Check and see if a zma is found from the filesystem
    ini_cnf_save_fs, ini_cnf_save_locs = filesys.build.cnf_fs_from_prefix(
        ini_ts_save_path, mod_ini_thy_info, cnf='min')
    if ini_cnf_save_locs:
        ini_zma_fs = autofile.fs.manager(
            ini_cnf_save_fs[-1].path(ini_cnf_save_locs), 'ZMATRIX')
        if ini_zma_fs[-1].file.zmatrix.exists(zma_locs):
            geo_path = ini_zma_fs[-1].file.zmatrix.exists(zma_locs)
            print(' - Z-Matrix found.')
            print(' - Reading Z-Matrix from path {}'.format(geo_path))
            guess_zmas.append(
                ini_zma_fs[-1].file.zmatrix.read(zma_locs))

    return guess_zmas


def check_filesys_for_guess2(ini_zma_fs, zma_locs=(0,)):
    """ Check if the filesystem for any TS structures at the input
        level of theory
    """

    guess_zmas = []
    print('ini_zma_fs test:', ini_zma_fs, zma_locs)
    if ini_zma_fs is not None:
        if ini_zma_fs[-1].file.zmatrix.exists(zma_locs):
            geo_path = ini_zma_fs[-1].file.zmatrix.exists(zma_locs)
            print(' - Z-Matrix found.')
            print(' - Reading Z-Matrix from path {}'.format(geo_path))
            guess_zmas.append(
                ini_zma_fs[-1].file.zmatrix.read(zma_locs))

    return guess_zmas


def scan_for_guess(rxn_typ, grid, dist_name, brk_name,
                   ts_zma, ts_info, mod_thy_info, thy_save_fs,
                   scn_run_fs, scn_save_fs, opt_script_str,
                   overwrite, update_guess, scn_typ='relaxed',
                   **opt_kwargs):
    """ saddle point scan code
    """

    # Build grid and names appropriate for reaction type
    if 'elimination' in rxn_typ:
        coord_grids = grid
        coord_names = [dist_name, brk_name]
    else:
        coord_grids = [grid]
        coord_names = [dist_name]

    scan.run_scan(
        zma=ts_zma,
        spc_info=ts_info,
        mod_thy_info=mod_thy_info,
        thy_save_fs=thy_save_fs,
        coord_names=coord_names,
        coord_grids=coord_grids,
        scn_run_fs=scn_run_fs,
        scn_save_fs=scn_save_fs,
        scn_typ=scn_typ,
        script_str=opt_script_str,
        overwrite=overwrite,
        update_guess=update_guess,
        reverse_sweep=False,
        saddle=False,
        constraint_dct=None,
        retryfail=False,
        chkstab=False,
        **opt_kwargs,
        )
    if 'elimination' in rxn_typ:
        coo_names = [dist_name, brk_name]
    else:
        coo_names = [dist_name]
    scan.save_scan(
        scn_run_fs=scn_run_fs,
        scn_save_fs=scn_save_fs,
        scn_typ=scn_typ,
        coo_names=coo_names,
        mod_thy_info=mod_thy_info,
        in_zma_fs=True
        )

    # Find the structure at the maximum on the grid opt scan
    if 'elimination' in rxn_typ:
        [grid1, grid2] = coord_grids
        max_zma = rxngrid.find_max_2d(
            grid1, grid2, dist_name, brk_name, scn_save_fs, mod_thy_info)
        guess_zmas = [max_zma]
    else:
        guess_zmas = rxngrid.find_max_1d(
            rxn_typ, grid, ts_zma, dist_name, scn_save_fs, mod_thy_info)

    return guess_zmas


def optimize_saddle_point(guess_zmas, ts_info, mod_thy_info,
                          run_fs, opt_script_str, overwrite,
                          **opt_kwargs):
    """ Optimize the transition state structure obtained from the grid search
    """

    if len(guess_zmas) == 1:
        print('\nThere is 1 guess Z-Matrix',
              'to attempt to find saddle point.')
    else:
        print('\nThere are {} guess Z-Matrices'.format(len(guess_zmas)),
              'to attempt to find saddle point.')

    # Loop over all the guess zmas to find a TS
    opt_ret = None
    for idx, zma in enumerate(guess_zmas):
        print('\nAttempting optimization of',
              'guess Z-Matrix {}...'.format(idx+1))

        # Run the transition state optimization
        es_runner.run_job(
            job='optimization',
            script_str=opt_script_str,
            run_fs=run_fs,
            geom=zma,
            spc_info=ts_info,
            thy_info=mod_thy_info,
            saddle=True,
            overwrite=overwrite,
            **opt_kwargs,
            )

        # Read the contents of the optimization
        opt_success, opt_ret = es_runner.read_job(
            job='optimization',
            run_fs=run_fs,
        )

        if opt_success:
            break

    return opt_ret


def saddle_point_hessian(opt_ret, ts_info, mod_thy_info,
                         run_fs, script_str, overwrite,
                         **opt_kwargs):
    """ run things for checking Hessian
    """

    # Obtain geometry from optimization
    opt_inf_obj, _, opt_out_str = opt_ret
    opt_prog = opt_inf_obj.prog
    geo = elstruct.reader.opt_geometry(opt_prog, opt_out_str)

    # Run a Hessian
    es_runner.run_job(
        job='hessian',
        script_str=script_str,
        run_fs=run_fs,
        geom=geo,
        spc_info=ts_info,
        thy_info=mod_thy_info,
        overwrite=overwrite,
        **opt_kwargs,
        )

    # Read the contents of the optimization
    hess_success, hess_ret = es_runner.read_job(
        job='hessian',
        run_fs=run_fs,
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

    saddle = True
    if len(imags) < 1:
        print('Warning: No imaginary modes for geometry')
        saddle = False
    elif len(imags) > 1:
        print('Warning: More than one imaginary mode for geometry')
        nbig_imag = 0
        for imag in imags:
            if imag >= 40.0:
                nbig_imag += 1
        if nbig_imag > 1:
            saddle = False

    for imag in imags:
        if imag <= 200.0:
            print('Imaginary mode {} is relatively low. Worth checking'.format(imag))

    return saddle


def save_saddle_point(
        opt_ret, hess_ret, freqs, imags,
        mod_thy_info,
        cnf_save_fs,
        ts_save_fs, ts_save_path,
        frm_bnd_keys, brk_bnd_keys, rcts_gra,
        zma_locs=(0,)):
    """ Optimize the transition state structure obtained from the grid search
    """

    # Read the geom, energy, and Hessian from output
    opt_inf_obj, opt_inp_str, opt_out_str = opt_ret
    opt_prog = opt_inf_obj.prog
    opt_method = opt_inf_obj.method
    ene = elstruct.reader.energy(opt_prog, opt_method, opt_out_str)
    geo = elstruct.reader.opt_geometry(opt_prog, opt_out_str)
    zma = elstruct.reader.opt_zmatrix(opt_prog, opt_out_str)

    # Build new zma using x2z and new torsion coordinates
    # zma = automol.geometry.zmatrix(
    #     geo, ts_bnd=(frm_bnd_keys, brk_bnd_keys))
    
    print(" - Reading hessian from output...")
    hess_inf_obj, hess_inp_str, hess_out_str = hess_ret
    hess_prog = hess_inf_obj.prog
    hess = elstruct.reader.hessian(hess_prog, hess_out_str)
    freqs = sorted([-1.0*val for val in imags] + freqs)

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
    zma_save_fs = fs.manager(cnf_save_path, 'ZMATRIX')
    zma_save_fs[-1].create(zma_locs)
    zma_save_fs[-1].file.geometry_info.write(opt_inf_obj, zma_locs)
    zma_save_fs[-1].file.geometry_input.write(opt_inp_str, zma_locs)
    zma_save_fs[-1].file.zmatrix.write(zma, zma_locs)

    # Save the form and break keys in the filesystem
    # shift_frm_bnd_keys = structure.geom.shift_vals_from_dummy(
    #     frm_bnd_keys, zma)
    # shift_brk_bnd_keys = structure.geom.shift_vals_from_dummy(
    #     brk_bnd_keys, zma)

    tra = (frozenset({frm_bnd_keys}),
           frozenset({brk_bnd_keys}))
    zma_save_fs[-1].file.transformation.write(tra, zma_locs)
    zma_save_fs[-1].file.reactant_graph.write(rcts_gra, zma_locs)

    # Save the energy in a single-point filesystem
    print(" - Saving energy...")
    sp_save_fs = autofile.fs.single_point(cnf_save_path)
    sp_save_fs[-1].create(mod_thy_info[1:4])
    sp_save_fs[-1].file.input.write(opt_inp_str, mod_thy_info[1:4])
    sp_save_fs[-1].file.info.write(opt_inf_obj, mod_thy_info[1:4])
    sp_save_fs[-1].file.energy.write(ene, mod_thy_info[1:4])


#####
# def check_filesys_for_guess(ini_ts_save_path, mod_ini_thy_info, zma_locs=(0,)):
#     """ Check if the filesystem for any TS structures at the input
#         level of theory
#     """
#     guess_zmas = []
#
#     # Check and see if a zma is found from the filesystem
#     ini_cnf_save_fs, ini_cnf_save_locs = filesys.build.cnf_fs_from_prefix(
#         ini_ts_save_path, mod_ini_thy_info, cnf='min')
#     if ini_cnf_save_locs:
#         ini_zma_fs = autofile.fs.manager(
#             ini_cnf_save_fs[-1].path(ini_cnf_save_locs), 'ZMATRIX')
#         if ini_zma_fs[-1].file.zmatrix.exists(zma_locs):
#             geo_path = ini_zma_fs[-1].file.zmatrix.exists(zma_locs)
#             print(' - Z-Matrix found.')
#             print(' - Reading Z-Matrix from path {}'.format(geo_path))
#             guess_zmas.append(
#                 ini_zma_fs[-1].file.zmatrix.read(zma_locs))
#
#     return guess_zmas

