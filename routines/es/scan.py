""" drivers for coordinate scans
"""
import numpy
import automol
import elstruct
from elstruct.reader._molpro2015.molecule import hess_geometry
import autofile

# New libs
from lib.phydat import phycon
from lib.runner import driver
from lib.runner import par as runpar
from lib.filesystem import minc as fsmin
from lib.filesystem import orb as fsorb
from routines.es import variational


def hr_prep(zma, geo, run_tors_names=(), scan_increment=30.0, ndim_tors='1dhr',
            saddle=False, frm_bnd_key=(), brk_bnd_key=()):
    """ set-up the hr for different rotor combinations
        tors_names = [ ['D1'], ['D2', 'D3'], ['D4'] ]
    """
    # Convert scan_increment to radians (was broken before)
    print('INCREMENT:')
    print(scan_increment)
    print(phycon.DEG2RAD)
    # scan_increment *= phycon.DEG2RAD

    # Get the tors names if thery have not already been supplied
    val_dct = automol.zmatrix.values(zma)
    if not run_tors_names:
        if not saddle:
            run_tors_names = [
                [name]
                for name in automol.geom.zmatrix_torsion_coordinate_names(geo)
            ]
            if ndim_tors == 'mdhr':
                run_tors_names = [[tors
                                   for rotor in run_tors_names
                                   for tors in rotor]]

    # Deal with the dimensionality of the rotors
    if ndim_tors == 'mdhr':
        run_tors_names = mdhr_prep(zma, run_tors_names)

    # Build the grids corresponding to the torsions
    run_tors_grids = []
    print('build linespace loop')
    for tors_names in run_tors_names:
        print('tors names', tors_names)
        tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
            zma, tors_names, scan_increment, frm_bnd_key=frm_bnd_key,
            brk_bnd_key=brk_bnd_key)
        print(tors_linspaces)
        run_tors_grids.append(
            [numpy.linspace(*linspace) + val_dct[name]
             for name, linspace in zip(tors_names, tors_linspaces)]
        )

    return run_tors_names, run_tors_grids


def mdhr_prep(zma, run_tors_names):
    """ Handle cases where the MDHR
    """

    # Figure out set of torsions are to be used: defined or AMech generated
    rotor_lst = run_tors_names
    print('rotor_lst', rotor_lst)

    # Check the dimensionality of each rotor to see if they are greater than 4
    # Call a function to reduce large rotors
    final_rotor_lst = []
    for rotor in rotor_lst:
        print('mdhr prep rotor', rotor)
        if len(rotor) > 4:
            print('LEN BAD')
            for reduced_rotor in reduce_rotor_dimensionality(zma, rotor):
                final_rotor_lst.append(reduced_rotor)
        else:
            print('LEN GOOD')
            final_rotor_lst.append(rotor)

    return final_rotor_lst


def reduce_rotor_dimensionality(zma, rotor):
    """ For rotors with a dimensionality greater than 4, try and take them out
    """

    # Find the methyl rotors for that are a part of the MDHR
    reduced_rotor_lst = []
    methyl_rotors = []
    for tors in rotor:
        # If a methyl rotor add to methyl rotor list
        if is_methyl_rotor():   # Add arguments when ID methyls
            methyl_rotors.append(tors)
        # Add to reduced rotor list
        else:
            reduced_rotor_lst.append(tors)

    # Add each of methyl rotors, if any exist
    if methyl_rotors:
        for methyl_rotor in methyl_rotors:
            reduced_rotor_lst.append(methyl_rotor)

    # Check new dimensionality of list; if still high, flatten to lst of 1DHRs
    if len(reduced_rotor_lst) > 4:
        reduced_rotor_lst = [tors
                             for rotor in reduced_rotor_lst
                             for tors in rotor]

    return reduced_rotor_lst


def is_methyl_rotor():
    """ Check if methyl rotor
    """
    return False


def hindered_rotor_scans(
        spc_info, thy_level, cnf_run_fs, cnf_save_fs, script_str, overwrite,
        scan_increment=30.0, saddle=False, run_tors_names=(), frm_bnd_key=(),
        brk_bnd_key=(), tors_model=('1dhr', False), **opt_kwargs):
    """ Perform 1d scans over each of the torsional coordinates
    """

    # Unpack tors model
    tors_model = ('mdhr', False)
    ndim_tors, freeze_all_tors = tors_model

    # Run with the old code
    min_cnf_locs = fsmin.min_energy_conformer_locators(cnf_save_fs)
    if min_cnf_locs:
        min_cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
        min_cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
        scn_run_fs = autofile.fs.scan(min_cnf_run_path)
        scn_save_fs = autofile.fs.scan(min_cnf_save_path)
        geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        zma = cnf_save_fs.leaf.file.zmatrix.read(min_cnf_locs)

        run_tors_names, run_tors_grids = hr_prep(
            zma, geo, run_tors_names=(),
            scan_increment=scan_increment, ndim_tors=ndim_tors,
            saddle=saddle, frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key)

        print('TEST: run_tors_names')
        print(run_tors_names)
        print(run_tors_grids)

        print('TEST: in loop')
        # for tors_name, tors_grid in zip(tors_names, tors_grids):
        for tors_names, tors_grids in zip(run_tors_names, run_tors_grids):

            # Get the dictionary for the torsional modes
            print(tors_names)
            print(tors_grids)
            if not tors_names:
                continue
            grid_dct = dict(zip(tors_names, tors_grids))

            # Get a list of the other tors coords to freeze
            print('freeze variable', freeze_all_tors)
            if freeze_all_tors:
                alt_constraints = [name
                                   for name_lst in run_tors_names
                                   for name in name_lst]
            else:
                alt_constraints = ()
            print('alt_constraints')
            print(alt_constraints)

            # Perform the scans
            save_scan(
                scn_run_fs=scn_run_fs,
                scn_save_fs=scn_save_fs,
                coo_names=tors_names,
            )

            run_scan(
                zma=zma,
                spc_info=spc_info,
                thy_level=thy_level,
                grid_dct=grid_dct,
                scn_run_fs=scn_run_fs,
                scn_save_fs=scn_save_fs,
                script_str=script_str,
                overwrite=overwrite,
                saddle=saddle,
                alt_constraints=alt_constraints,
                **opt_kwargs,
            )

            save_scan(
                scn_run_fs=scn_run_fs,
                scn_save_fs=scn_save_fs,
                coo_names=tors_names,
            )


def run_scan(
        zma, spc_info, thy_level, grid_dct, scn_run_fs, scn_save_fs,
        script_str, overwrite, update_guess=True,
        reverse_sweep=True, fix_failures=True, saddle=False,
        alt_constraints=(),
        **kwargs):
    """ run constrained optimization scan
    """

    vma = automol.zmatrix.var_(zma)
    if scn_save_fs.trunk.file.vmatrix.exists():
        existing_vma = scn_save_fs.trunk.file.vmatrix.read()
        assert vma == existing_vma

    coo_names = []
    grid_vals = []
    for item in grid_dct.items():
        (coo, coo_grid_vals) = item
        coo_names.append(coo)
        grid_vals.append(coo_grid_vals)

    print('grid_dct test:', grid_dct)
    print('grid_vals test:', grid_vals)
    # import sys
    # sys.exit()

    # for now, running only one-dimensional hindered rotor scans
    scn_save_fs.branch.create([coo_names])
    inf_obj = autofile.system.info.scan_branch(grid_dct)
    scn_save_fs.branch.file.info.write(inf_obj, [coo_names])
    npoint = 1
    for coo_grid_vals in grid_vals:
        npoint *= len(coo_grid_vals)
    grid_idxs = tuple(range(npoint))
    if len(grid_vals) == 1:
        for grid_val in grid_vals[0]:
            scn_run_fs.leaf.create([coo_names, [grid_val]])
        run_prefixes = tuple(scn_run_fs.leaf.path([coo_names, [grid_val]])
                             for grid_val in grid_vals[0])
        _run_1d_scan(
            script_str=script_str,
            run_prefixes=run_prefixes,
            scn_save_fs=scn_save_fs,
            guess_zma=zma,
            coo_name=coo_names[0],
            grid_idxs=grid_idxs,
            grid_vals=grid_vals[0],
            spc_info=spc_info,
            thy_level=thy_level,
            overwrite=overwrite,
            update_guess=update_guess,
            saddle=saddle,
            retry_failed=fix_failures,
            alt_constraints=alt_constraints,
            **kwargs
        )

        if reverse_sweep:
            _run_1d_scan(
                script_str=script_str,
                run_prefixes=list(reversed(run_prefixes)),
                scn_save_fs=scn_save_fs,
                guess_zma=zma,
                coo_name=coo_names[0],
                grid_idxs=list(reversed(grid_idxs)),
                grid_vals=list(reversed(grid_vals[0])),
                spc_info=spc_info,
                thy_level=thy_level,
                overwrite=overwrite,
                update_guess=update_guess,
                saddle=saddle,
                alt_constraints=alt_constraints,
                **kwargs
            )

    elif len(grid_vals) == 2:
        run_prefixes = []
        for grid_val_i in grid_vals[0]:
            for grid_val_j in grid_vals[1]:
                scn_run_fs.leaf.create([coo_names, [grid_val_i, grid_val_j]])
                run_prefixes.append(scn_run_fs.leaf.path(
                    [coo_names, [grid_val_i, grid_val_j]]))
        run_prefixes = tuple(run_prefixes)

        _run_2d_scan(
            script_str=script_str,
            run_prefixes=run_prefixes,
            scn_save_fs=scn_save_fs,
            guess_zma=zma,
            coo_names=coo_names,
            grid_idxs=grid_idxs,
            grid_vals=grid_vals,
            spc_info=spc_info,
            thy_level=thy_level,
            overwrite=overwrite,
            update_guess=update_guess,
            saddle=saddle,
            retry_failed=fix_failures,
            alt_constraints=alt_constraints,
            **kwargs
        )

        if reverse_sweep:
            run_prefixes = []
            for grid_val_i in grid_vals[0][::-1]:
                for grid_val_j in grid_vals[1][::-1]:
                    run_prefixes.append(scn_run_fs.leaf.path(
                        [coo_names, [grid_val_i, grid_val_j]]))
            run_prefixes = tuple(run_prefixes)
            _run_2d_scan(
                script_str=script_str,
                run_prefixes=run_prefixes,
                scn_save_fs=scn_save_fs,
                guess_zma=zma,
                coo_names=coo_names,
                grid_idxs=list(reversed(grid_idxs)),
                grid_vals=[list(reversed(grid_vals[0])),
                           list(reversed(grid_vals[1]))],
                spc_info=spc_info,
                thy_level=thy_level,
                overwrite=overwrite,
                update_guess=update_guess,
                saddle=saddle,
                alt_constraints=alt_constraints,
                **kwargs
            )

    elif len(grid_vals) == 3:
        run_prefixes = []
        for grid_val_i in grid_vals[0]:
            for grid_val_j in grid_vals[1]:
                for grid_val_k in grid_vals[2]:
                    scn_run_fs.leaf.create(
                        [coo_names, [grid_val_i, grid_val_j, grid_val_k]])
                    run_prefixes.append(scn_run_fs.leaf.path(
                        [coo_names, [grid_val_i, grid_val_j, grid_val_k]]))
        run_prefixes = tuple(run_prefixes)

        _run_3d_scan(
            script_str=script_str,
            run_prefixes=run_prefixes,
            scn_save_fs=scn_save_fs,
            guess_zma=zma,
            coo_names=coo_names,
            grid_idxs=grid_idxs,
            grid_vals=grid_vals,
            spc_info=spc_info,
            thy_level=thy_level,
            overwrite=overwrite,
            update_guess=update_guess,
            saddle=saddle,
            retry_failed=fix_failures,
            alt_constraints=alt_constraints,
            **kwargs
        )

        if reverse_sweep:
            run_prefixes = []
            for grid_val_i in grid_vals[0][::-1]:
                for grid_val_j in grid_vals[1][::-1]:
                    for grid_val_k in grid_vals[2][::-1]:
                        run_prefixes.append(scn_run_fs.leaf.path(
                            [coo_names, [grid_val_i, grid_val_j, grid_val_k]]))
            run_prefixes = tuple(run_prefixes)
            _run_3d_scan(
                script_str=script_str,
                run_prefixes=run_prefixes,
                scn_save_fs=scn_save_fs,
                guess_zma=zma,
                coo_names=coo_names,
                grid_idxs=list(reversed(grid_idxs)),
                grid_vals=[list(reversed(grid_vals[0])),
                           list(reversed(grid_vals[1])),
                           list(reversed(grid_vals[2]))],
                spc_info=spc_info,
                thy_level=thy_level,
                overwrite=overwrite,
                update_guess=update_guess,
                saddle=saddle,
                alt_constraints=alt_constraints,
                **kwargs
            )

    elif len(grid_vals) == 4:
        run_prefixes = []
        for grid_val_i in grid_vals[0]:
            for grid_val_j in grid_vals[1]:
                for grid_val_k in grid_vals[2]:
                    for grid_val_l in grid_vals[3]:
                        scn_run_fs.leaf.create(
                            [coo_names,
                             [grid_val_i, grid_val_j,
                              grid_val_k, grid_val_l]])
                        run_prefixes.append(scn_run_fs.leaf.path(
                            [coo_names,
                             [grid_val_i, grid_val_j,
                              grid_val_k, grid_val_l]]))
        run_prefixes = tuple(run_prefixes)

        _run_4d_scan(
            script_str=script_str,
            run_prefixes=run_prefixes,
            scn_save_fs=scn_save_fs,
            guess_zma=zma,
            coo_names=coo_names,
            grid_idxs=grid_idxs,
            grid_vals=grid_vals,
            spc_info=spc_info,
            thy_level=thy_level,
            overwrite=overwrite,
            update_guess=update_guess,
            saddle=saddle,
            retry_failed=fix_failures,
            alt_constraints=alt_constraints,
            **kwargs
        )

        if reverse_sweep:
            run_prefixes = []
            for grid_val_i in grid_vals[0][::-1]:
                for grid_val_j in grid_vals[1][::-1]:
                    for grid_val_k in grid_vals[2][::-1]:
                        for grid_val_l in grid_vals[3][::-1]:
                            run_prefixes.append(scn_run_fs.leaf.path(
                                [coo_names,
                                 [grid_val_i, grid_val_j,
                                  grid_val_k, grid_val_l]]))
            run_prefixes = tuple(run_prefixes)
            _run_4d_scan(
                script_str=script_str,
                run_prefixes=run_prefixes,
                scn_save_fs=scn_save_fs,
                guess_zma=zma,
                coo_names=coo_names,
                grid_idxs=list(reversed(grid_idxs)),
                grid_vals=[list(reversed(grid_vals[0])),
                           list(reversed(grid_vals[1])),
                           list(reversed(grid_vals[2])),
                           list(reversed(grid_vals[3]))],
                spc_info=spc_info,
                thy_level=thy_level,
                overwrite=overwrite,
                update_guess=update_guess,
                saddle=saddle,
                alt_constraints=alt_constraints,
                **kwargs
            )


def run_multiref_rscan(
        formula, high_mul, zma, spc_info, multi_level, dist_name, grid1, grid2,
        scn_run_fs, scn_save_fs, overwrite, update_guess=True,
        gradient=False, hessian=False, num_act_elc=None, num_act_orb=None,
        alt_constraints=()):
    """ run constrained optimization scan
    """

    vma = automol.zmatrix.var_(zma)
    if scn_save_fs.trunk.file.vmatrix.exists():
        existing_vma = scn_save_fs.trunk.file.vmatrix.read()
        assert vma == existing_vma

    grid = numpy.append(grid1, grid2)
    grid_dct = {dist_name: grid}
    if len(grid_dct) > 1:
        raise NotImplementedError

    coo_names = []
    grid_vals = []
    for item in grid_dct.items():
        (coo, coo_grid_vals) = item
        coo_names.append(coo)
        grid_vals.append(coo_grid_vals)

    scn_save_fs.branch.create([coo_names])
    inf_obj = autofile.system.info.scan_branch(grid_dct)
    scn_save_fs.branch.file.info.write(inf_obj, [coo_names])

    prog = multi_level[0]
    method = multi_level[1]

    _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(prog, method)

    ref_zma = automol.zmatrix.set_values(zma, {coo_names[0]: grid_vals[0][0]})

    # Build the elstruct CASSCF options list for multiref calcs
    cas_opt = []
    cas_opt.append(
        variational.wfn.cas_options(
            spc_info, formula, num_act_elc, num_act_orb,
            high_mul, add_two_closed=False))
    cas_opt.append(
        variational.wfn.cas_options(
            spc_info, formula, num_act_elc, num_act_orb,
            high_mul, add_two_closed=True))

    # Write the lines containing all the calcs for a guess wfn
    guess_str = variational.wfn.multiref_wavefunction_guess(
        high_mul, ref_zma, spc_info, multi_level, cas_opt)
    guess_lines = guess_str.splitlines()

    # Add the above-built objects to the elstruct opt_kwargs dct
    opt_kwargs['casscf_options'] = cas_opt[1]
    opt_kwargs['gen_lines'] = {1: guess_lines}

    # Add option to opt_kwargs dct to turn off symmetry
    opt_kwargs['mol_options'] = ['nosym']

    # Setup and run the first part of the scan to shorter distances
    coo_names = []
    grid1_vals = []
    grid1_dct = {dist_name: grid1}
    for item in grid1_dct.items():
        (coo, coo_grid1_vals) = item
        coo_names.append(coo)
        grid1_vals.append(coo_grid1_vals)

    npoint = 1
    for coo_grid1_vals in grid1_vals:
        npoint *= len(coo_grid1_vals)
    grid1_idxs = tuple(range(npoint))
    if len(grid1_vals) == 1:
        for grid1_val in grid1_vals[0]:
            scn_run_fs.leaf.create([coo_names, [grid1_val]])
        run_prefixes = tuple(scn_run_fs.leaf.path([coo_names, [grid1_val]])
                             for grid1_val in grid1_vals[0])

    _run_1d_scan(
        script_str=opt_script_str,
        run_prefixes=run_prefixes,
        scn_save_fs=scn_save_fs,
        guess_zma=zma,
        coo_name=coo_names[0],
        grid_idxs=grid1_idxs,
        grid_vals=grid1_vals[0],
        spc_info=spc_info,
        thy_level=multi_level,
        overwrite=overwrite,
        update_guess=update_guess,
        gradient=gradient,
        hessian=hessian,
        alt_constraints=alt_constraints,
        **opt_kwargs,
    )

    # Setup and run the second part of the scan to longer distances
    coo_names = []
    grid2_vals = []
    grid2_dct = {dist_name: grid2}
    for item in grid2_dct.items():
        (coo, coo_grid2_vals) = item
        coo_names.append(coo)
        grid2_vals.append(coo_grid2_vals)

    npoint = 1
    for coo_grid2_vals in grid2_vals:
        npoint *= len(coo_grid2_vals)
    grid2_idxs = tuple(range(npoint))
    if len(grid2_vals) == 1:
        for grid2_val in grid2_vals[0]:
            scn_run_fs.leaf.create([coo_names, [grid2_val]])
        run_prefixes = tuple(scn_run_fs.leaf.path([coo_names, [grid2_val]])
                             for grid2_val in grid2_vals[0])
    _run_1d_scan(
        script_str=opt_script_str,
        run_prefixes=run_prefixes,
        scn_save_fs=scn_save_fs,
        guess_zma=zma,
        coo_name=coo_names[0],
        grid_idxs=grid2_idxs,
        grid_vals=grid2_vals[0],
        spc_info=spc_info,
        thy_level=multi_level,
        overwrite=overwrite,
        update_guess=update_guess,
        gradient=gradient,
        hessian=hessian,
        alt_constraints=alt_constraints,
        **opt_kwargs,
    )


def _run_1d_scan(
        script_str, run_prefixes, scn_save_fs, guess_zma, coo_name,
        grid_idxs, grid_vals,
        spc_info, thy_level, overwrite, errors=(), options_mat=(),
        retry_failed=True, update_guess=True, saddle=False,
        gradient=False, hessian=False,
        alt_constraints=(),
        **kwargs):
    """ run 1 dimensional scan with constrained optimization
    """

    npoints = len(grid_idxs)
    assert len(grid_vals) == len(run_prefixes) == npoints
    grid_info = zip(grid_idxs, grid_vals, run_prefixes)
    for grid_idx, grid_val, run_prefix in grid_info:
        print("Point {}/{}".format(grid_idx+1, npoints))
        zma = automol.zmatrix.set_values(guess_zma, {coo_name: grid_val})
        run_fs = autofile.fs.run(run_prefix)

        frozen_coordinates = [coo_name] + list(alt_constraints)
        geo_exists = scn_save_fs.leaf.file.geometry.exists(
            [[coo_name], [grid_val]])
        if not geo_exists or overwrite:
            driver.run_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                run_fs=run_fs,
                geom=zma,
                spc_info=spc_info,
                thy_level=thy_level,
                overwrite=overwrite,
                frozen_coordinates=frozen_coordinates,
                errors=errors,
                options_mat=options_mat,
                retry_failed=retry_failed,
                saddle=saddle,
                **kwargs
            )

            ret = driver.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if ret is not None:
                inf_obj, _, out_str = ret
                prog = inf_obj.prog
                opt_zma = elstruct.reader.opt_zmatrix(prog, out_str)
                if update_guess:
                    guess_zma = opt_zma

                if gradient:
                    driver.run_job(
                        job=elstruct.Job.GRADIENT,
                        script_str=script_str,
                        run_fs=run_fs,
                        geom=opt_zma,
                        spc_info=spc_info,
                        thy_level=thy_level,
                        overwrite=overwrite,
                        frozen_coordinates=[coo_name],
                        errors=errors,
                        options_mat=options_mat,
                        retry_failed=retry_failed,
                        **kwargs
                    )

                    ret = driver.read_job(
                        job=elstruct.Job.GRADIENT, run_fs=run_fs)

                if hessian:
                    driver.run_job(
                        job=elstruct.Job.HESSIAN,
                        script_str=script_str,
                        run_fs=run_fs,
                        geom=opt_zma,
                        spc_info=spc_info,
                        thy_level=thy_level,
                        overwrite=overwrite,
                        frozen_coordinates=[coo_name],
                        errors=errors,
                        options_mat=options_mat,
                        retry_failed=retry_failed,
                        **kwargs
                    )

                    ret = driver.read_job(
                        job=elstruct.Job.HESSIAN, run_fs=run_fs)


def _run_2d_scan(
        script_str, run_prefixes, scn_save_fs, guess_zma, coo_names,
        grid_idxs, grid_vals,
        spc_info, thy_level, overwrite, errors=(),
        options_mat=(), retry_failed=True, update_guess=True, saddle=False,
        alt_constraints=(),
        **kwargs):
    """ run 2-dimensional scan with constrained optimization
     """

    npoints = len(grid_idxs)
    assert len(grid_vals[0])*len(grid_vals[1]) == len(run_prefixes) == npoints

    idx = 0
    for grid_val_i in grid_vals[0]:
        for grid_val_j in grid_vals[1]:
            grid_idx = grid_idxs[idx]
            run_prefix = run_prefixes[idx]
            print("Point {}/{}".format(grid_idx+1, npoints))
            zma = automol.zmatrix.set_values(
                guess_zma, {coo_names[0]: grid_val_i,
                            coo_names[1]: grid_val_j})
            run_fs = autofile.fs.run(run_prefix)
            idx += 1

            frozen_coordinates = coo_names + list(alt_constraints)
            if not scn_save_fs.leaf.file.geometry.exists(
                    [coo_names, [grid_val_i, grid_val_j]]) or overwrite:
                driver.run_job(
                    job=elstruct.Job.OPTIMIZATION,
                    script_str=script_str,
                    run_fs=run_fs,
                    geom=zma,
                    spc_info=spc_info,
                    thy_level=thy_level,
                    overwrite=overwrite,
                    frozen_coordinates=frozen_coordinates,
                    errors=errors,
                    options_mat=options_mat,
                    retry_failed=retry_failed,
                    saddle=saddle,
                    **kwargs
                )

                ret = driver.read_job(
                    job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
                if update_guess and ret is not None:
                    inf_obj, _, out_str = ret
                    prog = inf_obj.prog
                    guess_zma = elstruct.reader.opt_zmatrix(prog, out_str)


def _run_3d_scan(
        script_str, run_prefixes, scn_save_fs, guess_zma, coo_names,
        grid_idxs, grid_vals,
        spc_info, thy_level, overwrite, errors=(),
        options_mat=(), retry_failed=True, update_guess=True, saddle=False,
        alt_constraints=(),
        **kwargs):
    """ run 2-dimensional scan with constrained optimization
    """

    npoints = len(grid_idxs)
    assert len(grid_vals[0])*len(grid_vals[1])*len(grid_vals[2]) == len(run_prefixes)
    assert len(run_prefixes)== npoints

    idx = 0
    for grid_val_i in grid_vals[0]:
        for grid_val_j in grid_vals[1]:
            for grid_val_k in grid_vals[2]:
                grid_idx = grid_idxs[idx]
                run_prefix = run_prefixes[idx]
                print("Point {}/{}".format(grid_idx+1, npoints))
                zma = automol.zmatrix.set_values(
                    guess_zma, {coo_names[0]: grid_val_i,
                                coo_names[1]: grid_val_j,
                                coo_names[2]: grid_val_k})
                run_fs = autofile.fs.run(run_prefix)
                idx += 1

                frozen_coordinates = coo_names + list(alt_constraints)
                print('3d freeze test', coo_names, list(alt_constraints))
                exists = scn_save_fs.leaf.file.geometry.exists(
                    [coo_names, [grid_val_i, grid_val_j, grid_val_k]])
                if not exists or overwrite:
                    driver.run_job(
                        job=elstruct.Job.OPTIMIZATION,
                        script_str=script_str,
                        run_fs=run_fs,
                        geom=zma,
                        spc_info=spc_info,
                        thy_level=thy_level,
                        overwrite=overwrite,
                        frozen_coordinates=frozen_coordinates,
                        errors=errors,
                        options_mat=options_mat,
                        retry_failed=retry_failed,
                        saddle=saddle,
                        **kwargs
                    )

                    ret = driver.read_job(
                        job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
                    if update_guess and ret is not None:
                        inf_obj, _, out_str = ret
                        prog = inf_obj.prog
                        guess_zma = elstruct.reader.opt_zmatrix(prog, out_str)


def _run_4d_scan(
        script_str, run_prefixes, scn_save_fs, guess_zma, coo_names,
        grid_idxs, grid_vals,
        spc_info, thy_level, overwrite, errors=(),
        options_mat=(), retry_failed=True, update_guess=True, saddle=False,
        alt_constraints=(),
        **kwargs):
    """ run 2-dimensional scan with constrained optimization
    """

    npoints = len(grid_idxs)
    gmult = (len(grid_vals[0])*len(grid_vals[1]) *
             len(grid_vals[2])*len(grid_vals[3]))
    assert gmult == len(run_prefixes) == npoints

    idx = 0
    for grid_val_i in grid_vals[0]:
        for grid_val_j in grid_vals[1]:
            for grid_val_k in grid_vals[2]:
                for grid_val_l in grid_vals[3]:
                    grid_idx = grid_idxs[idx]
                    run_prefix = run_prefixes[idx]
                    print("Point {}/{}".format(grid_idx+1, npoints))
                    zma = automol.zmatrix.set_values(
                        guess_zma, {coo_names[0]: grid_val_i,
                                    coo_names[1]: grid_val_j,
                                    coo_names[2]: grid_val_k,
                                    coo_names[3]: grid_val_l})
                    run_fs = autofile.fs.run(run_prefix)
                    idx += 1

                    frozen_coordinates = coo_names + list(alt_constraints)
                    print('4d freeze test', coo_names, list(alt_constraints))
                    exists = scn_save_fs.leaf.file.geometry.exists(
                        [coo_names, [grid_val_i, grid_val_j, grid_val_k, grid_val_l]])
                    if not exists or overwrite:
                        driver.run_job(
                            job=elstruct.Job.OPTIMIZATION,
                            script_str=script_str,
                            run_fs=run_fs,
                            geom=zma,
                            spc_info=spc_info,
                            thy_level=thy_level,
                            overwrite=overwrite,
                            frozen_coordinates=frozen_coordinates,
                            errors=errors,
                            options_mat=options_mat,
                            retry_failed=retry_failed,
                            saddle=saddle,
                            **kwargs
                        )

                        ret = driver.read_job(
                            job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
                        if update_guess and ret is not None:
                            inf_obj, _, out_str = ret
                            prog = inf_obj.prog
                            guess_zma = elstruct.reader.opt_zmatrix(prog, out_str)


def save_scan(scn_run_fs, scn_save_fs, coo_names,
              gradient=False, hessian=False):
    """ save the scans that have been run so far
    """
    if not scn_run_fs.branch.exists([coo_names]):
        print("No scan to save. Skipping...")
    else:
        locs_lst = []
        for locs in scn_run_fs.leaf.existing([coo_names]):
            if not isinstance(locs[1][0], float):
                continue
            run_path = scn_run_fs.leaf.path(locs)
            run_fs = autofile.fs.run(run_path)
            print("Reading from scan run at {}".format(run_path))

            ret = driver.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if ret:
                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)
                geo = elstruct.reader.opt_geometry(prog, out_str)
                zma = elstruct.reader.opt_zmatrix(prog, out_str)

                save_path = scn_save_fs.leaf.path(locs)
                print(" - Saving...")
                print(" - Save path: {}".format(save_path))

                scn_save_fs.leaf.create(locs)
                scn_save_fs.leaf.file.geometry_info.write(inf_obj, locs)
                scn_save_fs.leaf.file.geometry_input.write(inp_str, locs)
                scn_save_fs.leaf.file.energy.write(ene, locs)
                scn_save_fs.leaf.file.geometry.write(geo, locs)
                scn_save_fs.leaf.file.zmatrix.write(zma, locs)

                locs_lst.append(locs)

                if gradient:
                    ret = driver.read_job(
                        job=elstruct.Job.GRADIENT, run_fs=run_fs)
                    if ret:
                        inf_obj, inp_str, out_str = ret
                        prog = inf_obj.prog
                        method = inf_obj.method
                        grad = elstruct.reader.gradient(prog, out_str)
                        scn_save_fs.leaf.file.gradient.write(grad, locs)

                if hessian:
                    ret = driver.read_job(
                        job=elstruct.Job.HESSIAN, run_fs=run_fs)
                    if ret:
                        inf_obj, inp_str, out_str = ret
                        prog = inf_obj.prog
                        method = inf_obj.method
                        hess = elstruct.reader.hessian(prog, out_str)
                        scn_save_fs.leaf.file.hessian.write(hess, locs)
                        if prog == 'molpro2015':
                            geo = hess_geometry(out_str)
                            scn_save_fs.leaf.file.geometry.write(geo, locs)

        if locs_lst:
            idxs_lst = [locs[-1] for locs in locs_lst]
            enes = [scn_save_fs.leaf.file.energy.read(locs)
                    for locs in locs_lst]
            geos = [scn_save_fs.leaf.file.geometry.read(locs)
                    for locs in locs_lst]

            traj = []
            for idxs, ene, geo in zip(idxs_lst, enes, geos):
                comment = 'energy: {:>15.10f}, grid idxs: {}'.format(ene, idxs)
                traj.append((comment, geo))

            traj_path = scn_save_fs.branch.file.trajectory.path([coo_names])
            print("Updating scan trajectory file at {}".format(traj_path))
            scn_save_fs.branch.file.trajectory.write(traj, [coo_names])


def infinite_separation_energy(
        spc_1_info, spc_2_info, ts_info, high_mul, ref_zma, ini_thy_info,
        thy_info,
        multi_info, run_prefix, save_prefix, scn_run_fs, scn_save_fs, locs,
        overwrite=False,
        num_act_elc=None, num_act_orb=None):
    """ Obtain the infinite separation energy from the multireference energy
        at a given reference point, the high-spin low-spin splitting at that
        reference point, and the high level energy for the high spin state
        at the reference geometry and for the fragments
    """

    # Initialize infinite sep energy
    inf_sep_ene = -1.0e12

    # set up all the file systems for the TS
    # start with the geo and reference theory info
    geo_run_path = scn_run_fs.leaf.path(locs)
    geo_save_path = scn_save_fs.leaf.path(locs)
    geo = scn_save_fs.leaf.file.geometry.read(locs)
    sp_run_fs = autofile.fs.single_point(geo_run_path)
    sp_save_fs = autofile.fs.single_point(geo_save_path)

    # get the multi reference ene for low spin state for the ref point on scan

    # file system for low spin multireference calculation

    multi_info[0] = 'molpro2015'
    multi_info[1] = 'caspt2'
    # ultimately the above should be properly passed
    prog = multi_info[0]
    method = multi_info[1]

    # get the multi reference energy for high spin state for ref point on scan
    hs_info = (ts_info[0], ts_info[1], high_mul)
    orb_restr = fsorb.orbital_restriction(hs_info, multi_info)
    multi_lvl = multi_info[0:3]
    multi_lvl.append(orb_restr)

    hs_run_fs = autofile.fs.high_spin(geo_run_path)
    hs_save_fs = autofile.fs.high_spin(geo_save_path)
    hs_run_fs.leaf.create(multi_lvl[1:4])
    hs_save_fs.leaf.create(multi_lvl[1:4])

    hs_mr_run_path = hs_run_fs.leaf.path(multi_lvl[1:4])
    hs_mr_save_path = hs_save_fs.leaf.path(multi_lvl[1:4])
    run_mr_fs = autofile.fs.run(hs_mr_run_path)

    mr_script_str, _, mr_kwargs, _ = runpar.run_qchem_par(prog, method)

    if num_act_elc is None and num_act_orb is None:
        num_act_elc = high_mul
        num_act_orb = num_act_elc
    ts_formula = automol.geom.formula(automol.zmatrix.geometry(ref_zma))

    cas_opt, _ = ts.cas_options_2(
        hs_info, ts_formula, num_act_elc, num_act_orb, high_mul)
    guess_str = ts.multiref_wavefunction_guess(
        high_mul, ref_zma, hs_info, multi_lvl, [cas_opt])
    guess_lines = guess_str.splitlines()

    mr_kwargs['casscf_options'] = cas_opt
    mr_kwargs['mol_options'] = ['nosym']
    mr_kwargs['gen_lines'] = {1: guess_lines}

    ret = driver.read_job(
        job='energy',
        run_fs=run_mr_fs,
    )
    if ret:
        print(" - Reading high spin multi reference energy from output...")
        inf_obj, inp_str, out_str = ret
        ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)
        hs_save_fs.leaf.file.energy.write(ene, multi_lvl[1:4])
        hs_save_fs.leaf.file.input.write(inp_str, multi_lvl[1:4])
        hs_save_fs.leaf.file.info.write(inf_obj, multi_lvl[1:4])

    if not hs_save_fs.leaf.file.energy.exists(multi_lvl[1:4]) or overwrite:
        print(" - Running high spin multi reference energy ...")
        driver.run_job(
            job='energy',
            script_str=mr_script_str,
            run_fs=run_mr_fs,
            geom=geo,
            spc_info=hs_info,
            thy_level=multi_lvl,
            overwrite=overwrite,
            **mr_kwargs,
        )

        ret = driver.read_job(
            job='energy',
            run_fs=run_mr_fs,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            print(" - Reading high spin multi reference energy from output...")
            hs_mr_ene = elstruct.reader.energy(
                inf_obj.prog, inf_obj.method, out_str)

            print(" - Saving high spin multi reference energy...")
            print(" - Save path: {}".format(hs_mr_save_path))
            hs_save_fs.leaf.file.energy.write(hs_mr_ene, multi_lvl[1:4])
            hs_save_fs.leaf.file.input.write(inp_str, multi_lvl[1:4])
            hs_save_fs.leaf.file.info.write(inf_obj, multi_lvl[1:4])
        else:
            print('ERROR: high spin multi reference energy job fails: ',
                  'Energy is needed to evaluate infinite separation energy')
            return inf_sep_ene

    else:
        hs_mr_ene = hs_save_fs.leaf.file.energy.read(multi_lvl[1:4])

    # file system for high spin single ireference calculation
    thy_info = ['molpro2015', 'ccsd(t)-f12', 'cc-pvdz-f12', 'RR']
    orb_restr = fsorb.orbital_restriction(hs_info, thy_info)
    thy_lvl = thy_info[0:3]
    thy_lvl.append(orb_restr)

    hs_run_fs.leaf.create(thy_lvl[1:4])
    hs_save_fs.leaf.create(thy_lvl[1:4])

    hs_sr_run_path = hs_run_fs.leaf.path(thy_lvl[1:4])
    hs_sr_save_path = hs_save_fs.leaf.path(thy_lvl[1:4])
    run_sr_fs = autofile.fs.run(hs_sr_run_path)

    sp_script_str, _, kwargs, _ = runpar.run_qchem_par(*thy_lvl[0:2])
    ret = driver.read_job(
        job='energy',
        run_fs=run_sr_fs,
    )
    if ret:
        print(" - Reading high spin single reference energy from output...")
        inf_obj, inp_str, out_str = ret
        ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)
        hs_save_fs.leaf.file.energy.write(ene, thy_lvl[1:4])
        hs_save_fs.leaf.file.input.write(inp_str, thy_lvl[1:4])
        hs_save_fs.leaf.file.info.write(inf_obj, thy_lvl[1:4])

    if not hs_save_fs.leaf.file.energy.exists(thy_lvl[1:4]) or overwrite:
        print(" - Running high spin single reference energy ...")

        errors, options_mat = runpar.set_molpro_options_mat(hs_info, geo)

        driver.run_job(
            job='energy',
            script_str=sp_script_str,
            run_fs=run_sr_fs,
            geom=geo,
            spc_info=hs_info,
            thy_level=thy_lvl,
            errors=errors,
            options_mat=options_mat,
            overwrite=overwrite,
            **kwargs,
        )

        ret = driver.read_job(
            job='energy',
            run_fs=run_sr_fs,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            print(" - Reading high spin single ref energy from output...")
            hs_sr_ene = elstruct.reader.energy(
                inf_obj.prog, inf_obj.method, out_str)

            print(" - Saving high spin single reference energy...")
            print(" - Save path: {}".format(hs_sr_save_path))
            hs_save_fs.leaf.file.energy.write(hs_sr_ene, thy_lvl[1:4])
            hs_save_fs.leaf.file.input.write(inp_str, thy_lvl[1:4])
            hs_save_fs.leaf.file.info.write(inf_obj, thy_lvl[1:4])
        else:
            print('ERROR: High spin single reference energy job fails: ',
                  'Energy is needed to evaluate infinite separation energy')
            return inf_sep_ene

    else:
        hs_sr_ene = hs_save_fs.leaf.file.energy.read(thy_lvl[1:4])

    # get the single reference energy for each of the reactant configurations
    spc_ene = []
    spc_infos = [spc_1_info, spc_2_info]
    for spc_info in spc_infos:
        # set up the file systems for the reactants one by one
        spc_run_fs = autofile.fs.species(run_prefix)
        spc_run_fs.leaf.create(spc_info)
        spc_run_path = spc_run_fs.leaf.path(spc_info)
        spc_save_fs = autofile.fs.species(save_prefix)
        spc_save_fs.leaf.create(spc_info)
        spc_save_path = spc_save_fs.leaf.path(spc_info)

        orb_restr = fsorb.orbital_restriction(spc_info, ini_thy_info)
        ini_thy_lvl = ini_thy_info[0:3]
        ini_thy_lvl.append(orb_restr)

        orb_restr = fsorb.orbital_restriction(spc_info, thy_info)
        thy_lvl = thy_info[0:3]
        thy_lvl.append(orb_restr)

        ini_thy_run_fs = autofile.fs.theory(spc_run_path)
        ini_thy_run_fs.leaf.create(ini_thy_lvl[1:4])
        ini_thy_run_path = ini_thy_run_fs.leaf.path(ini_thy_lvl[1:4])
        ini_thy_save_fs = autofile.fs.theory(spc_save_path)
        ini_thy_save_fs.leaf.create(ini_thy_lvl[1:4])
        ini_thy_save_path = ini_thy_save_fs.leaf.path(ini_thy_lvl[1:4])
        ini_cnf_run_fs = autofile.fs.conformer(ini_thy_run_path)
        ini_cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)
        min_cnf_locs = fsmin.min_energy_conformer_locators(
            ini_cnf_save_fs)
        min_cnf_run_path = ini_cnf_run_fs.leaf.path(min_cnf_locs)
        min_cnf_save_path = ini_cnf_save_fs.leaf.path(min_cnf_locs)

        thy_run_fs = autofile.fs.theory(spc_run_path)
        thy_run_fs.leaf.create(thy_lvl[1:4])
        thy_save_fs = autofile.fs.theory(spc_save_path)
        thy_save_fs.leaf.create(thy_lvl[1:4])

        geo = ini_cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)

        sp_run_fs = autofile.fs.single_point(min_cnf_run_path)
        sp_run_fs.leaf.create(thy_lvl[1:4])
        sp_save_fs = autofile.fs.single_point(min_cnf_save_path)
        sp_save_fs.leaf.create(thy_lvl[1:4])

        sp_sr_run_path = sp_run_fs.leaf.path(thy_lvl[1:4])
        sp_sr_save_path = sp_save_fs.leaf.path(thy_lvl[1:4])
        print('sp_sr_run_path')
        print(sp_sr_run_path)
        run_sr_fs = autofile.fs.run(sp_sr_run_path)

        ret = driver.read_job(
            job='energy',
            run_fs=run_sr_fs,
        )
        if ret:
            print(" - Reading single reference energy for",
                  "{} from output...".format(spc_info[0]))
            inf_obj, inp_str, out_str = ret
            ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)
            sp_save_fs.leaf.file.energy.write(ene, thy_lvl[1:4])
            sp_save_fs.leaf.file.input.write(inp_str, thy_lvl[1:4])
            sp_save_fs.leaf.file.info.write(inf_obj, thy_lvl[1:4])

        if not sp_save_fs.leaf.file.energy.exists(thy_lvl[1:4]) or overwrite:
            print(" - Running single reference energy for",
                  "{} from output...".format(spc_info[0]))
            driver.run_job(
                job='energy',
                script_str=sp_script_str,
                run_fs=run_sr_fs,
                geom=geo,
                spc_info=spc_info,
                thy_level=thy_lvl,
                overwrite=overwrite,
                **kwargs,
            )

            ret = driver.read_job(
                job='energy',
                run_fs=run_sr_fs,
            )

            if ret is not None:
                inf_obj, inp_str, out_str = ret

                print(" - Reading single reference energy for ",
                      "{} from output...".format(spc_info[0]))
                sp_sr_ene = elstruct.reader.energy(
                    inf_obj.prog, inf_obj.method, out_str)

                print(" - Saving single reference energy for ",
                      "{} from output...".format(spc_info[0]))
                print(" - Save path: {}".format(sp_sr_save_path))
                sp_save_fs.leaf.file.energy.write(sp_sr_ene, thy_lvl[1:4])
                sp_save_fs.leaf.file.input.write(inp_str, thy_lvl[1:4])
                sp_save_fs.leaf.file.info.write(inf_obj, thy_lvl[1:4])

            else:
                print('ERROR: Single reference energy job fails',
                      'for {}: '.format(spc_info[0]),
                      'Energy needed to evaluate infinite separation energy')
                return inf_sep_ene

        else:
            sp_sr_ene = sp_save_fs.leaf.file.energy.read(thy_lvl[1:4])

        spc_ene.append(sp_sr_ene)

    inf_sep_ene = spc_ene[0] + spc_ene[1] - hs_sr_ene + hs_mr_ene

    return inf_sep_ene
