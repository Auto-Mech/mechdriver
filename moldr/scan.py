""" drivers for coordinate scans
"""
import numpy
import automol
import elstruct
import autofile
import moldr
from elstruct.reader._molpro2015.molecule import hess_geometry


def hindered_rotor_scans(
        spc_info, thy_level, cnf_run_fs, cnf_save_fs, script_str, overwrite,
        scan_increment=30., saddle=False, tors_names='', frm_bnd_key=[],
        brk_bnd_key=[], **opt_kwargs):
    """ Perform 1d scans over each of the torsional coordinates
    """
    min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
    if min_cnf_locs:
        min_cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
        min_cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
        scn_run_fs = autofile.fs.scan(min_cnf_run_path)
        scn_save_fs = autofile.fs.scan(min_cnf_save_path)
        geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        zma = cnf_save_fs.leaf.file.zmatrix.read(min_cnf_locs)
        val_dct = automol.zmatrix.values(zma)
        if not saddle:
            tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
        if tors_names:
            tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
                zma, tors_names, scan_increment, frm_bnd_key=frm_bnd_key,
                brk_bnd_key=brk_bnd_key)
            tors_grids = [
                numpy.linspace(*linspace) + val_dct[name]
                for name, linspace in zip(tors_names, tors_linspaces)]

            for tors_name, tors_grid in zip(tors_names, tors_grids):
                save_scan(
                    scn_run_fs=scn_run_fs,
                    scn_save_fs=scn_save_fs,
                    coo_names=[tors_name],
                )

                run_scan(
                    zma=zma,
                    spc_info=spc_info,
                    thy_level=thy_level,
                    grid_dct={tors_name: tors_grid},
                    scn_run_fs=scn_run_fs,
                    scn_save_fs=scn_save_fs,
                    script_str=script_str,
                    overwrite=overwrite,
                    saddle=saddle,
                    **opt_kwargs,
                )

                save_scan(
                    scn_run_fs=scn_run_fs,
                    scn_save_fs=scn_save_fs,
                    coo_names=[tors_name],
                )


def run_scan(
        zma, spc_info, thy_level, grid_dct, scn_run_fs, scn_save_fs,
        script_str, overwrite, update_guess=True,
        reverse_sweep=True, fix_failures=True, saddle=False,
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
                **kwargs
            )

    elif len(grid_vals) == 2:
        run_prefixes = []
        for grid_val_i in grid_vals[0]:
            for grid_val_j in grid_vals[1]:
                scn_run_fs.leaf.create([coo_names, [grid_val_i, grid_val_j]])
                run_prefixes.append(scn_run_fs.leaf.path([coo_names, [grid_val_i, grid_val_j]]))
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
            **kwargs
        )

        if reverse_sweep:
            run_prefixes = []
            for grid_val_i in grid_vals[0][::-1]:
                for grid_val_j in grid_vals[1][::-1]:
                    run_prefixes.append(scn_run_fs.leaf.path([coo_names, [grid_val_i, grid_val_j]]))
            run_prefixes = tuple(run_prefixes)
            _run_2d_scan(
                script_str=script_str,
                run_prefixes=run_prefixes,
                scn_save_fs=scn_save_fs,
                guess_zma=zma,
                coo_names=coo_names,
                grid_idxs=list(reversed(grid_idxs)),
                grid_vals=[list(reversed(grid_vals[0])), list(reversed(grid_vals[1]))],
                spc_info=spc_info,
                thy_level=thy_level,
                overwrite=overwrite,
                update_guess=update_guess,
                saddle=saddle,
                **kwargs
            )


def run_multiref_rscan(
        formula, high_mul, zma, spc_info, multi_level, dist_name, grid1, grid2,
        scn_run_fs, scn_save_fs, script_str, overwrite, update_guess=True, gradient=False, hessian=False, num_act_elc=None, num_act_orb=None,
        **kwargs):
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

    _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(prog, method)

    if num_act_elc is None and num_act_orb is None:
        num_act_elc = high_mul - 1
        num_act_orb = num_act_elc

    ref_zma = automol.zmatrix.set_values(zma, {coo_names[0]: grid_vals[0][0]})
    cas_opt = ['', '']
    cas_opt[0], _ = moldr.ts.cas_options_1(spc_info, formula, num_act_elc, num_act_orb, high_mul)
    cas_opt[1], _ = moldr.ts.cas_options_2(spc_info, formula, num_act_elc, num_act_orb, high_mul)
    guess_str = moldr.ts.multiref_wavefunction_guess(
        high_mul, ref_zma, spc_info, multi_level, cas_opt)
    guess_lines = guess_str.splitlines()

    opt_kwargs['casscf_options'] = cas_opt[1]
    opt_kwargs['mol_options'] = ['nosym']
    opt_kwargs['gen_lines'] = {1: guess_lines}

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
        **opt_kwargs,
    )

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
        **opt_kwargs,
    )


def _run_1d_scan(
        script_str, run_prefixes, scn_save_fs, guess_zma, coo_name, grid_idxs, grid_vals,
        spc_info, thy_level, overwrite, errors=(), options_mat=(),
        retry_failed=True, update_guess=True, saddle=False, gradient=False, hessian=False,
        **kwargs):
    """ run 1 dimensional scan with constrained optimization
    """

    npoints = len(grid_idxs)
    assert len(grid_vals) == len(run_prefixes) == npoints
    for grid_idx, grid_val, run_prefix in zip(grid_idxs, grid_vals, run_prefixes):
        print("Point {}/{}".format(grid_idx+1, npoints))
        zma = automol.zmatrix.set_values(guess_zma, {coo_name: grid_val})
        run_fs = autofile.fs.run(run_prefix)

        if not scn_save_fs.leaf.file.geometry.exists([[coo_name], [grid_val]]) or overwrite:
            moldr.driver.run_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                run_fs=run_fs,
                geom=zma,
                spc_info=spc_info,
                thy_level=thy_level,
                overwrite=overwrite,
                frozen_coordinates=[coo_name],
                errors=errors,
                options_mat=options_mat,
                retry_failed=retry_failed,
                saddle=saddle,
                **kwargs
            )

            ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if ret is not None:
                inf_obj, _, out_str = ret
                prog = inf_obj.prog
                opt_zma = elstruct.reader.opt_zmatrix(prog, out_str)
                if update_guess:
                    guess_zma = opt_zma

                if gradient:
                    moldr.driver.run_job(
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

                    ret = moldr.driver.read_job(job=elstruct.Job.GRADIENT, run_fs=run_fs)

                if hessian:
                    moldr.driver.run_job(
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

                    ret = moldr.driver.read_job(job=elstruct.Job.HESSIAN, run_fs=run_fs)


def _run_2d_scan(
        script_str, run_prefixes, scn_save_fs, guess_zma, coo_names, grid_idxs, grid_vals,
        spc_info, thy_level, overwrite, errors=(),
        options_mat=(), retry_failed=True, update_guess=True, saddle=False, **kwargs):
    """ run 2-dimensional scan with constrained optimization
    """

    npoints = len(grid_idxs)
    assert len(grid_vals[0]) * len(grid_vals[1]) == len(run_prefixes) == npoints

    idx = 0
    for grid_val_i in grid_vals[0]:
        for grid_val_j in grid_vals[1]:
            grid_idx = grid_idxs[idx]
            run_prefix = run_prefixes[idx]
            print("Point {}/{}".format(grid_idx+1, npoints))
            zma = automol.zmatrix.set_values(
                guess_zma, {coo_names[0]: grid_val_i, coo_names[1]: grid_val_j})
            run_fs = autofile.fs.run(run_prefix)
            idx += 1

            if not scn_save_fs.leaf.file.geometry.exists(
                    [coo_names, [grid_val_i, grid_val_j]]) or overwrite:
                moldr.driver.run_job(
                    job=elstruct.Job.OPTIMIZATION,
                    script_str=script_str,
                    run_fs=run_fs,
                    geom=zma,
                    spc_info=spc_info,
                    thy_level=thy_level,
                    overwrite=overwrite,
                    frozen_coordinates=coo_names,
                    errors=errors,
                    options_mat=options_mat,
                    retry_failed=retry_failed,
                    saddle=saddle,
                    **kwargs
                )

                ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
                if update_guess and ret is not None:
                    inf_obj, _, out_str = ret
                    prog = inf_obj.prog
                    guess_zma = elstruct.reader.opt_zmatrix(prog, out_str)


def save_scan(scn_run_fs, scn_save_fs, coo_names, gradient=False, hessian=False):
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

            ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
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
                    ret = moldr.driver.read_job(job=elstruct.Job.GRADIENT, run_fs=run_fs)
                    if ret:
                        inf_obj, inp_str, out_str = ret
                        prog = inf_obj.prog
                        method = inf_obj.method
                        grad = elstruct.reader.gradient(prog, out_str)
                        scn_save_fs.leaf.file.gradient.write(grad, locs)

                if hessian:
                    ret = moldr.driver.read_job(job=elstruct.Job.HESSIAN, run_fs=run_fs)
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
        spc_1_info, spc_2_info, ts_info, high_mul, ref_zma, ini_thy_info, thy_info,
        multi_info, run_prefix, save_prefix, scn_run_fs, scn_save_fs, locs, overwrite=False,
        num_act_elc=None, num_act_orb=None):
    """ Obtain the infinite separation energy from the multireference energy at a given
    reference point, the high-spin low-spin splitting at that reference point, and the
    high level energy for the high spin state at the reference geometry and for the fragments
    """

    # set up all the file systems for the TS
    # start with the geo and reference theory info
    geo_run_path = scn_run_fs.leaf.path(locs)
    geo_save_path = scn_save_fs.leaf.path(locs)
    geo = scn_save_fs.leaf.file.geometry.read(locs)
    sp_run_fs = autofile.fs.single_point(geo_run_path)
    sp_save_fs = autofile.fs.single_point(geo_save_path)

    # get the multi reference energy for the low spin state for the reference point on the scan

    # file system for low spin multireference calculation

    multi_info[0] = 'molpro2015'
    multi_info[1] = 'caspt2'
    # ultimately the above should be properly passed
    prog = multi_info[0]
    method = multi_info[1]

#    orb_restr = moldr.util.orbital_restriction(ts_info, multi_info)
#    multi_lvl = multi_info[0:3]
#    multi_lvl.append(orb_restr)

#    sp_run_fs.leaf.create(multi_lvl[1:4])
#    sp_save_fs.leaf.create(multi_lvl[1:4])

#    sp_mr_run_path = sp_run_fs.leaf.path(multi_lvl[1:4])
#    sp_mr_save_path = sp_save_fs.leaf.path(multi_lvl[1:4])
#    run_mr_fs = autofile.fs.run(sp_mr_run_path)

#    mr_script_str, _, mr_kwargs, _ = moldr.util.run_qchem_par(prog, method)

#    num_act_elc = high_mul
#    num_act_orb = num_act_elc
#    ts_formula = automol.geom.formula(automol.zmatrix.geometry(ref_zma))

#    cas_opt, _ = moldr.ts.cas_options_2(ts_info, ts_formula, num_act_elc, num_act_orb, high_mul)
#    guess_str = moldr.ts.multiref_wavefunction_guess(high_mul, ref_zma, ts_info, multi_lvl, cas_opt)
#    guess_lines = guess_str.splitlines()

#    mr_kwargs['casscf_options'] = cas_opt
#    mr_kwargs['mol_options'] = ['nosym']
#    mr_kwargs['gen_lines'] = {1: guess_lines}

#    ret = moldr.driver.read_job(
#        job='energy',
#        run_fs=run_mr_fs,
#    )
#    if ret:
#        print(" - Reading low spin multi reference energy from output...")
#        inf_obj, inp_str, out_str = ret
#        ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)
#        sp_save_fs.leaf.file.energy.write(ene, multi_lvl[1:4])
#        sp_save_fs.leaf.file.input.write(inp_str, multi_lvl[1:4])
#        sp_save_fs.leaf.file.info.write(inf_obj, multi_lvl[1:4])
#
#    if not sp_save_fs.leaf.file.energy.exists(multi_lvl[1:4]) or overwrite:
#        print(" - Running low spin multi reference energy ...")
#        moldr.driver.run_job(
#            job='energy',
#            script_str=mr_script_str,
#            run_fs=run_mr_fs,
#            geom=geo,
#            spc_info=ts_info,
#            thy_level=multi_lvl,
#            overwrite=overwrite,
#            **mr_kwargs,
#        )
#
#        ret = moldr.driver.read_job(
#            job='energy',
#            run_fs=run_mr_fs,
#        )
#
#        if ret is not None:
#            inf_obj, inp_str, out_str = ret
#
#            print(" - Reading low spin multi reference energy from output...")
#            ls_mr_ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)

#            print(" - Saving low spin multi reference energy...")
#            print(" - Save path: {}".format(sp_mr_save_path))
#            sp_save_fs.leaf.file.energy.write(ls_mr_ene, multi_lvl[1:4])
#            sp_save_fs.leaf.file.input.write(inp_str, multi_lvl[1:4])
#            sp_save_fs.leaf.file.info.write(inf_obj, multi_lvl[1:4])
#        else:
#            print('ERROR: low spin multi reference energy job fails: ',
#                  'Energy is needed to evaluate infinite separation energy')
#            return

#    else:
#        ls_mr_ene = sp_save_fs.leaf.file.energy.read(multi_lvl[1:4])

#    print('low spin energy:', ls_mr_ene)
    # get the multi reference energy for the high spin state for the reference point on the scan

    hs_info = (ts_info[0], ts_info[1], high_mul)
    orb_restr = moldr.util.orbital_restriction(hs_info, multi_info)
    multi_lvl = multi_info[0:3]
    multi_lvl.append(orb_restr)

    hs_run_fs = autofile.fs.high_spin(geo_run_path)
    hs_save_fs = autofile.fs.high_spin(geo_save_path)
    hs_run_fs.leaf.create(multi_lvl[1:4])
    hs_save_fs.leaf.create(multi_lvl[1:4])

    hs_mr_run_path = hs_run_fs.leaf.path(multi_lvl[1:4])
    hs_mr_save_path = hs_save_fs.leaf.path(multi_lvl[1:4])
    run_mr_fs = autofile.fs.run(hs_mr_run_path)

    mr_script_str, _, mr_kwargs, _ = moldr.util.run_qchem_par(prog, method)

    if num_act_elc is None and num_act_orb is None:
        num_act_elc = high_mul
        num_act_orb = num_act_elc
    # num_act_elc = high_mul
    # num_act_orb = num_act_elc
    ts_formula = automol.geom.formula(automol.zmatrix.geometry(ref_zma))

    cas_opt, _ = moldr.ts.cas_options_2(hs_info, ts_formula, num_act_elc, num_act_orb, high_mul)
    guess_str = moldr.ts.multiref_wavefunction_guess(high_mul, ref_zma, hs_info, multi_lvl, [cas_opt])
    guess_lines = guess_str.splitlines()

    mr_kwargs['casscf_options'] = cas_opt
    mr_kwargs['mol_options'] = ['nosym']
    mr_kwargs['gen_lines'] = {1: guess_lines}

    ret = moldr.driver.read_job(
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
        moldr.driver.run_job(
            job='energy',
            script_str=mr_script_str,
            run_fs=run_mr_fs,
            geom=geo,
            spc_info=hs_info,
            thy_level=multi_lvl,
            overwrite=overwrite,
            **mr_kwargs,
        )

        ret = moldr.driver.read_job(
            job='energy',
            run_fs=run_mr_fs,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            print(" - Reading high spin multi reference energy from output...")
            hs_mr_ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)

            print(" - Saving high spin multi reference energy...")
            print(" - Save path: {}".format(hs_mr_save_path))
            hs_save_fs.leaf.file.energy.write(hs_mr_ene, multi_lvl[1:4])
            hs_save_fs.leaf.file.input.write(inp_str, multi_lvl[1:4])
            hs_save_fs.leaf.file.info.write(inf_obj, multi_lvl[1:4])
        else:
            print('ERROR: high spin multi reference energy job fails: ',
                  'Energy is needed to evaluate infinite separation energy')
            return

    else:
        hs_mr_ene = hs_save_fs.leaf.file.energy.read(multi_lvl[1:4])


    # get the single reference energy for the high spin state for the reference point on the scan
    # file system for high spin single ireference calculation
    thy_info = ['molpro2015', 'ccsd(t)-f12', 'cc-pvdz-f12', 'RR']
    orb_restr = moldr.util.orbital_restriction(hs_info, thy_info)
    thy_lvl = thy_info[0:3]
    thy_lvl.append(orb_restr)

    hs_run_fs.leaf.create(thy_lvl[1:4])
    hs_save_fs.leaf.create(thy_lvl[1:4])

    hs_sr_run_path = hs_run_fs.leaf.path(thy_lvl[1:4])
    hs_sr_save_path = hs_save_fs.leaf.path(thy_lvl[1:4])
    run_sr_fs = autofile.fs.run(hs_sr_run_path)

    sp_script_str, _, kwargs, _ = moldr.util.run_qchem_par(*thy_lvl[0:2])
    ret = moldr.driver.read_job(
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

        errors, options_mat = moldr.util.set_molpro_options_mat(hs_info, geo)

        moldr.driver.run_job(
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

        ret = moldr.driver.read_job(
            job='energy',
            run_fs=run_sr_fs,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            print(" - Reading high spin single reference energy from output...")
            hs_sr_ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)

            print(" - Saving high spin single reference energy...")
            print(" - Save path: {}".format(hs_sr_save_path))
            hs_save_fs.leaf.file.energy.write(hs_sr_ene, thy_lvl[1:4])
            hs_save_fs.leaf.file.input.write(inp_str, thy_lvl[1:4])
            hs_save_fs.leaf.file.info.write(inf_obj, thy_lvl[1:4])
        else:
            print('ERROR: High spin single reference energy job fails: ',
                  'Energy is needed to evaluate infinite separation energy')
            return

    else:
        hs_sr_ene = hs_save_fs.leaf.file.energy.read(thy_lvl[1:4])

    # get the single reference energy for each of the reactant configurations
    spc_ene = []
    spc_infos = [spc_1_info, spc_2_info]
    #for spc_info in zip(spc_1_info, spc_2_info):
    for spc_info in spc_infos:
        # set up the file systems for the reactants one by one
        spc_run_fs = autofile.fs.species(run_prefix)
        spc_run_fs.leaf.create(spc_info)
        spc_run_path = spc_run_fs.leaf.path(spc_info)
        spc_save_fs = autofile.fs.species(save_prefix)
        spc_save_fs.leaf.create(spc_info)
        spc_save_path = spc_save_fs.leaf.path(spc_info)

        orb_restr = moldr.util.orbital_restriction(spc_info, ini_thy_info)
        ini_thy_lvl = ini_thy_info[0:3]
        ini_thy_lvl.append(orb_restr)

        orb_restr = moldr.util.orbital_restriction(spc_info, thy_info)
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
        min_cnf_locs = moldr.util.min_energy_conformer_locators(ini_cnf_save_fs)
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

        # get the single reference energy for the high spin state for the reference point on the scan

        ret = moldr.driver.read_job(
            job='energy',
            run_fs=run_sr_fs,
        )
        if ret:
            print(" - Reading single reference energy for {} from output...".format(spc_info[0]))
            inf_obj, inp_str, out_str = ret
            ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)
            sp_save_fs.leaf.file.energy.write(ene, thy_lvl[1:4])
            sp_save_fs.leaf.file.input.write(inp_str, thy_lvl[1:4])
            sp_save_fs.leaf.file.info.write(inf_obj, thy_lvl[1:4])

        if not sp_save_fs.leaf.file.energy.exists(thy_lvl[1:4]) or overwrite:
            print(" - Running single reference energy for {} from output...".format(spc_info[0]))
            moldr.driver.run_job(
                job='energy',
                script_str=sp_script_str,
                run_fs=run_sr_fs,
                geom=geo,
                spc_info=spc_info,
                thy_level=thy_lvl,
                overwrite=overwrite,
                **kwargs,
            )

            ret = moldr.driver.read_job(
                job='energy',
                run_fs=run_sr_fs,
            )

            if ret is not None:
                inf_obj, inp_str, out_str = ret

                print(" - Reading single reference energy for {} from output...".format(spc_info[0]))
                sp_sr_ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)

                print(" - Saving single reference energy for {} from output...".format(spc_info[0]))
                print(" - Save path: {}".format(sp_sr_save_path))
                sp_save_fs.leaf.file.energy.write(sp_sr_ene, thy_lvl[1:4])
                sp_save_fs.leaf.file.input.write(inp_str, thy_lvl[1:4])
                sp_save_fs.leaf.file.info.write(inf_obj, thy_lvl[1:4])

            else:
                print('ERROR: Single reference energy job fails for {}: '.format(spc_info[0]),
                      'Energy is needed to evaluate infinite separation energy')
                return

        else:
            sp_sr_ene = sp_save_fs.leaf.file.energy.read(thy_lvl[1:4])

        spc_ene.append(sp_sr_ene)

    inf_sep_ene = spc_ene[0] + spc_ene[1] - hs_sr_ene + hs_mr_ene
    #inf_sep_ene = spc_ene[0] + spc_ene[1] - hs_sr_ene + hs_mr_ene - ls_mr_ene + ls_mr_ene
    # print('inf_sep_ene test:', inf_sep_ene)

    return inf_sep_ene
