""" es_runners for coordinate scans
"""

import automol
import elstruct
import autofile
from routines.es import runner as es_runner


def run_scan(
        zma, spc_info, thy_info, grid_dct, scn_run_fs, scn_save_fs,
        script_str, overwrite, update_guess=True,
        reverse_sweep=True, fix_failures=True, saddle=False,
        constraint_dct=None,
        **kwargs):
    """ run constrained optimization scan
    """

    vma = automol.zmatrix.var_(zma)
    if scn_save_fs[0].file.vmatrix.exists():
        existing_vma = scn_save_fs[0].file.vmatrix.read()
        assert vma == existing_vma

    coo_names = []
    grid_vals = []
    for item in grid_dct.items():
        (coo, coo_grid_vals) = item
        coo_names.append(coo)
        grid_vals.append(coo_grid_vals)

    # for now, running only one-dimensional hindered rotor scans
    scn_save_fs[1].create([coo_names])
    inf_obj = autofile.system.info.scan_branch(grid_dct)
    scn_save_fs[1].file.info.write(inf_obj, [coo_names])
    npoint = 1
    for coo_grid_vals in grid_vals:
        npoint *= len(coo_grid_vals)
    grid_idxs = tuple(range(npoint))
    if len(grid_vals) == 1:
        run_prefixes = []
        for grid_val in grid_vals[0]:
            if constraint_dct is not None:
                scn_run_fs[-1].create(
                    [coo_names, [grid_val], constraint_dct])
                run_prefixes.append(
                    scn_run_fs[-1].path(
                        [coo_names, [grid_val], constraint_dct]))
            else:
                scn_run_fs[-1].create([coo_names, [grid_val]])
                run_prefixes.append(
                    scn_run_fs[-1].path([coo_names, [grid_val]]))
        _run_1d_scan(
            script_str=script_str,
            run_prefixes=run_prefixes,
            scn_save_fs=scn_save_fs,
            guess_zma=zma,
            coo_name=coo_names[0],
            grid_idxs=grid_idxs,
            grid_vals=grid_vals[0],
            spc_info=spc_info,
            thy_info=thy_info,
            overwrite=overwrite,
            update_guess=update_guess,
            saddle=saddle,
            retry_failed=fix_failures,
            constraint_dct=constraint_dct,
            **kwargs
        )

        if reverse_sweep:
            print('\nDoing a reverse sweep of the HR scan to catch errors...')
            _run_1d_scan(
                script_str=script_str,
                run_prefixes=list(reversed(run_prefixes)),
                scn_save_fs=scn_save_fs,
                guess_zma=zma,
                coo_name=coo_names[0],
                grid_idxs=list(reversed(grid_idxs)),
                grid_vals=list(reversed(grid_vals[0])),
                spc_info=spc_info,
                thy_info=thy_info,
                overwrite=overwrite,
                update_guess=update_guess,
                saddle=saddle,
                constraint_dct=constraint_dct,
                **kwargs
            )

    elif len(grid_vals) == 2:
        run_prefixes = []
        for grid_val_i in grid_vals[0]:
            for grid_val_j in grid_vals[1]:
                if constraint_dct is not None:
                    run_prefixes.append(scn_run_fs[-1].path(
                        [coo_names,
                         [grid_val_i, grid_val_j],
                         constraint_dct]))
                else:
                    scn_run_fs[-1].create(
                        [coo_names, [grid_val_i, grid_val_j]])
                    run_prefixes.append(scn_run_fs[-1].path(
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
            thy_info=thy_info,
            overwrite=overwrite,
            update_guess=update_guess,
            saddle=saddle,
            retry_failed=fix_failures,
            constraint_dct=constraint_dct,
            **kwargs
        )

        if reverse_sweep:
            print('\nDoing a reverse sweep of the HR scans to catch errors...')
            run_prefixes = []
            for grid_val_i in grid_vals[0][::-1]:
                for grid_val_j in grid_vals[1][::-1]:
                    if constraint_dct is not None:
                        run_prefixes.append(scn_run_fs[-1].path(
                            [coo_names,
                             [grid_val_i, 2, grid_val_j, 2],
                             constraint_dct]))
                    else:
                        run_prefixes.append(scn_run_fs[-1].path(
                            [coo_names,
                             [grid_val_i, grid_val_j]]))
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
                thy_info=thy_info,
                overwrite=overwrite,
                update_guess=update_guess,
                saddle=saddle,
                constraint_dct=constraint_dct,
                **kwargs
            )

    elif len(grid_vals) == 3:
        run_prefixes = []
        for grid_val_i in grid_vals[0]:
            for grid_val_j in grid_vals[1]:
                for grid_val_k in grid_vals[2]:
                    if constraint_dct is not None:
                        scn_run_fs[-1].create(
                            [coo_names,
                             [grid_val_i, grid_val_j, grid_val_k],
                             constraint_dct])
                        run_prefixes.append(scn_run_fs[-1].path(
                            [coo_names,
                             [grid_val_i, grid_val_j, grid_val_k],
                             constraint_dct]))
                    else:
                        scn_run_fs[-1].create(
                            [coo_names, [grid_val_i, grid_val_j, grid_val_k]])
                        run_prefixes.append(scn_run_fs[-1].path(
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
            thy_info=thy_info,
            overwrite=overwrite,
            update_guess=update_guess,
            saddle=saddle,
            retry_failed=fix_failures,
            constraint_dct=constraint_dct,
            **kwargs
        )

        if reverse_sweep:
            print('\nDoing a reverse sweep of the HRs scan to catch errors...')
            run_prefixes = []
            for grid_val_i in grid_vals[0][::-1]:
                for grid_val_j in grid_vals[1][::-1]:
                    for grid_val_k in grid_vals[2][::-1]:
                        if constraint_dct is not None:
                            run_prefixes.append(scn_run_fs[-1].path(
                                [coo_names,
                                 [grid_val_i, grid_val_j, grid_val_k],
                                 constraint_dct]))
                        else:
                            run_prefixes.append(scn_run_fs[-1].path(
                                [coo_names,
                                 [grid_val_i, grid_val_j, grid_val_k],
                                 constraint_dct]))
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
                thy_info=thy_info,
                overwrite=overwrite,
                update_guess=update_guess,
                saddle=saddle,
                constraint_dct=constraint_dct,
                **kwargs
            )

    elif len(grid_vals) == 4:
        run_prefixes = []
        for grid_val_i in grid_vals[0]:
            for grid_val_j in grid_vals[1]:
                for grid_val_k in grid_vals[2]:
                    for grid_val_l in grid_vals[3]:
                        if constraint_dct is not None:
                            scn_run_fs[-1].create(
                                [coo_names,
                                 [grid_val_i, grid_val_j,
                                  grid_val_k, grid_val_l],
                                 constraint_dct])
                            run_prefixes.append(scn_run_fs[-1].path(
                                [coo_names,
                                 [grid_val_i, grid_val_j,
                                  grid_val_k, grid_val_l],
                                 constraint_dct]))
                        else:
                            scn_run_fs[-1].create(
                                [coo_names,
                                 [grid_val_i, grid_val_j,
                                  grid_val_k, grid_val_l]])
                            run_prefixes.append(scn_run_fs[-1].path(
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
            thy_info=thy_info,
            overwrite=overwrite,
            update_guess=update_guess,
            saddle=saddle,
            retry_failed=fix_failures,
            constraint_dct=constraint_dct,
            **kwargs
        )

        if reverse_sweep:
            print('\nDoing a reverse sweep of the HRs scan to catch errors...')
            run_prefixes = []
            for grid_val_i in grid_vals[0][::-1]:
                for grid_val_j in grid_vals[1][::-1]:
                    for grid_val_k in grid_vals[2][::-1]:
                        for grid_val_l in grid_vals[3][::-1]:
                            if constraint_dct is not None:
                                run_prefixes.append(scn_run_fs[-1].path(
                                    [coo_names,
                                     [grid_val_i, grid_val_j,
                                      grid_val_k, grid_val_l],
                                     constraint_dct]))
                            else:
                                run_prefixes.append(scn_run_fs[-1].path(
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
                thy_info=thy_info,
                overwrite=overwrite,
                update_guess=update_guess,
                saddle=saddle,
                constraint_dct=constraint_dct,
                **kwargs
            )


def _run_1d_scan(
        script_str, run_prefixes, scn_save_fs, guess_zma, coo_name,
        grid_idxs, grid_vals,
        spc_info, thy_info, overwrite, errors=(), options_mat=(),
        retry_failed=True, update_guess=True, saddle=False,
        constraint_dct=None,
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

        if constraint_dct is not None:
            geo_exists = scn_save_fs[-1].file.geometry.exists(
                [[coo_name], [grid_val], constraint_dct])
            frozen_coordinates = [coo_name] + list(constraint_dct)
        else:
            geo_exists = scn_save_fs[-1].file.geometry.exists(
                [[coo_name], [grid_val]])
            frozen_coordinates = [coo_name]
        if not geo_exists or overwrite:
            es_runner.run_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                run_fs=run_fs,
                geom=zma,
                spc_info=spc_info,
                thy_info=thy_info,
                overwrite=overwrite,
                frozen_coordinates=frozen_coordinates,
                errors=errors,
                options_mat=options_mat,
                retry_failed=retry_failed,
                saddle=saddle,
                **kwargs
            )

            ret = es_runner.read_job(
                job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if ret is not None:
                inf_obj, _, out_str = ret
                prog = inf_obj.prog
                opt_zma = elstruct.reader.opt_zmatrix(prog, out_str)
                if update_guess:
                    guess_zma = opt_zma


def _run_2d_scan(
        script_str, run_prefixes, scn_save_fs, guess_zma, coo_names,
        grid_idxs, grid_vals,
        spc_info, thy_info, overwrite, errors=(),
        options_mat=(), retry_failed=True, update_guess=True, saddle=False,
        constraint_dct=None,
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

            frozen_coordinates = coo_names + list(constraint_dct)
            if constraint_dct is not None:
                geo_exists = scn_save_fs[-1].file.geometry.exists(
                    [coo_names, [grid_val_i, grid_val_j], constraint_dct])
                frozen_coordinates = coo_names + list(constraint_dct)
            else:
                geo_exists = scn_save_fs[-1].file.geometry.exists(
                    [coo_names, [grid_val_i, grid_val_j]])
                frozen_coordinates = coo_names
            if not geo_exists or overwrite:
                es_runner.run_job(
                    job=elstruct.Job.OPTIMIZATION,
                    script_str=script_str,
                    run_fs=run_fs,
                    geom=zma,
                    spc_info=spc_info,
                    thy_info=thy_info,
                    overwrite=overwrite,
                    frozen_coordinates=frozen_coordinates,
                    errors=errors,
                    options_mat=options_mat,
                    retry_failed=retry_failed,
                    saddle=saddle,
                    **kwargs
                )

                ret = es_runner.read_job(
                    job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
                if ret is not None:
                    inf_obj, _, out_str = ret
                    prog = inf_obj.prog
                    opt_zma = elstruct.reader.opt_zmatrix(prog, out_str)
                    if update_guess:
                        guess_zma = opt_zma


def _run_3d_scan(
        script_str, run_prefixes, scn_save_fs, guess_zma, coo_names,
        grid_idxs, grid_vals,
        spc_info, thy_info, overwrite, errors=(),
        options_mat=(), retry_failed=True, update_guess=True, saddle=False,
        constraint_dct=None,
        **kwargs):
    """ run 2-dimensional scan with constrained optimization
    """

    npoints = len(grid_idxs)
    gmult = (len(grid_vals[0])*len(grid_vals[1])*len(grid_vals[2]))
    assert gmult == len(run_prefixes) == npoints

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

                if constraint_dct is not None:
                    geo_exists = scn_save_fs[-1].file.geometry.exists(
                        [coo_names, [grid_val_i, grid_val_j, grid_val_k],
                         constraint_dct])
                    frozen_coordinates = coo_names + list(constraint_dct)
                else:
                    geo_exists = scn_save_fs[-1].file.geometry.exists(
                        [coo_names, [grid_val_i, grid_val_j, grid_val_k]])
                    frozen_coordinates = coo_names
                if not geo_exists or overwrite:
                    es_runner.run_job(
                        job=elstruct.Job.OPTIMIZATION,
                        script_str=script_str,
                        run_fs=run_fs,
                        geom=zma,
                        spc_info=spc_info,
                        thy_info=thy_info,
                        overwrite=overwrite,
                        frozen_coordinates=frozen_coordinates,
                        errors=errors,
                        options_mat=options_mat,
                        retry_failed=retry_failed,
                        saddle=saddle,
                        **kwargs
                    )

                    ret = es_runner.read_job(
                        job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
                    if ret is not None:
                        inf_obj, _, out_str = ret
                        prog = inf_obj.prog
                        opt_zma = elstruct.reader.opt_zmatrix(prog, out_str)
                        if update_guess:
                            guess_zma = opt_zma


def _run_4d_scan(
        script_str, run_prefixes, scn_save_fs, guess_zma, coo_names,
        grid_idxs, grid_vals,
        spc_info, thy_info, overwrite, errors=(),
        options_mat=(), retry_failed=True, update_guess=True, saddle=False,
        constraint_dct=None,
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

                    if constraint_dct is not None:
                        geo_exists = scn_save_fs[-1].file.geometry.exists(
                            [coo_names,
                             [grid_val_i, grid_val_j, grid_val_k, grid_val_l],
                             constraint_dct])
                        frozen_coordinates = coo_names + list(constraint_dct)
                    else:
                        geo_exists = scn_save_fs[-1].file.geometry.exists(
                            [coo_names,
                             [grid_val_i, grid_val_j, grid_val_k, grid_val_l]])
                        frozen_coordinates = coo_names
                    if not geo_exists or overwrite:
                        es_runner.run_job(
                            job=elstruct.Job.OPTIMIZATION,
                            script_str=script_str,
                            run_fs=run_fs,
                            geom=zma,
                            spc_info=spc_info,
                            thy_info=thy_info,
                            overwrite=overwrite,
                            frozen_coordinates=frozen_coordinates,
                            errors=errors,
                            options_mat=options_mat,
                            retry_failed=retry_failed,
                            saddle=saddle,
                            **kwargs
                        )

                        ret = es_runner.read_job(
                            job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
                        if ret is not None:
                            inf_obj, _, out_str = ret
                            prog = inf_obj.prog
                            opt_zma = elstruct.reader.opt_zmatrix(
                                prog, out_str)
                            if update_guess:
                                guess_zma = opt_zma


def save_scan(scn_run_fs, scn_save_fs, coo_names, thy_info):
    """ save the scans that have been run so far
    """
    if not scn_run_fs[1].exists([coo_names]):
        print("No scan to save. Skipping...")
    else:
        locs_lst = []
        for locs in scn_run_fs[-1].existing([coo_names]):
            if not isinstance(locs[1][0], float):
                continue
            run_path = scn_run_fs[-1].path(locs)
            run_fs = autofile.fs.run(run_path)
            print("Reading from scan run at {}".format(run_path))

            ret = es_runner.read_job(
                job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if ret:
                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)
                geo = elstruct.reader.opt_geometry(prog, out_str)
                zma = elstruct.reader.opt_zmatrix(prog, out_str)

                save_path = scn_save_fs[-1].path(locs)
                print(" - Saving...")
                print(" - Save path: {}".format(save_path))

                scn_save_fs[-1].create(locs)
                scn_save_fs[-1].file.geometry_info.write(inf_obj, locs)
                scn_save_fs[-1].file.geometry_input.write(inp_str, locs)
                scn_save_fs[-1].file.energy.write(ene, locs)
                scn_save_fs[-1].file.geometry.write(geo, locs)
                scn_save_fs[-1].file.zmatrix.write(zma, locs)

                # Saving the energy to am SP filesys
                print(" - Saving energy...")
                sp_save_fs = autofile.fs.single_point(save_path)
                sp_save_fs[-1].create(thy_info[1:4])
                sp_save_fs[-1].file.input.write(inp_str, thy_info[1:4])
                sp_save_fs[-1].file.info.write(inf_obj, thy_info[1:4])
                sp_save_fs[-1].file.energy.write(ene, thy_info[1:4])

                locs_lst.append(locs)

        if locs_lst:
            idxs_lst = [locs[-1] for locs in locs_lst]
            enes = [scn_save_fs[-1].file.energy.read(locs)
                    for locs in locs_lst]
            geos = [scn_save_fs[-1].file.geometry.read(locs)
                    for locs in locs_lst]

            traj = []
            for idxs, ene, geo in zip(idxs_lst, enes, geos):
                comment = (
                    'energy: {:>15.10f}, '.format(ene) +
                    'grid idxs: {}'.format(idxs)
                )
                traj.append((comment, geo))

            traj_path = scn_save_fs[1].file.trajectory.path([coo_names])
            print("Updating scan trajectory file at {}".format(traj_path))
            scn_save_fs[1].file.trajectory.write(traj, [coo_names])


def save_cscan(cscn_run_fs, cscn_save_fs, coo_names):
    """ save the scans that have been run so far
    """

    print('cscn_path', cscn_run_fs[1].path([coo_names]))
    exists1 = cscn_run_fs[1].exists([coo_names])
    if exists1:
        locs_lst = []
        scn_locs1 = cscn_run_fs[2].existing([coo_names])
        for locs1 in scn_locs1:
            exists2 = cscn_run_fs[2].exists(locs1)
            print('cscn2_path', cscn_run_fs[2].path(locs1))
            if exists2:
                scn_locs2 = cscn_run_fs[3].existing(locs1)
                for locs2 in scn_locs2:
                    run_path = cscn_run_fs[-1].path(locs2)
                    print("Reading from scan run at {}".format(run_path))
                    run_fs = autofile.fs.run(run_path)
                    ret = es_runner.read_job(
                        job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
                    if ret:
                        inf_obj, inp_str, out_str = ret
                        prog = inf_obj.prog
                        method = inf_obj.method
                        ene = elstruct.reader.energy(prog, method, out_str)
                        geo = elstruct.reader.opt_geometry(prog, out_str)
                        zma = elstruct.reader.opt_zmatrix(prog, out_str)

                        save_path = cscn_save_fs[-1].path(locs2)
                        print(" - Saving...")
                        print(" - Save path: {}".format(save_path))

                        cscn_save_fs[-1].create(locs2)
                        cscn_save_fs[-1].file.geometry_info.write(
                            inf_obj, locs2)
                        cscn_save_fs[-1].file.geometry_input.write(
                            inp_str, locs2)
                        cscn_save_fs[-1].file.energy.write(ene, locs2)
                        cscn_save_fs[-1].file.geometry.write(geo, locs2)
                        cscn_save_fs[-1].file.zmatrix.write(zma, locs2)

                        locs_lst.append(locs2)

        if locs_lst:
            idxs_lst = [locs[-1] for locs in locs_lst]
            enes = [cscn_save_fs[-1].file.energy.read(locs)
                    for locs in locs_lst]
            geos = [cscn_save_fs[-1].file.geometry.read(locs)
                    for locs in locs_lst]

            traj = []
            for idxs, ene, geo in zip(idxs_lst, enes, geos):
                comment = (
                    'energy: {:>15.10f}, '.format(ene) +
                    'grid idxs: {}'.format(idxs)
                )
                traj.append((comment, geo))

            traj_path = cscn_save_fs[1].file.trajectory.path([coo_names])
            print("Updating scan trajectory file at {}".format(traj_path))
            cscn_save_fs[1].file.trajectory.write(traj, [coo_names])

    else:
        print("No cscan to save. (1) Skipping...")
