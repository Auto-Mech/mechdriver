""" drivers for coordinate scans
"""
import numpy
from qcelemental import constants as qcc
import automol
import elstruct
import autofile
import moldr

WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')


def hindered_rotor_scans(
        spc_info, thy_level, cnf_run_fs, cnf_save_fs,
        script_str, overwrite, scan_increment=30., saddle=False, tors_names='', new_grid=False, **opt_kwargs):
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
                zma, tors_names, scan_increment)
            tors_grids = [
                numpy.linspace(*linspace) + val_dct[name]
                for name, linspace in zip(tors_names, tors_linspaces)]

            for tors_name, tors_grid in zip(tors_names, tors_grids):
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
                    new_grid=new_grid,
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
        new_grid = True, **kwargs):
    """ run constrained optimization scan
    """
#    if len(grid_dct) > 1:
#        raise NotImplementedError

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
#    ((coo_name, grid_vals),) = grid_dct.items()
#    if scn_save_fs.branch.file.info.exists([[coo_name]]):
#        if not new_grid:
#            inf_obj = scn_save_fs.branch.file.info.read([[coo_name]])
#            existing_grid_dct = dict(inf_obj.grids)
#            existing_grid_vals = existing_grid_dct[coo_name]
#            print('grid vals test:', grid_vals, existing_grid_vals)
#            assert (numpy.shape(grid_vals) == numpy.shape(existing_grid_vals) and
#                    (numpy.allclose(grid_vals*180.0/numpy.pi, existing_grid_vals) or
#                     numpy.allclose(grid_vals, existing_grid_vals)))
#
    #inf_obj = autofile.system.info.scan_branch({coo_name: grid_vals})
    #scn_save_fs.branch.file.info.write(inf_obj, [[coo_name]])
    print('coo names test:', coo_names)
    scn_save_fs.branch.create([coo_names])
    inf_obj = autofile.system.info.scan_branch(grid_dct)
    scn_save_fs.branch.file.info.write(inf_obj, [coo_names])
    npoint = 1
    for coo_grid_vals in grid_vals:
        npoint *= len(coo_grid_vals)
    grid_idxs = tuple(range(npoint))
   # grid_names = ['{:.3f}'.format(val) for val in grid_vals]

    #for grid_idx in grid_idxs:
    #    scn_run_fs.leaf.create([[coo_name], [grid_idx]])
   # for grid_name in grid_names:
   #     scn_run_fs.leaf.create([[coo_name], [grid_name]])

   # prefixes = tuple(scn_run_fs.leaf.path([[coo_name], [grid_idx]])
   #                  for grid_idx in grid_idxs)
    #prefixes = tuple(scn_run_fs.leaf.path([[coo_name], [grid_name]])
    #                 for grid_name in grid_names)
    if len(grid_vals) == 1:
        for grid_val in grid_vals[0]:
            scn_run_fs.leaf.create([coo_names, [grid_val]])
        prefixes = tuple(scn_run_fs.leaf.path([coo_names, [grid_val]])
                        for grid_val in grid_vals[0])
        _run_1d_scan(
            script_str=script_str,
            prefixes=prefixes,
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
                prefixes=list(reversed(prefixes)),
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
        prefixes = []
        for grid_val_i in grid_vals[0]:
            for grid_val_j in grid_vals[1]:
                scn_run_fs.leaf.create([coo_names, [grid_val_i, grid_val_j]])
                prefixes.append(scn_run_fs.leaf.path([coo_names, [grid_val_i, grid_val_j]]))
        prefixes = tuple(prefixes)                

        _run_2d_scan(
            script_str=script_str,
            prefixes=prefixes,
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
            prefixes = []
            for grid_val_i in grid_vals[0][::-1]:
                for grid_val_j in grid_vals[1][::-1]:
                    prefixes.append(scn_run_fs.leaf.path([coo_names, [grid_val_i, grid_val_j]]))
            _run_2d_scan(
                script_str=script_str,
                prefixes=prefixes,
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


def infinite_separation_energy(
        rct_zmas, rcts_spc_info, 
        frag_geos, ref_sep_geo, ref_sep_ene, ts_mul, high_mul, frag_muls,
        thy_info, multi_thy_info, 
        scn_run_fs, scn_save_fs, script_str, overwrite, 
        **kwargs):
    """ Obtain the infinite separation energy from the multireference energy at a given 
    reference point, the high-spin low-spin splitting at that reference point, and the 
    high level energy for the high spin state at the reference geometry and for the fragments
    """

    delta_spin = multi_high_spin_energy - mult_low_spin_energy
    ene_thy_geo1 = 
    ene_thy_geo2 = 
    ene_thy_ref_high_spin
    ene_inf = ene_ref + delta_spin + ene_thy_geo1 + ene_thy_geo2 - ene_thy_ref_high_spin

#    for spc in enumerate(rcts_spc_info):

#        spc_run_fs = autofile.fs.species(run_prefix)
#        spc_run_fs.leaf.create(spc_info)
#        spc_run_path = spc_run_fs.leaf.path(spc_info)
#        spc_save_fs = autofile.fs.species(save_prefix)
#        spc_save_fs.leaf.create(spc_info)
#        spc_save_path = spc_save_fs.leaf.path(spc_info)
#
#        orb_restr = moldr.util.orbital_restriction(
#            spc_info, thy_info)
#        thy_lvl = thy_info[0:3]
#        thy_lvl.append(orb_restr)
#        thy_run_fs = autofile.fs.theory(spc_run_path)
#        thy_run_fs.leaf.create(thy_lvl)
#        thy_run_path = thy_run_fs.leaf.path(thy_lvl)
#        thy_save_fs = autofile.fs.theory(spc_save_path)
#        thy_save_fs.leaf.create(thy_lvl)
#        thy_save_path = thy_save_fs.leaf.path(thy_lvl)
#        cnf_run_fs = autofile.fs.conformer(thy_run_path)
#        cnf_save_fs = autofile.fs.conformer(thy_save_path)
#        min_cnf_locs = moldr.util.min_energy_conformer_locators(ini_cnf_save_fs)
#        min_cnf_run_path = ini_cnf_run_fs.leaf.path(min_cnf_locs)
#        min_cnf_save_path = ini_cnf_save_fs.leaf.path(min_cnf_locs)


#        for frag_geo in frag_geos
#            sp.run_energy(frag_spc_info, thy_level, frag_geo_run_fs, frag_geo_save_fs, locs, script_str, overwrite)
#            sp_save_fs = autofile.fs.single_point(geo_save_path)
#            sp_save_fs.leaf.create(thy_level[1:4])
#            sp_save_fs.leaf.file.energy.write(ene, thy_level[1:4])


def run_multiref_rscan(
        formula, high_mul, zma, spc_info, thy_level, dist_name, grid1, grid2,
        scn_run_fs, scn_save_fs, script_str, overwrite, update_guess=True,
        **kwargs):
    """ run constrained optimization scan
    """

    #print('grid1:', grid1)
    vma = automol.zmatrix.var_(zma)
    if scn_save_fs.trunk.file.vmatrix.exists():
        existing_vma = scn_save_fs.trunk.file.vmatrix.read()
        assert vma == existing_vma

    grid = numpy.append(grid1, grid2)
    grid_dct = {dist_name: grid}
    if len(grid_dct) > 1:
        raise NotImplementedError

    ((coo_name, grid_vals),) = grid_dct.items()
    scn_save_fs.branch.create([[coo_name]])
    print('scan_save_fs branch test:', scn_save_fs.trunk.path())
    if scn_save_fs.branch.file.info.exists([[coo_name]]):
        inf_obj = scn_save_fs.branch.file.info.read([[coo_name]])
        existing_grid_dct = dict(inf_obj.grids)
        existing_grid_vals = existing_grid_dct[coo_name]
        assert (numpy.shape(grid_vals) == numpy.shape(existing_grid_vals) and
                numpy.allclose(grid_vals, existing_grid_vals))

    inf_obj = autofile.system.info.scan_branch({coo_name: grid_vals})
    scn_save_fs.branch.file.info.write(inf_obj, [[coo_name]])
    charge = spc_info[1]
    mul = spc_info[2]
    basis = thy_level[2]
    prog = thy_level[0]
    #orb_restr = moldr.util.orbital_restriction(spc_info, thy_level)
    thy_level[0] = 'molpro2015'
    #print('thy_level test:', thy_level)
    thy_level[1] = 'caspt2'
    prog = 'molpro2015'
    method = 'caspt2'

    _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(prog, method)

    num_act_elc = 2
    num_act_orb = 2
    cas_opt, _ = moldr.ts.cas_options(spc_info, formula, num_act_elc, num_act_orb, high_mul)
    guess_str = moldr.ts.multiref_wavefunction_guess(
        formula, high_mul, zma, spc_info, thy_level, dist_name, coo_name, grid_vals, cas_opt)
    guess_lines = guess_str.splitlines()

    #print('guess_str test:', guess_str)
    opt_kwargs['casscf_options'] = cas_opt
    opt_kwargs['mol_options'] = ['nosym']
    opt_kwargs['gen_lines'] = guess_lines

    grid1_dct = {dist_name: grid1}
    ((_, grid1_vals),) = grid1_dct.items()
    #print('grid1_vals test:', grid1_vals)
    npoint1 = len(grid1_vals)
    grid1_idxs = tuple(range(npoint1))

    for grid_idx in grid1_idxs:
        scn_run_fs.leaf.create([[coo_name], [grid_idx]])
    prefixes = tuple(scn_run_fs.leaf.path([[coo_name], [grid_idx]])
                     for grid_idx in grid1_idxs)
    #print('theory level:', thy_level)
    _run_1d_scan(
        script_str=opt_script_str,
        prefixes=prefixes,
        guess_zma=zma,
        coo_name=coo_name,
        grid_idxs=grid1_idxs,
        grid_vals=grid1_vals,
        spc_info=spc_info,
        thy_level=thy_level,
        overwrite=overwrite,
        update_guess=update_guess,
        **opt_kwargs,
    )

    grid2_dct = {dist_name: grid2}
    ((_, grid2_vals),) = grid2_dct.items()
    npoint2 = len(grid2_vals)
    grid2_idxs = tuple(range(npoint2))
    grid2_idxs = tuple(x + npoint1 for x in grid2_idxs)

    for grid_idx in grid2_idxs:
        scn_run_fs.leaf.create([[coo_name], [grid_idx]])

    prefixes = tuple(scn_run_fs.leaf.path([[coo_name], [grid_idx]])
                     for grid_idx in grid2_idxs)
    _run_1d_scan(
        script_str=opt_script_str,
        prefixes=prefixes,
        guess_zma=zma,
        coo_name=coo_name,
        grid_idxs=grid2_idxs,
        grid_vals=grid2_vals,
        spc_info=spc_info,
        thy_level=thy_level,
        overwrite=overwrite,
        update_guess=update_guess,
        **opt_kwargs,
    )


def _run_1d_scan(
        script_str, prefixes, guess_zma, coo_name, grid_idxs, grid_vals, 
        spc_info, thy_level, overwrite, errors=(),
        options_mat=(), retry_failed=True, update_guess=True, saddle=False, **kwargs):

    npoints = len(grid_idxs)
    assert len(grid_vals) == len(prefixes) == npoints
    for grid_idx, grid_val, prefix in zip(grid_idxs, grid_vals, prefixes):
        print("Point {}/{}".format(grid_idx+1, npoints))
        zma = automol.zmatrix.set_values(guess_zma, {coo_name: grid_val})
        run_fs = autofile.fs.run(prefix)

        #print('kwargs test', kwargs)
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
        #print('update_guess and ret test:', update_guess, ret)
        if update_guess and ret is not None:
            inf_obj, _, out_str = ret
            prog = inf_obj.prog
            guess_zma = elstruct.reader.opt_zmatrix(prog, out_str)


def _run_2d_scan(
        script_str, prefixes, guess_zma, coo_names, grid_idxs, grid_vals, 
        spc_info, thy_level, overwrite, errors=(),
        options_mat=(), retry_failed=True, update_guess=True, saddle=False, **kwargs):

    npoints = len(grid_idxs)
    assert len(grid_vals[0]) * len(grid_vals[1]) == len(prefixes) == npoints

    idx = 0
    for grid_val_i in grid_vals[0]:
        for grid_val_j in grid_vals[1]:
            grid_idx = grid_idxs[idx]
            prefix = prefixes[idx]
            print("Point {}/{}".format(grid_idx+1, npoints))
            zma = automol.zmatrix.set_values(guess_zma, {coo_names[0]: grid_val_i, coo_names[1]: grid_val_j})
            run_fs = autofile.fs.run(prefix)
            idx += 1        

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

def save_scan(scn_run_fs, scn_save_fs, coo_names):
    """ save the scans that have been run so far
    """
    if not scn_run_fs.branch.exists([coo_names]):
        print("No scan to save. Skipping...")
    else:
        locs_lst = []
        for locs in scn_run_fs.leaf.existing([coo_names]):
            if '.' not in locs:
                continue
            print('locs test:', locs)
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


