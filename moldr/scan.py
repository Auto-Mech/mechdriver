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
        script_str, overwrite, scan_increment=30., **opt_kwargs):
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
        tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
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
                **opt_kwargs,
            )

            print('min_cnf_save_path in hindered_rotor_scan')
            print(min_cnf_save_path)

            save_scan(
                scn_run_fs=scn_run_fs,
                scn_save_fs=scn_save_fs,
                coo_names=[tors_name],
            )


def run_scan(
        zma, spc_info, thy_level, grid_dct, scn_run_fs, scn_save_fs,
        script_str, overwrite, update_guess=True,
        reverse_sweep=True, **kwargs):
    """ run constrained optimization scan
    """
    if len(grid_dct) > 1:
        raise NotImplementedError

#    scn_run_fs = autofile.fs.scan(run_prefix)
#    scn_save_fs = autofile.fs.scan(save_prefix)

    vma = automol.zmatrix.var_(zma)
    if scn_save_fs.trunk.file.vmatrix.exists():
        existing_vma = scn_save_fs.trunk.file.vmatrix.read()
        assert vma == existing_vma

    # for now, running only one-dimensional hindered rotor scans
    ((coo_name, grid_vals),) = grid_dct.items()
    scn_save_fs.branch.create([[coo_name]])

    if scn_save_fs.branch.file.info.exists([[coo_name]]):
        inf_obj = scn_save_fs.branch.file.info.read([[coo_name]])
        existing_grid_dct = dict(inf_obj.grids)
        existing_grid_vals = existing_grid_dct[coo_name]
#        print('grid test:', grid_vals, existing_grid_vals)
        assert (numpy.shape(grid_vals) == numpy.shape(existing_grid_vals) and
                (numpy.allclose(grid_vals*180.0/numpy.pi, existing_grid_vals) or
                 numpy.allclose(grid_vals, existing_grid_vals)))

    inf_obj = autofile.system.info.scan_branch({coo_name: grid_vals})
    scn_save_fs.branch.file.info.write(inf_obj, [[coo_name]])

    npoint = len(grid_vals)
    grid_idxs = tuple(range(npoint))

    for grid_idx in grid_idxs:
        scn_run_fs.leaf.create([[coo_name], [grid_idx]])

    prefixes = tuple(scn_run_fs.leaf.path([[coo_name], [grid_idx]])
                     for grid_idx in grid_idxs)

    _run_1d_scan(
        script_str=script_str,
        prefixes=prefixes,
        guess_zma=zma,
        coo_name=coo_name,
        grid_idxs=grid_idxs,
        grid_vals=grid_vals,
        spc_info=spc_info,
        thy_level=thy_level,
        overwrite=overwrite,
        update_guess=update_guess,
        **kwargs
    )

    if reverse_sweep:
        _run_1d_scan(
            script_str=script_str,
            prefixes=list(reversed(prefixes)),
            guess_zma=zma,
            coo_name=coo_name,
            grid_idxs=list(reversed(grid_idxs)),
            grid_vals=list(reversed(grid_vals)),
            spc_info=spc_info,
            thy_level=thy_level,
            overwrite=overwrite,
            update_guess=update_guess,
            **kwargs
        )


def run_multiref_rscan(
        formula, high_mul, zma, spc_info, thy_level, dist_name, grid1, grid2,
        run_prefix, save_prefix, script_str, overwrite, update_guess=True,
        **kwargs):
    """ run constrained optimization scan
    """

    electron_count = automol.formula._formula.electron_count(formula)
    # this is only for 2e,2o case
    closed_orb = electron_count//2 - 1
    occ_orb = electron_count//2 + 1
    # end of 2e,2o case
    two_spin = spc_info[2]-1
    chg = spc_info[1]
    cas_options = [
        elstruct.option.specify(elstruct.Option.Scf.MAXITER_, 40),
        elstruct.option.specify(elstruct.Option.Casscf.OCC_, occ_orb),
        elstruct.option.specify(elstruct.Option.Casscf.CLOSED_, closed_orb),
        elstruct.option.specify(elstruct.Option.Casscf.WFN_, electron_count, 1, two_spin, chg)
        ]

    scn_run_fs = autofile.fs.scan(run_prefix)
    scn_save_fs = autofile.fs.scan(save_prefix)

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
    if scn_save_fs.branch.file.info.exists([[coo_name]]):
        inf_obj = scn_save_fs.branch.file.info.read([[coo_name]])
        existing_grid_dct = dict(inf_obj.grids)
        existing_grid_vals = existing_grid_dct[coo_name]
#        print('grid_vals test:', grid_vals, existing_grid_vals)
        assert (numpy.shape(grid_vals) == numpy.shape(existing_grid_vals) and
                numpy.allclose(grid_vals, existing_grid_vals))

    inf_obj = autofile.system.info.scan_branch({coo_name: grid_vals})
    scn_save_fs.branch.file.info.write(inf_obj, [[coo_name]])
    charge = spc_info[1]
    mul = spc_info[2]
    basis = thy_level[2]
    prog = thy_level[0]
    orb_restr = moldr.util.orbital_restriction(spc_info, thy_level)
    #orb_restr = thy_level[3]
    thy_level[0] = 'molpro'
    prog = 'molpro'
    thy_level[1] = 'caspt2'
    _, script_str, _, kwargs = moldr.util.run_qchem_par(thy_level[0], thy_level[1])
#    print('orb_restr test:', orb_restr)

    guess_str1 = elstruct.writer.energy(
        geom=automol.zmatrix.set_values(zma, {coo_name: grid_vals[0]}),
        charge=charge,
        mult=high_mul,
        method='hf',
        basis=basis,
        prog=prog,
        orb_restricted=orb_restr,
        mol_options=['nosym'],
        )
    guess_str1 += '\n\n'
    guess_str1 = '\n'.join(guess_str1.splitlines()[2:])

    guess_str2 = elstruct.writer.energy(
        geom=automol.zmatrix.set_values(zma, {coo_name: grid_vals[0]}),
        charge=charge,
        mult=mul,
        method='casscf',
        basis=basis,
        prog=prog,
        orb_restricted=orb_restr,
        casscf_options=cas_options,
        mol_options=['nosym'],
        )
    guess_str2 += '\n\n'
    guess_str2 = '\n'.join(guess_str2.splitlines()[2:])

    guess_str = guess_str1 + guess_str2
    guess_lines = guess_str.splitlines()
    kwargs['casscf_options'] = cas_options
    kwargs['mol_options'] = ['nosym']
    kwargs['gen_lines'] = guess_lines

    grid1_dct = {dist_name: grid1}
    ((_, grid1_vals),) = grid1_dct.items()
    npoint1 = len(grid1_vals)
    grid1_idxs = tuple(range(npoint1))

    for grid_idx in grid1_idxs:
        scn_run_fs.leaf.create([[coo_name], [grid_idx]])
    prefixes = tuple(scn_run_fs.leaf.path([[coo_name], [grid_idx]])
                     for grid_idx in grid1_idxs)
#    print('update_guess0 test:', update_guess)
    _run_1d_scan(
        script_str=script_str,
        prefixes=prefixes,
        guess_zma=zma,
        coo_name=coo_name,
        grid_idxs=grid1_idxs,
        grid_vals=grid1_vals,
        spc_info=spc_info,
        thy_level=thy_level,
        overwrite=overwrite,
        update_guess=update_guess,
        **kwargs
    )

    # grid2 = numpy.append(grid1[0], grid2)
#    print('grid2 test in multi:', grid2)
    grid2_dct = {dist_name: grid2}
    ((_, grid2_vals),) = grid2_dct.items()
    npoint2 = len(grid2_vals)
    grid2_idxs = tuple(range(npoint2))
    grid2_idxs = tuple(x + npoint1 for x in grid2_idxs)
#    print('grid2_idxs test:', grid2_idxs)

    for grid_idx in grid2_idxs:
        scn_run_fs.leaf.create([[coo_name], [grid_idx]])

    prefixes = tuple(scn_run_fs.leaf.path([[coo_name], [grid_idx]])
                     for grid_idx in grid2_idxs)
#    print('update_guess0 test:', update_guess)
    _run_1d_scan(
        script_str=script_str,
        prefixes=prefixes,
        guess_zma=zma,
        coo_name=coo_name,
        grid_idxs=grid2_idxs,
        grid_vals=grid2_vals,
        spc_info=spc_info,
        thy_level=thy_level,
        overwrite=overwrite,
        update_guess=update_guess,
        **kwargs
    )


def _run_1d_scan(
        script_str, prefixes, guess_zma, coo_name, grid_idxs, grid_vals,
        spc_info, thy_level, overwrite, errors=(),
        options_mat=(), retry_failed=True, update_guess=True, **kwargs):

#    print('prefixes test:', prefixes)
    npoints = len(grid_idxs)
    assert len(grid_vals) == len(prefixes) == npoints
    for grid_idx, grid_val, prefix in zip(grid_idxs, grid_vals, prefixes):
        print("Point {}/{}".format(grid_idx+1, npoints))
        zma = automol.zmatrix.set_values(guess_zma, {coo_name: grid_val})
#        print('zma test:', zma)
        run_fs = autofile.fs.run(prefix)

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
            **kwargs
        )

        ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
#        print('prefix test:', prefix)
        # print('ret test:', ret)
        if update_guess and ret is not None:
            inf_obj, _, out_str = ret
            prog = inf_obj.prog
            guess_zma = elstruct.reader.opt_zmatrix(prog, out_str)
#            print('guess_zma test:', guess_zma)


def save_scan(scn_run_fs, scn_save_fs, coo_names):
    """ save the scans that have been run so far
    """
    print(coo_names)
    if len(coo_names) > 1:
        raise NotImplementedError

    if not scn_run_fs.branch.exists([coo_names]):
        print("No scan to save. Skipping...")
    else:
        locs_lst = []
        for locs in scn_run_fs.leaf.existing([coo_names]):
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


