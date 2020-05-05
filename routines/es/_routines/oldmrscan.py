""" es_runners for coordinate scans
"""

import automol
import elstruct
import autofile
from routines.es._routines import _scan as scan
from routines.es._routines import _wfn as wfn
from routines.es import runner as es_runner
from lib import filesys


def run_multiref_rscan(ts_zma, ts_info, ts_formula, high_mul,
                       grid1, grid2, dist_name,
                       num_act_orb, num_act_elc,
                       mod_var_scn_thy_info,
                       scn_run_fs, scn_save_fs,
                       overwrite, update_guess=True,
                       constraint_dct=None,
                       **opt_kwargs):
    """ run constrained optimization scan
    """

    # Build the elstruct CASSCF options list used to build the wfn guess
    # (1) Build wfn with active space
    # (2) Build wfn with active space + 2 closed orbitals for stability
    cas_opt = []
    cas_opt.append(
        wfn.cas_options(
            ts_info, ts_formula, num_act_elc, num_act_orb,
            add_two_closed=False))
    cas_opt.append(
        wfn.cas_options(
            ts_info, ts_formula, num_act_elc, num_act_orb,
            add_two_closed=True))

    # Write the string that has all the components for building the wfn guess
    ref_zma = automol.zmatrix.set_values(ts_zma, {dist_name: grid1[0]})
    guess_str = wfn.multiref_wavefunction_guess(
        high_mul, ref_zma, ts_info, mod_var_scn_thy_info, cas_opt)
    guess_lines = guess_str.splitlines()

    # Set the opt script string and build the opt_kwargs
    [prog, method, _, _] = mod_var_scn_thy_info
    _, opt_script_str, _, opt_kwargs = es_runner.par.run_qchem_par(
        prog, method)
    opt_kwargs['casscf_options'] = cas_opt[1]
    opt_kwargs['gen_lines'] = {1: guess_lines}
    opt_kwargs['mol_options'] = ['nosym']  # Turn off symmetry

    # Build the filesystem for the scan
    scn_save_fs[1].create([dist_name])
    inf_obj = autofile.system.info.scan_branch(grid_dct)
    scn_save_fs[1].file.info.write(inf_obj, [coo_names])

    scan.run_twoway_scan(
        ts_zma, ts_info, mod_var_scn_thy_info,
        grid1_dct, scn_run_fs, scn_save_fs,
        opt_script_str, overwrite,
        update_guess=True,
        reverse_sweep=True,
        fix_failures=True,
        saddle=False,
        constraint_dct=None,
        **opt_kwargs
    )


def infinite_separation_energy(
        spc_1_info, spc_2_info, ts_info, high_mul, ref_zma,
        mod_var_scn_thy_info,
        mod_var_sp1_thy_info, mod_var_sp2_thy_info,
        hs_var_scn_thy_info,
        hs_var_sp1_thy_info,
        hs_var_sp2_thy_info,
        mod_ini_thy_info,
        geo, geo_run_path, geo_save_path,
        run_prefix, save_prefix,
        overwrite=False,
        num_act_orb=None, num_act_elc=None):
    """ Obtain the infinite separation energy from the multireference energy
        at a given reference point, the high-spin low-spin splitting at that
        reference point, and the high level energy for the high spin state
        at the reference geometry and for the fragments
        scn = thy for optimizations
        sp1 = low-spin single points
        sp2 = high-spin single points for inf sep
    """

    # Initialize infinite sep energy
    inf_sep_ene = -1.0e12

    # Build filesys for high-spin multi-ref case
    hs_run_fs, hs_var_run_path = filesys.build.high_spin_from_prefix(
        geo_run_path, hs_var_sp1_thy_info)
    hs_save_fs, hs_var_save_path = filesys.build.high_spin_from_prefix(
        geo_save_path, hs_var_sp1_thy_info)

    # Calculate the high-spin case
    if not hs_save_fs[-1].file.energy.exists(multi_lvl[1:4]) or overwrite:
        print(" - Running high spin multi reference energy ...")
        opt_script_str, _, opt_kwargs, _ = es_runner.par.run_qchem_par(
            multi_info[0], multi_info[1])

        # Build wfn guess
        ts_formula = automol.geom.formula(automol.zmatrix.geometry(ref_zma))
        cas_opt = wfn.cas_options(
            hs_info, ts_formula, num_act_elc, num_act_orb, high_mul)
        guess_str = wfn.multiref_wavefunction_guess(
            high_mul, ref_zma, hs_info, multi_lvl, [cas_opt])
        guess_lines = guess_str.splitlines()
        opt_kwargs['casscf_options'] = cas_opt
        opt_kwargs['mol_options'] = ['nosym']
        opt_kwargs['gen_lines'] = {1: guess_lines}

        # Run the jobs
        es_runner.run_job(
            job='energy',
            script_str=var_script_str,
            run_fs=opt_var_fs,
            geom=geo,
            spc_info=hs_info,
            thy_level=multi_lvl,
            overwrite=overwrite,
            **opt_kwargs,
        )

        ret = es_runner.read_job(
            job='energy',
            run_fs=run_var_fs,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            print(" - Reading high spin multi reference",
                  "energy from output...")
            hs_var_ene = elstruct.reader.energy(
                inf_obj.prog, inf_obj.method, out_str)

            print(" - Saving high spin multi reference energy...")
            print(" - Save path: {}".format(hs_var_save_path))
            hs_save_fs[-1].file.energy.write(hs_var_ene, multi_lvl[1:4])
            hs_save_fs[-1].file.input.write(inp_str, multi_lvl[1:4])
            hs_save_fs[-1].file.info.write(inf_obj, multi_lvl[1:4])
        else:
            print('ERROR: high spin multi reference energy job fails: ',
                  'Energy is needed to evaluate infinite separation energy')
            return inf_sep_ene

    else:
        hs_var_ene = hs_save_fs[-1].file.energy.read(multi_lvl[1:4])

    # Build the high-spin, single-reference filesystem
    hs_run_fs, hs_sr_run_path = filesys.build.high_spin_from_prefix(
        geo_run_path, hs_var_sp2_thy_info)
    hs_save_fs, hs_sr_save_path = filesys.build.high_spin_from_prefix(
        geo_save_path, hs_var_sp2_thy_info)
    run_sr_fs = autofile.fs.run(hs_sr_run_path)

    # Calculate energy for high-spin, single-reference level
    exists = hs_save_fs[-1].file.energy.exists(mod_var_sp2_thy_info[1:4])
    if not exists or overwrite:
        print(" - Running high spin single reference energy ...")

        sp_script_str, _, kwargs, _ = es_runner.par.run_qchem_par(
            *mod_var_sp2_thy_info[0:2])
        errors, options_mat = es_runner.par.set_molpro_options_mat(
            hs_info, geo)

        es_runner.run_job(
            job='energy',
            script_str=sp_script_str,
            run_fs=run_sr_fs,
            geom=geo,
            spc_info=hs_info,
            thy_level=mod_var_sp2_thy_info,
            errors=errors,
            options_mat=options_mat,
            overwrite=overwrite,
            **kwargs,
        )

        ret = es_runner.read_job(
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
            hs_save_fs[-1].file.energy.write(
                hs_sr_ene, mod_var_sp2_thy_info[1:4])
            hs_save_fs[-1].file.input.write(inp_str, mod_var_sp2_thy_info[1:4])
            hs_save_fs[-1].file.info.write(inf_obj, mod_var_sp2_thy_info[1:4])
        else:
            print('ERROR: High spin single reference energy job fails: ',
                  'Energy is needed to evaluate infinite separation energy')
            return inf_sep_ene

    else:
        hs_sr_ene = hs_save_fs[-1].file.energy.read(thy_lvl[1:4])

    # get the single reference energy for each of the reactant configurations
    spc_ene = []
    for spc_info in (spc_1_info, spc_2_info):

        # Build filesys for ini thy info
        _, ini_thy_run_path = filesys.build.spc_thy_fs_from_root(
            run_prefix, spc_info, mod_ini_thy_info)
        _, ini_thy_save_path = filesys.build.spc_thy_fs_from_root(
            save_prefix, spc_info, mod_ini_thy_info)

        # Build conformer filesys
        ini_cnf_run_fs, _ = filesys.build.cnf_fs_from_prefix(
            ini_thy_run_path, cnf=None)
        ini_cnf_save_fs, ini_cnf_save_locs = filesys.build.cnf_fs_from_prefix(
            ini_thy_save_path, cnf='min')
        ini_cnf_save_path = filesys.build.cnf_paths_from_locs(
            ini_cnf_save_fs, ini_cnf_save_locs)[0]
        ini_cnf_run_path = filesys.build.cnf_paths_from_locs(
            ini_cnf_run_fs, ini_cnf_save_locs)[0]

        # Read the geometry
        geo = ini_cnf_save_fs[-1].file.geometry.read(ini_cnf_save_locs)

        # Build the single point filesys objects
        sp_run_fs, sp_sr_run_path = filesys.build.sp_from_prefix(
            ini_cnf_run_path, mod_var_sp2_thy_info[1:4])
        sp_save_fs, sp_sr_save_path = filesys.build.sp_from_prefix(
            ini_cnf_save_path, mod_var_sp2_thy_info[1:4])
        run_sr_fs = autofile.fs.run(sp_sr_run_path)

        # Calculate the single point energy
        if not sp_save_fs[-1].file.energy.exists(thy_lvl[1:4]) or overwrite:
            print(" - Running single reference energy for",
                  "{} from output...".format(spc_info[0]))
            es_runner.run_job(
                job='energy',
                script_str=sp_script_str,
                run_fs=run_sr_fs,
                geom=geo,
                spc_info=spc_info,
                thy_level=thy_lvl,
                overwrite=overwrite,
                **kwargs,
            )

            ret = es_runner.read_job(
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
                sp_save_fs[-1].file.energy.write(sp_sr_ene, thy_lvl[1:4])
                sp_save_fs[-1].file.input.write(inp_str, thy_lvl[1:4])
                sp_save_fs[-1].file.info.write(inf_obj, thy_lvl[1:4])

            else:
                print('ERROR: Single reference energy job fails',
                      'for {}: '.format(spc_info[0]),
                      'Energy needed to evaluate infinite separation energy')
                return inf_sep_ene

        else:
            sp_sr_ene = sp_save_fs[-1].file.energy.read(thy_lvl[1:4])

        spc_ene.append(sp_sr_ene)

    inf_sep_ene = spc_ene[0] + spc_ene[1] - hs_sr_ene + hs_var_ene

    return inf_sep_ene
