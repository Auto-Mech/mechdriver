""" Specialized functions that run and save scan calculations
    used for variational calculations
"""

import automol
import autofile
from routines.es._routines import sp
from routines.es._routines import _scan as scan
from routines.es._routines import _wfn as wfn
from routines.es.runner import qchem_params
from lib import filesys


def multiref_rscan(ts_zma, ts_info, ts_formula, high_mul,
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
    _, opt_script_str, _, opt_kwargs = qchem_params(
        prog, method)
    opt_kwargs['casscf_options'] = cas_opt[1]
    opt_kwargs['gen_lines'] = {1: guess_lines}
    opt_kwargs['mol_options'] = ['nosym']  # Turn off symmetry

    # Build the filesystem for the scan
    full_grid = grid1 + grid2  # wrong
    scn_save_fs[1].create([dist_name])
    inf_obj = autofile.schema.info_objects.scan_branch({dist_name: full_grid})
    scn_save_fs[1].file.info.write(inf_obj, [dist_name])

    # Run the scans
    run_two_way_scan(
        ts_zma, ts_info, mod_var_scn_thy_info,
        grid1, grid2, dist_name,
        scn_run_fs, scn_save_fs,
        opt_script_str, overwrite,
        update_guess=update_guess,
        reverse_sweep=False,
        fix_failures=True,
        saddle=False,
        constraint_dct=constraint_dct,
        **opt_kwargs
    )


# DERIVED FUNCTION THAT RUNS RUN_SCAN AND SAVE IN TWO DIRECTIONS #
def run_two_way_scan(ts_zma, ts_info, mod_var_scn_thy_info,
                     grid1, grid2, dist_name,
                     scn_run_fs, scn_save_fs,
                     opt_script_str, overwrite,
                     update_guess=True,
                     reverse_sweep=True,
                     fix_failures=True,
                     saddle=False,
                     constraint_dct=None,
                     **opt_kwargs):
    """ Run a two-part scan that goes into two directions, as for rxn path
    """
    # Setup and run the first part of the scan to shorter distances
    grid1_dct = {dist_name: grid1}
    if len(grid1_dct) > 1:
        raise NotImplementedError
    print('grid1', grid1_dct)

    scan.run_scan(
        ts_zma, ts_info, mod_var_scn_thy_info,
        grid1_dct, scn_run_fs, scn_save_fs,
        opt_script_str, overwrite,
        update_guess=update_guess,
        reverse_sweep=reverse_sweep,
        fix_failures=fix_failures,
        saddle=saddle,
        constraint_dct=constraint_dct,
        **opt_kwargs
    )

    # Setup and run the sectiond part of the scan to farther distances
    grid2_dct = {dist_name: grid2}
    if len(grid2_dct) > 1:
        raise NotImplementedError
    print('grid2', grid2_dct)

    scan.run_scan(
        ts_zma, ts_info, mod_var_scn_thy_info,
        grid2_dct, scn_run_fs, scn_save_fs,
        opt_script_str, overwrite,
        update_guess=update_guess,
        reverse_sweep=reverse_sweep,
        fix_failures=fix_failures,
        saddle=saddle,
        constraint_dct=constraint_dct,
        **opt_kwargs
    )

    # Save the scan
    scan.save_scan(
        scn_run_fs=scn_run_fs,
        scn_save_fs=scn_save_fs,
        coo_names=[dist_name],
        thy_info=mod_var_scn_thy_info
    )


# CALCULATE INFINITE SEPARATION ENERGY #
def molrad_inf_sep_ene(spc1_info, spc2_info,
                       run_prefix, save_prefix,
                       mod_thy_info, mod_ini_thy_info,
                       overwrite):
    """ Calculate the inf spe ene for a mol-rad ene
    """
    sp_script_str, _, kwargs, _ = qchem_params(
        *mod_thy_info[0:2])
    inf_sep_ene = reac_sep_ene(
        spc1_info, spc2_info,
        run_prefix, save_prefix,
        mod_thy_info, mod_ini_thy_info,
        overwrite, sp_script_str,
        **kwargs)

    return inf_sep_ene


# def radrad_inf_sep_ene(
#         spc1_info, spc2_info, ts_info, high_mul, ref_zma,
#         mod_var_scn_thy_info,
#         mod_var_sp1_thy_info, mod_var_sp2_thy_info,
#         hs_var_scn_thy_info,
#         hs_var_sp1_thy_info,
#         hs_var_sp2_thy_info,
#         mod_ini_thy_info,
#         geo, geo_run_path, geo_save_path,
#         run_prefix, save_prefix,
#         overwrite=False,
#         num_act_orb=None, num_act_elc=None):
#     """ Obtain the infinite separation energy from the multireference energy
#         at a given reference point, the high-spin low-spin splitting at that
#         reference point, and the high level energy for the high spin state
#         at the reference geometry and for the fragments
#         scn = thy for optimizations
#         sp1 = low-spin single points
#         sp2 = high-spin single points for inf sep
#     """
#
#     # Initialize infinite sep energy
#     inf_sep_ene = -1.0e12
#
#     # Prepare filesys and guesses for the multi reference calc
#     hs_run_fs, hs_var_run_path = filesys.build.high_spin_from_prefix(
#         geo_run_path, hs_var_sp1_thy_info)
#     hs_save_fs, hs_var_save_path = filesys.build.high_spin_from_prefix(
#         geo_save_path, hs_var_sp1_thy_info)
#
#     opt_script_str, _, opt_kwargs, _ = qchem_params(
#         multi_info[0], multi_info[1])
#     ts_formula = automol.geom.formula(automol.zmatrix.geometry(ref_zma))
#     cas_opt = wfn.cas_options(
#         hs_info, ts_formula, num_act_elc, num_act_orb, high_mul)
#     guess_str = wfn.multiref_wavefunction_guess(
#         high_mul, ref_zma, hs_info, multi_lvl, [cas_opt])
#     guess_lines = guess_str.splitlines()
#     opt_kwargs['casscf_options'] = cas_opt
#     opt_kwargs['mol_options'] = ['nosym']
#     opt_kwargs['gen_lines'] = {1: guess_lines}
#
#     opt_script_str, _, opt_kwargs, _ = qchem_params(
#         multi_info[0], multi_info[1])
#
#     # Prepare filesys and guesses for the single reference calc
#     hs_run_fs, hs_sr_run_path = filesys.build.high_spin_from_prefix(
#         geo_run_path, hs_var_sp2_thy_info)
#     hs_save_fs, hs_sr_save_path = filesys.build.high_spin_from_prefix(
#         geo_save_path, hs_var_sp2_thy_info)
#     run_sr_fs = autofile.fs.run(hs_sr_run_path)
#
#     sp_script_str, _, kwargs, _ = qchem_params(
#         *mod_var_sp2_thy_info[0:2])
#     errors, options_mat = es_runner.par.set_molpro_options_mat(
#         hs_info, geo)
#
#     # Calculate the energies for the two cases
#     for x in s:
#        if not hs_save_fs[-1].file.energy.exists(multi_lvl[1:4]) or overwrite:
#             print(" - Running high spin multi reference energy ...")
#             opt_script_str, _, opt_kwargs, _ = qchem_params(
#                 multi_info[0], multi_info[1])
#
#             # Calculate the single point energy
#          if not sp_save_fs[-1].file.energy.exists(thy_lvl[1:4]) or overwrite:
#                 sp.run_energy(zma, geo, spc_info, thy_info,
#                               geo_save_fs, geo_run_path, geo_save_path, locs,
#                               script_str, overwrite, **kwargs)
#                 if not sp_save_fs[-1].file.energy.exists(thy_lvl[1:4]):
#               print('ERROR: High spin single reference energy job fails: ',
#                   'Energy is needed to evaluate infinite separation energy')
#                     hs_sr_ene = None
#                 else:
#               print(" - Reading high spin single ref energy from output...")
#                     hs_sr_ene = sp_save_fs[-1].file.energy.read(thy_lvl[1:4])
#             else:
#                 hs_sr_ene = sp_save_fs[-1].file.energy.read(thy_lvl[1:4])
#
#         else:
#             hs_var_ene = hs_save_fs[-1].file.energy.read(multi_lvl[1:4])
#
#     # get the single reference energy for each of the reactant configurations
#     reac1_ene, reac2_ene = reac_sep_ene(
#         spc1_info, spc2_info,
#         run_prefix, save_prefix,
#         mod_var_sp2_thy_info, mod_ini_thy_info,
#         overwrite, sp_script_str,
#         **kwargs)
#
#     # Calculate the infinite seperation energy
#     all_enes = (reac1_ene, reac2_ene, hs_sr_ene, hs_var_ene)
#     if all(ene is not None for ene in all_enes):
#         inf_sep_ene = reac1_ene + reac2_ene - hs_sr_ene + hs_var_ene
#     else:
#         inf_sep_ene = None
#
#     return inf_sep_ene


def reac_sep_ene(spc1_info, spc2_info,
                 run_prefix, save_prefix,
                 thy_info, ini_thy_info,
                 overwrite, sp_script_str,
                 **kwargs):
    """ Calculate the sum of two reactants
    """

    # get the single reference energy for each of the reactant configurations
    spc_enes = []
    for spc_info in (spc1_info, spc2_info):
        print('spc_info', spc_info)

        # Modify thy info for species
        mod_thy_info = filesys.inf.modify_orb_restrict(spc_info, thy_info)
        mod_ini_thy_info = filesys.inf.modify_orb_restrict(
            spc_info, ini_thy_info)

        # Build filesys for ini thy info
        _, ini_thy_run_path = filesys.build.spc_thy_fs_from_root(
            run_prefix, spc_info, mod_ini_thy_info)
        _, ini_thy_save_path = filesys.build.spc_thy_fs_from_root(
            save_prefix, spc_info, mod_ini_thy_info)

        # Build conformer filesys
        ini_cnf_run_fs, _ = filesys.build.cnf_fs_from_prefix(
            ini_thy_run_path, mod_ini_thy_info, cnf=None)
        ini_cnf_save_fs, ini_cnf_save_locs = filesys.build.cnf_fs_from_prefix(
            ini_thy_save_path, mod_ini_thy_info, cnf='min')
        ini_cnf_save_paths = filesys.build.cnf_paths_from_locs(
            ini_cnf_save_fs, ini_cnf_save_locs)

        # Read the geometry and set paths
        zma = ini_cnf_save_fs[-1].file.zmatrix.read(ini_cnf_save_locs)
        geo = ini_cnf_save_fs[-1].file.geometry.read(ini_cnf_save_locs)
        geo_run_path = ini_cnf_run_fs[-1].path(ini_cnf_save_locs)
        geo_save_path = ini_cnf_save_fs[-1].path(ini_cnf_save_locs)

        # Build the single point filesys objects
        sp_save_fs, _ = filesys.build.sp_from_prefix(
            ini_cnf_save_paths[0], mod_thy_info)

        # Calculate the save single point energy
        sp.run_energy(zma, geo, spc_info, mod_thy_info,
                      ini_cnf_save_fs, geo_run_path, geo_save_path,
                      ini_cnf_save_locs,
                      sp_script_str, overwrite, **kwargs)
        exists = sp_save_fs[-1].file.energy.exists(mod_thy_info[1:4])
        if not exists:
            print('No ene found')
            ene = None
        else:
            print('Reading energy')
            ene = sp_save_fs[-1].file.energy.read(mod_thy_info[1:4])

        # Append ene to list
        spc_enes.append(ene)

    # Analyze the energies in the list
    inf_ene = 0.0
    for i, ene in enumerate(spc_enes):
        if ene is not None:
            inf_ene += ene
        else:
            print('ERROR: Single reference energy job fails',
                  'for {}: '.format(spc_info[i]),
                  'Energy needed to evaluate infinite separation energy')
            inf_ene = None

    return inf_ene
