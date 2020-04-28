""" es_runners for coordinate scans
"""

# import numpy
# import automol
# import elstruct
# import autofile
# from routines.es import _run_1d_scan
# from routines.es import _wfn as wfn
# from routines.es import runner as es_runner
# from lib import filesys
#
#
# def run_multiref_rscan(
#       formula, high_mul, zma, spc_info, multi_level, dist_name, grid1, grid2,
#         scn_run_fs, scn_save_fs, overwrite, update_guess=True,
#         num_act_elc=None, num_act_orb=None,
#         constraint_dct=None):
#     """ run constrained optimization scan
#     """
#
#     vma = automol.zmatrix.var_(zma)
#     if scn_save_fs[0].file.vmatrix.exists():
#         existing_vma = scn_save_fs[0].file.vmatrix.read()
#         assert vma == existing_vma
#
#     grid = numpy.append(grid1, grid2)
#     grid_dct = {dist_name: grid}
#     if len(grid_dct) > 1:
#         raise NotImplementedError
#
#     coo_names = []
#     grid_vals = []
#     for item in grid_dct.items():
#         (coo, coo_grid_vals) = item
#         coo_names.append(coo)
#         grid_vals.append(coo_grid_vals)
#
#     scn_save_fs[1].create([coo_names])
#     inf_obj = autofile.system.info.scan_branch(grid_dct)
#     scn_save_fs[1].file.info.write(inf_obj, [coo_names])
#
#     prog = multi_level[0]
#     method = multi_level[1]
#
#     _, opt_script_str, _, opt_kwargs = es_runner.par.run_qchem_par(
#         prog, method)
#
#    ref_zma = automol.zmatrix.set_values(zma, {coo_names[0]: grid_vals[0][0]})
#
#     # Build the elstruct CASSCF options list for multiref calcs
#     cas_opt = []
#     cas_opt.append(
#         wfn.cas_options(
#             spc_info, formula, num_act_elc, num_act_orb,
#             high_mul, add_two_closed=False))
#     cas_opt.append(
#         wfn.cas_options(
#             spc_info, formula, num_act_elc, num_act_orb,
#             high_mul, add_two_closed=True))
#
#     # Write the lines containing all the calcs for a guess wfn
#     guess_str = wfn.multiref_wavefunction_guess(
#         high_mul, ref_zma, spc_info, multi_level, cas_opt)
#     guess_lines = guess_str.splitlines()
#
#     # Add the above-built objects to the elstruct opt_kwargs dct
#     opt_kwargs['casscf_options'] = cas_opt[1]
#     opt_kwargs['gen_lines'] = {1: guess_lines}
#
#     # Add option to opt_kwargs dct to turn off symmetry
#     opt_kwargs['mol_options'] = ['nosym']
#
#     # Setup and run the first part of the scan to shorter distances
#     coo_names = []
#     grid1_vals = []
#     grid1_dct = {dist_name: grid1}
#     for item in grid1_dct.items():
#         (coo, coo_grid1_vals) = item
#         coo_names.append(coo)
#         grid1_vals.append(coo_grid1_vals)
#
#     npoint = 1
#     for coo_grid1_vals in grid1_vals:
#         npoint *= len(coo_grid1_vals)
#     grid1_idxs = tuple(range(npoint))
#     if len(grid1_vals) == 1:
#         if constraint_dct is not None:
#             for grid1_val in grid1_vals[0]:
#              scn_run_fs[-1].create([coo_names, [grid1_val], constraint_dct])
#             run_prefixes = tuple(
#                 scn_run_fs[-1].path([coo_names, [grid1_val], constraint_dct])
#                 for grid1_val in grid1_vals[0])
#         else:
#             for grid1_val in grid1_vals[0]:
#                 scn_run_fs[-1].create([coo_names, [grid1_val]])
#             run_prefixes = tuple(
#                 scn_run_fs[-1].path([coo_names, [grid1_val]])
#                 for grid1_val in grid1_vals[0])
#     _run_1d_scan(
#         script_str=opt_script_str,
#         run_prefixes=run_prefixes,
#         scn_save_fs=scn_save_fs,
#         guess_zma=zma,
#         coo_name=coo_names[0],
#         grid_idxs=grid1_idxs,
#         grid_vals=grid1_vals[0],
#         spc_info=spc_info,
#         thy_level=multi_level,
#         overwrite=overwrite,
#         update_guess=update_guess,
#         constraint_dct=constraint_dct,
#         **opt_kwargs,
#     )
#
#     # Setup and run the second part of the scan to longer distances
#     coo_names = []
#     grid2_vals = []
#     grid2_dct = {dist_name: grid2}
#     for item in grid2_dct.items():
#         (coo, coo_grid2_vals) = item
#         coo_names.append(coo)
#         grid2_vals.append(coo_grid2_vals)
#
#     npoint = 1
#     for coo_grid2_vals in grid2_vals:
#         npoint *= len(coo_grid2_vals)
#     grid2_idxs = tuple(range(npoint))
#     if len(grid2_vals) == 1:
#         if constraint_dct is not None:
#             for grid2_val in grid2_vals[0]:
#                 scn_run_fs[-1].create(
#                     [coo_names, [grid2_val], constraint_dct])
#             run_prefixes = tuple(
#                 scn_run_fs[-1].path([coo_names, [grid2_val], constraint_dct])
#                 for grid2_val in grid2_vals[0])
#         else:
#             for grid2_val in grid2_vals[0]:
#                 scn_run_fs[-1].create(
#                     [coo_names, [grid2_val]])
#             run_prefixes = tuple(
#                 scn_run_fs[-1].path([coo_names, [grid2_val]])
#                 for grid2_val in grid2_vals[0])
#     _run_1d_scan(
#         script_str=opt_script_str,
#         run_prefixes=run_prefixes,
#         scn_save_fs=scn_save_fs,
#         guess_zma=zma,
#         coo_name=coo_names[0],
#         grid_idxs=grid2_idxs,
#         grid_vals=grid2_vals[0],
#         spc_info=spc_info,
#         thy_level=multi_level,
#         overwrite=overwrite,
#         update_guess=update_guess,
#         constraint_dct=constraint_dct,
#         **opt_kwargs,
#     )
#
#
# def infinite_separation_energy(
#         spc_1_info, spc_2_info, ts_info, high_mul, ref_zma, ini_thy_info,
#         thy_info,
#         multi_info, run_prefix, save_prefix, scn_run_fs, scn_save_fs, locs,
#         overwrite=False,
#         num_act_elc=None, num_act_orb=None):
#     """ Obtain the infinite separation energy from the multireference energy
#         at a given reference point, the high-spin low-spin splitting at that
#         reference point, and the high level energy for the high spin state
#         at the reference geometry and for the fragments
#     """
#
#     # Initialize infinite sep energy
#     inf_sep_ene = -1.0e12
#
#     # set up all the file systems for the TS
#     # start with the geo and reference theory info
#     geo_run_path = scn_run_fs[-1].path(locs)
#     geo_save_path = scn_save_fs[-1].path(locs)
#     geo = scn_save_fs[-1].file.geometry.read(locs)
#     sp_run_fs = autofile.fs.single_point(geo_run_path)
#     sp_save_fs = autofile.fs.single_point(geo_save_path)
#
#    # get the multi reference energy for high spin state for ref point on scan
#     hs_info = (ts_info[0], ts_info[1], high_mul)
#     orb_restr = filesys.inf.orbital_restriction(hs_info, multi_info)
#     multi_lvl = multi_info[0:3]
#     multi_lvl.append(orb_restr)
#
#     hs_run_fs = autofile.fs.high_spin(geo_run_path)
#     hs_save_fs = autofile.fs.high_spin(geo_save_path)
#     hs_run_fs[-1].create(multi_lvl[1:4])
#     hs_save_fs[-1].create(multi_lvl[1:4])
#
#     hs_mr_run_path = hs_run_fs[-1].path(multi_lvl[1:4])
#     hs_mr_save_path = hs_save_fs[-1].path(multi_lvl[1:4])
#     run_mr_fs = autofile.fs.run(hs_mr_run_path)
#
#     mr_script_str, _, mr_kwargs, _ = es_runner.par.run_qchem_par(
#         multi_info[0], multi_info[1])
#
#     if num_act_elc is None and num_act_orb is None:
#         num_act_elc = high_mul
#         num_act_orb = num_act_elc
#     ts_formula = automol.geom.formula(automol.zmatrix.geometry(ref_zma))
#
#     cas_opt = wfn.cas_options(
#         hs_info, ts_formula, num_act_elc, num_act_orb, high_mul)
#     guess_str = wfn.multiref_wavefunction_guess(
#         high_mul, ref_zma, hs_info, multi_lvl, [cas_opt])
#     guess_lines = guess_str.splitlines()
#
#     mr_kwargs['casscf_options'] = cas_opt
#     mr_kwargs['mol_options'] = ['nosym']
#     mr_kwargs['gen_lines'] = {1: guess_lines}
#
#     ret = es_runner.read_job(
#         job='energy',
#         run_fs=run_mr_fs,
#     )
#     if ret:
#         print(" - Reading high spin multi reference energy from output...")
#         inf_obj, inp_str, out_str = ret
#         ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)
#         hs_save_fs[-1].file.energy.write(ene, multi_lvl[1:4])
#         hs_save_fs[-1].file.input.write(inp_str, multi_lvl[1:4])
#         hs_save_fs[-1].file.info.write(inf_obj, multi_lvl[1:4])
#
#     if not hs_save_fs[-1].file.energy.exists(multi_lvl[1:4]) or overwrite:
#         print(" - Running high spin multi reference energy ...")
#         es_runner.run_job(
#             job='energy',
#             script_str=mr_script_str,
#             run_fs=run_mr_fs,
#             geom=geo,
#             spc_info=hs_info,
#             thy_level=multi_lvl,
#             overwrite=overwrite,
#             **mr_kwargs,
#         )
#
#         ret = es_runner.read_job(
#             job='energy',
#             run_fs=run_mr_fs,
#         )
#
#         if ret is not None:
#             inf_obj, inp_str, out_str = ret
#
#           print(" - Reading high spin multi reference energy from output...")
#             hs_mr_ene = elstruct.reader.energy(
#                 inf_obj.prog, inf_obj.method, out_str)
#
#             print(" - Saving high spin multi reference energy...")
#             print(" - Save path: {}".format(hs_mr_save_path))
#             hs_save_fs[-1].file.energy.write(hs_mr_ene, multi_lvl[1:4])
#             hs_save_fs[-1].file.input.write(inp_str, multi_lvl[1:4])
#             hs_save_fs[-1].file.info.write(inf_obj, multi_lvl[1:4])
#         else:
#             print('ERROR: high spin multi reference energy job fails: ',
#                   'Energy is needed to evaluate infinite separation energy')
#             return inf_sep_ene
#
#     else:
#         hs_mr_ene = hs_save_fs[-1].file.energy.read(multi_lvl[1:4])
#
#     # file system for high spin single ireference calculation
#     thy_info = ['molpro2015', 'ccsd(t)-f12', 'cc-pvdz-f12', 'RR']
#     orb_restr = filesys.inf.orbital_restriction(hs_info, thy_info)
#     thy_lvl = thy_info[0:3]
#     thy_lvl.append(orb_restr)
#
#     hs_run_fs[-1].create(thy_lvl[1:4])
#     hs_save_fs[-1].create(thy_lvl[1:4])
#
#     hs_sr_run_path = hs_run_fs[-1].path(thy_lvl[1:4])
#     hs_sr_save_path = hs_save_fs[-1].path(thy_lvl[1:4])
#     run_sr_fs = autofile.fs.run(hs_sr_run_path)
#
#     sp_script_str, _, kwargs, _ = es_runner.par.run_qchem_par(*thy_lvl[0:2])
#     ret = es_runner.read_job(
#         job='energy',
#         run_fs=run_sr_fs,
#     )
#     if ret:
#         print(" - Reading high spin single reference energy from output...")
#         inf_obj, inp_str, out_str = ret
#         ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)
#         hs_save_fs[-1].file.energy.write(ene, thy_lvl[1:4])
#         hs_save_fs[-1].file.input.write(inp_str, thy_lvl[1:4])
#         hs_save_fs[-1].file.info.write(inf_obj, thy_lvl[1:4])
#
#     if not hs_save_fs[-1].file.energy.exists(thy_lvl[1:4]) or overwrite:
#         print(" - Running high spin single reference energy ...")
#
#         errors, options_mat = es_runner.par.set_molpro_options_mat(
#             hs_info, geo)
#
#         es_runner.run_job(
#             job='energy',
#             script_str=sp_script_str,
#             run_fs=run_sr_fs,
#             geom=geo,
#             spc_info=hs_info,
#             thy_level=thy_lvl,
#             errors=errors,
#             options_mat=options_mat,
#             overwrite=overwrite,
#             **kwargs,
#         )
#
#         ret = es_runner.read_job(
#             job='energy',
#             run_fs=run_sr_fs,
#         )
#
#         if ret is not None:
#             inf_obj, inp_str, out_str = ret
#
#             print(" - Reading high spin single ref energy from output...")
#             hs_sr_ene = elstruct.reader.energy(
#                 inf_obj.prog, inf_obj.method, out_str)
#
#             print(" - Saving high spin single reference energy...")
#             print(" - Save path: {}".format(hs_sr_save_path))
#             hs_save_fs[-1].file.energy.write(hs_sr_ene, thy_lvl[1:4])
#             hs_save_fs[-1].file.input.write(inp_str, thy_lvl[1:4])
#             hs_save_fs[-1].file.info.write(inf_obj, thy_lvl[1:4])
#         else:
#             print('ERROR: High spin single reference energy job fails: ',
#                   'Energy is needed to evaluate infinite separation energy')
#             return inf_sep_ene
#
#     else:
#         hs_sr_ene = hs_save_fs[-1].file.energy.read(thy_lvl[1:4])
#
#     # get the single reference energy for each of the reactant configurations
#     spc_ene = []
#     spc_infos = [spc_1_info, spc_2_info]
#     for spc_info in spc_infos:
#         # set up the file systems for the reactants one by one
#         spc_run_fs = autofile.fs.species(run_prefix)
#         spc_run_fs[-1].create(spc_info)
#         spc_run_path = spc_run_fs[-1].path(spc_info)
#         spc_save_fs = autofile.fs.species(save_prefix)
#         spc_save_fs[-1].create(spc_info)
#         spc_save_path = spc_save_fs[-1].path(spc_info)
#
#         orb_restr = filesys.inf.orbital_restriction(spc_info, ini_thy_info)
#         ini_thy_lvl = ini_thy_info[0:3]
#         ini_thy_lvl.append(orb_restr)
#
#         orb_restr = filesys.inf.orbital_restriction(spc_info, thy_info)
#         thy_lvl = thy_info[0:3]
#         thy_lvl.append(orb_restr)
#
#         ini_thy_run_fs = autofile.fs.theory(spc_run_path)
#         ini_thy_run_fs[-1].create(ini_thy_lvl[1:4])
#         ini_thy_run_path = ini_thy_run_fs[-1].path(ini_thy_lvl[1:4])
#         ini_thy_save_fs = autofile.fs.theory(spc_save_path)
#         ini_thy_save_fs[-1].create(ini_thy_lvl[1:4])
#         ini_thy_save_path = ini_thy_save_fs[-1].path(ini_thy_lvl[1:4])
#         ini_cnf_run_fs = autofile.fs.conformer(ini_thy_run_path)
#         ini_cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)
#         min_cnf_locs = filesys.minc.min_energy_conformer_locators(
#             ini_cnf_save_fs)
#         min_cnf_run_path = ini_cnf_run_fs[-1].path(min_cnf_locs)
#         min_cnf_save_path = ini_cnf_save_fs[-1].path(min_cnf_locs)
#
#         thy_run_fs = autofile.fs.theory(spc_run_path)
#         thy_run_fs[-1].create(thy_lvl[1:4])
#         thy_save_fs = autofile.fs.theory(spc_save_path)
#         thy_save_fs[-1].create(thy_lvl[1:4])
#
#         geo = ini_cnf_save_fs[-1].file.geometry.read(min_cnf_locs)
#
#         sp_run_fs = autofile.fs.single_point(min_cnf_run_path)
#         sp_run_fs[-1].create(thy_lvl[1:4])
#         sp_save_fs = autofile.fs.single_point(min_cnf_save_path)
#         sp_save_fs[-1].create(thy_lvl[1:4])
#
#         sp_sr_run_path = sp_run_fs[-1].path(thy_lvl[1:4])
#         sp_sr_save_path = sp_save_fs[-1].path(thy_lvl[1:4])
#         run_sr_fs = autofile.fs.run(sp_sr_run_path)
#
#         ret = es_runner.read_job(
#             job='energy',
#             run_fs=run_sr_fs,
#         )
#         if ret:
#             print(" - Reading single reference energy for",
#                   "{} from output...".format(spc_info[0]))
#             inf_obj, inp_str, out_str = ret
#           ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)
#             sp_save_fs[-1].file.energy.write(ene, thy_lvl[1:4])
#             sp_save_fs[-1].file.input.write(inp_str, thy_lvl[1:4])
#             sp_save_fs[-1].file.info.write(inf_obj, thy_lvl[1:4])
#
#         if not sp_save_fs[-1].file.energy.exists(thy_lvl[1:4]) or overwrite:
#             print(" - Running single reference energy for",
#                   "{} from output...".format(spc_info[0]))
#             es_runner.run_job(
#                 job='energy',
#                 script_str=sp_script_str,
#                 run_fs=run_sr_fs,
#                 geom=geo,
#                 spc_info=spc_info,
#                 thy_level=thy_lvl,
#                 overwrite=overwrite,
#                 **kwargs,
#             )
#
#             ret = es_runner.read_job(
#                 job='energy',
#                 run_fs=run_sr_fs,
#             )
#
#             if ret is not None:
#                 inf_obj, inp_str, out_str = ret
#
#                 print(" - Reading single reference energy for ",
#                       "{} from output...".format(spc_info[0]))
#                 sp_sr_ene = elstruct.reader.energy(
#                     inf_obj.prog, inf_obj.method, out_str)
#
#                 print(" - Saving single reference energy for ",
#                       "{} from output...".format(spc_info[0]))
#                 print(" - Save path: {}".format(sp_sr_save_path))
#                 sp_save_fs[-1].file.energy.write(sp_sr_ene, thy_lvl[1:4])
#                 sp_save_fs[-1].file.input.write(inp_str, thy_lvl[1:4])
#                 sp_save_fs[-1].file.info.write(inf_obj, thy_lvl[1:4])
#
#             else:
#                 print('ERROR: Single reference energy job fails',
#                       'for {}: '.format(spc_info[0]),
#                       'Energy needed to evaluate infinite separation energy')
#                 return inf_sep_ene
#
#         else:
#             sp_sr_ene = sp_save_fs[-1].file.energy.read(thy_lvl[1:4])
#
#         spc_ene.append(sp_sr_ene)
#
#     inf_sep_ene = spc_ene[0] + spc_ene[1] - hs_sr_ene + hs_mr_ene
#
#     return inf_sep_ene
