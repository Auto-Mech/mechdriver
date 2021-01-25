"""
Do scan stuff
"""
# def run_scan():
#     """ run and save a the scan
#     """
#     if 'elimination' in typ:
#         grid1, grid2 = grid
#         grid_dct = {dist_name: grid1, brk_name: grid2}
#     else:
#         grid_dct = {dist_name: grid}
#     print('grid_dct')
#     print(grid_dct)
#     moldr.scan.run_scan(
#         zma=ts_zma,
#         spc_info=ts_info,
#         thy_level=ref_level,
#         grid_dct=grid_dct,
#         scn_run_fs=scn_run_fs,
#         scn_save_fs=scn_save_fs,
#         script_str=opt_script_str,
#         saddle=False,
#         overwrite=overwrite,
#         update_guess=update_guess,
#         reverse_sweep=False,
#         fix_failures=False,
#         **opt_kwargs,
#         )
#     if 'elimination' in typ:
#         moldr.scan.save_scan(
#             scn_run_fs=scn_run_fs,
#             scn_save_fs=scn_save_fs,
#             coo_names=[dist_name, brk_name],
#             )
#     else:
#         moldr.scan.save_scan(
#             scn_run_fs=scn_run_fs,
#             scn_save_fs=scn_save_fs,
#             coo_names=[dist_name],
#             )
#
#     return None
