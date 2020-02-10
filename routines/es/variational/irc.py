"""
IRC calcs
"""

def _run_irc(
        ts_info, thy_level, ts_run_fs, ts_save_fs, locs,
        overwrite, new_grid=False, **opt_kwargs):
    """ Run the IRC
    """

    # Obtain saddle-point minimmum-energy conformer from filesystem
    ts_run_path = ts_run_fs.leaf.path(locs)
    # ts_save_path = ts_save_fs.leaf.path(locs)
    geo = ts_save_fs.leaf.file.geometry.read(locs)
    
    # Check if IRC run to desired specs
    # If Not run the IRC calculation
    for grid_idx, grid_val, run_prefix in zip(grid_idxs, grid_vals, run_prefixes):
        if not scn_save_fs.leaf.file.geometry.exists([['RX'], [grid_val]]) or overwrite:
            run_irc = True

    if run_irc:
        # Set up the IRC run filesystem
        run_fs = autofile.fs.run(ts_run_path)

        # Run the IRC in the forward and reverse direction
        for irc_direction in ('forward', 'reverse'):
            moldr.driver.run_job(
                job='irc',
                script_str=irc_script_str,
                run_fs=run_fs,
                geom=zma,
                spc_info=ts_info,
                thy_level=ref_level,
                overwrite=overwrite,
                irc_direction=irc_direction,
                **opt_kwargs,
                )

    # Read the IRC from the file system
    for irc_direction in ('forward', 'reverse'):
        opt_ret = moldr.driver.read_job(
            job=elstruct.Job.IRC,
            run_fs=run_fs,
        )

        if opt_ret is not None:
            inf_obj, _, out_str = opt_ret
            prog = inf_obj.prog
            geos, gras, hessians = elstruct.reader.irc_points(prog, out_str)
            enes = elstruct.reader.irc_energies(prog, out_str)
            coords = elstruct.reader.irc_coordinates(prog, out_str)

            print(" - Saving...")
            print(" - Save path: {}".format())

            dist_name = 'RX'
            for idx, coord in enumerate(coords):
                locs = [[dist_name], [coord]]
                # save_fs.leaf.file.energy.write(enes[idx], locs)
                # save_fs.leaf.file.geometry.write(geos[idx], locs)
                # save_fs.leaf.file.gradient.write(gras[idx], locs)
                # save_fs.leaf.file.hessian.write(hessians[idx], locs)
