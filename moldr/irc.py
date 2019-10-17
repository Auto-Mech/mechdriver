"""
toy irc code
"""

def _run_irc(
        spc_info, thy_level, cnf_run_fs, cnf_save_fs,
        irc_script_str, overwrite, new_grid=False, **opt_kwargs):
    """ Run the IRC 
    """
    
    # Obtain saddle-point minimmum-energy conformer from filesystem or argument?

    # Set up the IRC run filesystem
    run_fs = autofile.fs.run()

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

        opt_ret = moldr.driver.read_job(
            job='irc',
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

            dist_name = '??????'
            for idx, coord in enumerate(coords):
                locs = [[dist_name], [coord]]
                # save_fs.leaf.file.energy.write(enes[idx], locs)
                # save_fs.leaf.file.geometry.write(geos[idx], locs)
                # save_fs.leaf.file.gradient.write(gras[idx], locs)
                # save_fs.leaf.file.hessian.write(hessians[idx], locs)
