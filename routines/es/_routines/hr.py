""" es_runners for coordinate scans
"""

from routines.es._routines import _scan as scan


def hindered_rotor_scans(
        zma, spc_info, thy_info, scn_run_fs, scn_save_fs,
        run_tors_names, run_tors_grids,
        script_str, overwrite,
        saddle=False, constraint_dct=None,
        retryfail=True, **opt_kwargs):
    """ Perform scans over each of the torsional coordinates
    """

    print('\nRunning hindered rotor scans for the following rotors...')
    for names in run_tors_names:
        print(names)
    if constraint_dct is not None:
        print('\nUser requested that all torsions of system will be fixed.')

    # for tors_name, tors_grid in zip(tors_names, tors_grids):
    for tors_names, tors_grids in zip(run_tors_names, run_tors_grids):

        print('\nRunning Rotor: {}...'.format(tors_names))

        # Get the dictionary for the torsional modes
        if not tors_names:
            continue
        grid_dct = dict(zip(tors_names, tors_grids))

        print('\nSaving any HR in run filesys...')
        if constraint_dct is None:
            scan.save_scan(
                scn_run_fs=scn_run_fs,
                scn_save_fs=scn_save_fs,
                coo_names=tors_names,
                thy_info=thy_info)
        else:
            scan.save_cscan(
                cscn_run_fs=scn_run_fs,
                cscn_save_fs=scn_save_fs,
                coo_names=tors_names,
                thy_info=thy_info)

        print('\nRunning any HR Scans if needed...')
        scan.run_scan(
            zma=zma,
            spc_info=spc_info,
            thy_info=thy_info,
            grid_dct=grid_dct,
            scn_run_fs=scn_run_fs,
            scn_save_fs=scn_save_fs,
            script_str=script_str,
            overwrite=overwrite,
            saddle=saddle,
            retryfail=retryfail,
            constraint_dct=constraint_dct,
            **opt_kwargs,
        )

        print('\nSaving any newly run HR scans in run filesys...')
        if constraint_dct is None:
            scan.save_scan(
                scn_run_fs=scn_run_fs,
                scn_save_fs=scn_save_fs,
                coo_names=tors_names,
                thy_info=thy_info
            )
        else:
            scan.save_cscan(
                cscn_run_fs=scn_run_fs,
                cscn_save_fs=scn_save_fs,
                coo_names=tors_names,
                thy_info=thy_info
            )
