""" es_runners for coordinate scans
"""

import itertools
from routines.es._routines import _scan as scan
from lib import filesys
from lib.structure import tors as torsprep


def hindered_rotor_scans(
        zma, spc_info, mod_thy_info, thy_save_fs,
        zma_run_path, zma_save_path,
        run_tors_names, run_tors_grids,
        script_str, overwrite,
        scn_typ='relaxed',
        saddle=False, const_names=None,
        retryfail=True, chkstab=None, **opt_kwargs):
    """ Perform scans over each of the torsional coordinates
    """

    # Set appropriate value for check stability
    # If not set, don't check if saddle=True
    if chkstab is None:
        chkstab = bool(not saddle)

    print('\nRunning hindered rotor scans for the following rotors...')
    for names in run_tors_names:
        print(names)
    # print(*run_tors_names)
    # print(list(itertools.chain(*run_tors_names)))
    # print(const_names)
    if const_names is not None:
        if set(list(itertools.chain(*run_tors_names))) == set(const_names):
            print('\nUser requested all torsions of system will be fixed.')

    # for tors_name, tors_grid in zip(tors_names, tors_grids):
    for tors_names, tors_grids in zip(run_tors_names, run_tors_grids):

        print('\nRunning Rotor: {}...'.format(tors_names))

        # Setting the constraints
        constraint_dct = torsprep.build_constraint_dct(
            zma, const_names, tors_names)

        # Setting the filesystem
        print('hr constraint dct', constraint_dct)
        scn_run_fs = filesys.build.scn_fs_from_cnf(
            zma_run_path, constraint_dct=constraint_dct)
        scn_save_fs = filesys.build.scn_fs_from_cnf(
            zma_save_path, constraint_dct=constraint_dct)

        print('\nSaving any HR in run filesys...')
        if constraint_dct is None:
            scan.save_scan(
                scn_run_fs=scn_run_fs,
                scn_save_fs=scn_save_fs,
                scn_typ=scn_typ,
                coo_names=tors_names,
                mod_thy_info=mod_thy_info,
                in_zma_fs=True)
        else:
            scan.save_cscan(
                cscn_run_fs=scn_run_fs,
                cscn_save_fs=scn_save_fs,
                scn_typ=scn_typ,
                coo_names=tors_names,
                constraint_dct=constraint_dct,
                mod_thy_info=mod_thy_info,
                in_zma_fs=True)

        print('\nRunning any HR Scans if needed...')
        scan.run_scan(
            zma=zma,
            spc_info=spc_info,
            mod_thy_info=mod_thy_info,
            thy_save_fs=thy_save_fs,
            coord_names=tors_names,
            coord_grids=tors_grids,
            scn_run_fs=scn_run_fs,
            scn_save_fs=scn_save_fs,
            scn_typ=scn_typ,
            script_str=script_str,
            overwrite=overwrite,
            update_guess=True,
            reverse_sweep=True,
            saddle=saddle,
            constraint_dct=constraint_dct,
            retryfail=retryfail,
            chkstab=chkstab,
            **opt_kwargs
        )

        print('\nSaving any newly run HR scans in run filesys...')
        if constraint_dct is None:
            scan.save_scan(
                scn_run_fs=scn_run_fs,
                scn_save_fs=scn_save_fs,
                scn_typ=scn_typ,
                coo_names=tors_names,
                mod_thy_info=mod_thy_info,
                in_zma_fs=True)
        else:
            scan.save_cscan(
                cscn_run_fs=scn_run_fs,
                cscn_save_fs=scn_save_fs,
                scn_typ=scn_typ,
                coo_names=tors_names,
                constraint_dct=constraint_dct,
                mod_thy_info=mod_thy_info,
                in_zma_fs=True)
