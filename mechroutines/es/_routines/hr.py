""" es_runners for coordinate scans
"""

import automol
import elstruct
from mechlib.amech_io import printer as ioprinter
from mechroutines.es.runner import scan, qchem_params


def hindered_rotor_scans(
        zma, spc_info, mod_thy_info,
        scn_run_fs, scn_save_fs,
        rotors, tors_model, method_dct,
        overwrite,
        zrxn=None,
        saddle=False,
        increment=0.5235987756,
        retryfail=True):
    """ Perform scans over each of the torsional coordinates
    """

    if tors_model != '1dhrfa':
        script_str, kwargs = qchem_params(
            method_dct, spc_info=spc_info,
            geo=automol.zmat.geometry(zma),
            job=elstruct.Job.OPTIMIZATION)
        scn_typ = 'relaxed'
        update_guess = True
        backstep = True
        reverse_sweep = True
    else:
        script_str, kwargs = qchem_params(
            method_dct, spc_info=spc_info,
            geo=automol.zmat.geometry(zma),
            job=elstruct.Job.ENERGY)
        scn_typ = 'rigid'
        update_guess = False
        backstep = False
        reverse_sweep = False

    # backstep = False
    run_tors_names = automol.data.rotor.rotors_torsion_names(rotors)
    run_tors_grids = automol.data.rotor.rotors_torsion_grids(rotors, increment=increment)

    # Set constraints
    const_names = automol.zmat.set_constraint_names(
        zma, run_tors_names, tors_model)

    ioprinter.run_rotors(run_tors_names, const_names)
    for tors_names, tors_grids in zip(run_tors_names, run_tors_grids):

        ioprinter.info_message(
            f'Running Rotor: {"-".join(tors_names)}', newline=1)

        # Setting the constraints
        constraint_dct = automol.zmat.constraint_dict(
            zma, const_names, tors_names)
        print('in hr', tors_names)
        scan.execute_scan(
            zma=zma,
            spc_info=spc_info,
            mod_thy_info=mod_thy_info,
            coord_names=tors_names,
            coord_grids=tors_grids,
            scn_run_fs=scn_run_fs,
            scn_save_fs=scn_save_fs,
            scn_typ=scn_typ,
            script_str=script_str,
            overwrite=overwrite,
            zrxn=zrxn,
            update_guess=update_guess,
            reverse_sweep=reverse_sweep,
            saddle=saddle,
            constraint_dct=constraint_dct,
            retryfail=retryfail,
            **kwargs,
        )
        if backstep:
            ioprinter.info_message('Attempting backstep routine')
            scan.run_backsteps(
                zma=zma,
                spc_info=spc_info,
                mod_thy_info=mod_thy_info,
                coord_names=tors_names,
                coord_grids=tors_grids,
                scn_run_fs=scn_run_fs,
                scn_save_fs=scn_save_fs,
                scn_typ=scn_typ,
                script_str=script_str,
                overwrite=overwrite,
                saddle=saddle,
                constraint_dct=constraint_dct,
                retryfail=retryfail,
                **kwargs,
            )


def check_hr_pot(tors_pots, tors_zmas, tors_paths, emax=-0.5, emin=-10.0):
    """ Check hr pot to see if a new mimnimum is needed
    """

    new_min_zma = None
    new_emin = emax
    new_path = None

    print('\nAssessing the HR potential...')
    for name in tors_pots:

        print(f'- Rotor {name}')
        pots = tors_pots[name].values()
        zmas = tors_zmas[name].values()
        paths = tors_paths[name].values()
        for pot, zma, path in zip(pots, zmas, paths):
            if emin < pot < emax:
                if pot < new_emin:
                    new_min_zma = zma
                    new_path = path
                    new_emin = pot
    if new_min_zma is not None:
        print(' - New minimmum energy ZMA found for torsion')
        print(f' - Ene = {new_emin}')
        print(f' - Found at path: {new_path}')
        print(automol.zmat.string(new_min_zma))

    return new_min_zma
