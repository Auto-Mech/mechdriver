""" es_runners for coordinate scans
"""

import automol
import elstruct
from mechroutines.es.runner import scan
from mechroutines.es.runner import qchem_params
from mechlib.amech_io import printer as ioprinter
from mechlib import structure
from phydat import phycon


def hindered_rotor_scans(
        zma, spc_info, mod_thy_info, thy_save_fs,
        scn_run_fs, scn_save_fs,
        rotors, tors_model, method_dct,
        overwrite,
        saddle=False,
        increment=30.0*phycon.DEG2RAD,
        retryfail=True, chkstab=None):
    """ Perform scans over each of the torsional coordinates
    """

    if tors_model != '1dhrfa':
        script_str, kwargs = qchem_params(
            method_dct, job=elstruct.Job.OPTIMIZATION)
        scn_typ = 'relaxed'
    else:
        script_str, kwargs = qchem_params(
            method_dct, job=elstruct.Job.ENERGY)
        scn_typ = 'rigid'

    run_tors_names = automol.rotor.names(rotors)
    run_tors_grids = automol.rotor.grids(rotors, increment=increment)

    # Set constraints
    const_names = structure.tors.set_constraint_names(
        zma, run_tors_names, tors_model)

    # Set appropriate value for check stability
    # If not set, don't check if saddle=True
    if chkstab is None:
        chkstab = bool(not saddle)

    ioprinter.run_rotors(run_tors_names, const_names)

    # for tors_name, tors_grid in zip(tors_names, tors_grids):
    for tors_names, tors_grids in zip(run_tors_names, run_tors_grids):

        ioprinter.info_message(
            'Running Rotor: {}...'.format(tors_names),
            newline=1)

        # Setting the constraints
        constraint_dct = automol.zmat.constraint_dct(
            zma, const_names, tors_names)

        scan.execute_scan(
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
            chkstab=False,
            **kwargs,
        )


def check_hr_pot(tors_pots, tors_zmas, tors_paths, emax=-0.5, emin=-10.0):
    """ Check hr pot to see if a new mimnimum is needed
    """

    new_min_zma = None

    print('\nAssessing the HR potential...')
    for name in tors_pots:

        print('- Rotor {}'.format(name))
        pots = tors_pots[name].values()
        zmas = tors_zmas[name].values()
        paths = tors_paths[name].values()
        for pot, zma, path in zip(pots, zmas, paths):
            if emin < pot < emax:
                new_min_zma = zma
                emin = pot
                print(' - New minimmum energy ZMA found for torsion')
                print(' - Ene = {}'.format(pot))
                print(' - Found at path: {}'.format(path))
                print(automol.zmat.string(zma))

    return new_min_zma


# Read and print the potential
# sp_fs = autofile.fs.single_point(ini_cnf_save_path)
# ref_ene = sp_fs[-1].file.energy.read(mod_ini_thy_info[1:4])
# ref_ene = ini_cnf_save_fs[-1].file.energy.read(ini_min_cnf_locs)
# tors_pots, tors_zmas = {}, {}
# for tors_names, tors_grids in zip(run_tors_names, run_tors_grids):
#     constraint_dct = structure.tors.build_constraint_dct(
#         zma, const_names, tors_names)
#     pot, _, _, _, zmas, _ = structure.tors.read_hr_pot(
#         tors_names, tors_grids,
#         ini_cnf_save_path,
#         mod_ini_thy_info, ref_ene,
#         constraint_dct,
#         read_zma=True)
#     tors_pots[tors_names] = pot
#     tors_zmas[tors_names] = zmas

# # Print potential
# structure.tors.print_hr_pot(tors_pots)
