"""
 functions for reading and writing the filesystem
"""

import automol
import elstruct
import autofile
from mechroutines.es import runner as es_runner


def structure(run_fs, save_fs, locs, job, mod_thy_info,
              zma_locs=(0,), in_zma_fs=False, cart_to_zma=False):
    """ Save a geometry and associated information from some
        electronic structure routine into the filesystem.
    """

    success, ret = es_runner.read_job(job=job, run_fs=run_fs)

    if success:

        # Get the geo, zma, and ene based on job type
        ene, zma, geo = _read(run_fs, job, cart_to_zma=cart_to_zma)

        # Obtain inf obj and inp str to write in filesys
        inf_obj, inp_str, _ = ret

        # Set and print the save path information
        save_path = save_fs[-1].path(locs)
        print(" - Saving...")
        print(" - Save path: {}".format(save_path))

        # Save the geometry information
        save_fs[-1].create(locs)
        save_fs[-1].file.geometry_info.write(inf_obj, locs)
        save_fs[-1].file.geometry_input.write(inp_str, locs)
        save_fs[-1].file.geometry.write(geo, locs)
        save_fs[-1].file.energy.write(ene, locs)

        # Save zma information seperately, if required
        if not in_zma_fs:
            zma_save_fs = autofile.fs.zmatrix(save_path)
            zma_save_fs[-1].create(zma_locs)
            zma_save_fs[-1].file.geometry_info.write(inf_obj, zma_locs)
            zma_save_fs[-1].file.geometry_input.write(inp_str, zma_locs)
            zma_save_fs[-1].file.zmatrix.write(zma, zma_locs)
        elif zma:
            save_fs[-1].file.zmatrix.write(zma, locs)

        # Saving the energy to an SP filesys
        print(" - Saving energy...")
        sp_save_fs = autofile.fs.single_point(save_path)
        sp_save_fs[-1].create(mod_thy_info[1:4])
        sp_save_fs[-1].file.input.write(inp_str, mod_thy_info[1:4])
        sp_save_fs[-1].file.info.write(inf_obj, mod_thy_info[1:4])
        sp_save_fs[-1].file.energy.write(ene, mod_thy_info[1:4])

        saved = True

    else:
        saved = False

    return saved


def _read(run_fs, job, cart_to_zma=False):
    """ Read the output
    """

    assert job in (elstruct.Job.OPTIMIZATION, elstruct.Job.ENERGY), (
        'Job not opt or energy'
    )
    success, ret = es_runner.read_job(
        job=job, run_fs=run_fs)
    if ret and success:
        inf_obj, _, out_str = ret
        prog = inf_obj.prog
        method = inf_obj.method

        try:
            ene = elstruct.reader.energy(prog, method, out_str)
        except TypeError:
            ene = None
            print('No ENE')
        except AttributeError:
            ene = None
            print('No ENE')

        if job == elstruct.Job.OPTIMIZATION:
            # Read the optimized structs from the output file
            try:
                geo = elstruct.reader.opt_geometry(prog, out_str)
            except TypeError:
                geo = None
                print('No GEOM')

            # Read the optimized zma
            if cart_to_zma:
                zma = automol.geom.zmatrix(geo)
            elif 'psi' not in prog:
                try:
                    zma = elstruct.reader.opt_zmatrix(prog, out_str)
                except TypeError:
                    zma = None
                    print('No ZMA')
            else:
                zma = None
                print('No ZMA')
                

        elif job == elstruct.Job.ENERGY:
            # Read the initial structs stored in the run filesys
            zma = run_fs[-1].file.zmatrix.read([job])
            geo = run_fs[-1].file.geometry.read([job])
    else:
        ene, zma, geo = None, None, None

    return ene, zma, geo


def instability(conn_zma, disconn_zmas,
                instab_save_fs, cnf_save_fs,
                zma_locs=(0,),
                save_cnf=False):
    """ write the instability files
    """

    # Get a connected geometry
    conn_geo = automol.zmat.geometry(conn_zma)

    # Save the geometry information
    instab_save_fs[-1].create()
    instab_save_fs[-1].file.geometry.write(conn_geo)
    instab_save_path = instab_save_fs[-1].path()

    # Grab the zma and instability transformation
    zrxn, conn_zma = automol.reac.instability_transformation(
        conn_zma, disconn_zmas)

    # Save zma information seperately, if required
    zma_save_fs = autofile.fs.zmatrix(instab_save_path)
    zma_save_fs[-1].create(zma_locs)
    zma_save_fs[-1].file.zmatrix.write(conn_zma, zma_locs)

    # Write the files into the filesystem
    zma_save_fs[-1].file.reaction.write(zrxn, zma_locs)

    if save_cnf:

        # Save the geometry information
        cnf_locs = [autofile.schema.generate_new_conformer_id()]
        cnf_save_fs[-1].create(cnf_locs)
        cnf_save_fs[-1].file.geometry.write(conn_geo, cnf_locs)
        cnf_save_path = cnf_save_fs[-1].path(cnf_locs)

        # Save zma information seperately, if required
        zma_save_fs = autofile.fs.zmatrix(cnf_save_path)
        zma_save_fs[-1].create(zma_locs)
        zma_save_fs[-1].file.zmatrix.write(conn_zma, zma_locs)

    # Set and print the save path information
    print(" - Saving...")
    print(" - Save path: {}".format(instab_save_path))
    if save_cnf:
        print(" - Save path: {}".format(cnf_save_path))


def flux(vrc_ret, ts_run_fs, ts_save_fs, ts_locs=(0,), vrc_locs=(0,)):
    """ Save the VaReCoF flux and input
    """

    # Unpack the ret
    inf_obj, inp_strs, out_str = vrc_ret
    tst_str, divsur_str, molpro_str, tml_str, struct_str, pot_str = inp_strs

    # Get the flux string (somehow)
    flux_str = ''

    # Save the files
    ts_save_path = ts_save_fs[-1].path(ts_locs)

    vrc_fs = autofile.fs.vrctst(ts_save_path)
    vrc_fs[-1].create(vrc_locs)
    vrc_fs[-1].file.vrctst_tst.write(tst_str, vrc_locs)
    vrc_fs[-1].file.vrctst_divsur.write(divsur_str, vrc_locs)
    vrc_fs[-1].file.vrctst_molpro.write(molpro_str, vrc_locs)
    vrc_fs[-1].file.vrctst_tml.write(tml_str, vrc_locs)
    vrc_fs[-1].file.vrctst_struct.write(struct_str, vrc_locs)
    vrc_fs[-1].file.vrctst_pot.write(pot_str, vrc_locs)
    vrc_fs[-1].file.vrctst_flux.write(flux_str, vrc_locs)
