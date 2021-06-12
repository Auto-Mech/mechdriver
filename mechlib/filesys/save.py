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
        cnf_locs = (
            autofile.schema.generate_new_ring_id(),
            autofile.schema.generate_new_conformer_id()
        )
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


def energy_transfer(etrans_save_fs, etrans_locs,
                    run_epsilons, run_sigmas, run_geos, run_ranseeds,
                    version, input_str, els_str):
    """ Save lennard-jones energy transfer params
    """

    # Read any epsilons and sigma currently in the filesystem
    ioprinter.info_message(
        'Reading Lennard-Jones parameters and Geoms from filesystem...',
        newline=1)
    if etrans_save_fs[-1].file.epsilon.exists(etrans_locs):
        fs_eps = etrans_save_fs[-1].file.epsilon.read(etrans_locs)
    else:
        fs_eps = ()
    if etrans_save_fs[-1].file.sigma.exists(etrans_locs):
        fs_sig = etrans_save_fs[-1].file.sigma.read(etrans_locs)
    else:
        fs_sig = ()
    if etrans_save_fs[-1].file.trajectory.exists(etrans_locs):
        traj = etrans_save_fs[-1].file.trajectory.read(etrans_locs)
        fs_geos = tuple(x[0] for x in traj)
        fs_ranseeds = tuple(x[1] for x in traj)
    else:
        fs_geos, fs_ranseeds = (), ()

    # Set nsamp numbers
    ini_nsampd = len(run_geos)
    cur_nsampd = len(run_sigmas)
    full_nsampd = ini_nsampd + cur_nsampd

    # Add the lists from the two together
    geoms = fs_geos + run_geos
    sigmas = (fs_sig,) + run_sigmas
    epsilons = (fs_eps,) + run_epsilons
    ranseeds = fs_ranseeds + run_ranseeds

    assert len(geoms) == len(sigmas) == len(epsilons), (
        'Number of geoms, sigmas, and epsilons not the same'
    )

    avg_sigma = sigmas / full_nsampd
    avg_epsilon = epsilons / full_nsampd
    ioprinter.info_message(
        'Averaging new vals with one in filesys')
    ioprinter.info_message(
        'Average Sigma to save [unit]:', avg_sigma, newline=1)
    ioprinter.info_message('Number of values = ', full_nsampd)

    # Update the trajectory file
    traj = []
    for geo, eps, sig, ranseed in zip(geoms, epsilons, sigmas, ranseeds):
        comment = 'Epsilon: {} cm-1   Sigma: {} Ang   RandSeed: {}'.format(
            eps, sig, ranseed)
        traj.append((comment, geo))

    # Write the info obj add ranseeds to this I think
    inf_obj = autofile.schema.info_objects.lennard_jones(
        full_nsampd, program='OneDMin', version=version)

    # Write the params to the save file system
    etrans_save_fs[-1].file.lj_input.write(input_str, etrans_locs)
    etrans_save_fs[-1].file.info.write(inf_obj, etrans_locs)
    etrans_save_fs[-1].file.molpro_inp_file.write(els_str, etrans_locs)
    etrans_save_fs[-1].file.epsilon.write(avg_epsilon, etrans_locs)
    etrans_save_fs[-1].file.sigma.write(avg_sigma, etrans_locs)


# Job readers
# Constituent functions for saving various bits of information
def read_job_zma(ret, init_zma=None, rebuild=False):
    """ trys to read a zma from a job and does processing on it as needed

        (1) try to read opt zma from output
        (2) update supplied init zma using opt geo (in outpt)
        (3) update init zma (in output) using opt geo (in outpt)
        (4) convert opt geo into a zma [dies for disconn zmas]
    """

    # if zma is None:
    #     if zrxn is None:
    #         zma = automol.geom.zmatrix(geo)
    #     else:
    #         zma = automol.reac.ts_zmatrix(zrxn, geo)

    _, _, out_str, prog, _ = _unpack_ret(ret)
    zma = elstruct.reader.opt_zmatrix(prog, out_str)
    if zma is None or rebuild:
        print('Getting ZMA from a geometry...')
        geo = elstruct.reader.opt_geometry(prog, out_str)
        if init_zma is not None:
            print('Resetting ZMA coords using opt geoms...')
            zma = rebuild_zma_from_opt_geo(init_zma, geo)
        else:
            init_zma = elstruct.reader.inp_zmatrix(prog, out_str)
            if init_zma is not None:
                print('Resetting ZMA coords using opt geoms...')
                zma = rebuild_zma_from_opt_geo(init_zma, geo)
            else:
                zma = automol.geom.zmatrix(geo)

    return zma


def rebuild_zma_from_opt_geo(init_zma, opt_geo):
    """ build zma from opt
    """

    dummy_key_dct = automol.zmat.dummy_key_dictionary(init_zma)
    geo_wdummy = automol.geom.insert_dummies(opt_geo, dummy_key_dct)
    zma = automol.zmat.from_geometry(init_zma, geo_wdummy)

    return zma


def _save_geom(ret, cnf_fs, cnf_locs):
    """ Saving a Cartesian geometry. could be for cnf fs or scn fs
    """

    print(" - Reading geometry from output...")
    inf_obj, inp_str, out_str, prog, _ = _unpack_ret(ret)

    geo = elstruct.reader.opt_geometry(prog, out_str)

    cnf_fs[-1].create(cnf_locs)
    cnf_path = cnf_fs[-1].path(cnf_locs)
    print(" - Saving at {}".format(cnf_path))

    cnf_fs[-1].file.geometry_info.write(inf_obj, cnf_locs)
    cnf_fs[-1].file.geometry_input.write(inp_str, cnf_locs)
    cnf_fs[-1].file.geometry.write(geo, cnf_locs)


def _save_grad(ret, cnf_fs, cnf_locs):
    """ Saving a geometry
    """

    print(" - Reading gradient from output...")
    inf_obj, inp_str, out_str, prog, _ = _unpack_ret(ret)

    grad = elstruct.reader.gradient(prog, out_str)

    cnf_fs[-1].create(cnf_locs)
    cnf_path = cnf_fs[-1].path(cnf_locs)
    print(" - Saving at {}".format(cnf_path))

    cnf_fs[-1].file.gradient_info.write(inf_obj, cnf_locs)
    cnf_fs[-1].file.gradient_input.write(inp_str, cnf_locs)
    cnf_fs[-1].file.gradient.write(grad, cnf_locs)


def _save_zmatrix(ret, zma_fs, zma_locs, init_zma=None):
    """ Read a zma from output. If one is not found,
        then attempt to convert it from geometry.

        might not work if redundants print like zma
    """

    print(" - Reading Z-Matrix from output...")
    inf_obj, inp_str, _, _, _ = _unpack_ret(ret)

    zma = read_job_zma(ret, init_zma=init_zma, rebuild=False)

    zma_fs[-1].create(zma_locs)
    zma_path = zma_fs[-1].path(zma_locs)
    print(" - Saving at {}".format(zma_path))

    zma_fs[-1].file.geometry_info.write(inf_obj, zma_locs)
    zma_fs[-1].file.geometry_input.write(inp_str, zma_locs)
    zma_fs[-1].file.zmatrix.write(zma, zma_locs)


def _save_energy(ret, sp_fs, sp_locs):
    """ Sving a energy
    """

    print(" - Reading energy from output...")
    inf_obj, inp_str, out_str, prog, method = _unpack_ret(ret)

    ene = elstruct.reader.energy(prog, method, out_str)

    sp_fs[-1].create(sp_locs)
    sp_path = sp_fs[-1].path(sp_locs)
    print(" - Saving at {}".format(sp_path))

    sp_fs[-1].file.input.write(inp_str, sp_locs)
    sp_fs[-1].file.info.write(inf_obj, sp_locs)
    sp_fs[-1].file.energy.write(ene, sp_locs)


def _save_hessian(ret, cnf_fs, cnf_locs):
    """ saved the hessian
    """

    print(" - Reading hessian and harmonic frequencies from output...")
    inf_obj, inp_str, out_str, prog, _ = _unpack_ret(ret)

    hess = elstruct.reader.hessian(prog, out_str)
    freqs = elstruct.reader.harmonic_frequencies(prog, out_str)

    cnf_fs[-1].create(cnf_locs)
    cnf_path = cnf_fs[-1].path(cnf_locs)
    print(" - Saving at {}".format(cnf_path))

    cnf_fs[-1].file.hessian_info.write(inf_obj, cnf_locs)
    cnf_fs[-1].file.hessian_input.write(inp_str, cnf_locs)
    cnf_fs[-1].file.hessian.write(hess, cnf_locs)
    cnf_fs[-1].file.harmonic_frequencies.write(freqs, cnf_locs)


def _save_rotors(zma_fs, zma_locs, zrxn=None):
    """ Save the rotors

        Reading the ZMA from the filesystem to build rotors
    """

    zma = zma_fs[-1].file.zmatrix.read(zma_locs)
    rotors = automol.rotor.from_zmatrix(zma, zrxn=zrxn)
    if any(rotors):
        zma_path = zma_fs[-1].path(zma_locs)
        print(" - Rotors identified from Z-Matrix at {}".format(zma_path))
        print(" - Saving rotor information at same location")
        zma_fs[-1].file.torsions.write(rotors, zma_locs)


def _save_rings(zma_fs, zma_locs, zrxn=None):
    """ Save rings
    """

    zma = zma_fs[-1].file.zmatrix.read(zma_locs)
    rings_atoms = automol.zmat.all_rings_atoms(zma, zrxn=zrxn)
    ring_dct = automol.zmat.all_rings_dct(zma, rings_atoms)
    if ring_dct:
        zma_fs[-1].file.ring_torsions.write(ring_dct, zma_locs)


def _save_reaction(zma_fs, zma_locs, zrxn=None):
    """ Save the reaction object
    """

    if zrxn is not None:
        print(" - Saving reaction")
        zma_fs[-1].file.reaction.write(zrxn, zma_locs)


def _save_instab(ret, instab_fs):
    """ save the instability filesys
    """

    print(" - Reading hessian from output...")
    inf_obj, inp_str, _, _, _ = _unpack_ret(ret)

    geo = ()

    # Save the geometry information
    instab_fs[-1].create()
    instab_path = instab_fs[-1].path()
    print(" - Saving at {}".format(instab_path))
    instab_fs[-1].file.geometry_info.write(inf_obj)
    instab_fs[-1].file.geometry_input.write(inp_str)
    instab_fs[-1].file.geometry.write(geo)


# Helper functions for reading jobs
def _generate_locs(rng_locs=None, tors_locs=None, zma_locs=None):
    """ Determine locators
    """

    if rng_locs is None:
        rng_locs = [autofile.schema.generate_new_ring_id()]
    if tors_locs is None:
        tors_locs = [autofile.schema.generate_new_conformer_id()]
    if zma_locs is None:
        zma_locs = (0,)

    cnf_locs = rng_locs + tors_locs

    return cnf_locs, zma_locs


def _unpack_ret(ret):
    """ Unpack ret object and get other commonly useful data
    """

    inf_obj, inp_str, out_str = ret
    prog = inf_obj.prog
    method = inf_obj.method

    return inf_obj, inp_str, out_str, prog, method
