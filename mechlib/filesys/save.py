"""
 functions for reading and writing the filesystem
"""

import automol
import elstruct
import autofile
from mechlib.amech_io import printer as ioprinter


def atom(sp_ret, cnf_fs, thy_locs, zma,
         rng_locs=None, tors_locs=None, zma_locs=(0,)):
    """ Save an all of the information for an atom into the
        CONF layer of the save filesystem: SPC/THY/CONFS/Z
        that can be parsed from a single-point energy job.

        :param sp_ret:
        :type sp_ret:
        :param cnf_fs:
    """

    # Build filesystem locs and objects
    cnf_locs, zma_locs = _generate_locs(
        rng_locs=rng_locs, tors_locs=tors_locs, zma_locs=zma_locs)

    # Save an arbitrary geom and zma
    cnf_fs[-1].create(cnf_locs)
    zma_fs = autofile.fs.zmatrix(cnf_fs[-1].path(cnf_locs))
    zma_fs[-1].create(zma_locs)

    geo = automol.zmat.geometry(zma)
    cnf_fs[-1].file.geometry.write(geo, cnf_locs)
    zma_fs[-1].file.zmatrix.write(zma, zma_locs)

    # Save data from energy job
    sp_fs = autofile.fs.single_point(cnf_fs[-1].path(cnf_locs))
    _save_energy(sp_ret, sp_fs, thy_locs)


def _conformer_setup(cnf_fs, rng_locs, tors_locs, zma_locs):
    """ sets up the save locs and fs for conformer
    """
    # Build filesystem locs and objects
    cnf_locs, zma_locs = _generate_locs(
        rng_locs=rng_locs, tors_locs=tors_locs, zma_locs=zma_locs)

    zma_fs = autofile.fs.zmatrix(cnf_fs[-1].path(cnf_locs))
    sp_fs = autofile.fs.single_point(cnf_fs[-1].path(cnf_locs))
    return cnf_locs, zma_locs, zma_fs, sp_fs


def _conformer_aux_info(zma_fs, zma_locs, zrxn=None):
    """ Save rings, rotors, and reaction objects
        into the ZMA fs
    """
    _save_rotors(zma_fs, zma_locs, zrxn=zrxn)
    _save_rings(zma_fs, zma_locs, zrxn=zrxn)
    _save_reaction(zma_fs, zma_locs, zrxn=zrxn)


def conformer(opt_ret, hess_ret, cnf_fs, thy_locs,
              zrxn=None, init_zma=None,
              rng_locs=None, tors_locs=None, zma_locs=None):
    """ Save a conformer and relevant information.

        thy_fs = (_fs, _locs)
    """
    cnf_locs, zma_locs, zma_fs, sp_fs = _conformer_setup(
        cnf_fs, rng_locs, tors_locs, zma_locs)

    # Save data from optimization and hessian jobs
    _save_geom(opt_ret, cnf_fs, cnf_locs)
    _save_zmatrix(opt_ret, zma_fs, zma_locs, init_zma=init_zma)
    _save_energy(opt_ret, sp_fs, thy_locs)

    if hess_ret is not None:
        _save_hessian(hess_ret, cnf_fs, cnf_locs)

    # Save ring and cnf samp files, if needed
    init_cnf_samp(cnf_fs, cnf_locs)

    # Save auxiliary information for the structure, if needed
    _conformer_aux_info(zma_fs, zma_locs, zrxn=zrxn)


def parsed_conformer(
        save_info, cnf_fs, thy_locs, rng_locs,
        zma_locs=(0,), tors_locs=None, zrxn=None, hess_ret=None):
    """ Save a conformer and relavent information without an output ret
    """
    cnf_locs, zma_locs, zma_fs, sp_fs = _conformer_setup(
        cnf_fs, rng_locs, tors_locs, zma_locs)

    # Save data from optimization and hessian jobs
    geo, zma, ene, inf_obj, inp_str = save_info
    _save_geom_parsed(geo, inf_obj, inp_str, cnf_fs, cnf_locs)
    _save_zmatrix_parsed(zma, inf_obj, inp_str, zma_fs, zma_locs)
    _save_energy_parsed(ene, inf_obj, inp_str, sp_fs, thy_locs)
    if hess_ret is not None:
        _save_hessian(hess_ret, cnf_fs, cnf_locs)

    # Save ring and cnf samp files, if needed
    init_cnf_samp(cnf_fs, cnf_locs)

    # Save auxiliary information for the structure, if needed
    _conformer_aux_info(zma_fs, zma_locs, zrxn=zrxn)


def sym_indistinct_conformer(geo, cnf_fs, cnf_tosave_locs, cnf_saved_locs):
    """ Save conformer that is symmetryically similar to another conformer
        that is in the filesystem.
    """

    # Set the path to the conf filesys with sym similar
    cnf_save_path = cnf_fs[-1].path(cnf_saved_locs)

    # Build the sym file sys
    sym_save_fs = autofile.fs.symmetry(cnf_save_path)
    sym_save_path = cnf_fs[-1].path(cnf_saved_locs)
    ioprinter.save_symmetry(sym_save_path)
    sym_save_fs[-1].create([cnf_tosave_locs[-1]])
    sym_save_fs[-1].file.geometry.write(geo, [cnf_tosave_locs[-1]])


def scan_point_structure(opt_ret, scn_fs, scn_locs, thy_locs, job,
                         init_zma=None, init_geo=None):
    """ save info for the hindered rotor
    """

    if job == elstruct.Job.OPTIMIZATION:
        _save_geom(opt_ret, scn_fs, scn_locs)
        _save_zmatrix(opt_ret, scn_fs, scn_locs, init_zma=init_zma)
    elif job == elstruct.Job.ENERGY:
        scn_fs[-1].create(scn_locs)
        scn_fs[-1].file.zmatrix.write(init_zma, scn_locs)
        if init_geo is None:
            init_geo = automol.zmat.geometry(init_zma)
        scn_fs[-1].file.geometry.write(init_geo, scn_locs)

    sp_fs = autofile.fs.single_point(scn_fs[-1].path(scn_locs))
    _save_energy(opt_ret, sp_fs, thy_locs)


def init_cnf_samp(cnf_fs, cnf_locs):
    """ init cnf samp
    """

    inf_obj = autofile.schema.info_objects.conformer_trunk(0)
    inf_obj.nsamp = 1

    if not cnf_fs[0].file.info.exists():
        rinf_obj = autofile.schema.info_objects.conformer_trunk(0)
        rinf_obj.nsamp = 1
        cnf_fs[0].file.info.write(rinf_obj)
    if not cnf_fs[1].file.info.exists([cnf_locs[0]]):
        cinf_obj = autofile.schema.info_objects.conformer_branch(0)
        cinf_obj.nsamp = 1
        cnf_fs[1].file.info.write(cinf_obj, [cnf_locs[0]])


# Saving other
def instability(conn_zma, disconn_zmas, cnf_save_fs,
                rng_locs=None, tors_locs=None, zma_locs=(0,)):
    """ write the instability files
    """

    cnf_locs, zma_locs = _generate_locs(
        rng_locs=rng_locs, tors_locs=tors_locs, zma_locs=zma_locs)

    # Get structural and instability info for storage
    conn_geo = automol.zmat.geometry(conn_zma)

    zrxn, conn_zma = automol.reac.instability_transformation(
        conn_zma, disconn_zmas)

    # Save structure and stability info in CONF and Z filesys
    cnf_save_fs[-1].create(cnf_locs)
    cnf_save_fs[-1].file.geometry.write(conn_geo, cnf_locs)
    cnf_save_path = cnf_save_fs[-1].path(cnf_locs)

    # Save zma information seperately, if required
    zma_save_fs = autofile.fs.zmatrix(cnf_save_path)
    zma_save_fs[-1].create(zma_locs)
    zma_save_fs[-1].file.zmatrix.write(conn_zma, zma_locs)
    zma_save_fs[-1].file.instability.write(zrxn, zma_locs)

    # Set and print the save path information
    print(" - Saving...")
    print(f" - Save path: {cnf_save_path}")


def flux(flux_str, corr_pot_str, inp_strs, vrc_fs, vrc_locs=(0,)):
    """ Save the VaReCoF flux and input
    """

    # Unpack the ret
    divsur_str = inp_strs['divsur.inp']
    tst_str = inp_strs['tst.inp']
    molpro_str = inp_strs['molpro.inp']
    struct_str = inp_strs['structure.inp']
    tml_str = inp_strs['run.tml']

    # Save the files
    vrc_fs[-1].create(vrc_locs)
    vrc_fs[-1].file.vrctst_tst.write(tst_str, vrc_locs)
    vrc_fs[-1].file.vrctst_divsur.write(divsur_str, vrc_locs)
    vrc_fs[-1].file.vrctst_molpro.write(molpro_str, vrc_locs)
    vrc_fs[-1].file.vrctst_tml.write(tml_str, vrc_locs)
    vrc_fs[-1].file.vrctst_struct.write(struct_str, vrc_locs)
    vrc_fs[-1].file.vrctst_pot.write(corr_pot_str, vrc_locs)
    vrc_fs[-1].file.vrctst_flux.write(flux_str, vrc_locs)

    # Set and print the save path information
    vrc_path = vrc_fs[-1].path(vrc_locs)
    print(" - Saving VRCTST input files and final flux file...")
    print(f" - Save path: {vrc_path}")


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
        comment = f'Epsilon: {eps} cm-1  Sigma: {sig} Ang  RandSeed: {ranseed}'
        traj.append((comment, geo))

    # Write the info obj add ranseeds to this I think
    inf_obj = autofile.schema.info_objects.lennard_jones(
        full_nsampd, program='OneDMin', version=version)

    # Write the params to the save file system
    etrans_save_fs[-1].file.lj_input.write(input_str, etrans_locs)
    etrans_save_fs[-1].file.lennard_jones_info.write(inf_obj, etrans_locs)
    etrans_save_fs[-1].file.molpro_inp_file.write(els_str, etrans_locs)
    etrans_save_fs[-1].file.epsilon.write(avg_epsilon, etrans_locs)
    etrans_save_fs[-1].file.sigma.write(avg_sigma, etrans_locs)


# Job readers
# Constituent functions for saving various bits of information
def read_zma_from_geo(init_zma, geo):
    """ trys to read a zma from a job and does processing on it as needed
        (1) update init zma (in output) using opt geo (in outpt)
    """
    zma = None
    if init_zma is not None:
        print('Resetting ZMA coords using opt geoms...')
        zma = rebuild_zma_from_opt_geo(init_zma, geo)
    return zma


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

    dtt = automol.zmat.conversion_info(init_zma)
    geo_wdummy = automol.geom.apply_zmatrix_conversion(opt_geo, dtt)
    zma = automol.zmat.from_geometry(init_zma, geo_wdummy)

    return zma


def _save_geom_parsed(geo, inf_obj, inp_str, cnf_fs, cnf_locs):
    """ Saving a Cartesian geometry. could be for cnf fs or scn fs
    """

    cnf_fs[-1].create(cnf_locs)
    cnf_path = cnf_fs[-1].path(cnf_locs)
    print(f" - Saving at {cnf_path}")

    cnf_fs[-1].file.geometry_info.write(inf_obj, cnf_locs)
    cnf_fs[-1].file.geometry_input.write(inp_str, cnf_locs)
    cnf_fs[-1].file.geometry.write(geo, cnf_locs)


def _save_geom(ret, cnf_fs, cnf_locs):
    """ Saving a Cartesian geometry. could be for cnf fs or scn fs
    """

    print(" - Reading geometry from output...")
    inf_obj, inp_str, out_str, prog, _ = _unpack_ret(ret)
    geo = elstruct.reader.opt_geometry(prog, out_str)
    _save_geom_parsed(geo, inf_obj, inp_str, cnf_fs, cnf_locs)


def _save_grad(ret, cnf_fs, cnf_locs):
    """ Saving a geometry
    """

    print(" - Reading gradient from output...")
    inf_obj, inp_str, out_str, prog, _ = _unpack_ret(ret)

    grad = elstruct.reader.gradient(prog, out_str)

    cnf_fs[-1].create(cnf_locs)
    cnf_path = cnf_fs[-1].path(cnf_locs)
    print(f" - Saving at {cnf_path}")

    cnf_fs[-1].file.gradient_info.write(inf_obj, cnf_locs)
    cnf_fs[-1].file.gradient_input.write(inp_str, cnf_locs)
    cnf_fs[-1].file.gradient.write(grad, cnf_locs)


def _save_zmatrix_parsed(zma, inf_obj, inp_str, zma_fs, zma_locs):
    """ save zmatrix from zmat object
    """
    zma_fs[-1].create(zma_locs)
    zma_path = zma_fs[-1].path(zma_locs)
    print(f" - Saving at {zma_path}")

    zma_fs[-1].file.geometry_info.write(inf_obj, zma_locs)
    zma_fs[-1].file.geometry_input.write(inp_str, zma_locs)
    zma_fs[-1].file.zmatrix.write(zma, zma_locs)


def _save_zmatrix(ret, zma_fs, zma_locs, init_zma=None):
    """ Read a zma from output. If one is not found,
        then attempt to convert it from geometry.

        might not work if redundants print like zma
    """

    print(" - Reading Z-Matrix from output...")
    inf_obj, inp_str, out_str, prog, _ = _unpack_ret(ret)
    zma = None
    if init_zma is not None:
        print('using opt geo fro zma')
        geo = elstruct.reader.opt_geometry(prog, out_str)
        zma = read_zma_from_geo(init_zma, geo)
    if zma is None:
        zma = read_job_zma(ret, init_zma=init_zma)
    _save_zmatrix_parsed(zma, inf_obj, inp_str, zma_fs, zma_locs)


def _save_energy_parsed(ene, inf_obj, inp_str, sp_fs, sp_locs):
    """ Save an energy that has been parsed
    """
    sp_fs[-1].create(sp_locs)
    sp_path = sp_fs[-1].path(sp_locs)
    print(f" - Saving at {sp_path}")

    sp_fs[-1].file.input.write(inp_str, sp_locs)
    sp_fs[-1].file.info.write(inf_obj, sp_locs)
    sp_fs[-1].file.energy.write(ene, sp_locs)


def _save_energy(ret, sp_fs, sp_locs):
    """ Sving a energy
    """

    print(" - Reading energy from output...")
    inf_obj, inp_str, out_str, prog, method = _unpack_ret(ret)

    ene = elstruct.reader.energy(prog, method, out_str)
    _save_energy_parsed(ene, inf_obj, inp_str, sp_fs, sp_locs)


def _save_hessian_parsed(hess, freqs, inf_obj, inp_str, cnf_fs, cnf_locs):
    """ save a hessian that has been parsed
    """
    cnf_fs[-1].create(cnf_locs)
    cnf_path = cnf_fs[-1].path(cnf_locs)
    print(f" - Saving at {cnf_path}")

    cnf_fs[-1].file.hessian_info.write(inf_obj, cnf_locs)
    cnf_fs[-1].file.hessian_input.write(inp_str, cnf_locs)
    cnf_fs[-1].file.hessian.write(hess, cnf_locs)
    cnf_fs[-1].file.harmonic_frequencies.write(freqs, cnf_locs)


def _save_hessian(ret, cnf_fs, cnf_locs):
    """ saved the hessian
    """

    print(" - Reading hessian and harmonic frequencies from output...")
    inf_obj, inp_str, out_str, prog, _ = _unpack_ret(ret)

    hess = elstruct.reader.hessian(prog, out_str)
    freqs = elstruct.reader.harmonic_frequencies(prog, out_str)
    _save_hessian_parsed(hess, freqs, inf_obj, inp_str, cnf_fs, cnf_locs)


def _save_rotors(zma_fs, zma_locs, zrxn=None):
    """ Save the rotors

        Reading the ZMA from the filesystem to build rotors
    """

    zma = zma_fs[-1].file.zmatrix.read(zma_locs)
    # print('zma')
    # print(automol.zmat.string(zma))
    # print('\nzrxn')
    # print(zrxn)
    gra = None if zrxn is None else automol.reac.ts_graph(zrxn)
    rotors = automol.data.rotor.rotors_from_zmatrix(zma, gra=gra)
    if rotors:
        tors_lst = automol.data.rotor.rotors_torsions(rotors, sort=True)
        zma_path = zma_fs[-1].path(zma_locs)
        print(f" - Rotors identified from Z-Matrix at {zma_path}")
        print(" - Saving rotor information at same location.")
        zma_fs[-1].file.torsions.write(tors_lst, zma_locs)


def _save_rings(zma_fs, zma_locs, zrxn=None):
    """ Save rings
    """

    zma = zma_fs[-1].file.zmatrix.read(zma_locs)
    tsg = None if zrxn is None else automol.reac.ts_graph(zrxn)
    rings_atoms = automol.zmat.all_rings_atoms(zma, tsg=tsg)
    ring_dct = automol.zmat.all_rings_dct(zma, rings_atoms)
    if ring_dct:
        zma_path = zma_fs[-1].path(zma_locs)
        print(f" - Ring torsions identified from Z-Matrix at {zma_path}")
        print(" - Saving ring torsions information at same location.")
        zma_fs[-1].file.ring_torsions.write(ring_dct, zma_locs)


def _save_reaction(zma_fs, zma_locs, zrxn=None):
    """ Save the reaction object
    """

    if zrxn is not None:
        zma_path = zma_fs[-1].path(zma_locs)
        print(f" - Saving Reaction Class+Graph object at {zma_path}")
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
    print(f" - Saving at {instab_path}")
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
