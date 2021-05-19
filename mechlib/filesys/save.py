"""
 functions for reading and writing the filesystem
"""

import automol
import elstruct
import autofile
from mechroutines.es import runner as es_runner


def conformer(opt_ret, hess_ret, cnf_fs, thy_locs,
              zrxn=None, init_zma=None,
              rng_locs=None, tors_locs=None, zma_locs=None):
    """ Save a conformer and relevant information.

        thy_fs = (_fs, _locs)
    """

    # Build filesystem locs and objects
    cnf_locs, zma_locs = _generate_locs(
        rng_locs=rng_locs, tors_locs=tors_locs, zma_locs=zma_locs)

    # cnf_fs = autofile.fs.conformer(PREFIX)
    # make paths and pass?
    zma_fs = autofile.fs.zmatrix(cnf_fs[-1].path(cnf_locs))
    sp_fs = autofile.fs.single_point(cnf_fs[-1].path(cnf_locs))

    # Save data
    _save_geom(opt_ret, cnf_fs, cnf_locs)
    _save_zmatrix(opt_ret, zma_fs, zma_locs)
    _save_rotors(_, zma_fs, zma_locs)  # Need a zma, prob just read from filesys?
    _save_reaction(zma_fs, zma_locs, zrxn=zrxn)
    _save_energy(opt_ret, sp_fs, thy_locs)

    if hess_ret is not None:
        _save_hessian(hess_ret, cnf_fs, cnf_locs)


def sym_indistinct_conformer(geo, cnf_fs, cnf_tosave_locs, cnf_saved__locs):
    """
    """

    # Set the path to the conf filesys with sym similar
    cnf_save_path = cnf_save_fs[-1].path(cnf_saved_locs)

    # Build the sym file sys
    sym_save_fs = autofile.fs.symmetry(cnf_save_path)
    sym_save_path = cnf_save_fs[-1].path(cnf_saved_locs)
    ioprinter.save_symmetry(sym_save_path)
    sym_save_fs[-1].create([cnf_tosave_locs[-1]])
    sym_save_fs[-1].file.geometry.write(geo, [cnf_tosave_locs[-1]])


def hindered_rotor_point(opt_ret, scn_fs, scn_locs, job):
    """ save info for the hindered rotor
    """

    _save_geom(opt_ret, scn_fs, scn_locs)
    _save_zmatrix(opt_ret, scn_fs, scn_locs)


def init_cnf_samp(cnf_save_fs, cnf_locs):
    """ init cnf samp
    """

    inf_obj = autofile.schema.info_objects.conformer_trunk(0)
    inf_obj.nsamp = 1

    if not cnf_save_fs[0].info.exists():
        cnf_save_fs[0].file.info.write(inf_obj)
    if not cnf_save_fs[1].file.info.exists([locs[0]]):
        cnf_save_fs[1].file.info.write(inf_obj, [locs[0]])


# Saving other
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


# Job readers
# Constituent functions for saving various bits of information
def read_job_zma(ret, init_zma=None, rebuild=False):
    """ trys to read a zma from a job and does processing on it as needed
    """

    _, _, out_str, prog, _ = _unpack_ret(ret)
    zma = elstruct.reader.opt_zmatrix(prog, out_str)
    if zma is None or rebuild:
        print('Getting ZMA from a geometry...')
        geo = elstruct.reader.opt_geometry(prog, out_str)
        if init_zma is not None:
            zma = rebuild_zma_from_opt_geo(init_zma, geo)
        else:
            zma = automol.geom.zmatrix(geo)

    return zma


def rebuild_zma_from_opt_geo(init_zma, opt_geo):
    """ build zma from opt
    """
    dummy_key_dct = automol.zmat.dummy_key_dictionary(init_zma)
    geo_wdummy = automol.geom.insert_dummies(geo, dummy_key_dct)
    zma = automol.zmat.from_geometry(zma, geo_wdummy)

    return zma


def _save_geom(ret, cnf_fs, cnf_locs):
    """ Saving a Cartesian geometry. could be for cnf fs or scn fs
    """

    print(" - Reading geometry from output...")
    inf_obj, inp_str, out_str, prog, _ = _unpack_ret(ret)

    geo = elstruct.reader.opt_geometry(prog, out_str)

    cnf_fs[-1].create(cnf_locs)
    cnf_save_path = cnf_fs[-1].path(cnf_locs)
    print(" - Saving at {}".format(cnf_save_path))
    cnf_fs[-1].file.geometry_info.write(inf_obj, cnf_locs)
    cnf_fs[-1].file.geometry_input.write(inp_str, cnf_locs)
    cnf_fs[-1].file.geometry.write(geo, cnf_locs)


def _save_grad(ret, cnf_fs, cnf_locs):
    """ Saving a geometry
    """

    print(" - Reading geometry from output...")
    inf_obj, inp_str, out_str, prog, _ = _unpack_ret(ret)

    grad = elstruct.reader.gradient(prog, out_str)

    cnf_fs[-1].create(cnf_locs)
    cnf_save_path = cnf_fs[-1].path(cnf_locs)
    print(" - Saving at {}".format(cnf_save_path))
    cnf_fs[-1].file.gradient_info.write(inf_obj, cnf_locs)
    cnf_fs[-1].file.gradient_input.write(inp_str, cnf_locs)
    cnf_fs[-1].file.gradient.write(grad, cnf_locs)


def _save_zmatrix(ret, zma_fs, zma_locs, init_zma=None):
    """ Read a zma from output. If one is not found,
        then attempt to convert it from geometry.

        might not work if redundants print like zma
    """

    print(" - Reading geometry from output...")
    inf_obj, inp_str, out_str, prog, _ = _unpack_ret(ret)

    zma = read_job_zma(ret, init_zma=None, rebuild=False)

    zma_fs[-1].create(zma_locs)
    zma_fs[-1].file.geometry_info.write(inf_obj, zma_locs)
    zma_fs[-1].file.geometry_input.write(inp_str, zma_locs)
    zma_fs[-1].file.zmatrix.write(zma, zma_locs)


def _save_energy(ret, sp_fs, sp_locs):
    """ Sving a energy
    """

    print(" - Reading energy from output...")
    inf_obj, inp_str, out_str, prog, method = _unpack_ret(ret)

    ene = elstruct.reader.energy(prog, method, out_str)

    print(" - Saving energy of unique geometry...")
    sp_fs[-1].create(sp_locs)
    sp_fs[-1].file.input.write(inp_str, sp_locs)
    sp_fs[-1].file.info.write(inf_obj, sp_locs)
    sp_fs[-1].file.energy.write(ene, sp_locs)


def _save_hessian(ret, cnf_fs, cnf_locs):
    """ saved the hessian
    """

    print(" - Reading hessian from output...")
    inf_obj, inp_str, out_str, prog, _ = _unpack_ret(ret)

    hess = elstruct.reader.hessian(prog, out_str)

    cnf_fs[-1].file.harmonic_frequencies.write(freqs, cnf_locs)
    # neg, pos = automol.util.separate_negatives(lst)
    # freqs = sorted([-1.0*val for val in imags] + freqs)
    # print('TS freqs: {}'.format(' '.join(str(freq) for freq in freqs)))

    cnf_fs[-1].create(cnf_locs)
    cnf_save_path = cnf_fs[-1].path(cnf_locs)
    print(" - Saving at {}".format(cnf_save_path))
    cnf_fs[-1].file.hessian_info.write(inf_obj, cnf_locs)
    cnf_fs[-1].file.hessian_input.write(inp_str, cnf_locs)
    cnf_fs[-1].file.hessian.write(hess, cnf_locs)
    cnf_fs[-1].file.harmonic_frequencies.write(freqs, cnf_locs)


def _save_rotors(zma, zma_fs, zma_locs):
    """ Save the rotors
    """

    rotors = automol.rotor.from_zmatrix(zma)
    if any(rotors):
        print(" - Species has rotors, saving them...")
        zma_fs[-1].file.torsions.write(rotors, zma_locs)


def _save_ring_tors(zma, zma_fs, zma_locs):
    """ Save the rotors
    """

    rings_atoms =  _get_ring_atoms(zma, zrxn)
    tors_dct = {}
    for ring_atoms in rings_atoms:
        dct_label = '-'.join(str(atm+1) for atm in ring_atoms)
        samp_range_dct = _get_ring_samp_ranges(zma, ring_atoms)
        tors_dct[dct_label] = samp_range_dct
    if tors_dct:
        print(" - Species has ring torsions, saving them...")
        zma_save_fs[-1].file.ring_torsions.write(tors_dct, zma_locs)


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
    inf_obj, inp_str, out_str, prog, _ = _unpack_ret(ret)

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
