""" drivers for coordinate scans
"""

import numpy

import automol
import autofile
from phydat import phycon
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo
from mechanalyzer.inf import rxn as rinfo
from mechlib.filesys._build import build_fs
from mechlib.filesys.mincnf import min_energy_conformer_locators


def potential(names, grid_vals, cnf_save_path,
              mod_tors_ene_info, ref_ene,
              constraint_dct,
              read_geom=False, read_grad=False,
              read_hess=False, read_zma=False):
    """ Get the potential for a hindered rotor
    """

    # Build initial lists for storing potential energies and Hessians
    # grid_points = automol.pot.points(grid_vals)
    grid_coords = automol.pot.coords(grid_vals)
    back_coords = tuple([tuple([
        val + 4*numpy.pi for val in grid]) for grid in grid_coords])
    pot, geoms, grads, hessians, zmas, paths = {}, {}, {}, {}, {}, {}

    # Set up filesystem information
    zma_fs = autofile.fs.zmatrix(cnf_save_path)
    zma_path = zma_fs[-1].path([0])
    if constraint_dct is None:
        scn_fs = autofile.fs.scan(zma_path)
    else:
        scn_fs = autofile.fs.cscan(zma_path)

    # Read the filesystem
    for idx, vals in enumerate(grid_coords):

        # Get angles in degrees for potential for now
        vals_conv = tuple(val*phycon.RAD2DEG for val in vals)

        # Get locs for reading filesysten
        locs = [names, vals]
        back_locs = [names, back_coords[idx]]
        if constraint_dct is not None:
            locs = [constraint_dct] + locs
            back_locs = [constraint_dct] + back_locs

        # Read values of interest
        ene = energy(scn_fs, locs, mod_tors_ene_info)
        back_ene = energy(scn_fs, back_locs, mod_tors_ene_info)
        step_ene = None
        if ene is not None:
            if back_ene is not None:
                step_ene = min(ene, back_ene)
            else:
                step_ene = ene
        elif back_ene is not None:
            step_ene = back_ene

        if step_ene is not None:
            pot[vals_conv] = (step_ene - ref_ene) * phycon.EH2KCAL
        else:
            pot[vals_conv] = -10.0

        if read_geom:
            if scn_fs[-1].file.geometry.exists(locs):
                geoms[vals_conv] = scn_fs[-1].file.geometry.read(locs)
            else:
                geoms[vals_conv] = None

        if read_grad:
            if scn_fs[-1].file.gradient.exists(locs):
                grads[vals_conv] = scn_fs[-1].file.gradient.read(locs)
            else:
                grads[vals_conv] = None

        if read_hess:
            if scn_fs[-1].file.hessian.exists(locs):
                hessians[vals_conv] = scn_fs[-1].file.hessian.read(locs)
            else:
                hessians[vals_conv] = None

        if read_zma:
            if scn_fs[-1].file.zmatrix.exists(locs):
                zmas[vals_conv] = scn_fs[-1].file.zmatrix.read(locs)
            else:
                zmas[vals_conv] = None

        paths[vals] = scn_fs[-1].path(locs)

    return pot, geoms, grads, hessians, zmas, paths


# Single data point readers
def geometry(cnf_save_fs, mod_thy_info, conf='sphere'):
    """ get the geometry
    """

    assert conf in ('minimum', 'sphere')

    # Read the file system
    if conf == 'minimum':
        geom = _min_energy_conformer(cnf_save_fs, mod_thy_info)
    elif conf == 'sphere':
        geom = _spherical_conformer(cnf_save_fs)

    return geom


def _min_energy_conformer(cnf_save_fs, mod_thy_info):
    """ Reads the minimum-energy conformer from the save FileSystem
    """

    ini_loc_info = min_energy_conformer_locators(
        cnf_save_fs, mod_thy_info)
    locs, path = ini_loc_info
    if path:
        min_conf = cnf_save_fs[-1].file.geometry.read(locs)
    else:
        min_conf = None

    return min_conf


def _spherical_conformer(cnf_save_fs):
    """ Reads the conformer from the Save FileSystem that is the most
        spherical.
    """

    cnf_locs_lst = cnf_save_fs[-1].existing()
    if cnf_locs_lst:
        cnf_geoms = [cnf_save_fs[-1].file.geometry.read(locs)
                     for locs in cnf_locs_lst]
        round_geom = automol.geom.minimum_volume_geometry(cnf_geoms)
    else:
        round_geom = None

    return round_geom


def energy(filesys, locs, mod_tors_ene_info):
    """ Read the energy from an SP filesystem that is located in some
        root 'filesys object'
    """

    if filesys[-1].exists(locs):
        path = filesys[-1].path(locs)
        sp_fs = autofile.fs.single_point(path)
        if sp_fs[-1].file.energy.exists(mod_tors_ene_info[1:4]):
            ene = sp_fs[-1].file.energy.read(mod_tors_ene_info[1:4])
        else:
            ene = None
    else:
        ene = None

    return ene


def reaction(rxn_info, ini_thy_info, zma_locs, save_prefix, ts_locs=(0,)):
    """ Check if reaction exists in the filesystem and has been identified
    """

    zrxn, zma = None, None

    sort_rxn_info = rinfo.sort(rxn_info, scheme='autofile')
    ts_info = rinfo.ts_info(rxn_info)
    mod_ini_thy_info = tinfo.modify_orb_label(ini_thy_info, ts_info)

    rxn_fs = autofile.fs.reaction(save_prefix)
    if rxn_fs[-1].exists(sort_rxn_info):
        _, cnf_save_fs = build_fs(
            save_prefix, save_prefix, 'CONFORMER',
            rxn_locs=sort_rxn_info,
            thy_locs=mod_ini_thy_info[1:],
            # this needs to be fixed for any case with more than one TS
            ts_locs=ts_locs)

        _, ini_min_cnf_path = min_energy_conformer_locators(
            cnf_save_fs, mod_ini_thy_info)
        if ini_min_cnf_path:
            zma_fs = autofile.fs.zmatrix(ini_min_cnf_path)
            if zma_fs[-1].file.reaction.exists(zma_locs):
                zrxn = zma_fs[-1].file.reaction.read(zma_locs)
                zma = zma_fs[-1].file.zmatrix.read(zma_locs)

        if zrxn is None:
            _, zma_fs = build_fs(
                '', save_prefix, 'ZMATRIX',
                rxn_locs=sort_rxn_info, ts_locs=ts_locs,
                thy_locs=mod_ini_thy_info[1:])

            if zma_fs[-1].file.reaction.exists(zma_locs):
                zrxn = zma_fs[-1].file.reaction.read(zma_locs)
                zma = zma_fs[-1].file.zmatrix.read(zma_locs)

    return zrxn, zma


def instability_transformation(spc_dct, spc_name, thy_info, save_prefix,
                               zma_locs=(0,)):
    """ see if a species and unstable and handle task management
    """

    spc_info = sinfo.from_dct(spc_dct[spc_name])
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)

    _, cnf_save_fs = build_fs(
        '', save_prefix, 'CONFORMER',
        spc_locs=spc_info,
        thy_locs=mod_thy_info[1:])

    # Check if any locs exist first?
    ini_loc_info = min_energy_conformer_locators(
        cnf_save_fs, mod_thy_info)
    _, min_cnf_path = ini_loc_info

    zma_save_fs = autofile.fs.zmatrix(min_cnf_path)

    # Check if the instability files exist
    if zma_save_fs[-1].file.instability.exists(zma_locs):
        instab_trans = zma_save_fs[-1].file.instability.read(zma_locs)
        zma = zma_save_fs[-1].file.zmatrix.read(zma_locs)
        _instab = (instab_trans, zma)
        path = zma_save_fs[-1].file.zmatrix.path(zma_locs)
    else:
        _instab = None
        path = None

    return _instab, path


def energy_trans(etrans_save_fs, etrans_locs):
    """ Read out the the enery transfer parameters in the filesys
    """

    nsamp = etrans_save_fs[-1].file.read(etrans_locs)
    epsilon = etrans_save_fs[-1].file.read.epsilon(etrans_locs)
    sigma = etrans_save_fs[-1].file.read.sigma(etrans_locs)
    # alpha = etrans_fs[-1].file.alpha.read(etrans_locs)

    min_geo_traj = etrans_save_fs[-1].file.read.min_geos(etrans_locs)
    min_geos = automol.geom.from_xyz_trajectory_string(min_geo_traj)

    return nsamp, epsilon, sigma, min_geos
