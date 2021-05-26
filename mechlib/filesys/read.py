""" drivers for coordinate scans
"""

import automol
import autofile
from phydat import phycon
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo
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
    grid_points = automol.pot.points(grid_vals)
    grid_coords = automol.pot.coords(grid_vals)
    pot, geoms, grads, hessians, zmas, paths = {}, {}, {}, {}, {}, {}

    # Set up filesystem information
    zma_fs = autofile.fs.zmatrix(cnf_save_path)
    zma_path = zma_fs[-1].path([0])
    if constraint_dct is None:
        scn_fs = autofile.fs.scan(zma_path)
    else:
        scn_fs = autofile.fs.cscan(zma_path)

    # Read the energies and Hessians from the filesystem
    for point, vals in zip(grid_points, grid_coords):

        locs = [names, vals]
        if constraint_dct is not None:
            locs = [constraint_dct] + locs

        ene = energy(scn_fs, locs, mod_tors_ene_info)
        if ene is not None:
            pot[point] = (ene - ref_ene) * phycon.EH2KCAL
        else:
            pot[point] = -10.0

        if read_geom:
            if scn_fs[-1].file.geometry.exists(locs):
                geoms[point] = scn_fs[-1].file.geometry.read(locs)
            else:
                geoms[point] = None

        if read_grad:
            if scn_fs[-1].file.gradient.exists(locs):
                grads[point] = scn_fs[-1].file.gradient.read(locs)
            else:
                grads[point] = None

        if read_hess:
            if scn_fs[-1].file.hessian.exists(locs):
                hessians[point] = scn_fs[-1].file.hessian.read(locs)
            else:
                hessians[point] = None

        if read_zma:
            if scn_fs[-1].file.zmatrix.exists(locs):
                zmas[point] = scn_fs[-1].file.zmatrix.read(locs)
            else:
                zmas[point] = None

        paths[point] = scn_fs[-1].path(locs)

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
    if zma_save_fs[-1].file.reaction.exists(zma_locs):
        zrxn = zma_save_fs[-1].file.reaction.read(zma_locs)
        zma = zma_save_fs[-1].file.zmatrix.read(zma_locs)
        _instab = (zrxn, zma)
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

    min_geo_traj = etrans_save_fs[-1].file.read.min_geos(etrans_locs)
    min_geos = automol.geom.from_xyz_trajectory_string(min_geo_traj)

    return nsamp, epsilon, sigma, min_geos
