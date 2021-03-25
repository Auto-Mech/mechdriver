""" drivers for coordinate scans
"""

import itertools
import automol
import autorun
import autofile
from autofile import fs
from phydat import phycon


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
    zma_fs = fs.zmatrix(cnf_save_path)
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
    _, zma_save_fs = build_fs(
        '', save_prefix, 'ZMATRIX',
        spc_locs=spc_info,
        thy_locs=mod_thy_info[1:],
        instab_locs=())

    # Check if the instability files exist
     if zma_save_fs[-1].file.reaction.exists(zma_locs):
         zrxn = zma_save_fs[-1].file.reaction.exists(zma_locs))
         zma = zma_save_fs[-1].file.zmatrix.exists(zma_locs))
         _instab = (zrxn, zma)
         print('- Found files denoting species instability at path')
         print('    {}'.format(zma_fs[-1].path(zma_locs)))
     else:
         _instab = None
         print('- No files denoting instability were found at path')
         print('    {}'.format(zma_fs[-1].path(zma_locs)))

    return _instab
