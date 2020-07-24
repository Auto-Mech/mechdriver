"""
  Functions to read the filesystem and pull objects from it
"""

import sys
import automol
import autofile


def get_zma_geo(filesys, locs):
    """ Get the geometry and zmatrix from a filesystem
    """
    if filesys[-1].file.zmatrix.exists(locs):
        zma = filesys[-1].file.zmatrix.read(locs)
    else:
        zma = None

    if filesys[-1].file.geometry.read(locs):
        geo = filesys[-1].file.geometry.read(locs)
    else:
        geo = None

    # Check
    if zma is None and geo is None:
        print('*ERROR: Neither a Z-Matrix or a Cartesian Geom exists level')
        sys.exit()

    return zma, geo


def min_energy_conformer_locators(cnf_save_fs, mod_thy_info,
                                  zpe_corrd=False):
    """ locators for minimum energy conformer
    """

    cnf_locs_lst = cnf_save_fs[-1].existing()
    if cnf_locs_lst:
        cnf_enes = []
        for locs in cnf_locs_lst:
            cnf_path = cnf_save_fs[-1].path(locs)
            sp_fs = autofile.fs.single_point(cnf_path)
            cnf_enes.append(
                sp_fs[-1].file.energy.read(mod_thy_info[1:4]))
        if zpe_corrd:
            pass
            # # How to calculate the ZPE (just use harm zpe?A
            # # Ignore confs without zpe? or break the run
            # cnf_hessians = [cnf_save_fs[-1].file.hessian.read(locs)
            #             for locs in cnf_locs_lst]
            # cnf_corrs = [ene+zpe for ene, zpe
            # in zip(cnf_enes, cnf_zpes)]
        else:
            min_cnf_locs = cnf_locs_lst[cnf_enes.index(min(cnf_enes))]
    else:
        min_cnf_locs = []
    return min_cnf_locs


def min_dist_conformer_zma(dist_name, cnf_save_fs):
    """ locators for minimum energy conformer """
    cnf_locs_lst = cnf_save_fs[-1].existing()
    cnf_zmas = []
    for locs in cnf_locs_lst:
        zma_fs = autofile.fs.manager(cnf_save_fs[-1].path(locs), 'ZMATRIX')
        cnf_zmas.append(zma_fs[-1].file.zmatrix.read([0]))
    min_dist = 100.
    min_zma = []
    for zma in cnf_zmas:
        dist = automol.zmatrix.values(zma)[dist_name]
        if dist < min_dist:
            min_dist = dist
            min_zma = zma
    min_zma = [min_zma]
    return min_zma


def min_dist_conformer_zma_geo(dist_coords, cnf_save_fs):
    """ locators for minimum energy conformer """
    cnf_locs_lst = cnf_save_fs[-1].existing()
    cnf_zmas = []
    for locs in cnf_locs_lst:
        zma_fs = autofile.fs.manager(cnf_save_fs[-1].path(locs), 'ZMATRIX')
        cnf_zmas.append(zma_fs[-1].file.zmatrix.read([0]))
    min_dist = 100.
    min_zma = []
    for zma in cnf_zmas:
        zmas, _ = automol.zmatrix.shifted_standard_zmas_graphs([zma])
        zma = zmas[0]
        geo = automol.zmatrix.geometry(zma)
        dist = automol.geom.distance(geo, *list(dist_coords))
        if dist < min_dist:
            min_dist = dist
            min_zma = zma
    min_zma = [min_zma]
    return min_zma


def locs_sort(save_fs):
    """ sort trajectory file according to energies
    """
    locs_lst = save_fs[-1].existing()
    if locs_lst:
        enes = [save_fs[-1].file.energy.read(locs)
                for locs in locs_lst]
        sorted_locs = []
        for _, loc in sorted(zip(enes, locs_lst), key=lambda x: x[0]):
            sorted_locs.append(loc)
    return sorted_locs


def traj_sort(save_fs, mod_thy_info):
    """ sort trajectory file according to energies
    """
    locs_lst = save_fs[-1].existing()
    if locs_lst:
        enes = []
        for locs in locs_lst:
            cnf_path = save_fs[-1].path(locs)
            sp_fs = autofile.fs.single_point(cnf_path)
            enes.append(
                sp_fs[-1].file.energy.read(mod_thy_info[1:4]))
        # enes = [save_fs[-1].file.energy.read(locs)
        #         for locs in locs_lst]
        geos = [save_fs[-1].file.geometry.read(locs)
                for locs in locs_lst]
        traj = []
        traj_sort_data = sorted(zip(enes, geos, locs_lst), key=lambda x: x[0])
        for ene, geo, locs in traj_sort_data:
            comment = 'energy: {0:>15.10f} \t {1}'.format(ene, locs[0])
            traj.append((comment, geo))
        traj_path = save_fs[0].file.trajectory.path()
        print("Updating trajectory file at {}".format(traj_path))
        save_fs[0].file.trajectory.write(traj)
