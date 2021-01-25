"""
  Functions to read the filesystem and pull objects from it
"""

import sys
import automol
import autofile
from phydat import phycon


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


def min_energy_conformer_locators(cnf_save_fs, mod_thy_info):
    """ wrapper for minimum locs
    """
    locs, paths = conformer_locators(
        cnf_save_fs, mod_thy_info, cnf_range='min')
    if locs and paths:
        ret = locs[0], paths[0]
    else:
        ret = [], []
    return ret


def conformer_locators(cnf_save_fs, mod_thy_info, cnf_range='min'):
    """ locators for minimum energy conformer
    """

    cnf_locs_lst = cnf_save_fs[-1].existing()
    fin_locs_lst, fin_paths_lst = [], []

    if cnf_locs_lst:

        cnf_locs_lst, cnf_enes_lst = _sorted_cnf_lsts(
            cnf_locs_lst, cnf_save_fs, mod_thy_info)

        if cnf_range == 'min':
            fin_locs_lst = [cnf_locs_lst[0]]
        elif cnf_range == 'all':
            fin_locs_lst = cnf_locs_lst
        elif 'e' in cnf_range:
            fin_locs_lst = _erange_locs(cnf_locs_lst, cnf_enes_lst, cnf_range)
        elif 'n' in cnf_range:
            fin_locs_lst = _nrange_locs(cnf_locs_lst, cnf_range)

    else:
        print('No conformers located in {}'.format(cnf_save_fs[-1].root.path()))
    for locs in fin_locs_lst:
        fin_paths_lst.append(cnf_save_fs[-1].path(locs))

    return fin_locs_lst, fin_paths_lst


def _sorted_cnf_lsts(cnf_locs_lst, cnf_save_fs, mod_thy_info):
    """ Get a list of conformer locs and energies, sorted by energies
    """

    cnf_enes_lst = []
    if len(cnf_locs_lst) == 1:
        cnf_enes_lst = [10]
    else:
        for locs in cnf_locs_lst:
            cnf_path = cnf_save_fs[-1].path(locs)
            sp_fs = autofile.fs.single_point(cnf_path)
            cnf_enes_lst.append(sp_fs[-1].file.energy.read(mod_thy_info[1:4]))

    # Sort the cnf locs and cnf enes
    cnf_enes_lst, cnf_locs_lst = zip(*sorted(zip(cnf_enes_lst, cnf_locs_lst)))

    return cnf_locs_lst, cnf_enes_lst


def _erange_locs(cnf_locs, cnf_enes, ethresh):
    """ Get a range of e values
    """

    thresh = float(ethresh.split('e')[1])

    min_cnf_locs = []
    min_ene = cnf_enes[0]
    for locs, ene in zip(cnf_locs, cnf_enes):
        rel_ene = (ene - min_ene) * phycon.EH2KCAL
        print('rel_ene:', rel_ene, thresh)
        if rel_ene <= thresh:
            min_cnf_locs.append(locs)

    return min_cnf_locs


def _nrange_locs(cnf_locs, nthresh):
    """ Get a range of n values
    """

    thresh = int(nthresh.split('n')[1])

    min_cnf_locs = []
    for idx, locs in enumerate(cnf_locs):
        if idx+1 <= thresh:
            min_cnf_locs.append(locs)

    return min_cnf_locs


def min_dist_conformer_zma(dist_name, cnf_save_fs):
    """ locators for minimum energy conformer """
    cnf_locs_lst = cnf_save_fs[-1].existing()
    cnf_zmas = []
    for locs in cnf_locs_lst:
        zma_fs = autofile.fs.zmatrix(cnf_save_fs[-1].path(locs))
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
        zma_fs = autofile.fs.zmatrix(cnf_save_fs[-1].path(locs))
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
