"""
  Functions to read the filesystem and pull objects from it
"""

import sys
import automol
import autofile
from phydat import phycon
from mechlib.filesys import build_fs

# def min_energy_ring_conformer_locators(cnf_save_fs, mod_thy_info):
#     """ wrapper for minimum locs
#     """
#     rng_locs_lst = cnf_save_fs[-1].existing()
#     rng_cnf_locs_lst = []
#     rng_cnf_enes_lst = []
#     for rng_locs in rng_locs_lst:
#         rng_save_path = rng_save_fs[-1].path(rng_locs)
#         _, rng_cnf_save_fs = build_fs(
#             rng_save_path, rng_save_path, 'CONFORMER')
#         locs, paths = conformer_locators(
#             rng_cnf_save_fs, mod_thy_info, cnf_range='min')
#         if locs and paths:
#             cnf_locs_lst, cnf_enes_lst = _sorted_cnf_lsts(
#                 locs, rng_cnf_save_fs, mod_thy_info)
#             rng_cnf_locs_lst.append((rng_locs, cnf_locs_lst[0],))
#             rng_cnf_enes_lst.append(cnf_enes_lst[0])
#     if rng_cnf_locs_lst:
#         rng_cnf_enes_lst, rng_cnf_locs_lst = zip(*sorted(zip(rng_cnf_enes_lst, rng_cnf_locs_lst)))
#         locs = rng_cnf_locs_lst[0]
#         rng_save_path = rng_save_fs[-1].path(locs[0])
#         _, rng_cnf_save_fs = build_fs(
#             rng_save_path, rng_save_path, 'CONFORMER')
#         rng_cnf_save_path = rng_cnf_save_fs[-1].path(locs[1])
#         paths = (rng_save_path, rng_cnf_save_path)
#         ret = locs, paths
#     else:
#         ret = ('',''), ('', '')
#     return ret


# def ring_conformer_locators(cnf_save_fs, mod_thy_info, cnf_range='min'):
#    """ wrapper for minimum locs
#    """
#    locs_lst = cnf_save_fs[-1].existing()
#    rng_cnf_locs_lst = []
#    rng_cnf_enes_lst = []
#    if locs_lst:
#        rng_cnf_locs_lst, rng_cnf_enes_lst = _sorted_cnf_lsts(
#            locs_lst, cnf_save_fs, mod_thy_info)
#    if rng_cnf_locs_lst:
#        rng_cnf_enes_lst, rng_cnf_locs_lst = zip(*sorted(zip(rng_cnf_enes_lst, rng_cnf_locs_lst)))
#        if cnf_range == 'min':
#            fin_locs_lst = [rng_cnf_locs_lst[0]]
#        elif cnf_range == 'all':
#            fin_locs_lst = rng_cnf_locs_lst
#        elif 'e' in cnf_range:
#            fin_locs_lst = _erange_locs(rng_cnf_locs_lst, rng_cnf_enes_lst, cnf_range)
#        elif 'n' in cnf_range:
#            fin_locs_lst = _nrange_locs(rng_cnf_locs_lst, cnf_range)
#        fin_paths = []
#        for locs in fin_locs_lst:
#            rng_locs, cnf_locs = locs
#            rng_save_path = rng_save_fs[-1].path(rng_locs)
#            _, rng_cnf_save_fs = build_fs(
#                rng_save_path, rng_save_path, 'CONFORMER')
#            rng_cnf_save_path = rng_cnf_save_fs[-1].path(cnf_locs)
#            fin_paths.append((rng_save_path, rng_cnf_save_path))
#        ret = fin_locs_lst, fin_paths
#    else:
#        ret = [([], [])], [('', '')]
#    return ret


def min_energy_conformer_locators(cnf_save_fs, mod_thy_info):
    """ wrapper for minimum locs
    """
    locs, paths = conformer_locators(
        cnf_save_fs, mod_thy_info, cnf_range='min')
    if locs and paths:
        ret = locs[0], paths[0]
    else:
        ret = ['',''], ''
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
        elif 'r' in cnf_range:
            fin_locs_lst = _rrange_locs(cnf_locs_lst, cnf_range)

    else:
        print('No conformers located in {}'.format(
            cnf_save_fs[0].path()))
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


def _rrange_locs(cnf_locs, nthresh):
    """ Get a range of n values
    """

    thresh = int(nthresh.split('r')[1])

    min_cnf_locs = []
    used_rids = []
    for idx, locs in enumerate(cnf_locs):
        if not list(locs)[0] in used_rids:
            if idx+1 <= thresh:
                min_cnf_locs.append(locs)
                used_rids.append(list(locs)[0])
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


def traj_sort(save_fs, mod_thy_info, rid=None):
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
        geos = [save_fs[-1].file.geometry.read(locs)
                for locs in locs_lst]
        traj = []
        traj_sort_data = sorted(zip(enes, geos, locs_lst), key=lambda x: x[0])
        for ene, geo, locs in traj_sort_data:
            comment = 'energy: {0:>15.10f} \t {1}'.format(ene, locs[0])
            traj.append((geo, comment))
        traj_path = save_fs[0].file.trajectory.path()
        print("Updating trajectory file at {}".format(traj_path))
        save_fs[0].file.trajectory.write(traj)

        if rid is not None:
            locs_lst = save_fs[-1].existing()
            if locs_lst:
                enes = []
                geos = []
                for locs in locs_lst:
                    trid, _ = locs
                    if trid == rid:
                        cnf_path = save_fs[-1].path(locs)
                        sp_fs = autofile.fs.single_point(cnf_path)
                        enes.append(
                            sp_fs[-1].file.energy.read(mod_thy_info[1:4]))
                        geos.append(save_fs[-1].file.geometry.read(locs))
                traj = []
                traj_sort_data = sorted(zip(enes, geos, locs_lst), key=lambda x: x[0])
                for ene, geo, locs in traj_sort_data:
                    comment = 'energy: {0:>15.10f} \t {1} {2}'.format(ene, locs[0], locs[1])
                    traj.append((geo, comment))
                traj_path = save_fs[1].file.trajectory.path([rid])
                print("Updating trajectory file at {}".format(traj_path))
                save_fs[1].file.trajectory.write(traj, [rid])
