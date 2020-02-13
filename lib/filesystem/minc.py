""" Dealing withe the minimum-energy conf and trajectories
"""


def min_energy_conformer_locators(cnf_save_fs):
    """ locators for minimum energy conformer """
    cnf_locs_lst = cnf_save_fs.leaf.existing()
    if cnf_locs_lst:
        cnf_enes = [cnf_save_fs.leaf.file.energy.read(locs)
                    for locs in cnf_locs_lst]
        min_cnf_locs = cnf_locs_lst[cnf_enes.index(min(cnf_enes))]
    else:
        min_cnf_locs = None
    return min_cnf_locs


def locs_sort(save_fs):
    """ sort trajectory file according to energies
    """
    locs_lst = save_fs.leaf.existing()
    if locs_lst:
        enes = [save_fs.leaf.file.energy.read(locs)
                for locs in locs_lst]
        sorted_locs = []
        for _, loc in sorted(zip(enes, locs_lst), key=lambda x: x[0]):
            sorted_locs.append(loc)
    return sorted_locs


def traj_sort(save_fs):
    """ sort trajectory file according to energies
    """
    locs_lst = save_fs.leaf.existing()
    if locs_lst:
        enes = [save_fs.leaf.file.energy.read(locs)
                for locs in locs_lst]
        geos = [save_fs.leaf.file.geometry.read(locs)
                for locs in locs_lst]
        traj = []
        traj_sort_data = sorted(zip(enes, geos, locs_lst), key=lambda x: x[0])
        for ene, geo, locs in traj_sort_data:
            comment = 'energy: {0:>15.10f} \t {1}'.format(ene, locs[0])
            traj.append((comment, geo))
        traj_path = save_fs.trunk.file.trajectory.path()
        print("Updating trajectory file at {}".format(traj_path))
        save_fs.trunk.file.trajectory.write(traj)
