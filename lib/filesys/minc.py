""" Dealing withe the minimum-energy conf and trajectories
"""


def min_energy_conformer_locators(cnf_save_fs, zpe_corrd=False):
    """ locators for minimum energy conformer """
    cnf_locs_lst = cnf_save_fs[-1].existing()
    if cnf_locs_lst:
        cnf_enes = [cnf_save_fs[-1].file.energy.read(locs)
                    for locs in cnf_locs_lst]
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


def traj_sort(save_fs):
    """ sort trajectory file according to energies
    """
    locs_lst = save_fs[-1].existing()
    if locs_lst:
        enes = [save_fs[-1].file.energy.read(locs)
                for locs in locs_lst]
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
