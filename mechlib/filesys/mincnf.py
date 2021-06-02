"""
  Functions to read the filesystem and pull objects from it
"""

import time
import autofile
from phydat import phycon
from mechlib.amech_io import printer as ioprinter


def min_energy_conformer_locators(cnf_save_fs, mod_thy_info):
    """ Obtain the (ring-id, tors-id) filesystem locator pair and
        path for the conformer of a species with the lowest energy
        for the specified electronic structure method that currently
        exists in the save filesystem.

        :param cnf_save_fs: CONF object with save filesys prefix
        :type cnf_save_fs: autofile.fs.conformer obj
        :param mod_thy_info: ???
        :type mod_thy_info: ???
        :rtype: (tuple(str, str), str)
    """

    locs, paths = conformer_locators(
        cnf_save_fs, mod_thy_info, cnf_range='min')
    if locs and paths:
        ret = locs[0], paths[0]
    else:
        ret = ('', ''), ''

    return ret


def conformer_locators(cnf_save_fs, mod_thy_info, cnf_range='min'):
    """ Obtain the (ring-id, tors-id) filesystem locator pair and
        path for all conformers meeting

        --the desired criteria
        for the specified electronic structure method that currently
        exists in the save filesystem--

        The conformer locators are filtered using various criteria

        All conformer locators are supported by energy in ascending order.
        Using the energies at the SP layer:
        save/SPC/THY/CONFS/SP/THY
        save/RXN/THY/TS/CONFS/SP/THY

        :param cnf_save_fs: CONF object with save filesys prefix
        :type cnf_save_fs: autofile.fs.conformer obj
        :param mod_thy_info: ???
        :type mod_thy_info: ???
        :param cnf_range:
        :rtype: (tuple(str, str), str)
    """

    fin_locs_lst, fin_paths_lst = (), ()

    cnf_locs_lst = cnf_save_fs[-1].existing()
    if cnf_locs_lst:
        cnf_locs_lst, cnf_enes_lst = _sorted_cnf_lsts(
            cnf_locs_lst, cnf_save_fs, mod_thy_info)
        if cnf_range == 'min':
            fin_locs_lst = (cnf_locs_lst[0],)
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
        fin_paths_lst += (cnf_save_fs[-1].path(locs),)

    return fin_locs_lst, fin_paths_lst


def _sorted_cnf_lsts(cnf_locs_lst, cnf_save_fs, mod_thy_info):
    """ Sort the list of conformer locators in the save filesystem
        using the energies from the specified electronic structure method.
        The conformers are sorted such that the energies are sorted
        in ascedning order.

        :param cnf_locs_lst:
        :type cnf_locs_lst: tuple(tuple(tuple(str),tuple(str)))
        :param cnf_save_fs: CONF object with save filesys prefix
        :type cnf_save_fs: autofile.fs.conformer obj
        :param mod_thy_info: ???
        :type mod_thy_info: ???
        :rtype (tuple(tuple(tuple(str),tuple(str))), tuple(float))
    """

    fnd_cnf_enes_lst = []
    fnd_cnf_locs_lst = []
    if len(cnf_locs_lst) == 1:
        fnd_cnf_enes_lst = [10]
        fnd_cnf_locs_lst = cnf_locs_lst
    else:
        for idx, locs in enumerate(cnf_locs_lst):
            cnf_path = cnf_save_fs[-1].path(locs)
            sp_fs = autofile.fs.single_point(cnf_path)

            if sp_fs[-1].file.energy.exists(mod_thy_info[1:4]):
                fnd_cnf_enes_lst.append(sp_fs[-1].file.energy.read(
                    mod_thy_info[1:4]))
                fnd_cnf_locs_lst.append(cnf_locs_lst[idx])
            else:
                ioprinter.info_message(
                    'No energy saved in single point directory for {}'
                    .format(cnf_path))
                geo_inf_obj = cnf_save_fs[-1].file.geometry_info.read(
                    mod_thy_info[1:4])
                geo_end_time = geo_inf_obj.utc_end_time
                current_time = autofile.schema.utc_time()
                if (current_time - geo_end_time).total_seconds() < 120:
                    wait_time = (
                        120 - (current_time - geo_end_time).total_seconds()
                    )
                    ioprinter.info_message(
                        'Geo was saved in the last ' +
                        '{:3.2f} seconds, waiting for {:3.2f} seconds'.format(
                            (current_time - geo_end_time).total_seconds(),
                            wait_time))
                    time.sleep(wait_time)
                    if sp_fs[-1].file.energy.exists(mod_thy_info[1:4]):
                        fnd_cnf_enes_lst.append(sp_fs[-1].file.energy.read(
                            mod_thy_info[1:4]))
                        fnd_cnf_locs_lst.append(cnf_locs_lst[idx])
                        ioprinter.info_message('the energy is now found')
                    else:
                        ioprinter.info_message('waiting helped nothing')

    # Sort the cnf locs and cnf enes
    cnf_enes_lst, cnf_locs_lst = zip(
        *sorted(zip(fnd_cnf_enes_lst, fnd_cnf_locs_lst)))

    return cnf_locs_lst, cnf_enes_lst


def _erange_locs(cnf_locs, cnf_enes, ethresh):
    """ Obtain a list of conformer locators that includes the
        minimum-energy conformer and all conformers that lie
        above it that are within some energy threshold.

        :param cnf_locs: list of conformer filesys locators
        :type cnf_locs:
        :param cnf_enes: list of conformer energies
        :type cnf_enes:
        :param ethresh: energy threshold (in kcal)
        :type ethres: float
        :rtype: tuple()
    """

    thresh = float(ethresh.split('e')[1])

    min_cnf_locs = []
    min_ene = cnf_enes[0]
    for locs, ene in zip(cnf_locs, cnf_enes):
        rel_ene = (ene - min_ene) * phycon.EH2KCAL
        if rel_ene <= thresh:
            min_cnf_locs.append(locs)

    return min_cnf_locs


def _rrange_locs(cnf_locs, nthresh):
    """ Obtain a list of conformer locators that includes the
        minimum-energy conformer and
        above it that are within some energy threshold.
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
    """ Obtain a list of up-to N conformer locators where N is the
        number requested by the user. The count (N=1) begins from the
        minimum-energy conformer.

        :param cnf_locs: (ring-id, tors-id) CONF filesys locators
        :type cnf_locs: tuple(tuple(str), tuple(str))
        :rtype: tuple(str)
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
    """ Reads all geometries and energies which exist at the
        lowest sub-layer of some specified layer in the save
        filesystem.

        The energies are obtained from the `SP`  sub-layer that
        corresponds to the specified electronic structure method.
        The geometries and energies are then sorted together such
        that the energies are in ascending order.

        These are then written to a new .xyz trajectory file in
        save filesystem layer. This occurs in both the ring-id
        and tors-id CONF sublayers.

        :param save_fs:
        :type save_fs: autofile.fs object
        :param mod_thy_info: ???
        :type mod_thy_info: ???
        :param rid: ring-id locator for CONF filesystem
        :type rid: str
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
                enes, geos = [], []
                for locs in locs_lst:
                    trid, _ = locs
                    if trid == rid:
                        cnf_path = save_fs[-1].path(locs)
                        sp_fs = autofile.fs.single_point(cnf_path)
                        enes.append(
                            sp_fs[-1].file.energy.read(mod_thy_info[1:4]))
                        geos.append(save_fs[-1].file.geometry.read(locs))
                traj = []
                traj_sort_data = sorted(
                    zip(enes, geos, locs_lst), key=lambda x: x[0])
                for ene, geo, locs in traj_sort_data:
                    comment = (
                        'energy: {0:>15.10f} \t {1} {2}'
                    ).format(ene, locs[0], locs[1])
                    traj.append((geo, comment))
                traj_path = save_fs[1].file.trajectory.path([rid])
                print("Updating trajectory file at {}".format(traj_path))
                save_fs[1].file.trajectory.write(traj, [rid])
