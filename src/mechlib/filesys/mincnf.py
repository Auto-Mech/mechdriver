"""
  Functions to read the filesystem and pull objects from it
"""

import os
import time
import numpy
import autofile
import elstruct
import thermfit
from phydat import phycon
import automol.zmat
import automol.geom
from automol.reac import with_structures
from automol.geom import hydrogen_bonded_structure
from autorun import execute_function_in_parallel
from mechanalyzer.inf import thy as tinfo
from mechlib.amech_io import printer as ioprinter


def min_energy_conformer_locators(
        cnf_save_fs, mod_thy_info, hbond_cutoffs=None, nprocs=1):
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

    if nprocs is None:
        nprocs = 1

    locs, paths = conformer_locators(
        cnf_save_fs, mod_thy_info,
        cnf_range='min', hbond_cutoffs=hbond_cutoffs, nprocs=nprocs)
    if locs and paths:
        ret = locs[0], paths[0]
    else:
        ret = ('', ''), ''

    return ret


def conformer_locators(
        cnf_save_fs, mod_thy_info,
        cnf_range='min', sort_info_lst=None, print_enes=False,
        hbond_cutoffs=None, nprocs=1):
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
        :param cnf_range: the range of conformers to grab
        :type cnf_range: str
        :param sort_info_lst: level info to include sp or zpe in sorting
        :type sort_info_lst: list of tuples or Nones
        :rtype: (tuple(str, str), str)
    """

    def _conformer_locators(
            cnf_save_fs, mod_thy_info, cnf_range='min',
            only_hbnds=True, only_nonhbnds=True,
            freq_info=None, sp_info=None, sort_prop_dct=None,
            print_enes=False, already_counted_locs_lst=(), hbond_cutoffs=None,
            nprocs=1):

        fin_locs_lst, fin_paths_lst = (), ()
        cnf_locs_lst = cnf_save_fs[-1].existing()
        if cnf_locs_lst:
            cnf_locs_lst, cnf_enes_lst = _sorted_cnf_lsts(
                cnf_locs_lst, cnf_save_fs, mod_thy_info,
                freq_info=freq_info, sp_info=sp_info,
                sort_prop_dct=sort_prop_dct, nprocs=nprocs)
            if only_hbnds:
                cnf_locs_lst, cnf_enes_lst = _remove_nonhbonded_structures(
                    cnf_save_fs, cnf_locs_lst, cnf_enes_lst,
                    hbond_cutoffs=hbond_cutoffs)
            elif only_nonhbnds:
                cnf_locs_lst, cnf_enes_lst = _remove_hbonded_structures(
                    cnf_save_fs, cnf_locs_lst, cnf_enes_lst,
                    hbond_cutoffs=hbond_cutoffs)

            if cnf_locs_lst:
                if cnf_range == 'min':
                    fin_locs_lst = (cnf_locs_lst[0],)
                elif cnf_range == 'all':
                    fin_locs_lst = tuple(cnf_locs_lst)
                elif 'e' in cnf_range:
                    fin_locs_lst = tuple(_erange_locs(
                        cnf_locs_lst, cnf_enes_lst, cnf_range,
                        already_counted_locs_lst))
                elif 'n' in cnf_range:
                    fin_locs_lst = tuple(_nrange_locs(
                        cnf_locs_lst, cnf_range, already_counted_locs_lst))
                elif 'r' in cnf_range:
                    fin_locs_lst = tuple(_rrange_locs(
                        cnf_locs_lst, cnf_range, already_counted_locs_lst))

            for idx, locs in enumerate(fin_locs_lst):
                fin_paths_lst += (cnf_save_fs[-1].path(locs),)

            if print_enes:
                header = '\nConformer Ordering'
                if only_hbnds:
                    header += ' for only hydrogen bonded conformers'
                if only_nonhbnds:
                    header += ' for only non-hydrogen bonded conformers'
                elif not only_hbnds:
                    header += ' for all conformers'
                print(header)
                print(f'{"rid":<16}{"cid":<16}{"energy[kcal/mol]":<16}')
                print(f'{"-------":<16}{"-------":<16}{"-------":<16}')
                for idx, locs in enumerate(cnf_locs_lst):
                    _ene = (cnf_enes_lst[idx] - cnf_enes_lst[0]) * phycon.EH2KCAL
                    if locs in fin_locs_lst:
                        mark = '*'
                    else:
                        mark = ''
                    print(f'{locs[0]:<16}{locs[1]:<16}{_ene:<6.2f}{mark:<3}')
                print()
        else:
            print(f'No conformers located in {cnf_save_fs[0].path()}')

        return fin_locs_lst, fin_paths_lst

    if nprocs is None:
        nprocs = 1

    cnf_range_nohb, cnf_range_hb, cnf_range_any = _process_cnf_range(
        cnf_range)
    freq_info, sp_info, sort_prop_dct = _process_sort_info(sort_info_lst)
    union_locs_lst, union_paths_lst = (), ()
    if cnf_range_hb is not None:
        tmp_locs_lst, tmp_paths_lst = _conformer_locators(
            cnf_save_fs, mod_thy_info, cnf_range=cnf_range_hb,
            only_hbnds=True, only_nonhbnds=False,
            freq_info=freq_info, sp_info=sp_info,
            sort_prop_dct=sort_prop_dct,
            print_enes=print_enes,
            already_counted_locs_lst=union_locs_lst,
            hbond_cutoffs=hbond_cutoffs,
            nprocs=nprocs)
        union_locs_lst += tmp_locs_lst
        union_paths_lst += tmp_paths_lst
    if cnf_range_nohb is not None:
        tmp_locs_lst, tmp_paths_lst = _conformer_locators(
            cnf_save_fs, mod_thy_info, cnf_range=cnf_range_nohb,
            only_hbnds=False, only_nonhbnds=True,
            freq_info=freq_info, sp_info=sp_info,
            sort_prop_dct=sort_prop_dct,
            print_enes=print_enes,
            already_counted_locs_lst=union_locs_lst,
            hbond_cutoffs=hbond_cutoffs,
            nprocs=nprocs)
        union_locs_lst += tmp_locs_lst
        union_paths_lst += tmp_paths_lst

    if cnf_range_any is not None:
        tmp_locs_lst, tmp_paths_lst = _conformer_locators(
            cnf_save_fs, mod_thy_info, cnf_range=cnf_range_any,
            only_hbnds=False, only_nonhbnds=False,
            freq_info=freq_info, sp_info=sp_info,
            sort_prop_dct=sort_prop_dct,
            print_enes=print_enes,
            already_counted_locs_lst=union_locs_lst,
            hbond_cutoffs=hbond_cutoffs,
            nprocs=nprocs)
        union_locs_lst += tmp_locs_lst
        union_paths_lst += tmp_paths_lst

    return tuple(union_locs_lst), tuple(union_paths_lst)


def _sorted_cnf_lsts(
        cnf_locs_lst, cnf_save_fs, mod_thy_info,
        freq_info=None, sp_info=None, sort_prop_dct=None, nprocs=1):
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

    def _parallel_get_sort_energy_parameters(
            cnf_save_fs, mod_thy_info, freq_info,
            sp_info, sort_prop_dct, cnf_locs_lst,
            output_queue=None):
        locs_enes_dct = {}
        first_enes = None
        # print('locs on proc', os.getpid(),  cnf_locs_lst)
        for locs in cnf_locs_lst:
            sort_ene, first_enes = _sort_energy_parameter(
                locs, cnf_save_fs, mod_thy_info, freq_info,
                sp_info, sort_prop_dct, first_enes=first_enes)
            if first_enes is not None:
                locs_enes_dct[tuple(locs)] = (sort_ene, sum(first_enes))
            else:
                locs_enes_dct[tuple(locs)] = (sort_ene, first_enes)
        output_queue.put((locs_enes_dct,))

    fnd_cnf_enes_lst = []
    fnd_cnf_locs_lst = []
    if len(cnf_locs_lst) == 1:
        fnd_cnf_enes_lst = [10]
        fnd_cnf_locs_lst = cnf_locs_lst
    else:
        args = (
                cnf_save_fs, mod_thy_info, freq_info,
                sp_info, sort_prop_dct
                )
        locs_enes_dct_lst = execute_function_in_parallel(
            _parallel_get_sort_energy_parameters, cnf_locs_lst,
            args, nprocs=nprocs)
        first_ene = None
        for locs_enes_dct in locs_enes_dct_lst:
            for locs in locs_enes_dct:
                _, locs_first_ene = locs_enes_dct[locs]
                if locs_first_ene is not None:
                    first_ene = locs_first_ene
                break

        for locs in cnf_locs_lst:
            for locs_enes_dct in locs_enes_dct_lst:
                if tuple(locs) in locs_enes_dct:
                    sort_ene, tmp_first_ene = locs_enes_dct[tuple(locs)]
                    if sort_ene is not None:
                        if first_ene is not None:
                            sort_ene = sort_ene + (
                                (tmp_first_ene - first_ene) / phycon.EH2KCAL)
                        fnd_cnf_enes_lst.append(sort_ene)
                        fnd_cnf_locs_lst.append(locs)
            # commenting out from merge conflict
            # elif cnf_save_fs[-1].file.geometry_info.exists(locs):
            #     ioprinter.info_message(
            #         'No energy saved in single point directory for '
            #         f'{cnf_path}')
            #     geo_inf_obj = cnf_save_fs[-1].file.geometry_info.read(
            #         locs)
            #     geo_end_time = geo_inf_obj.utc_end_time
            #     current_time = autofile.schema.utc_time()
            #     save_time = (current_time - geo_end_time).total_seconds()
            #     if save_time < 120:
            #         wait_time = 120 - save_time
            #         ioprinter.info_message(
            #             f'Geo saved in the last {save_time:3.2f} seconds,'
            #             f' waiting for {wait_time:3.2f} seconds')
            #         time.sleep(wait_time)
            #         if sp_fs[-1].file.energy.exists(mod_thy_info[1:4]):
            #             fnd_cnf_enes_lst.append(sp_fs[-1].file.energy.read(
            #                 mod_thy_info[1:4]))
            #             fnd_cnf_locs_lst.append(cnf_locs_lst[idx])
            #             ioprinter.info_message('the energy is now found')
            #         else:
            #             ioprinter.info_message('waiting helped nothing')
            # print('found', fnd_cnf_enes_lst, fnd_cnf_locs_lst)

    # Sort the cnf locs and cnf enes
    if fnd_cnf_locs_lst:
        cnf_enes_lst, cnf_locs_lst = zip(
            *sorted(zip(fnd_cnf_enes_lst, fnd_cnf_locs_lst)))
    else:
        cnf_locs_lst, cnf_enes_lst = [], []

    return cnf_locs_lst, cnf_enes_lst


def _wait_for_energy_to_be_saved(cnf_save_fs, locs, sp_fs, sp_info):
    """ in case a geo was just written and its about to write and ene
    """
    ene = None
    # cnf_path = cnf_save_fs[-1].path(locs)
    if cnf_save_fs[-1].file.geometry_info.exists(locs):
        # ioprinter.info_message(
        #     f'No energy saved in single point directory for {cnf_path}')
        geo_inf_obj = cnf_save_fs[-1].file.geometry_info.read(
            locs)
        geo_path = cnf_save_fs[-1].file.geometry_info.path(locs)
        geo_end_time = geo_inf_obj.utc_end_time
        current_time = autofile.schema.utc_time()
        _time = (current_time - geo_end_time).total_seconds()
        if _time < 120:
            wait_time = 120 - _time
            ioprinter.info_message(
                f'Geo was saved in the last {_time:3.2f} seconds, '
                f'waiting for {wait_time:3.2f} seconds')
            time.sleep(wait_time)
            if sp_fs[-1].file.energy.exists(sp_info[1:4]):
                ene = sp_fs[-1].file.energy.read(sp_info[1:4])
    return ene


def _erange_locs(cnf_locs, cnf_enes, ethresh, ignore_locs_lst=()):
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
            if locs not in ignore_locs_lst:
                min_cnf_locs.append(locs)

    return min_cnf_locs


def _rrange_locs(cnf_locs, nthresh, ignore_locs_lst=()):
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
                if locs not in ignore_locs_lst:
                    min_cnf_locs.append(locs)
                    used_rids.append(list(locs)[0])

    return min_cnf_locs


def _nrange_locs(cnf_locs, nthresh, ignore_locs_lst=()):
    """ Obtain a list of up-to N conformer locators where N is the
        number requested by the user. The count (N=1) begins from the
        minimum-energy conformer.

        :param cnf_locs: (ring-id, tors-id) CONF filesys locators
        :type cnf_locs: tuple(tuple(str), tuple(str))
        :rtype: tuple(str)
    """

    thresh = int(nthresh.split('n')[1])

    min_cnf_locs = []
    idx = 0
    for locs in cnf_locs:
        if locs not in ignore_locs_lst:
            if idx+1 <= thresh:
                min_cnf_locs.append(locs)
            idx += 1

    return min_cnf_locs


def sort_info_lst(sort_str, thy_dct):
    """ Return the levels to sort conformers by if zpve or sp
        levels were assigned in input

        if we ask for zpe(lvl_wbs),sp(lvl_b2t),gibbs(700)
        out sort_info_lst will be [('gaussian', 'wb97xd', '6-31*', 'RU'),
        ('gaussian', 'b2plypd3', 'cc-pvtz', 'RU'), None, None, 700.]
    """
    sort_lvls = [None, None, None, None, None]
    sort_typ_lst = ['freqs', 'sp', 'enthalpy', 'entropy', 'gibbs']
    if sort_str is not None:
        for sort_param in sort_str.split(','):
            idx = None
            for typ_idx, typ_str in enumerate(sort_typ_lst):
                if typ_str in sort_param:
                    lvl_key = sort_str.split(typ_str + '(')[1].split(')')[0]
                    idx = typ_idx
            if idx is not None:
                if idx < 2:
                    method_dct = thy_dct.get(lvl_key)
                    if method_dct is None:
                        ioprinter.warning_message(
                            f'no {lvl_key} in theory.dat, '
                            f'not using {sort_typ_lst[idx]} in sorting')
                        continue
                    thy_info = tinfo.from_dct(method_dct)
                    sort_lvls[idx] = thy_info
                else:
                    sort_lvls[idx] = float(lvl_key)
    return sort_lvls


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
        # Update the trajectory file in the CONFS/rid level for rings
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
            comment = f'energy: {ene:<15.10f} \t {locs[0]}'
            traj.append((geo, comment))
        traj_path = save_fs[0].file.trajectory.path()

        print(f"Updating ring-torsion trajectory file at {traj_path}")
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
                        f'energy: {ene:>15.10f} \t {locs[0]} {locs[1]}')
                    traj.append((geo, comment))
                traj_path = save_fs[1].file.trajectory.path([rid])
                print(f"Updating torsion trajectory file at {traj_path}")
                save_fs[1].file.trajectory.write(traj, [rid])


def _remove_hbonded_structures(
        cnf_save_fs, cnf_locs_lst, cnf_enes_lst, hbond_cutoffs=None):
    """ Remove conformer locations for conformers with hydrogen bonds
    """
    fin_locs_lst = ()
    fin_enes_lst = ()
    for locs, enes in zip(cnf_locs_lst, cnf_enes_lst):
        if cnf_save_fs[-1].file.geometry.exists(locs):
            geo = cnf_save_fs[-1].file.geometry.read(locs)

            # Try and get a reaction object for transition state
            zma_fs = autofile.fs.zmatrix(cnf_save_fs[-1].path(locs))
            if zma_fs[-1].file.reaction.exists((0,)):
                zrxn = zma_fs[-1].file.reaction.read((0,))
                grxn = with_structures(zrxn, "geom")
                tsg = automol.reac.ts_graph(grxn)
            else:
                tsg = None

            if hbond_cutoffs is not None:
                hydrogen_bonded_structure_ = hydrogen_bonded_structure(
                    geo, *hbond_cutoffs, tsg=tsg)
            else:
                hydrogen_bonded_structure_ = hydrogen_bonded_structure(
                    geo, tsg=tsg)
            if not hydrogen_bonded_structure_:
                fin_locs_lst += (locs,)
                fin_enes_lst += (enes,)
            else:
                print(f'Removing {locs} from list because its hbonded. '
                      f'Cutoffs are: {hbond_cutoffs}')
    return fin_locs_lst, fin_enes_lst


def _remove_nonhbonded_structures(
        cnf_save_fs, cnf_locs_lst, cnf_enes_lst, hbond_cutoffs=None):
    """ Remove conformer locations for conformers without hydrogen bonds
    """
    fin_locs_lst = ()
    fin_enes_lst = ()
    for locs, enes in zip(cnf_locs_lst, cnf_enes_lst):
        if cnf_save_fs[-1].file.geometry.exists(locs):
            geo = cnf_save_fs[-1].file.geometry.read(locs)

            # Try and get a reaction object for transition state
            zma_fs = autofile.fs.zmatrix(cnf_save_fs[-1].path(locs))
            if zma_fs[-1].file.reaction.exists((0,)):
                zrxn = zma_fs[-1].file.reaction.read((0,))
                grxn = with_structures(zrxn, "geom")
                tsg = automol.reac.ts_graph(grxn)
            else:
                tsg = None

            if hbond_cutoffs is not None:
                hydrogen_bonded_structure_ = hydrogen_bonded_structure(
                    geo, *hbond_cutoffs, tsg=tsg)
            else:
                hydrogen_bonded_structure_ = hydrogen_bonded_structure(
                    geo, tsg=tsg)
            if hydrogen_bonded_structure_:
                fin_locs_lst += (locs,)
                fin_enes_lst += (enes,)
            else:
                print(
                    'Removing ', locs, ' from list because its not hbonded.',
                    'Cutoffs are', hbond_cutoffs)
    return fin_locs_lst, fin_enes_lst


def _process_cnf_range(cnf_range):
    """ Split specific info out of cnf_range
    """
    cnf_range_nohb, cnf_range_hb, cnf_range_any = None, None, None
    if '_noHB' in cnf_range:
        cnf_range_nohb = cnf_range.replace('_noHB', '')
    elif '_HB' in cnf_range:
        cnf_range_hb = cnf_range.replace('_HB', '')
    elif 'union' in cnf_range:
        cnf_range = cnf_range.replace('union', '')
        cnf_range = cnf_range.replace('(', '').replace(')', '')
        cnf_range = cnf_range.replace(' ', '')
        cnf_range = cnf_range.split(',')
        if len(cnf_range) == 3:
            cnf_range_nohb, cnf_range_hb, cnf_range_any = cnf_range
        if len(cnf_range) == 2:
            cnf_range_nohb, cnf_range_hb = cnf_range
    else:
        cnf_range_any = cnf_range
    return cnf_range_nohb, cnf_range_hb, cnf_range_any


def _process_sort_info(sort_info_lst):
    """ Split out the zpe and sp sort info
    """
    freq_info = None
    sp_info = None
    sort_prop = {}
    if sort_info_lst is not None:
        freq_info = sort_info_lst[0]
        sp_info = sort_info_lst[1]
        if sort_info_lst[2] is not None:
            sort_prop['enthalpy'] = sort_info_lst[2]
        elif sort_info_lst[3] is not None:
            sort_prop['entropy'] = sort_info_lst[3]
        elif sort_info_lst[4] is not None:
            sort_prop['gibbs'] = sort_info_lst[4]
        else:
            sort_prop['electronic'] = True
    return freq_info, sp_info, sort_prop


def zpe_from_harmonic_frequencies(
        cnf_fs, locs, mod_thy_info, freq_info):
    """ gets zpe from the harmonic frequencies
        that are saved in the filesystem
    """

    freqs = None
    if freq_info is not None:
        if mod_thy_info != freq_info:
            print(f'geometry level {mod_thy_info} '
                  'does not match requested zpe level zpe_info')
            print('Will read zpe from geometry level instead')
        if cnf_fs[-1].file.harmonic_frequencies.exists(locs):
            freqs = cnf_fs[-1].file.harmonic_frequencies.read(locs)
            zpe = 0.5 * sum(freqs) * phycon.WAVEN2EH
        else:
            print(f'No harmonic frequencies at: {cnf_fs[-1].path(locs)}')
            zpe = None
    return zpe


#def this_conformer_was_run_in_run(zma, cnf_fs):
def this_conformer_was_run_in_run(zma, cnf_fs, save_cnf_fs, thy_info):
    """ Assess if this conformer was run in RUN.
    """
    match_locs = None
    job = elstruct.Job.OPTIMIZATION
    sym_locs = []
    run_locs_lst = cnf_fs[-1].existing(ignore_bad_formats=True)
    save_locs_lst = save_cnf_fs[-1].existing(ignore_bad_formats=True)
    # This is to check if it was not found because the geometry 
    # changed so much during the optimization
    for locs in run_locs_lst:
        cnf_path = cnf_fs[-1].path(locs)
        run_fs = autofile.fs.run(cnf_path)
        run_path = run_fs[-1].path([job])
        if run_fs[-1].file.info.exists([job]):
            inf_obj = run_fs[-1].file.info.read([job])
            status = inf_obj.status
            if status == autofile.schema.RunStatus.SUCCESS:
                inp_str = run_fs[-1].file.input.read([job])
                inp_str = inp_str.replace('=', '')
                prog = inf_obj.prog
                method = inf_obj.method
                out_str = run_fs[-1].file.output.read([job])
                ran_ene = elstruct.reader.energy(prog, method, out_str)
                ran_geo = elstruct.reader.opt_geometry(prog, out_str)
                # try:
                inp_zma = elstruct.reader.inp_zmatrix(prog, inp_str)
                if inp_zma is not None:
                    if automol.zmat.almost_equal(inp_zma, zma,
                                                 dist_rtol=0.18, ang_atol=.2):
                        ioprinter.info_message(
                            'This conformer was already run ' +
                            f'in {run_path}.')
                        match_locs = locs
                    # Except:
                    #     ioprinter.info_message(
                    #        f'Program {prog} lacks inp ZMA reader for check')
                    if match_locs is not None:
                        break
    # This is to find if it was not saved becaue its equivalent
    # to other conformers
    if match_locs is not None:
        out_enes = []
        out_geos = []
        #for idx, locs in enumerate(run_locs_lst):
        for locs in save_locs_lst:
            #cnf_path = cnf_fs[-1].path(locs)
            #run_fs = autofile.fs.run(cnf_path)
            #run_path = run_fs[-1].path([job])
            #if run_fs[-1].file.info.exists([job]):
            #    inf_obj = run_fs[-1].file.info.read([job])
            #    status = inf_obj.status
            #    if status == autofile.schema.RunStatus.SUCCESS:
            #        method = inf_obj.method
            #        prog = inf_obj.prog
            #        out_str = run_fs[-1].file.output.read([job])
            #        idx_ene = elstruct.reader.energy(prog, method, out_str)
            #        idx_geo = elstruct.reader.opt_geometry(prog, out_str)
            #        if idx == locs_idx:
            #            # out_enes.append(10000)
            #            # out_geos.append(None)
            #            ran_ene = idx_ene
            #            ran_geo = idx_geo
            #        # else:
            #        out_enes.append(idx_ene)
            #        out_geos.append(idx_geo)
            cnf_path = save_cnf_fs[-1].path(locs)
            sp_fs = autofile.fs.single_point(cnf_path)
            if save_cnf_fs[-1].file.geometry.exists(locs):
                idx_ene = sp_fs[-1].file.energy.read(thy_info[1:4])
                idx_geo = save_cnf_fs[-1].file.geometry.read(locs)
                if not ran_geo == idx_geo:
                    out_enes.append(idx_ene)
                    out_geos.append(idx_geo)
            else:
                out_enes.append(10000)
                out_geos.append(None)
        for idx, _ in enumerate(out_enes):
            sym_idx = _sym_unique(
                ran_geo, ran_ene,
                [out_geos[idx]], [out_enes[idx]], ethresh=1.0e-5)
            if sym_idx is not None:
                sym_locs.append(save_locs_lst[idx])
                 
    return match_locs is not None, sym_locs


def _sym_unique(geo, ene, saved_geos, saved_enes, ethresh=1.0e-5):
    """ Check if a conformer is symmetrically distinct from the
        existing conformers in the filesystem
    """

    sym_idx = None
    new_saved_geos = []
    idx_dct = {}
    for i, (sene, sgeo) in enumerate(zip(saved_enes, saved_geos)):
        if abs(ene - sene) < ethresh:
            idx_dct[len(new_saved_geos)] = i
            new_saved_geos.append(sgeo)
    if new_saved_geos:
        _, sym_idx = automol.geom.is_unique(
            geo, new_saved_geos, check_dct={'coulomb': 1e-2})

    if sym_idx is not None:
        print(' - Structure is not symmetrically unique.')
        sym_idx = idx_dct[sym_idx]

    return sym_idx


def collect_rrho_params(
        cnf_save_fs, locs, sp_info, freq_info, mod_thy_info):
    """ get geo, freqs, and elec. ene from filesystem
    """
    geo = None
    freqs = None
    ene = None
    if cnf_save_fs[-1].file.geometry.exists(locs):
        geo = cnf_save_fs[-1].file.geometry.read(locs)
        if freq_info is None or freq_info == mod_thy_info:
            freq_fs = cnf_save_fs
            freq_locs = locs
        else:
            freq_fs, freq_locs = get_freq_location(cnf_save_fs, geo, freq_info[1:4], locs)
        if freq_locs is not None:
            if freq_fs[-1].file.harmonic_frequencies.exists(freq_locs):
                freqs = freq_fs[-1].file.harmonic_frequencies.read(freq_locs)

            cnf_path = freq_fs[-1].path(freq_locs)
            sp_fs = autofile.fs.single_point(cnf_path)
            if sp_info is not None:
                sp_thy_info = sp_info[1:4]
            else:
                sp_thy_info = mod_thy_info[1:4]
            if sp_fs[-1].file.energy.exists(sp_thy_info):
                ene = sp_fs[-1].file.energy.read(sp_thy_info)
            else:
                ene = _wait_for_energy_to_be_saved(
                    cnf_save_fs, locs, sp_fs, sp_info)

    if freqs is not None:
        freqs = [freq for freq in freqs if freq > 0.]

    return geo, freqs, ene


def get_freq_location(cnf_fs, geo, freq_thy_locs, cnf_locs):
    """ find the frequencies for a conformer at a different level of theory
    """
    path_prefix = autofile.fs.path_prefix(
        cnf_fs[-1].path(cnf_locs), ['THEORY', 'CONFORMER'])
    freq_cnf_fs = autofile.fs.manager(
        path_prefix, [['THEORY', freq_thy_locs]], 'CONFORMER')
    freq_locs = []
    for freq_cnf_locs in freq_cnf_fs[-1].existing():
        if freq_cnf_fs[-1].file.hessian.exists(freq_cnf_locs):
            freq_locs.append(freq_cnf_locs)
    match_dct = fs_confs_dict(
        freq_cnf_fs, freq_locs, cnf_fs, [cnf_locs])
    if match_dct[tuple(cnf_locs)] is None:
        match_freqs_locs = None
        # Check TS filesystem
        ts_path_prefix = autofile.fs.path_prefix(
            cnf_fs[-1].path(cnf_locs), [
                'THEORY', 'TRANSITION STATE', 'CONFORMER'])
        ts_cnf_fs = autofile.fs.manager(
            ts_path_prefix, [['THEORY', freq_thy_locs]], 'TRANSITION STATE')
        for ts_locs in ts_cnf_fs[-1].existing():
            ts_path = ts_cnf_fs[-1].path(ts_locs)
            ts_freq_cnf_fs = autofile.fs.conformer(ts_path)
            freq_locs = []
            for freq_cnf_locs in ts_freq_cnf_fs[-1].existing():
                if ts_freq_cnf_fs[-1].file.hessian.exists(freq_cnf_locs):
                    freq_locs.append(freq_cnf_locs)
            match_dct = fs_confs_dict(
                ts_freq_cnf_fs, freq_locs, cnf_fs, [cnf_locs], saddle=True)
            if match_dct[tuple(cnf_locs)] is not None:
                match_freqs_locs = tuple(match_dct[tuple(cnf_locs)])
                freq_cnf_fs = ts_freq_cnf_fs
                break
        if match_freqs_locs is None:
            print(
                'No freqs found that match ', cnf_fs[-1].path(cnf_locs),
                'in the freq location', freq_cnf_fs[0].path())
    else:
        match_freqs_locs = tuple(match_dct[tuple(cnf_locs)])
    return freq_cnf_fs, match_freqs_locs


def _check_prop_requirements(sort_prop_dct, geo, freqs, sp_ene, locs):
    """ Check that there is the correct data available for the desired
        conformer sorting method.
    """
    sort_prop = 'electronic'
    if 'enthalpy' in sort_prop_dct:
        if sort_prop_dct['enthalpy'] == 0:
            sort_prop = 'ground'
    elif 'entropy' in sort_prop_dct:
        sort_prop = 'entropy'
    elif 'gibbs' in sort_prop_dct:
        sort_prop = 'gibbs'
    if sort_prop in ['enthalpy', 'entropy', 'gibbs']:
        if geo is None:
            # ioprinter.warning_message('No geo found for ', locs)
            sort_prop = None
        if freqs is None:
            # ioprinter.warning_message('No freqs found for ', locs)
            sort_prop = None
    if sp_ene is None:
        # ioprinter.warning_message('No energy found for ', locs)
        sort_prop = None
    return sort_prop


def _sort_energy_parameter(
            locs, cnf_save_fs, mod_thy_info, freq_info,
            sp_info, sort_prop_dct, first_enes=None):
    """ find the correct energy (gibbs, entropy, enthalpy)
        at the zpe and sp or inp lvls of theory
    """
    sort_ene = None
    geo, freqs, sp_ene = collect_rrho_params(
        cnf_save_fs, locs, sp_info, freq_info, mod_thy_info)
    sort_prop = _check_prop_requirements(
        sort_prop_dct, geo, freqs, sp_ene, locs)
    if sort_prop in ['electronic', 'ground']:
        # ioprinter.debug_message('sorting by electronic energy')
        sort_ene = sp_ene
    if sort_prop == 'ground':
        zpe = 0.5 * sum(freqs) * phycon.WAVEN2EH
        sort_ene += zpe
    elif sort_prop == 'enthalpy':
        ioprinter.info_message('Sorting by Enthalpy(T) not implemented yet')
        sort_ene = None
    elif sort_prop == 'entropy':
        sort_ene = thermfit.rrho_entropy(
            geo, freqs, sort_prop_dct[sort_prop])
        sort_ene = sort_ene / phycon.EH2KCAL
    elif sort_prop == 'gibbs':
        # ioprinter.debug_message('sorting by Gibbs')
        zpe = 0.5 * sum(freqs) * phycon.WAVEN2EH
        zpe = (zpe) * phycon.EH2KCAL
        spe = sp_ene * phycon.EH2KCAL
        if first_enes is None:
            first_enes = (spe, zpe)
            rel_zero_ene = 0.
            rel_sp_ene = 0.
        else:
            rel_sp_ene = spe - first_enes[0]
            rel_zero_ene = zpe - first_enes[1]
        sort_ene = thermfit.pf.rrho_gibbs_factor(
            geo, freqs, rel_zero_ene, sort_prop_dct[sort_prop])
        sort_ene += rel_sp_ene
        sort_ene = sort_ene / phycon.EH2KCAL
    return sort_ene, first_enes


def fs_confs_dict(cnf_save_fs, cnf_save_locs_lst,
                  ini_cnf_save_fs, ini_cnf_save_locs_lst, saddle=False):
    """ Assess which structures from the cnf_save_fs currently exist
        within the ini_cnf_save_fs. Generate a dictionary to connect
        the two
    """
    match_dct = {}
    dist_tol = 0.1
    if saddle:
        dist_tol = 0.25
    for ini_locs in ini_cnf_save_locs_lst:

        match_dct[tuple(ini_locs)] = None
        # Loop over structs in cnf_save, see if they match the current struct
        # inigeo = ini_cnf_save_fs[-1].file.geometry.read(ini_locs)
        # inizma = automol.geom.zmatrix(inigeo)
        # inizma =  ini_cnf_save_fs[-1].file.zmatrix.read(ini_locs)
        ini_cnf_save_path = ini_cnf_save_fs[-1].path(ini_locs)
        # ioprinter.checking('structures', ini_cnf_save_path)
        ini_zma_save_fs = autofile.fs.zmatrix(ini_cnf_save_path)
        inizmas = [ini_zma_save_fs[-1].file.zmatrix.read((0,))]
        ini_sym_fs = autofile.fs.symmetry(ini_cnf_save_path)
        dtt = automol.zmat.conversion_info(inizmas[0])
        for sym_locs in ini_sym_fs[-1].existing():
            geo = ini_sym_fs[-1].file.geometry.read(sym_locs)
            geo_wdummy = automol.geom.apply_zmatrix_conversion(geo, dtt)
            try:
                inizmas.append(automol.zmat.from_geometry(inizmas[0], geo_wdummy))
            except:
                print('some structures have a different zmatrix')
        for inizma in inizmas:
            for locs in cnf_save_locs_lst:
                # geo = cnf_save_fs[-1].file.geometry.read(locs)
                # zma = automol.geom.zmatrix(geo)
                zma_save_fs = autofile.fs.zmatrix(cnf_save_fs[-1].path(locs))
                zma = zma_save_fs[-1].file.zmatrix.read((0,))
                if automol.zmat.almost_equal(
                        inizma, zma,
                        dist_rtol=dist_tol, ang_atol=.4):
                    # cnf_save_path = cnf_save_fs[-1].path(locs)
                    # ioprinter.info_message(
                    #     f'- Similar structure found at {cnf_save_path}')
                    match_dct[tuple(ini_locs)] = tuple(locs)
                    break
                else:
                    sym_fs = autofile.fs.symmetry(cnf_save_fs[-1].path(locs))
                    dtt = automol.zmat.conversion_info(zma)
                    for sym_locs in sym_fs[-1].existing():
                        geo = sym_fs[-1].file.geometry.read(sym_locs)
                        geo_wdummy = automol.geom.apply_zmatrix_conversion(geo, dtt)
                        try:
                            sym_zma = automol.zmat.from_geometry(zma, geo_wdummy)
                            if automol.zmat.almost_equal(
                                    inizma, sym_zma, dist_rtol=dist_tol, ang_atol=.4):
                                match_dct[tuple(ini_locs)] = tuple(locs)
                                break
                        except:
                            print('some symmetry structures have a different zmatrix')
    return match_dct
