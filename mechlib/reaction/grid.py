""" Build the grid for a transition state search
"""

import automol
import autofile


# Functions for locating maxima
def find_max_1d(typ, grid, ts_zma, scan_name, scn_save_fs,
                mod_thy_info, constraint_dct):
    """ Find the maxmimum of a one-dimensional grid of points
        calculated along some coordinate.

        :param typ: reaction class type
        :type typ: str
        :param grid: set of points that comprise the grid
        :type grid: tuple(numpy.ndarray)
        :param ts_zma: initial ts_zma from ???
        :type ts_zma:
        :param scan_name: name of coordinate in zma along which scan conducted
        :type scan_name: str
        :param scn_save_fs: SCAN/CSCAN object with save filesys prefix
        :type scn_save_fs: autofile.fs.scan or autofile.fs.cscan object
        :param mod_thy_info: ???
        :type mod_thy_info: ???
        :param constraint_dct: values of coordinates to constrain during scan
        :type constraint_dct: dict[str: float]
    """

    # Get the locs and energies along the grid
    locs_lst, enes_lst = _grid_vals(
        grid, scan_name, scn_save_fs,
        mod_thy_info, constraint_dct)

    # Get the max zma
    max_idx = automol.pot.find_max1d(enes_lst)

    # Build lst of guess zmas
    guess_zmas = []

    # Get zma at maximum
    max_locs = locs_lst[max_idx]
    max_zma = scn_save_fs[-1].file.zmatrix.read(max_locs)
    guess_zmas.append(max_zma)

    # # Add second guess zma for migrations
    if 'migration' in typ:
        max_grid_val = grid[max_idx]
        mig_zma = automol.zmat.set_values_by_name(
            ts_zma, {scan_name: max_grid_val})
        guess_zmas.append(mig_zma)

    return guess_zmas


def find_max_2d(grid1, grid2, scan_name1, scan_name2, scn_save_fs,
                mod_thy_info, constraint_dct):
    """ Find the maxmimum of a two-dimensional grid of points
        calculated along two coordinates.

        :param grid1: set of points that comprise first dimension of grid
        :type: tuple(numpy.ndarray)
        :param grid2: set of points that comprise second dimension of grid
        :type: tuple(numpy.ndarray)
        :param scan_name1: name of coordinate in zma along dirst dimension
        :type scan_name1: str
        :param scan_name2: name of coordinate in zma along second dimension
        :type scan_name2: str
        :param scn_save_fs: SCAN/CSCAN object with save filesys prefix
        :type scn_save_fs: autofile.fs.scan or autofile.fs.cscan object
        :param mod_thy_info: ???
        :type mod_thy_info: ???
        :param constraint_dct: values of coordinates to constrain during scan
        :type constraint_dct: dict[str: float]
    """

    enes_lst = []
    locs_lst_lst = []
    for grid_val_j in grid2:
        locs_list = []
        for grid_val_i in grid1:
            if constraint_dct is None:
                locs_list.append([[scan_name1, scan_name2],
                                  [grid_val_i, grid_val_j]])
            else:
                locs_list.append([constraint_dct, [scan_name1, scan_name2],
                                  [grid_val_i, grid_val_j]])
        enes = []
        locs_lst = []
        for locs in locs_list:
            if scn_save_fs[-1].exists(locs):
                scn_path = scn_save_fs[-1].path(locs)
                sp_save_fs = autofile.fs.single_point(scn_path)
                enes.append(sp_save_fs[-1].file.energy.read(mod_thy_info[1:4]))
                locs_lst.append(locs)
        locs_lst_lst.append(locs_lst)
        if enes:
            enes_lst.append(enes)
    max_enes = []
    max_locs = []
    for idx_j, enes in enumerate(enes_lst):
        max_ene = -10000.
        max_loc = ''
        for idx_i, ene in enumerate(enes):
            if ene > max_ene:
                max_ene = ene
                max_loc = locs_lst_lst[idx_j][idx_i]
        max_enes.append(max_ene)
        max_locs.append(max_loc)
    min_ene = 10000.
    locs = []
    for idx_j, ene in enumerate(max_enes):
        if ene < min_ene:
            min_ene = ene
            locs = max_locs[idx_j]
    max_locs = locs
    max_ene = min_ene
    max_zma = scn_save_fs[-1].file.zmatrix.read(max_locs)

    return max_zma


def _grid_vals(grid, scan_name, scn_save_fs,
               mod_thy_info, constraint_dct):
    """ Build a list of the filesystem locators for each optimized
        structure as well as a list of their corresponding energies.

        :param grid: set of points that comprise the grid
        :type grid: tuple(numpy.ndarray)
        :param scan_name: name of coordinate in zma along which scan conducted
        :type scan_name: str
        :param scn_save_fs: SCAN/CSCAN object with save filesys prefix
        :type scn_save_fs: autofile.fs.scan or autofile.fs.cscan object
        :param mod_thy_info: ???
        :type mod_thy_info: ???
        :param constraint_dct: values of coordinates to constrain during scan
        :type constraint_dct: dict[str: float]
    """

    # Initialize the lists
    locs_lst = []
    enes_lst = []

    # Build the lists of all the locs for the grid
    grid_locs = []
    for grid_val_i in grid:
        if constraint_dct is None:
            grid_locs.append([[scan_name], [grid_val_i]])
        else:
            grid_locs.append([constraint_dct, [scan_name], [grid_val_i]])

    # Get the energies along the grid
    for locs in grid_locs:
        if scn_save_fs[-1].exists(locs):
            scn_path = scn_save_fs[-1].path(locs)
            sp_save_fs = autofile.fs.single_point(scn_path)
            enes_lst.append(sp_save_fs[-1].file.energy.read(mod_thy_info[1:4]))
            locs_lst.append(locs)

    return locs_lst, enes_lst
