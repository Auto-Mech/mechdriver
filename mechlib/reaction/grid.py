""" Build the grid for a transition state search
"""

import automol
import autofile


# Functions for locating maxima
def find_max_1d(typ, grid, ts_zma, dist_name, scn_save_fs,
                mod_thy_info, constraint_dct):
    """ Find the maxmimum of the grid along one dimension
    """

    # Get the locs and energies along the grid
    locs_lst, enes_lst = _grid_vals(
        grid, dist_name, scn_save_fs,
        mod_thy_info, constraint_dct)
    print('lst', locs_lst, enes_lst)

    # Get the max zma
    max_idx = automol.pot.find_max1d(enes_lst)
    print('max', max_idx)
    
    # Build lst of guess zmas
    guess_zmas = []

    # Get zma at maximum
    max_locs = locs_lst[max_idx]
    try:
        max_zma = scn_save_fs[-1].file.zmatrix.read(max_locs)
    except:
        max_geo = scn_save_fs[-1].file.geometry.read(max_locs)
        max_zma = automol.geom.zmatrix(max_geo)

    guess_zmas.append(max_zma)

    # # Add second guess zma for migrations
    if 'migration' in typ:
        max_grid_val = grid[max_idx]
        mig_zma = automol.zmat.set_values(
            ts_zma, {dist_name: max_grid_val})
        guess_zmas.append(mig_zma)

    return guess_zmas


def find_max_2d(grid1, grid2, dist_name, brk_name, scn_save_fs,
                mod_thy_info, constraint_dct):
    """ Find the maxmimum of the grid along two dimensions
    """
    enes_lst = []
    locs_lst_lst = []
    for grid_val_j in grid2:
        locs_list = []
        for grid_val_i in grid1:
            if constraint_dct is None:
                locs_list.append([[dist_name, brk_name],
                                  [grid_val_i, grid_val_j]])
            else:
                locs_list.append([constraint_dct, [dist_name, brk_name],
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
            print('ene max_ene', ene, max_ene)
            if ene > max_ene:
                max_ene = ene
                max_loc = locs_lst_lst[idx_j][idx_i]
                print('new max', max_ene, max_loc)
        max_enes.append(max_ene)
        max_locs.append(max_loc)
    print('max enes', max_enes)
    min_ene = 10000.
    locs = []
    for idx_j, ene in enumerate(max_enes):
        print('ene min_ene', ene, min_ene)
        if ene < min_ene:
            min_ene = ene
            locs = max_locs[idx_j]
            print('locs', locs)
    max_locs = locs
    max_ene = min_ene
    print('min max loc', max_ene, max_locs)
    print('min max loc', scn_save_fs[-1].path(max_locs))
    max_zma = scn_save_fs[-1].file.zmatrix.read(max_locs)

    # print('geometry for maximum along scan:', max_zma)
    # print('energy for maximum along scan:', max_ene)

    return max_zma


def _grid_vals(grid, dist_name, scn_save_fs,
               mod_thy_info, constraint_dct):
    """ efef
    """

    # Initialize the lists
    locs_lst = []
    enes_lst = []

    # Build the lists of all the locs for the grid
    grid_locs = []
    for grid_val_i in grid:
        if constraint_dct is None:
            grid_locs.append([[dist_name], [grid_val_i]])
        else:
            grid_locs.append([constraint_dct, [dist_name], [grid_val_i]])

    # Get the energies along the grid
    for locs in grid_locs:
        print('locs', locs)
        print(scn_save_fs[-1].path(locs))
        if scn_save_fs[-1].exists(locs):
            print('exists')
            scn_path = scn_save_fs[-1].path(locs)
            sp_save_fs = autofile.fs.single_point(scn_path)
            enes_lst.append(sp_save_fs[-1].file.energy.read(mod_thy_info[1:4]))
            locs_lst.append(locs)

    return locs_lst, enes_lst
