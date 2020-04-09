""" Build the grid for a transition state search
"""

import numpy
import automol
from lib.phydat import phycon
from lib.phydat import bnd


def find_max_1d(typ, grid, ts_zma, dist_name, scn_save_fs):
    """ Find the maxmimum of the grid along one dimension
    """
    locs_list = []
    locs_lst = []
    enes = []
    for grid_val_i in grid:
        locs_list.append([[dist_name], [grid_val_i]])
    for locs in locs_list:
        if scn_save_fs[-1].exists(locs):
            enes.append(scn_save_fs[-1].file.energy.read(locs))
            locs_lst.append(locs)
    max_ene = max(enes)
    max_idx = enes.index(max_ene)
    if 'migration' in typ:
        max_grid_val = grid[max_idx]
        max_zma = automol.zmatrix.set_values(
            ts_zma, {dist_name: max_grid_val})
    else:
        max_locs = locs_lst[max_idx]
        max_zma = scn_save_fs[-1].file.zmatrix.read(max_locs)

    return max_zma, max_ene


def find_max_2d(grid1, grid2, dist_name, brk_name, scn_save_fs):
    """ Find the maxmimum of the grid along two dimensions
    """
    enes_lst = []
    locs_lst_lst = []
    for grid_val_j in grid2:
        locs_list = []
        for grid_val_i in grid1:
            locs_list.append([[dist_name, brk_name], [grid_val_i, grid_val_j]])
        enes = []
        locs_lst = []
        for locs in locs_list:
            if scn_save_fs[-1].exists(locs):
                enes.append(scn_save_fs[-1].file.energy.read(locs))
                locs_lst.append(locs)
        locs_lst_lst.append(locs_lst)
        enes_lst.append(enes)
        print('enes_lst', enes_lst)
    max_enes = []
    max_locs = []
    for idx_j, enes in enumerate(enes_lst):
        max_ene = -10000.
        max_loc = ''
        for idx_i, ene in enumerate(enes):
            if ene > max_ene:
                max_ene = ene
                max_loc = locs_lst_lst[idx_j][idx_i]
                print('new max', max_ene, max_loc)
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
    print('min max loc', max_ene, max_locs)
    print('min max loc', scn_save_fs[-1].path(max_locs))
    max_zma = scn_save_fs[-1].file.zmatrix.read(max_locs)

    return max_zma, max_ene


def build_grid(rtype, rbktype, ts_bnd_len, ts_zma, dist_name, npoints=None):
    """ Set the grid for a transition state search
    """

    # Set up the backup type
    if 'beta_scission' in rbktype:
        bkp_grid, bkp_update_guess = beta_scission_bkp_grid(
            npoints, ts_bnd_len)
    elif 'addition' in rtype:
        bkp_grid, bkp_update_guess = addition_bkp_grid(
            npoints, ts_bnd_len)
    else:
        bkp_grid = []
        bkp_update_guess = False

    # Set the main type
    # if spin == 'high':
    if 'beta scission' in rtype:
        grid, update_guess = beta_scission_grid(npoints, ts_bnd_len)
    elif 'addition' in rtype and 'rad' not in rtype:
        grid, update_guess = addition_grid(npoints, ts_bnd_len)
    elif 'hydrogen migration' in rtype and 'rad' not in rtype:
        grid, update_guess = hydrogen_migration_grid(npoints, ts_bnd_len,
                                                     ts_zma, dist_name)
    elif 'unimolecular elimination' in rtype:
        grid, update_guess = unimolecular_elimination_grid(npoints, ts_bnd_len,
                                                           ts_zma, dist_name)
    elif 'hydrogen abstraction' in rtype:
        grid, update_guess = hydrogen_abstraction(npoints, ts_bnd_len)
    elif 'substitution' in rtype:
        grid, update_guess = substitution(npoints, ts_bnd_len)
    elif 'insertion' in rtype:
        grid, update_guess = insertion(npoints, ts_bnd_len)
    # elif spin == 'low':
    elif 'radical radical' in rtype and 'addition' in rtype:
        grid, update_guess = radrad_addition_grid(npoints, ts_bnd_len)
    elif 'radical radical' in rtype and 'hydrogen abstraction' in rtype:
        grid, update_guess = radrad_hydrogen_abstraction(npoints, ts_bnd_len)
    else:
        raise NotImplementedError

    return grid, update_guess, bkp_grid, bkp_update_guess


# Tight TS grid

def beta_scission_grid(npoints, ts_bnd_len):
    """ Build forward 1D grid for a beta scission reaction
    """
    npoints = 8 if npoints is not None else npoints
    rmin = 1.4 * phycon.ANG2BOHR
    rmax = 2.0 * phycon.ANG2BOHR
    if ts_bnd_len in bnd.LEN_DCT:
        npoints = 14
        bnd_len = bnd.LEN_DCT[ts_bnd_len]
        rmin = bnd_len + 0.1 * phycon.ANG2BOHR
        rmax = bnd_len + 0.8 * phycon.ANG2BOHR
    grid = numpy.linspace(rmin, rmax, npoints)
    update_guess = False

    return grid, update_guess


def addition_grid(npoints, ts_bnd_len):
    """ Build forward 1D grid for addition reaction
    """

    npoints = 14
    rmin = 1.6 * phycon.ANG2BOHR
    rmax = 2.8 * phycon.ANG2BOHR
    if ts_bnd_len in bnd.LEN_DCT:
        bnd_len = bnd.LEN_DCT[ts_bnd_len]
        rmin = bnd_len + 0.1 * phycon.ANG2BOHR
        rmax = bnd_len + 1.2 * phycon.ANG2BOHR
    grid = _geometric_progression(
        rmin, rmax, npoints, gfact=1.1, rstp=0.05)
    update_guess = False

    return grid, update_guess


def hydrogen_migration_grid(npoints, ts_bnd_len, ts_zma, dist_name):
    """ Build forward 1D grid for addition reaction
    """
    interval = 0.3*phycon.ANG2BOHR
    # get rmax from ts_zma
    rmax = automol.zmatrix.values(ts_zma)[dist_name]
    rmin1 = 2.0*phycon.ANG2BOHR
    rmin2 = 1.3*phycon.ANG2BOHR
    if ts_bnd_len in bnd.LEN_DCT:
        bnd_len = bnd.LEN_DCT[ts_bnd_len]
        rmin2 = bnd_len + 0.05 * phycon.ANG2BOHR
    if rmax > rmin1:
        npoints = (rmax-rmin1)/interval
        if npoints < 1:
            grid1 = []
        else:
            grid1 = numpy.linspace(rmax, rmin1, npoints)
    else:
        grid1 = []
    grid2 = numpy.linspace(rmin1, rmin2, 18)
    grid = numpy.concatenate((grid1, grid2), axis=None)
    update_guess = True

    return grid, update_guess


def unimolecular_elimination_grid(npoints, ts_bnd_len, ts_zma, brk_name):
    """ Build forward 2D grid for elimination reaction
    """
    brk_coo, = automol.zmatrix.coordinates(ts_zma)[brk_name]
    brk_len_key = tuple(sorted(map(syms.__getitem__, brk_coo)))

    interval = 0.2 * phycon.ANG2BOHR
    rmin = 1.4 * phycon.ANG2BOHR
    rmax = 2.8 * phycon.ANG2BOHR
    if ts_bnd_len in bnd.LEN_DCT:
        bnd_len = bnd.LEN_DCT[ts_bnd_len]
        brk_len = bnd.LEN_DCT[brk_len_key]
        r1min = bnd_len + 0.2 * phycon.ANG2BOHR
        r1max = bnd_len + 1.4 * phycon.ANG2BOHR
        r2min = brk_len + 0.2 * phycon.ANG2BOHR
        r2max = brk_len + 0.8 * phycon.ANG2BOHR
        grid1 = numpy.linspace(r1min, r1max, 8)
        grid2 = numpy.linspace(r2min, r2max, 4)
        grid = [grid1, grid2]
        update_guess = False

    return grid, update_guess


def hydrogen_abstraction(npoints, ts_bnd_len):
    """ Build forward 1D grid for hydrogen abstraction reaction
    """
    npoints = 16
    rmin = 0.7 * phycon.ANG2BOHR
    rmax = 2.2 * phycon.ANG2BOHR
    if ts_bnd_len in bnd.LEN_DCT:
        bnd_len = bnd.LEN_DCT[ts_bnd_len]
        rmin = bnd_len
        rmax = bnd_len + 1.0 * phycon.ANG2BOHR
    grid = numpy.linspace(rmin, rmax, npoints)
    update_guess = False

    return grid, update_guess


def substitution(npoints, ts_bnd_len):
    """ Build forward 1D grid for substitution reaction
    """
    npoints = 14
    rmin = 0.7 * phycon.ANG2BOHR
    rmax = 2.4 * phycon.ANG2BOHR
    if ts_bnd_len in bnd.LEN_DCT:
        bnd_len = bnd.LEN_DCT[ts_bnd_len]
        rmin = bnd_len
        rmax = bnd_len + 1.4 * phycon.ANG2BOHR
    grid = numpy.linspace(rmin, rmax, npoints)
    update_guess = False

    return grid, update_guess


def insertion(npoints, ts_bnd_len):
    """ Build forward 1D grid for insertion reaction
    """
    npoints = 16
    rmin = 1.4 * phycon.ANG2BOHR
    rmax = 2.4 * phycon.ANG2BOHR
    if ts_bnd_len in bnd.LEN_DCT:
        bnd_len = bnd.LEN_DCT[ts_bnd_len]
        rmin = bnd_len
        rmax = bnd_len + 1.4 * phycon.ANG2BOHR
    grid = numpy.linspace(rmin, rmax, npoints)
    update_guess = False

    return grid, update_guess


# Barrierless TS grid

def radrad_addition_grid(npoints, ts_bnd_len):
    """ Build forward 1D grid for a beta scission reaction
    """

    npoints1 = 4
    npoints2 = 4
    rstart = 2.4 * phycon.ANG2BOHR
    rend1 = 1.8 * phycon.ANG2BOHR
    rend2 = 3.0 * phycon.ANG2BOHR
    grid1 = numpy.linspace(rstart, rend1, npoints1)
    grid2 = numpy.linspace(rstart, rend2, npoints2)
    grid2 = numpy.delete(grid2, 0)
    grid = [grid1, grid2]
    update_guess = True

    return grid, update_guess


def radrad_hydrogen_abstraction(npoints, ts_bnd_len):
    """ Build forward 1D grid for elimination reaction
    """
    rstart = 2.4 * phycon.ANG2BOHR
    rend1 = 1.4 * phycon.ANG2BOHR
    rend2 = 3.0 * phycon.ANG2BOHR
    npoints1 = 8
    npoints2 = 4
    grid1 = numpy.linspace(rstart, rend1, npoints1)
    grid2 = numpy.linspace(rstart, rend2, npoints2)
    grid2 = numpy.delete(grid2, 0)
    grid = [grid1, grid2]
    update_guess = True

    return grid, update_guess


# Backup Grid
def beta_scission_bkp_grid(npoints, ts_bnd_len):
    """ Build backward 1D grid for a beta scission reaction
    """
    rmin = 1.4 * phycon.ANG2BOHR
    rmax = 2.0 * phycon.ANG2BOHR
    if ts_bnd_len in bnd.LEN_DCT:
        bnd_len = bnd.LEN_DCT[ts_bnd_len]
        npoints = 14
        rmin = bnd_len + 0.1 * phycon.ANG2BOHR
        rmax = bnd_len + 0.8 * phycon.ANG2BOHR
    bkp_grid = numpy.linspace(rmin, rmax, npoints)
    bkp_update_guess = False

    return bkp_grid, bkp_update_guess


def addition_bkp_grid(npoints, ts_bnd_len):
    """ Build backward 1D grid for a beta scission reaction
    """
    rmin = 1.6 * phycon.ANG2BOHR
    rmax = 2.8 * phycon.ANG2BOHR
    if ts_bnd_len in bnd.LEN_DCT:
        npoints = 14
        bnd_len = bnd.LEN_DCT[ts_bnd_len]
        rmin = bnd_len + 0.1 * phycon.ANG2BOHR
        rmax = bnd_len + 1.4 * phycon.ANG2BOHR
    bkp_grid = numpy.linspace(rmin, rmax, npoints)
    bkp_update_guess = False

    return bkp_grid, bkp_update_guess


# Special grid progressions
def _geometric_progression(rmin, rmax, npoints, gfact=1.1, rstp=0.05):
    """ Build a grid using a geometric progresion
    """
    grid = [rmin]
    rgrid = rmin
    for _ in range(npoints):
        rgrid += rstp
        if rgrid == rmax:
            break
        grid.append(rgrid)
        rstp = rstp * gfact
    grid = numpy.array(grid)

    return grid
