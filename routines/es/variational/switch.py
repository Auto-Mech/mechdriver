"""
 Handle a switch from a barrierless reaction to a fixed sadpt
"""

def assess_mep_for_sadpt(typ, grid, ts_zma, dist_name, scn_save_fs):
    """ Find the maxmimum of the grid along one dimension
    """
    # Read the energies along the scan for the filesystem
    locs_list = []
    locs_lst = []
    enes = []
    for grid_val_i in grid:
        locs_list.append([[dist_name], [grid_val_i]])
    for locs in locs_list:
        if scn_save_fs.leaf.exists(locs):
            enes.append(scn_save_fs.leaf.file.energy.read(locs))
            locs_lst.append(locs)
    max_ene = max(enes)
    max_idx = enes.index(max_ene)

    # Analyze the energies to see if a maxima has been located
    saddle_found = False

    # If a true saddle point has been found return a max zma
    if saddle_found:
        max_locs = locs_lst[max_idx]
        max_zma = scn_save_fs.leaf.file.zmatrix.read(max_locs)
    else:
        max_zma = ((), {})

    return max_zma
