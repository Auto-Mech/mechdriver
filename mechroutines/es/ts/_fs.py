""" rpath task function
"""

from mechlib import filesys
from mechlib.filesys import build_fs
from mechlib.filesys import root_locs


def rpath_fs(ts_dct, tsname,
             mod_ini_thy_info,
             es_keyword_dct,
             run_prefix, save_prefix):
    """ reaction path filesystem
    """

    # Set up coordinate name
    rxn_coord = es_keyword_dct.get('rxncoord')

    # Get the zma and ts locs
    zma_locs = (ts_dct['zma_idx'],)
    ts_locs = (int(tsname.split('_')[-1]),)

    # Build filesys object down to TS FS
    ts_fs = build_fs(
        run_prefix, save_prefix, 'TRANSITION STATE',
        thy_locs=mod_ini_thy_info[1:],
        **root_locs(ts_dct, saddle=True))
    ini_ts_run_fs, ini_ts_save_fs = ts_fs

    # generate fs
    if rxn_coord == 'irc':
        # Try and locate a minimum-energy conformer
        cnf_fs = build_fs(
            run_prefix, save_prefix, 'CONFORMER',
            thy_locs=mod_ini_thy_info[1:],
            **root_locs(ts_dct, saddle=True, name=tsname))
        ini_cnf_run_fs, ini_cnf_save_fs = cnf_fs

        ini_loc_info = filesys.mincnf.min_energy_conformer_locators(
            ini_cnf_save_fs, mod_ini_thy_info)
        ini_min_locs, ini_pfx_save_path = ini_loc_info

        if any(ini_min_locs):
            # Run IRC from saddle point minimum-energy conformer
            ini_pfx_run_path = ini_cnf_run_fs[-1].path(ini_min_locs)
            ini_pfx_save_path = ini_cnf_save_fs[-1].path(ini_min_locs)
            scn_alg = 'irc-sadpt'
        else:
            # Run IRC from series of points {Rmax, Rmax-1, ...}
            ini_pfx_run_path = ini_ts_run_fs[-1].path(ts_locs)
            ini_pfx_save_path = ini_ts_save_fs[-1].path(ts_locs)
            scn_alg = 'irc-rmax'
    else:
        # Run a scan along the requested reaction coordinates
        # Have an auto option that just selects the coordinate?
        # THIS IS THE DEFAULT AND ISNT FINISHED
        ini_pfx_run_path = ini_ts_run_fs[-1].path(ts_locs)
        ini_pfx_save_path = ini_ts_save_fs[-1].path(ts_locs)
        scn_alg = 'drp'

    # Set up the scan filesystem objects using the predefined prefix
    scn_fs = build_fs(
        ini_pfx_run_path, ini_pfx_save_path, 'SCAN',
        zma_locs=zma_locs)
    cscn_fs = build_fs(
        ini_pfx_run_path, ini_pfx_save_path, 'CSCAN',
        zma_locs=zma_locs)

    return scn_alg, scn_fs, cscn_fs, cnf_fs, ini_min_locs
