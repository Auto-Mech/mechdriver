""" rpath task function

- build to TS/00/ (has Z and CONFS)
- if scan
    - if irc
        - if sadpt exists
            use sadpt
        - else:
            use rmax-1 series
    - else
        - take rn scan
- else:
    - if irc
        - build min cnfs, if cnfs get scan fs
    - else
        - build z fs
"""

from mechlib import filesys
from mechlib.filesys import build_fs
from mechlib.filesys import root_locs


def rpath_fs(spc_dct_i, spc_name,
             mod_ini_thy_info, ts_info,
             es_keyword_dct,
             run_prefix, save_prefix):
    """ reaction path filesystem
    """

    # Set up coordinate name
    rxn_coord = es_keyword_dct.get('rxn_coord')

    # Build filesys object down to TS FS
    _root = root_locs(spc_dct_i, saddle=True, name=spc_name)
    ts_fs = build_fs(
        run_prefix, save_prefix, 'TRANSITION STATE',
        thy_locs=mod_ini_thy_info[1:],
        **_root)
    ini_ts_run_fs, ini_ts_save_fs = ts_fs

    # generate fs
    if rxn_coord == 'irc':
        # Try and locate a minimum-energy conformer
        cnf_fs = build_fs(
            run_prefix, save_prefix, 'CONFORMER',
            thy_locs=mod_ini_thy_info[1:],
            **_root)
        ini_cnf_run_fs, ini_cnf_save_fs = cnf_fs

        ini_loc_info = filesys.mincnf.min_energy_conformer_locators(
            ini_cnf_save_fs, mod_ini_thy_info)
        ini_min_locs, ini_pfx_save_path = ini_loc_info

        if ini_min_locs:
            # Run IRC from saddle point minimum-energy conformer
            ini_pfx_run_path = ini_cnf_run_fs[-1].path(ini_min_locs)
            ini_pfx_save_path = ini_cnf_save_fs[-1].path(ini_min_locs)
            scn_alg = 'irc-sadpt'
        else:
            # Run IRC from series of points {Rmax, Rmax-1, ...}
            ini_pfx_run_path = ini_ts_run_fs[-1].path(ts_info)
            ini_pfx_save_path = ini_ts_save_fs[-1].path(ts_info)
            scn_alg = 'irc-rmax'
    else:
        # Run a scan along the requested reaction coordinates
        ini_pfx_run_path = ini_ts_run_fs[-1].path(ts_info)
        ini_pfx_save_path = ini_ts_save_fs[-1].path(ts_info)
        scn_alg = 'drp'

    # Set up the scan filesystem objects using the predefined prefix
    scn_fs = build_fs(
        ini_pfx_run_path, ini_pfx_save_path, 'SCAN',
        zma_locs=(0,))

    return scn_alg, scn_fs, cnf_fs, ini_min_locs
