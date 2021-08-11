"""
  Build filesystem objects
"""

import autofile
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo
from mechanalyzer.inf import rxn as rinfo
from mechlib.filesys.mincnf import min_energy_conformer_locators
from mechlib.filesys.mincnf import conformer_locators
from mechlib.filesys._build import build_fs
from mechlib.filesys._build import root_locs


def pf_rngs_filesys(spc_dct_i, spc_model_dct_i,
                    run_prefix, save_prefix, saddle, name=None):
    """ Create various filesystems needed
    """

    pf_filesystems = {}
    pf_filesystems['harm'] = set_model_filesys(
        spc_dct_i, spc_model_dct_i['vib']['geolvl'][1][1],
        run_prefix, save_prefix, saddle, name=name, cnf_range='r100')
    if spc_model_dct_i['symm']['mod'] == 'sampling':
        pf_filesystems['symm'] = set_model_filesys(
            spc_dct_i, spc_model_dct_i['symm']['geolvl'][1][1],
            run_prefix, save_prefix, saddle, name=name, cnf_range='r100')
    else:
        pf_filesystems['symm'] = None
    if spc_model_dct_i['tors']['mod'] != 'rigid':
        pf_filesystems['tors'] = set_model_filesys(
            spc_dct_i, spc_model_dct_i['tors']['geolvl'][1][1],
            run_prefix, save_prefix, saddle, name=name, cnf_range='r100')
    else:
        pf_filesystems['tors'] = None
    if spc_model_dct_i['vib']['mod'] == 'vpt2':
        pf_filesystems['vpt2'] = set_model_filesys(
            spc_dct_i, spc_model_dct_i['vib']['vpt2lvl'][1][1],
            run_prefix, save_prefix, saddle, name=name, cnf_range='r100')
    else:
        pf_filesystems['vpt2'] = None

    # Add the prefixes for now
    pf_filesystems['run_prefix'] = run_prefix
    pf_filesystems['save_prefix'] = save_prefix

    return pf_filesystems


def pf_filesys(spc_dct_i, spc_model_dct_i,
               run_prefix, save_prefix, saddle, name=None, spc_locs=None):
    """ Create various filesystems needed
    """

    pf_filesystems = {}
    cnf_range='min'
    if spc_locs is not None:
        cnf_range='specified'
    pf_filesystems['harm'] = set_model_filesys(
        spc_dct_i, spc_model_dct_i['vib']['geolvl'][1][1],
        run_prefix, save_prefix, saddle, name=name,
        cnf_range=cnf_range, spc_locs=spc_locs)
    if 'mod' in spc_model_dct_i['symm']:
        pf_filesystems['symm'] = set_model_filesys(
            spc_dct_i, spc_model_dct_i['symm']['geolvl'][1][1],
            run_prefix, save_prefix, saddle, name=name,
            cnf_range=cnf_range, spc_locs=spc_locs)
    else:
        pf_filesystems['symm'] = None
    if spc_model_dct_i['tors']['mod'] != 'rigid':
        pf_filesystems['tors'] = set_model_filesys(
            spc_dct_i, spc_model_dct_i['tors']['geolvl'][1][1],
            run_prefix, save_prefix, saddle, name=name,
            cnf_range=cnf_range, spc_locs=spc_locs)
    else:
        pf_filesystems['tors'] = None
    if spc_model_dct_i['vib']['mod'] == 'vpt2':
        pf_filesystems['vpt2'] = set_model_filesys(
            spc_dct_i, spc_model_dct_i['vib']['vpt2lvl'][1][1],
            run_prefix, save_prefix, saddle, name=name,
            cnf_range=cnf_range, spc_locs=spc_locs)
    else:
        pf_filesystems['vpt2'] = None

    # Add the prefixes for now
    pf_filesystems['run_prefix'] = run_prefix
    pf_filesystems['save_prefix'] = save_prefix

    return pf_filesystems


def set_model_filesys(spc_dct_i, level,
                      run_prefix, save_prefix, saddle, name=None,
                      cnf_range='min', spc_locs=None):
    """ Gets filesystem objects for reading many calculations
    """

    # Set the spc_info
    if saddle:
        rxn_info = spc_dct_i['rxn_info']
        spc_info = rinfo.ts_info(rxn_info)
    else:
        spc_info = sinfo.from_dct(spc_dct_i)

    levelp = tinfo.modify_orb_label(level, spc_info)

    _root = root_locs(spc_dct_i, saddle=saddle, name=name)
    cnf_run_fs, cnf_save_fs = build_fs(
        run_prefix, save_prefix, 'CONFORMER',
        thy_locs=levelp[1:],
        **_root)

    if cnf_range == 'specified':
        min_rngs_locs = spc_locs
        min_rngs_path = cnf_save_fs[-1].path(min_rngs_locs)
        cnf_run_fs[-1].create(min_rngs_locs)
    elif cnf_range == 'min':
        min_rngs_locs, min_rngs_path = min_energy_conformer_locators(
            cnf_save_fs, levelp)
        cnf_run_fs[-1].create(min_rngs_locs)
    else:
        min_rngs_locs, min_rngs_path = conformer_locators(
            cnf_save_fs, levelp, cnf_range=cnf_range)
        for min_locs in min_rngs_locs:
            cnf_run_fs[-1].create(min_locs)
    # Create run fs if that directory has been deleted to run the jobs

    return [cnf_save_fs, min_rngs_path, min_rngs_locs, '', cnf_run_fs]
    # return [cnf_save_fs, cnf_save_path, min_cnf_locs, save_path, cnf_run_fs]


def set_rpath_filesys(ts_dct, level):
    """ Gets filesystem objects for reading many calculations
    """

    # Set the spc_info
    spc_info = sinfo.from_dct(ts_dct)

    # Set some path stuff
    save_path = ts_dct['rxn_fs'][3]
    run_path = ts_dct['rxn_fs'][2]

    # Set theory filesystem used throughout
    thy_save_fs = autofile.fs.theory(save_path)
    thy_run_fs = autofile.fs.theory(run_path)

    levelp = tinfo.modify_orb_label(level[1], spc_info)

    # Get the save fileystem path
    save_path = thy_save_fs[-1].path(levelp[1:4])
    run_path = thy_run_fs[-1].path(levelp[1:4])

    thy_save_fs[-1].create(levelp[1:4])
    thy_run_fs[-1].create(levelp[1:4])

    thy_save_path = thy_save_fs[-1].path(levelp[1:4])
    thy_run_path = thy_run_fs[-1].path(levelp[1:4])

    ts_save_fs = autofile.fs.transition_state(thy_save_path)
    ts_save_fs[0].create()
    ts_save_path = ts_save_fs[0].path()
    ts_run_fs = autofile.fs.transition_state(thy_run_path)
    ts_run_fs[0].create()
    ts_run_path = ts_run_fs[0].path()

    return ts_run_path, ts_save_path, thy_run_path, thy_save_path


def get_rxn_scn_coords(ts_path, coord_name, zma_locs=(0,)):
    """ Get the values along the reaction coordinate
    """

    # Build ZMA filesys
    zma_fs = autofile.fs.zmatrix(ts_path)
    zma_path = zma_fs[-1].path(zma_locs)

    # Read the values of the reaction coord
    scn_save_fs = autofile.fs.scan(zma_path)
    scn_locs = scn_save_fs[-1].existing([[coord_name]])
    scn_grids = [locs[1][0] for locs in scn_locs
                 if locs[1][0] != 1000.0]

    return scn_grids


def make_run_path(pf_filesystems, choice):
    """ Make a run path from pf filesystems
    """
    [_, _, min_locs, _, run_fs] = pf_filesystems[choice]
    run_fs[-1].create(min_locs)
    run_path = run_fs[-1].path(min_locs)

    return run_path


def get_spc_locs_lst(
        spc_dct_i, spc_model_dct_i,
        run_prefix, save_prefix, saddle,
        cnf_range='min', sort_info_lst=None, name=None):
    """ return the locations for a pf level
    """
    # Set the spc_info
    if saddle:
        rxn_info = spc_dct_i['rxn_info']
        spc_info = rinfo.ts_info(rxn_info)
    else:
        spc_info = sinfo.from_dct(spc_dct_i)

    level = spc_model_dct_i['vib']['geolvl'][1][1]
    levelp = tinfo.modify_orb_label(level, spc_info)
    mod_info_lst = []
    if sort_info_lst is not None:
        for idx, info in enumerate(sort_info_lst):
            if info is not None:
                mod_info_lst.append(tinfo.modify_orb_label(info, spc_info))
            else:
                mod_info_lst.append(None)
    _root = root_locs(spc_dct_i, saddle=saddle, name=name)
    cnf_run_fs, cnf_save_fs = build_fs(
        run_prefix, save_prefix, 'CONFORMER',
        thy_locs=levelp[1:],
        **_root)

    if cnf_range == 'min':
        min_locs, _ = min_energy_conformer_locators(
            cnf_save_fs, levelp)
        cnf_run_fs[-1].create(min_locs)
        min_locs_lst = [min_locs]
    else:
        min_locs_lst, _ = conformer_locators(
            cnf_save_fs, levelp, cnf_range=cnf_range,
            sort_info_lst=mod_info_lst, print_enes=True)
        for min_locs in min_locs_lst:
            cnf_run_fs[-1].create(min_locs)

    return min_locs_lst
