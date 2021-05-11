"""
  Build filesystem objects
"""

import autofile
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo
from mechanalyzer.inf import rxn as rinfo
from mechlib.filesys import mincnf
from mechlib.filesys import build_fs
from mechlib.filesys import root_locs


def pf_rngs_filesys(spc_dct_i, spc_model_dct_i,
                    run_prefix, save_prefix, saddle, name=None):
    """ Create various filesystems needed
    """

    pf_filesystems = {}
    pf_filesystems['harm'] = set_model_filesys(
        spc_dct_i, spc_model_dct_i['vib']['geolvl'][1][1],
        run_prefix, save_prefix, saddle, name=name)
    if spc_model_dct_i['symm']['mod'] == 'sampling':
        pf_filesystems['symm'] = set_model_filesys(
            spc_dct_i, spc_model_dct_i['symm']['geolvl'][1][1],
            run_prefix, save_prefix, saddle, name=name)
    if spc_model_dct_i['tors']['mod'] != 'rigid':
        pf_filesystems['tors'] = set_model_filesys(
            spc_dct_i, spc_model_dct_i['tors']['geolvl'][1][1],
            run_prefix, save_prefix, saddle, name=name)
    if spc_model_dct_i['vib']['mod'] == 'vpt2':
        pf_filesystems['vpt2'] = set_model_filesys(
            spc_dct_i, spc_model_dct_i['vib']['vpt2lvl'][1][1],
            run_prefix, save_prefix, saddle, name=name)

    # Add the prefixes for now
    pf_filesystems['run_prefix'] = run_prefix
    pf_filesystems['save_prefix'] = save_prefix

    return pf_filesystems


def pf_filesys(spc_dct_i, spc_model_dct_i,
               run_prefix, save_prefix, saddle, name=None):
    """ Create various filesystems needed
    """

    pf_filesystems = {}
    print('harm')
    pf_filesystems['harm'] = set_model_filesys(
        spc_dct_i, spc_model_dct_i['vib']['geolvl'][1][1],
        run_prefix, save_prefix, saddle, name=name, rings='min')
    if spc_model_dct_i['symm']['mod'] == 'sampling':
        print('symm')
        pf_filesystems['symm'] = set_model_filesys(
            spc_dct_i, spc_model_dct_i['symm']['geolvl'][1][1],
            run_prefix, save_prefix, saddle, name=name, rings='min')
    if spc_model_dct_i['tors']['mod'] != 'rigid':
        print('tors')
        pf_filesystems['tors'] = set_model_filesys(
            spc_dct_i, spc_model_dct_i['tors']['geolvl'][1][1],
            run_prefix, save_prefix, saddle, name=name, rings='min')
    if spc_model_dct_i['vib']['mod'] == 'vpt2':
        print('vpt2')
        if spc_model_dct_i['vib']['mod'] == 'vpt2':
            pf_filesystems['vpt2'] = set_model_filesys(
                spc_dct_i, spc_model_dct_i['vib']['vpt2lvl'][1][1],
                run_prefix, save_prefix, saddle, name=name, rings='min')

    # Add the prefixes for now
    pf_filesystems['run_prefix'] = run_prefix
    pf_filesystems['save_prefix'] = save_prefix

    return pf_filesystems


def set_model_filesys(spc_dct_i, level,
                      run_prefix, save_prefix, saddle, name=None, rings='all'):
    """ Gets filesystem objects for reading many calculations
    """

    # Set the spc_info
    if saddle:
        rxn_info = spc_dct_i['rxn_info']
        spc_info = rinfo.ts_info(rxn_info)
    else:
        spc_info = sinfo.from_dct(spc_dct_i)

    print('level', level)
    levelp = tinfo.modify_orb_label(level, spc_info)

    _root = root_locs(spc_dct_i, saddle=saddle, name=name)
    cnf_run_fs, cnf_save_fs = build_fs(
        run_prefix, save_prefix, 'CONFORMER',
        thy_locs=levelp[1:],
        **_root)

    if rings == 'min':
        min_rngs_locs, min_rngs_path = mincnf.min_energy_conformer_locators(
            cnf_save_fs, levelp)
        cnf_run_fs[-1].create(min_rngs_locs)
    else:
        min_rngs_locs, min_rngs_path = mincnf.conformer_locators(
             cnf_save_fs, levelp, cnf_range='r100')
        for min_locs in min_rngs_locs:
            cnf_run_fs[-1].create(min_locs)

    print('model filesys', min_rngs_locs, min_rngs_path)
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
    print('level', levelp)
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


def set_etrans_fs(pf_filesystems, lj_info, thy_info):
    """ build etrans filesystems
    """

    # Get the harmonic filesys information
    [cnf_fs, harm_path, min_cnf_locs, _, run_path] = pf_filesystems['harm']

    # Build the energy transfer filesys
    lj_mod_thy_info = filesys.inf.modify_orb_restrict(
        lj_info, run_thy_info)
    etrans_save_fs, etrans_locs = filesys.build.etrans_fs_from_prefix(
        tgt_cnf_save_path, bath_info, lj_mod_thy_info)

    return etrans_save_fs, etrans_locs


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
