"""
  Build filesystem objects
"""

import autofile
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo
from mechlib.filesys import mincnf


def pf_filesys(spc_dct_i, pf_levels,
               run_prefix, save_prefix, saddle):
    """ Create various filesystems needed
    """

    pf_filesystems = {}
    pf_filesystems['harm'] = set_model_filesys(
        spc_dct_i, pf_levels['harm'][1], run_prefix, save_prefix, saddle)
    if pf_levels['sym']:
        pf_filesystems['sym'] = set_model_filesys(
            spc_dct_i, pf_levels['sym'][1], run_prefix, save_prefix, saddle)
    if pf_levels['tors']:
        pf_filesystems['tors'] = set_model_filesys(
            spc_dct_i, pf_levels['tors'][1][0],
            run_prefix, save_prefix, saddle)
    if pf_levels['vpt2']:
        pf_filesystems['vpt2'] = set_model_filesys(
            spc_dct_i, pf_levels['vpt2'][1], run_prefix, save_prefix, saddle)

    return pf_filesystems


def set_model_filesys(spc_dct_i, level, run_prefix, save_prefix, saddle):
    """ Gets filesystem objects for reading many calculations
    """

    # Set the spc_info
    spc_info = sinfo.from_dct(spc_dct_i)
    # Set some path stuff
    if saddle:
        save_path = spc_dct_i['rxn_fs'][3]
        run_path = spc_dct_i['rxn_fs'][2]
    else:
        spc_save_fs = autofile.fs.species(save_prefix)
        spc_save_fs[-1].create(spc_info)
        save_path = spc_save_fs[-1].path(spc_info)
        spc_run_fs = autofile.fs.species(run_prefix)
        spc_run_fs[-1].create(spc_info)
        run_path = spc_run_fs[-1].path(spc_info)

    # Set theory filesystem used throughout
    thy_save_fs = autofile.fs.theory(save_path)
    thy_run_fs = autofile.fs.theory(run_path)

    # Set the level for the model
    levelp = tinfo.modify_orb_label(level, spc_info)

    # Get the save fileystem path
    save_path = thy_save_fs[-1].path(levelp[1:4])
    run_path = thy_run_fs[-1].path(levelp[1:4])
    thy_run_fs[-1].create(levelp[1:4])
    if saddle:
        save_fs = autofile.fs.transition_state(save_path)
        save_fs[0].create()
        save_path = save_fs[0].path()
        run_fs = autofile.fs.transition_state(run_path)
        run_fs[0].create()
        run_path = run_fs[0].path()

    # Get the fs object and the locs
    cnf_run_fs = autofile.fs.conformer(run_path)
    cnf_save_fs = autofile.fs.conformer(save_path)
    min_cnf_locs, cnf_save_path = mincnf.min_energy_conformer_locators(
        cnf_save_fs, levelp)

    # Get the save path for the conformers
    # if min_cnf_locs:
    #     cnf_save_path = cnf_save_fs[-1].path(min_cnf_locs)
    # else:
    #     cnf_save_path = ''

    return [cnf_save_fs, cnf_save_path, min_cnf_locs, save_path, cnf_run_fs]


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
