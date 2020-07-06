"""
  Build filesystem objects
"""

import autofile
from lib.filesys import inf as finf
from lib.filesys import mincnf


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
    """ Gets filesystem objects for torsional calculations
    """

    # Set the spc_info
    spc_info = finf.get_spc_info(spc_dct_i)

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
    levelp = finf.modify_orb_restrict(spc_info, level)

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
    min_cnf_locs = mincnf.min_energy_conformer_locators(cnf_save_fs)

    # Get the save path for the conformers
    if min_cnf_locs:
        cnf_save_path = cnf_save_fs[-1].path(min_cnf_locs)
    else:
        cnf_save_path = ''

    return [cnf_save_fs, cnf_save_path, min_cnf_locs, save_path, cnf_run_fs]


def make_run_path(pf_filesystems, choice):
    """ Make a run path from pf filesystems
    """
    [_, _, min_locs, _, run_fs] = pf_filesystems[choice]
    run_fs[-1].create(min_locs)
    run_path = run_fs[-1].path(min_locs)

    return run_path
