"""
temp home for filesys
"""

def set_model_filesys(thy_save_fs, spc_info, level, saddle=False):
    """ Gets filesystem objects for torsional calculations
    """
    # Set the level for the model
    levelp = filesys.inf.modify_orb_restrict(spc_info, level)

    # Get the save fileystem path
    save_path = thy_save_fs[-1].path(levelp[1:4])
    if saddle:
        save_fs = autofile.fs.transition_state(save_path)
        save_fs[0].create()
        save_path = save_fs[0].path()

    # Get the fs object and the locs
    cnf_save_fs = autofile.fs.conformer(save_path)
    min_cnf_locs = filesys.mincnf.min_energy_conformer_locators(cnf_save_fs)

    # Get the save path for the conformers
    if min_cnf_locs:
        cnf_save_path = cnf_save_fs[-1].path(min_cnf_locs)
    else:
        cnf_save_path = ''

    return cnf_save_fs, cnf_save_path, min_cnf_locs, save_path





def _cnf_filesys(spc_dct, rxn, pf_levels, save_prefix,
                 saddle=False, level='harm'):
    """ Set needed conformer filesys objects
    """

    if level == 'harm':
        thy_info = pf_levels[2]
    elif level == 'vpt2':
        thy_info = pf_levels[2]

    # Set the filesystem objects
    spc_info = (spc_dct['ich'], spc_dct['chg'], spc_dct['mul'])
    mod_thy_info = fsorb.mod_orb_restrict(spc_info, thy_info)
    if not saddle:
        _, thy_save_path = fbuild.spc_thy_fs_from_root(
            save_prefix, spc_info, mod_thy_info)
    else:
        rxn_info = finf.rxn_info(
            rxn['reacs'], rxn['prods'], spc_dct)
        _, thy_save_path = fbuild.rxn_thy_fs_from_root(
            save_prefix, rxn_info, mod_thy_info)
        _, thy_save_path = fbuild.ts_fs_from_thy(thy_save_path)

    # Get the cnf filesys needed everything based off geo+freq (also vpt2)
    cnf_save_fs, cnf_save_locs = fbuild.cnf_fs_from_prefix(
        thy_save_path, cnf='min')
    cnf_save_paths = fbuild.cnf_paths_from_locs(
        cnf_save_fs, cnf_save_locs)

    return cnf_save_fs, cnf_save_paths[0], cnf_save_locs, thy_save_path


def _scn_filesys(cnf_save_path, run_tors_names):
    """ Set needed conformer filesys objects
    """

    # Get a list of the other tors coords to freeze and set the filesystem
    frz_all_tors = es_keyword_dct['frz_all_tors']
    if frz_all_tors:
        scn_save_fs = autofile.fs.cscan(cnf_save_path)
        scn_locs = fbuild.cscn_locs_from_fs(scn_save_fs, run_tors_names)
    else:
        scn_save_fs = autofile.fs.scan(cnf_save_path)
        scn_locs = fbuild.scn_locs_from_fs(scn_save_fs, run_tors_names)

    return scn_save_fs, scn_locs


def _ts_filesys(spc_dct, rxn, pf_levels, save_prefix, level='harm'):
    """ Set needed conformer filesys objects
    """

    if level == 'harm':
        thy_info = pf_levels[2]
    elif level == 'vpt2':
        thy_info = pf_levels[2]

    # Set the filesystem objects
    spc_info = (spc_dct['ich'], spc_dct['chg'], spc_dct['mul'])
    mod_thy_info = fsorb.mod_orb_restrict(spc_info, thy_info)
    rxn_info = finf.rxn_info(
        rxn['reacs'], rxn['prods'], spc_dct)
    _, thy_save_path = fbuild.rxn_thy_fs_from_root(
        save_prefix, rxn_info, mod_thy_info)
    ts_save_fs, ts_save_path = fbuild.ts_fs_from_thy(thy_save_path)

    return ts_save_fs, ts_save_path

