"""
Build paths for spc and thy file systems
"""

import autofile
from lib.filesystem import inf as finf
from lib.filesystem import minc as fsmin
from lib.filesystem import orb as fsorb


def set_fs(spc_dct, spc, thy_info,
           run_prefix, save_prefix,
           setfs_chk=True, ini_fs=False):
    """ set up filesystem """

    # Build the species filesys objs
    print('spc info')
    print(spc_dct[spc])
    [spc_info,
     spc_run_fs, spc_save_fs,
     spc_run_path, spc_save_path] = set_spc_fs(
         spc_dct, spc, run_prefix, save_prefix)

    # Initialize the filesystems
    thy_run_fs = None
    thy_run_path = None
    thy_save_fs = None
    thy_save_path = None
    cnf_run_fs = None
    cnf_save_fs = None
    tau_run_fs = None
    tau_save_fs = None
    thy_level = thy_info
    scn_run_fs = None
    scn_save_fs = None

    # Build the theory filesys obj
    thy_level = thy_info
    orb_restr = fsorb.orbital_restriction(
        spc_info, thy_info)
    thy_level = thy_info[0:3]
    thy_level.append(orb_restr)
    thy_run_fs = autofile.fs.theory(spc_run_path)
    thy_save_fs = autofile.fs.theory(spc_save_path)

    if setfs_chk:
        if 'ts_' in spc:
            thy_run_fs[-1].create(thy_level[1:4])
            thy_run_path = thy_run_fs[-1].path(thy_level[1:4])
            thy_save_fs[-1].create(thy_level[1:4])
            thy_save_path = thy_save_fs[-1].path(thy_level[1:4])

            thy_run_fs = autofile.fs.ts(thy_run_path)
            thy_run_fs[0].create()
            thy_run_path = thy_run_fs[0].path()

            thy_save_fs = autofile.fs.ts(thy_save_path)
            thy_save_fs[0].create()
            thy_save_path = thy_save_fs[0].path()

        else:
            thy_run_fs[-1].create(thy_level[1:4])
            thy_run_path = thy_run_fs[-1].path(thy_level[1:4])
            thy_save_fs[-1].create(thy_level[1:4])
            thy_save_path = thy_save_fs[-1].path(thy_level[1:4])

        cnf_run_fs = autofile.fs.conformer(thy_run_path)
        cnf_save_fs = autofile.fs.conformer(thy_save_path)
        tau_run_fs = autofile.fs.tau(thy_run_path)
        tau_save_fs = autofile.fs.tau(thy_save_path)
        min_cnf_locs = fsmin.min_energy_conformer_locators(
            cnf_save_fs)
        if min_cnf_locs:
            min_cnf_run_path = cnf_run_fs[-1].path(min_cnf_locs)
            min_cnf_save_path = cnf_save_fs[-1].path(min_cnf_locs)
            scn_run_fs = autofile.fs.conformer(min_cnf_run_path)
            scn_save_fs = autofile.fs.conformer(min_cnf_save_path)

    # filesys = [spc_run_fs, spc_save_fs, thy_run_fs, thy_save_fs,
    #            cnf_run_fs, cnf_save_fs, tau_run_fs, tau_save_fs,
    #            scn_run_fs, scn_save_fs]

    # Add run fs if needed
    if not ini_fs:
        filesys = [spc_run_fs, spc_save_fs, 
                   thy_run_fs, thy_save_fs,
                   cnf_run_fs, cnf_save_fs,
                   tau_run_fs, tau_save_fs,
                   scn_run_fs, scn_save_fs]
        run_fs = autofile.fs.run(thy_run_path)
        run_fs[0].create()
        filesys.append(run_fs)
    else:
        filesys = [thy_run_fs, thy_save_fs,
                   cnf_run_fs, cnf_save_fs,
                   tau_run_fs, tau_save_fs,
                   scn_run_fs, scn_save_fs]

    return filesys, thy_level


def set_spc_fs(spc_dct, spc, run_prefix, save_prefix):
    """ set the species filesys objects
    """
    if 'ts_' in spc:  # and rad_rad_ts != 'pst':
        run_fs, save_fs, run_path, save_path = spc_dct[spc]['rxn_fs']
        info = finf.get_spc_info(spc_dct[spc])
    else:
        info = finf.get_spc_info(spc_dct[spc])
        run_fs = autofile.fs.species(run_prefix)
        run_fs[-1].create(info)
        run_path = run_fs[-1].path(info)

        save_fs = autofile.fs.species(save_prefix)
        save_fs[-1].create(info)
        save_path = save_fs[-1].path(info)

    return info, run_fs, save_fs, run_path, save_path


def spc_fs_from_root(root_prefix, spc_info):
    """ Create species filesystem object and path
    """
    spc_fs = autofile.fs.species(root_prefix)
    spc_fs[-1].create(spc_info)
    spc_path = spc_run_fs[-1].path(spc_info)

    return spc_fs, spc_path


def rxn_fs_from_root(root_prefix, rxn_info):
    """ get filesystems for a reaction
        rxn_info = [rxn_ichs, rxn_chgs, rxn_muls, rxn_mul]
    """
    rxn_fs = autofile.fs.reaction(root_prefix)
    rxn_fs[-1].create(rxn_info)
    rxn_path = rxn_run_fs[-1].path(rxn_info)

    return rxn_fs, rxn_path


def thy_fs_from_root(root_prefix, spc_info, thy_info):
    """ create theory run path
    """
    # Build the species filesystem
    _, spc_path = spc_fs_from_root(root_prefix, spc_info)
    
    # Set the thy level needed for the filesys
    orb_restr = fsorb.orbital_restriction(
        spc_info, thy_info)
    thy_lvl = thy_info[0:3]
    thy_lvl.append(orb_restr)

    # Build the theory filesystem
    thy_fs = autofile.fs.theory(spc_path)
    thy_fs[-1].create(thy_lvl)
    thy_path = thy_run_fs[-1].path(thy_lvl)

    return thy_fs, thy_path


def thy_fs_from_spc(spc_prefix, thy_info):
    """ create theory run path
    """
    # Set the thy level needed for the filesys
    orb_restr = fsorb.orbital_restriction(
        spc_info, thy_info)
    thy_lvl = thy_info[0:3]
    thy_lvl.append(orb_restr)

    # Build the theory filesystem
    thy_fs = autofile.fs.theory(spc_prefix)
    thy_fs[-1].create(thy_lvl)
    thy_path = thy_run_fs[-1].path(thy_lvl)

    return thy_fs, thy_path


def cnf_fs_from_root(root_prefix, spc_info, thy_info, cnf='min'):
    """ create theory run path
    """
    assert cnf in ('min', 'all')

    # Build the theory filesystem
    _, thy_path = thy_fs_from_root(root_prefix, spc_info, thy_info)

    # Conformer filesys using theory
    cnf_fs = autofile.fs.conformer(thy_path)
    
    # Set the locs and the path
    if cnf == 'min':
        cnf_locs = fsmin.min_energy_conformer_locators(cnf_fs)
    elif cnf == 'all':
        cnf_locs = thy_fs.leaf.existing()
    cnf_paths = []
    if cnf_locs:
        for locs in cnf_locs:
            cnf_paths.append(cnf_run_fs[-1].path(cnf_locs))

    return cnf_fs, cnf_locs


def cnf_fs_from_thy(thy_prefix, spc_info, thy_info):
    """ create theory run path
    """

    # Conformer filesys using theory
    cnf_fs = autofile.fs.conformer(thy_prefix)

    # Set the locs and the path
    if cnf == 'min':
        cnf_locs = fsmin.min_energy_conformer_locators(cnf_fs)
    elif cnf == 'all':
        cnf_locs = thy_fs.leaf.existing()
    cnf_paths = []
    if cnf_locs:
        for locs in cnf_locs:
            cnf_paths.append(cnf_fs[-1].path(cnf_locs))

    return cnf_fs, cnf_locs, cnf_paths


def tau_fs_from_root(root_prefix, spc_info, thy_info, tau='all'):
    """ create theory run path
    """
    assert tau in ('all')

    # Build the theory filesystem
    _, thy_path = thy_fs_from_root(root_prefix, spc_info, thy_info)

    # Conformer filesys
    tau_fs = autofile.fs.tau(thy_path)
    
    # Set the locs and the path
    if tau == 'all':
        tau_locs = thy_fs.leaf.existing()
    tau_paths = []
    if tau_locs:
        for locs in tau_locs:
            tau_paths.append(tau_fs[-1].path(tau_locs))

    return tau_fs, tau_locs


def tau_fs_from_thy(thy_prefix, tau='all'):
    """ create theory run path
    """
    assert tau in ('all')

    # Conformer filesys
    tau_fs = autofile.fs.tau(thy_prefix)
    
    # Set the locs and the path
    if tau == 'all':
        tau_locs = thy_fs.leaf.existing()
    tau_paths = []
    if tau_locs:
        for locs in tau_locs:
            tau_paths.append(tau_fs[-1].path(tau_locs))

    return tau_fs, tau_locs


def scn_fs_from_root(root_prefix, spc_info, thy_info, cnf='min'):
    """ create theory run path
    """

    # Build the conformer filesystem
    _, _, cnf_paths = cnf_fs_from_root(root_prefix, spc_info, thy_info, cnf=cnf)
    
    # Set scan filesys
    scn_fs = autofile.fs.scan(cnf_paths[0])
    scn_locs = scn_fs[-1].existing([coo_names]):

    return scn_fs, scn_locs


def scn_fs_from_cnf(cnf_prefix, spc_info, thy_info, cnf='min'):
    """ create theory run path
    """

    # Set scan filesys
    scn_fs = autofile.fs.scan(cnf_prefix)
    scn_locs = scn_fs[-1].existing([coo_names]):

    return scn_fs, scn_locs


def cscn_fs_from_root(root_prefix, spc_info, thy_info, cnf='min'):
    """ create theory run path
    """

    # Build the conformer filesystem
    _, _, cnf_paths = cnf_fs_from_root(root_prefix, spc_info, thy_info, cnf=cnf)
    
    # Set scan filesys
    cscn_fs = autofile.fs.cscan(cnf_paths[0])
    cscn_locs = cscn_fs[-1].existing([coo_names]):

    return cscn_fs, cscn_locs


def cscn_fs_from_cnf(cnf_prefix):
    """ create theory run path
    """

    # Set scan filesys
    cscn_fs = autofile.fs.cscan(cnf_prefix)
    cscn_locs, cscn_paths = [], []
    if cscn_run_fs[1].exists([coo_names]):
        scn_locs1 = cscn_run_fs[2].existing([coo_names])
        for locs1 in scn_locs1:
            if cscn_run_fs[2].exists(locs1):
                scn_locs2 = cscn_run_fs[3].existing(locs1)
                for locs2 in scn_locs2:
                    cscn_path = cscn_run_fs[-1].path(locs2)

    return cscn_fs, cscn_locs


def cscn_fs_from_ts(ts_prefix):
    """ create theory run path
    """
    return cscn_fs, cscn_locs


def ts_fs_from_root(root_prefix, spc_info, thy_info):
    """ Create species filesystem object and path
    """
    # Build the theory filesystem
    _, thy_path = thy_fs_from_root(root_prefix, spc_info, thy_info)
    
    ts_fs = autofile.fs.ts(thy_path)
    ts_fs[0].create()
    ts_path = thy_fs[0].path()

    return ts_fs, ts_path


def ts_fs_from_thy(thy_prefix, spc_info, thy_info):
    """ Create species filesystem object and path
    """
    
    ts_fs = autofile.fs.ts(thy_prefix)
    ts_fs[0].create()
    ts_path = thy_fs[0].path()

    return ts_fs, ts_path



def get_ts_fs(rxn_run_path, rxn_save_path, ref_level):
    """ Build a transition state filesystem
    """
    thy_run_fs = autofile.fs.theory(rxn_run_path)
    thy_run_fs[-1].create(ref_level[1:4])
    thy_run_path = thy_run_fs[-1].path(ref_level[1:4])

    thy_save_fs = autofile.fs.theory(rxn_save_path)
    thy_save_fs[-1].create(ref_level[1:4])
    thy_save_path = thy_save_fs[-1].path(ref_level[1:4])

    scn_run_fs = autofile.fs.scan(thy_run_path)
    scn_save_fs = autofile.fs.scan(thy_save_path)

    ts_run_fs = autofile.fs.ts(thy_run_path)
    ts_run_fs[0].create()
    ts_run_path = ts_run_fs[0].path()
    run_fs = autofile.fs.run(ts_run_path)

    ts_save_fs = autofile.fs.ts(thy_save_path)
    ts_save_fs[0].create()
    ts_save_path = ts_save_fs[0].path()

    # cnf_run_fs = autofile.fs.conformer(ts_run_path)
    # cnf_save_fs = autofile.fs.conformer(ts_save_path)
    # cnf_save_fs[0].create()

    return ts_run_fs, ts_save_fs, ts_run_path, ts_save_path


def get_rxn_fs(run_prefix, save_prefix,
               rxn_ichs, rxn_chgs, rxn_muls, ts_mul):
    """ get filesystems for a reaction
    """
    rxn_run_fs = autofile.fs.reaction(run_prefix)
    rxn_run_fs[-1].create(
        [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
    rxn_run_path = rxn_run_fs[-1].path(
        [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

    rxn_ichs = tuple(map(tuple, rxn_ichs))
    rxn_chgs = tuple(map(tuple, rxn_chgs))
    rxn_muls = tuple(map(tuple, rxn_muls))
    rxn_save_fs = autofile.fs.reaction(save_prefix)
    rxn_save_fs[-1].create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
    rxn_save_path = rxn_save_fs[-1].path(
        [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

    return rxn_run_fs, rxn_save_fs, rxn_run_path, rxn_save_path
