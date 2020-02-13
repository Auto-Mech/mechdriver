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
            thy_run_fs.leaf.create(thy_level[1:4])
            thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])
            thy_save_fs.leaf.create(thy_level[1:4])
            thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])

            thy_run_fs = autofile.fs.ts(thy_run_path)
            thy_run_fs.trunk.create()
            thy_run_path = thy_run_fs.trunk.path()

            thy_save_fs = autofile.fs.ts(thy_save_path)
            thy_save_fs.trunk.create()
            thy_save_path = thy_save_fs.trunk.path()

        else:
            thy_run_fs.leaf.create(thy_level[1:4])
            thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])
            thy_save_fs.leaf.create(thy_level[1:4])
            thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])

        cnf_run_fs = autofile.fs.conformer(thy_run_path)
        cnf_save_fs = autofile.fs.conformer(thy_save_path)
        tau_run_fs = autofile.fs.tau(thy_run_path)
        tau_save_fs = autofile.fs.tau(thy_save_path)
        min_cnf_locs = fsmin.min_energy_conformer_locators(
            cnf_save_fs)
        if min_cnf_locs:
            min_cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
            min_cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
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
        run_fs.trunk.create()
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
        run_fs.leaf.create(info)
        run_path = run_fs.leaf.path(info)

        save_fs = autofile.fs.species(save_prefix)
        save_fs.leaf.create(info)
        save_path = save_fs.leaf.path(info)

    return info, run_fs, save_fs, run_path, save_path


def get_spc_fs(prefix, spc_info):
    """ Create species filesystem object and path
    """
    spc_fs = autofile.fs.species(prefix)
    spc_fs.leaf.create(spc_info)
    spc_path = spc_run_fs.leaf.path(spc_info)

    return spc_fs, spc_path


def get_thy_fs(prefix, spc_info, thy_info):
    """ create theory run path
    """
    # Set the thy level needed for the filesys
    orb_restr = fsorb.orbital_restriction(
        spc_info, thy_info)
    thy_lvl = thy_info[0:3]
    thy_lvl.append(orb_restr)

    # Build the species filesystem
    _, spc_path = get_spc_run_path(prefix, spc_info)

    # Build the theory filesystem
    thy_fs = autofile.fs.theory(spc_path)
    thy_fs.leaf.create(thy_lvl)
    thy_path = thy_run_fs.leaf.path(thy_lvl)

    return thy_fs, thy_path


def get_cnf_fs(thy_prefix, spc_info, thy_info):
    """ create theory run path
    """
    cnf_fs = autofile.fs.conformer(thy_run_path)
    return cnf_fs


def get_min_cnf_fs(thy_prefix, spc_info, thy_info):
    """ create theory run path
    """
    cnf_fs = autofile.fs.conformer(thy_run_path)
    min_cnf_locs = fsmin.min_energy_conformer_locators(
        cnf_save_fs)
    if min_cnf_locs:
        min_cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
    return cnf_fs


def get_scn_fs(thy_prefix, spc_info, thy_info):
    """ create theory run path
    """
    cnf_fs = autofile.fs.conformer(thy_run_path)
    min_cnf_locs = fsmin.min_energy_conformer_locators(
        cnf_save_fs)
    if min_cnf_locs:
        scn_run_fs = autofile.fs.conformer(min_cnf_run_path)
        min_cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
    return cnf_fs



def get_rxn_fs(run_prefix, save_prefix,
               rxn_ichs, rxn_chgs, rxn_muls, ts_mul):
    """ get filesystems for a reaction
    """
    rxn_run_fs = autofile.fs.reaction(run_prefix)
    rxn_run_fs.leaf.create(
        [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
    rxn_run_path = rxn_run_fs.leaf.path(
        [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

    rxn_ichs = tuple(map(tuple, rxn_ichs))
    rxn_chgs = tuple(map(tuple, rxn_chgs))
    rxn_muls = tuple(map(tuple, rxn_muls))
    rxn_save_fs = autofile.fs.reaction(save_prefix)
    rxn_save_fs.leaf.create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
    rxn_save_path = rxn_save_fs.leaf.path(
        [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

    return rxn_run_fs, rxn_save_fs, rxn_run_path, rxn_save_path
