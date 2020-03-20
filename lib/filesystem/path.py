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
