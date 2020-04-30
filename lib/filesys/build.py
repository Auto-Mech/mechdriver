"""
Library of functions to interact with the filesystem
"""

import os
import sys
import autofile
from lib.filesys.mincnf import min_energy_conformer_locators


# BUILD HIGH LEVEL SPC AND RXN FILESYS #
def prefix_fs(prefix):
    """ Physically make the run/save filesys root given a prefix path
        :param str prefix: file path - /path/to/root/run
    """
    if not os.path.exists(prefix):
        try:
            os.mkdir(prefix)
        except FileNotFoundError:
            print('Cannot make directory at path specified in run.dat.')
            print('Path: {}'.format(prefix))
            sys.exit()


# BUILD HIGH LEVEL SPC AND RXN FILESYS #
def spc_fs_from_root(root_prefix, spc_info):
    """ Create root species filesystem:
            /root_prefix/SPC
        :param root_prefix: path where SPC filesys will be created
        :param list spc_info: [InChI String, Charge, Multiplicity]
        :return spc_fs: SPC filesystem object
        :rtype: autofile.DataSeries object
        :return spc_path: filesystem path
        :rtype: string
    """
    spc_fs = autofile.fs.species(root_prefix)
    spc_fs[-1].create(spc_info)
    spc_path = spc_fs[-1].path(spc_info)

    return spc_fs, spc_path


def rxn_fs_from_root(root_prefix, rxn_info):
    """ Create root reaction filesystem:
            /root_prefix/RXN
        :param root_prefix: path where RXN filesys will be created
        :param list Rxn_info: [RP InChIs, RP Charges, RP Mult, RXN Mult]
                              where RP=Reacs and Prods, RXN=TS
        :return rxn_fs: RXN filesystem object
        :rtype: autofile.DataSeries object
        :return rxn_path: filesystem path
        :rtype: string
    """
    rxn_fs = autofile.fs.reaction(root_prefix)
    rxn_fs[-1].create(rxn_info)
    rxn_path = rxn_fs[-1].path(rxn_info)

    return rxn_fs, rxn_path


# BUILD THEORY FILESYS THAT ARE UNDERNEATH SPC AND RXN FILESYS #
def spc_thy_fs_from_root(root_prefix, spc_info, mod_thy_info):
    """ Create ES Theory filesystem within a SPC filesytem:
            /root_prefix/SPC/THY
        :param root_prefix: path where SPC filesys will be created
        :param list spc_info: [InChI String, Charge, Multiplicity]
        :param list mod_thy_info: [Program, Method, Basis, OrbLetter(R/U)]
        :return thy_fs: THY filesystem object
        :rtype: autofile.DataSeries object
        :return thy_path: filesystem path
        :rtype: string
    """
    # Build the species filesystem
    _, spc_path = spc_fs_from_root(root_prefix, spc_info)

    # Build the theory filesystem
    thy_fs = autofile.fs.theory(spc_path)
    thy_fs[-1].create(mod_thy_info[1:4])
    thy_path = thy_fs[-1].path(mod_thy_info[1:4])

    return thy_fs, thy_path


def rxn_thy_fs_from_root(root_prefix, rxn_info, mod_thy_info):
    """ Create ES Theory filesystem within a RXN filesystem:
            /root_prefix/RXN/THY
        :param root_prefix: path where RXN filesys will be created
        :param list rxn_info: [RP InChIs, RP Charges, RP Mult, RXN Mult]
                              where RP=Reacs and Prods, RXN=TS
        :param list mod_thy_info: [Program, Method, Basis, OrbLetter(R/U)]
        :return thy_fs: THY filesystem object
        :rtype: autofile.DataSeries object
        :return thy_path: filesystem path
        :rtype: string
    """
    # Build the species filesystem
    _, rxn_path = rxn_fs_from_root(root_prefix, rxn_info)

    # Build the theory filesystem
    thy_fs = autofile.fs.theory(rxn_path)
    thy_fs[-1].create(mod_thy_info[1:4])
    thy_path = thy_fs[-1].path(mod_thy_info[1:4])

    return thy_fs, thy_path


def thy_fs_from_prefix(prefix, mod_thy_info):
    """ Create root reaction filesystem using generic prefix path:
            /prefix/THY
        :param prefix: path where RXN filesys will be created
        :param list mod_thy_info: [Program, Method, Basis, OrbLetter(R/U)]
        :return thy_fs: THY filesystem object
        :rtype: autofile.DataSeries object
        :return thy_path: filesystem path
        :rtype: string
    """
    thy_fs = autofile.fs.theory(prefix)
    thy_fs[-1].create(mod_thy_info[1:4])
    thy_path = thy_fs[-1].path(mod_thy_info[1:4])

    return thy_fs, thy_path


# BUILD CONFORMER FILESYS THAT ARE UNDERNEATH SPC(RXN)/THY FILESYS #
def spc_cnf_fs_from_root(root_prefix, spc_info, mod_thy_info, saddle=False):
    """ create theory run path
    """
    # Build the theory filesystem
    _, thy_path = spc_thy_fs_from_root(root_prefix, spc_info, mod_thy_info)

    # Build an intermediate TS filesystem if needed
    if saddle:
        _, cnf_prefix = ts_fs_from_thy(thy_path)
    else:
        cnf_prefix = thy_path

    # Conformer filesys using theory
    cnf_fs = autofile.fs.conformer(cnf_prefix)

    return cnf_fs


def rxn_cnf_fs_from_root(root_prefix, spc_info, mod_thy_info, saddle=False):
    """ create theory run path
    """
    # Build the theory filesystem
    _, thy_path = spc_thy_fs_from_root(root_prefix, spc_info, mod_thy_info)

    # Build an intermediate TS filesystem if needed
    if saddle:
        _, cnf_prefix = ts_fs_from_thy(thy_path)
    else:
        cnf_prefix = thy_path

    # Conformer filesys using theory
    cnf_fs = autofile.fs.conformer(cnf_prefix)

    return cnf_fs


def cnf_fs_from_thy(thy_prefix, cnf=None, saddle=False):
    """ create theory run path
    """
    # Build an intermediate TS filesystem if needed
    if saddle:
        ts_fs = autofile.fs.ts(thy_prefix)
        ts_fs[0].create()
        cnf_prefix = ts_fs[0].path()
    else:
        cnf_prefix = thy_prefix

    # Conformer filesys using theory
    cnf_fs = autofile.fs.conformer(cnf_prefix)

    # Set the locs and the path
    cnf_locs = []
    if cnf is not None:
        if cnf == 'min':
            cnf_locs = min_energy_conformer_locators(cnf_fs)
        elif cnf == 'all':
            cnf_locs = cnf_fs[1].existing()

    # elif selection == 'subset':
    #     min_locs = min_energy_conformer_locators(save_dir)
    #     min_ene = energy.read(min_lcs)
    #     ini_locs_lst = save_dir[-1].existing()
    #     locs_lst = []
    #     for locs in ini_locs_lst:
    #         ene = energy.read(locs) - min_ene
    #         if ene <= ene_cut_off:
    #             locs_lst.append(locs)

    return cnf_fs, cnf_locs


def cnf_fs_from_prefix(cnf_prefix, cnf=None):
    """ create theory run path
    """
    # Conformer filesys using theory
    cnf_fs = autofile.fs.conformer(cnf_prefix)

    # Set the locs and the path
    cnf_locs = []
    if cnf is not None:
        if cnf == 'min':
            cnf_locs = min_energy_conformer_locators(cnf_fs)
        elif cnf == 'all':
            cnf_locs = cnf_fs[1].existing()

    return cnf_fs, cnf_locs


def cnf_paths_from_locs(cnf_fs, cnf_locs):
    """ Get paths
    """
    cnf_paths = []
    if cnf_locs:
        for locs in cnf_locs:
            cnf_paths.append(cnf_fs[-1].path([locs]))

    return cnf_paths


def cnf_create(cnf_fs, cids):
    """ create a new cnf filesys using list of cids
    """
    for cid in cids:
        cnf_fs[-1].create([cid])


def tau_fs_from_root(root_prefix, spc_info, mod_thy_info, tau='all'):
    """ create theory run path
    """
    assert tau == 'all'

    # Build the theory filesystem
    _, thy_path = spc_thy_fs_from_root(root_prefix, spc_info, mod_thy_info)

    # Conformer filesys
    tau_fs = autofile.fs.tau(thy_path)

    # Set the locs and the path
    if tau == 'all':
        tau_locs = tau_fs[1].existing()
    tau_paths = []
    if tau_locs:
        for locs in tau_locs:
            tau_paths.append(tau_fs[-1].path(locs))

    return tau_fs, tau_locs


def tau_fs_from_thy(thy_prefix, tau='all'):
    """ create theory run path
    """
    assert tau == 'all'

    # Conformer filesys
    tau_fs = autofile.fs.tau(thy_prefix)

    # Set the locs and the path
    if tau == 'all':
        tau_locs = tau_fs[1].existing()
    tau_paths = []
    if tau_locs:
        for locs in tau_locs:
            tau_paths.append(tau_fs[-1].path(locs))

    return tau_fs, tau_locs


def scn_fs_from_cnf(cnf_prefix, constraint_dct=None):
    """ build either a SCAN or CSCAN FS
    """
    if constraint_dct is None:
        scn_fs = autofile.fs.scan(cnf_prefix)
    else:
        scn_fs = autofile.fs.cscan(cnf_prefix)

    return scn_fs


def scn_locs_from_fs(scn_fs, coo_names, constraint_dct=None):
    """ get locs for a SCAN or CSCAN fs
    """
    if constraint_dct is not None:
        scn_locs = scn_fs[-1].existing([coo_names])
    else:
        # scn_locs, scn_paths = [], []
        scn_locs = []
        if scn_fs[1].exists([coo_names]):
            scn_locs1 = scn_fs[2].existing([coo_names])
            for locs1 in scn_locs1:
                if scn_fs[2].exists(locs1):
                    scn_locs2 = scn_fs[3].existing(locs1)
                    for locs2 in scn_locs2:
                        scn_locs.append(locs2)
                        # scn_paths.append(scn_fs[-1].path(locs2))

    return scn_locs


def cscn_fs_from_ts(ts_prefix, coo_names):
    """ create theory run path
    """
    cscn_fs = autofile.fs.cscan(ts_prefix)
    cscn_locs, cscn_paths = [], []
    if cscn_fs[1].exists([coo_names]):
        scn_locs1 = cscn_fs[2].existing([coo_names])
        for locs1 in scn_locs1:
            if cscn_fs[2].exists(locs1):
                scn_locs2 = cscn_fs[3].existing(locs1)
                for locs2 in scn_locs2:
                    cscn_paths.append(cscn_fs[-1].path(locs2))

    return cscn_fs, cscn_locs


def sp_from_prefix(prefix, thy_info):
    """ Build single point fs from prefix
    """

    sp_fs = autofile.fs.single_point(prefix)
    sp_fs[-1].create(thy_info[1:4])
    sp_fs_path = sp_fs[-1].path(thy_info[1:4])

    return sp_fs, sp_fs_path


def high_spin_from_prefix(prefix, thy_info):
    """
    """

    hs_fs = autofile.fs.high_spin(prefix)
    hs_fs[-1].create(thy_info[1:4])
    hs_fs_path = hs_fs[-1].path(thy_info[1:4])

    return hs_fs, hs_fs_path


def ts_fs_from_root(root_prefix, spc_info, thy_info):
    """ Create species filesystem object and path
    """
    # Build the theory filesystem
    _, thy_path = rxn_thy_fs_from_root(root_prefix, spc_info, thy_info)

    ts_fs = autofile.fs.ts(thy_path)
    ts_fs[0].create()
    ts_path = ts_fs[0].path()

    return ts_fs, ts_path


def ts_fs_from_thy(thy_prefix):
    """ Create species filesystem object and path
    """

    ts_fs = autofile.fs.ts(thy_prefix)
    ts_fs[0].create()
    ts_path = ts_fs[0].path()

    return ts_fs, ts_path


def run_fs_from_prefix(prefix):
    """ Build a run filesystem object forom some generic prefix
    """
    run_fs = autofile.fs.run(prefix)
    run_fs[0].create()

    return run_fs

# Old function that I need to get rid of
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
