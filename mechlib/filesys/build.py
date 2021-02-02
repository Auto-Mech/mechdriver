"""
Library of functions to interact with the filesystem
"""

import os
import sys
import autofile
from autofile import fs
from mechlib.filesys.mincnf import conformer_locators


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


def cnf_fs_from_thy(thy_prefix, mod_thy_info, cnf=None, saddle=False):
    """ create theory run path
    """
    # Build an intermediate TS filesystem if needed
    if saddle:
        ts_fs = autofile.fs.transition_state(thy_prefix)
        ts_fs[0].create()
        cnf_prefix = ts_fs[0].path()
    else:
        cnf_prefix = thy_prefix

    # Conformer filesys using theory
    cnf_fs = autofile.fs.conformer(cnf_prefix)

    # Set the locs and the path
    cnf_locs = []
    if cnf is not None:
        if cnf == 'all':
            cnf_locs = cnf_fs[1].existing()
        else:
            cnf_locs = conformer_locators(
                cnf_fs, mod_thy_info, cnf_range=cnf)

    return cnf_fs, cnf_locs


def cnf_fs_from_prefix(cnf_prefix, mod_thy_info, cnf=None):
    """ create theory run path
    """
    # Conformer filesys using theory
    cnf_fs = autofile.fs.conformer(cnf_prefix)

    # Set the locs and the path
    cnf_locs = []
    if cnf is not None:
        if cnf == 'all':
            cnf_locs = cnf_fs[1].existing()
        else:
            cnf_locs = conformer_locators(
                cnf_fs, mod_thy_info, cnf_range=cnf)

    return cnf_fs, cnf_locs


def cnf_paths_from_locs(cnf_fs, cnf_locs):
    """ Get paths
    """
    cnf_paths = []
    if cnf_locs:
        for locs in cnf_locs:
            cnf_paths.append(cnf_fs[-1].path(locs))

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
    if constraint_dct is None:
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
    """ Build high-spin filesystem from prefix
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

    ts_fs = autofile.fs.transition_state(thy_path)
    ts_fs[0].create()
    ts_path = ts_fs[0].path()

    return ts_fs, ts_path


def ts_fs_from_thy(thy_prefix):
    """ Create species filesystem object and path
    """

    ts_fs = autofile.fs.transition_state(thy_prefix)
    ts_fs[0].create()
    ts_path = ts_fs[0].path()

    return ts_fs, ts_path


def vrc_fs_from_thy(thy_prefix):
    """ Create species filesystem object and path
    """

    vrc_fs = autofile.fs.vrctst(thy_prefix)
    vrc_fs[0].create()
    vrc_path = vrc_fs[0].path()

    return vrc_fs, vrc_path


def run_fs_from_prefix(prefix):
    """ Build a run filesystem object forom some generic prefix
    """
    run_fs = autofile.fs.run(prefix)
    run_fs[0].create()

    return run_fs


def zma_fs_from_prefix(prefix, zma_idxs=(0,)):
    """ Build a zma filesys object
    """
    zma_fs = fs.zmatrix(prefix)
    zma_fs[-1].create(zma_idxs)
    zma_path = zma_fs[-1].path(zma_idxs)

    return zma_fs, zma_path


def zma_from_prefix(prefix, spc_dct_i, rxn_ichs=None, dirn='forw'):
    """ get locs
    """

    # Build the zma fs
    zma_fs = fs.zmatrix(prefix)

    # Get the locs
    # zma_locs = get_zma_locs(zma_fs, spc_dct_i, rxn_ichs=rxn_ichs, wanted_dirn=dirn)
    zma_locs = (0,)

    # print(zma_locs)
    zma_fs[-1].create(zma_locs)
    zma_path = zma_fs[-1].path(zma_locs)

    return zma_fs, zma_path
    # return zma_fs, zma_locs


def get_zma_locs(zma_fs, spc_dct_i, rxn_ichs=None, wanted_dirn=('forw',)):
    """ get locs
        rxn = None, or (rct_name/ich, prd_name/ich)
    """
    
    assert set(wanted_dirn) <= {'forw', 'rev'}

    # Grab the user-input zma idx from the dct if they exist
    dct_zma_idx = spc_dct_i.get('zma_idx', None)

    if dct_zma_idx is None:
        if rxn_ichs is None:
            # Right now we assume that we will just grab 0th idx species
            zma_locs = (0,)
        else:
            # Grab the lowest idx associated with the direction
            zma_locs_lst = rxn_zma_locs_lst(zma_fs, rxn_ichs)
            if zma_locs_lst:
                want_lst = tuple(locs for locs, dirn in zma_locs_lst
                                 if dirn in wanted_dirn)
                zma_locs = (want_lst[0],)
            else:
                zma_locs = None
    else:
        zma_locs = (dct_zma_idx,)

    return zma_locs


def spc_zma_locs_lst(zma_fs):
    """ Get the zma locs for the species
    """
    locs_lst = tuple()
    for locs in zma_fs[-1].existing():
        locs_lst += ((locs, None),)

    return locs_lst


def rxn_zma_locs_lst(zma_fs, rxn_ichs):
    """ Get the zma locs that are needed
    """

    assert wanted_dirn in ('forw', 'rev', 'any')

    locs_lst = tuple()
    for locs in zma_fs[-1].existing():
        # Not placed in the scans filesystem
        # maybe place the vmatrix and transformation file?
        if zma_fs[-1].file.zmatrix.exists(locs):
            ts_zma = zma_fs[-1].file.zmatrix.read(locs)
            frm_bnd_keys, brk_bnd_keys = zma_fs[-1].file.transformation.read(locs)
            dirn = _chk_direction(rxn_ichs, ts_zma,
                                  frm_bnd_keys, brk_bnd_keys)
            locs_lst += ((locs, dirn),)

    return zma_locs


def _chk_direction(rxn_ichs, ts_zma,
                   frm_bnd_keys, brk_bnd_keys):
    """ Get the zma locs that are needed
    """

    fs_rct_ichs = automol.zmatrix.ts.zmatrix_reactant_inchis(
        ts_zma, frm_bnd_keys, brk_bnd_keys, remove_stereo=False)
    fs_rct_ichs = set(fs_rct_ichs)

    fs_prd_ichs = automol.zmatrix.ts.zmatrix_product_inchis(
        ts_zma, frm_bnd_keys, brk_bnd_keys, remove_stereo=False)
    fs_prd_ichs = set(fs_prd_ichs)

    rct_ichs, prd_ichs = rxn_ichs[0], rxn_ichs[1]
    print('rxn ich lst check')
    print('rct ich', set(rct_ichs))
    print('fs rct ich', fs_rct_ichs)
    print('prd ich', set(prd_ichs))
    print('fs prd ich', fs_prd_ichs)
    if set(rct_ichs) == fs_rct_ichs and set(prd_ichs) == fs_prd_ichs:
        dirn = 'forw'
    elif set(rct_ichs) == fs_prd_ichs and set(prd_ichs) == fs_rct_ichs:
        dirn = 'rev'
    else:
        dirn = None

    return dirn


def etrans_fs_from_prefix(prefix, bath_info, thy_info):
    """ Build an energy transfer filesystem obj
    """

    # Build the energy transfer filesys object
    print('prefix', prefix)
    etrans_fs = autofile.fs.energy_transfer(prefix)

    # Build the energy transfer locs object
    etrans_locs = bath_info + thy_info[1:4]
   
    # Create directory
    print('etrans locs', etrans_locs)
    etrans_fs[-1].create(etrans_locs)

    return etrans_fs, etrans_locs


# Old function that I need to get rid of




def etrans_fs_from_prefix(prefix, bath_info, thy_info):
    """ Build an energy transfer filesystem obj
    """

    # Build the energy transfer filesys object
    etrans_fs = autofile.fs.energy_transfer(prefix)

    # Build the energy transfer locs object
    etrans_locs = bath_info + thy_info[1:4]

    # Create directory
    print('etrans locs', etrans_locs)
    etrans_fs[-1].create(etrans_locs)

    return etrans_fs, etrans_locs


def build_fs(prefix, run, locs_idx=None):
    """ Create a build filesys object for some procedure
    """

    # Initialize the build object
    bld_run_fs = autofile.fs.build(prefix)

    # Set the build locs
    if locs_idx is not None:
        assert isinstance(locs_idx, int), (
            'locs idx {} is not an integer'.format(locs_idx)
        )
    else:
        bld_run_fs[-1].create([run, 0])
        existing_locs = bld_run_fs[-1].existing()
        current_idxs = tuple(idx for [name, idx] in existing_locs
                             if name == run)
        locs_idx = max(current_idxs) + 1

    bld_locs = [run, locs_idx]

    return bld_run_fs, bld_locs


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

