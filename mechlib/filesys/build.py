"""
Library of functions to interact with the filesystem
"""

import os
import sys
import autofile
from autofile import fs
from mechlib.filesys.mincnf import conformer_locators


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

