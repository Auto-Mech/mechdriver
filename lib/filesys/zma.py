"""
  Toy functions for building the ZMA locs
"""

import automol


def get_zma_locs(spc_dct, zma_fs, rxn_info=None, dirn='forw'):
    """ get locs

        rxn = None, or (rct_name/ich, prd_name/ich)
    """

    # Grab the user-input zma idx from the dct if they exist
    dct_zma_idx = spc_dct.get('zma_idx', None)

    if dct_zma_idx is None:
        # Get idx 0 for species; get idx for correct dirn for rxn
        if rxn is None:
            zma_locs = (0,)
        else:
            zma_locs = rxn_zma_locs(zma_fs, rxn_info, dirn)
    else:
        zma_locs = (dct_zma_idx,)

    return zma_locs


def rxn_zma_locs(zma_fs, rxn_info, wanted_dirn):
    """ Get the zma locs that are needed
    """

    assert dirn in ('forw', 'any')

    locs_dirn_lst = tuple()
    for locs in zma_fs[-1].existing():
        ts_zma = zma_fs[-1].file.zmatrix.read(locs)
        tra = zma_fs[-1].file.transformation.read(locs)
        dirn = chk_direction(rxn_info, ts_zma,
                             frm_bnd_keys, brk_bnd_keys)
        locs_dirn_lst += ((locs, dirn),)

    zma_locs = None
    for (locs, dirn) in locs_dirn_lst:
        if dirn == 'forw' and wanted_dirn == 'forw':
            zma_locs = locs
            break
        elif dirn == 'forw' and wanted_dirn == 'any':
            zma_locs = locs
            break

    return zma_locs


def chk_direction(rxn_info, ts_zma,
                  frm_bnd_keys, brk_bnd_keys):
    """ Get the zma locs that are needed
    """

    fs_rct_ichs = automol.zmatrix.ts.zmatrix_product_inchis(
        ts_zma, frm_bnd_keys, brk_bnd_keys, remove_stereo=False)
    fs_rct_ichs = set(fs_rct_ichs)

    fs_prd_ichs = automol.zmatrix.ts.zmatrix_product_inchis(
        ts_zma, frm_bnd_keys, brk_bnd_keys, remove_stereo=False)
    fs_prd_ichs = set(fs_prd_ichs)

    rct_ichs, prd_ichs = rxn_info[0][0], rxn_info[0][1]
    if set(rct_ichs) == fs_rct_ichs and set(prd_ichs) == fs_prd_ichs:
        dirn = 'forw'
    elif set(rct_ichs) == fs_prd_ichs and set(prd_ichs) == fs_rct_ichs:
        dirn = 'rev'
    else:
        dirn = None

    return dirn
