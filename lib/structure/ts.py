"""
  Some TS stuff
"""

from autofile import fs


def rxn_bnd_keys(cnf_fs, cnf_locs, zma_locs=(0)):
    """ get bond broken and formed keys for a transition state
    """

    cnf_path = cnf_fs[-1].path(cnf_locs)
    zma_fs = fs.manager(cnf_path, 'ZMATRIX')
    tra = zma_fs[-1].file.transformation.read(zma_locs)
    frm_bnd_keys, brk_bnd_keys = tra
    frm_bnd_keys = next(iter(frm_bnd_keys))
    brk_bnd_keys = next(iter(brk_bnd_keys))

    return frm_bnd_keys, brk_bnd_keys
