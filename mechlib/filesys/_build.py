""" New autofile build interface
"""

import os
import autofile
from mechanalyzer.inf import rxn as rinfo
from mechanalyzer.inf import spc as sinfo


def root_locs(spc_dct_i, saddle=False, name=None):
    """ Set the root locatores for the species and TS
    """

    if not saddle:
        spc_info = sinfo.from_dct(spc_dct_i)
        rxn_info = None
        ts_info = None
    else:
        spc_info = None
        rxn_info = rinfo.sort(spc_dct_i['rxn_info'])
        ts_num = int(name.split('_')[-1])
        ts_info = (ts_num,)
        # print('TEST ts locs:', ts_info)
        # ts_info = (0,)  # may be more complicated
        # ts_info = ()

    return {'spc_locs': spc_info, 'rxn_locs': rxn_info, 'ts_locs': ts_info}


def build_fs(run_prefix, save_prefix, end,
             rxn_locs=None, spc_locs=None,
             thy_locs=None, ts_locs=None,
             cnf_locs=None, tau_locs=None,
             zma_locs=None,
             scn_locs=None, cscn_locs=None):
    """ Build the filesystems
    """

    _fs = []
    for prefix in (run_prefix, save_prefix):
        _fs.append(
            _build_fs(
                prefix, end,
                rxn_locs=rxn_locs, spc_locs=spc_locs,
                thy_locs=thy_locs, ts_locs=ts_locs,
                cnf_locs=cnf_locs, tau_locs=tau_locs,
                zma_locs=zma_locs,
                scn_locs=scn_locs, cscn_locs=cscn_locs)
        )

    return _fs[0], _fs[1]


def _build_fs(prefix, end,
              rxn_locs=None, spc_locs=None,
              thy_locs=None, ts_locs=None,
              cnf_locs=None, tau_locs=None,
              zma_locs=None,
              scn_locs=None, cscn_locs=None):
    """ Build the filesystems
    """

    prefix_locs = []
    if rxn_locs is not None:
        prefix_locs.append(('REACTION', rxn_locs))
    if spc_locs is not None:
        prefix_locs.append(('SPECIES', spc_locs))
    if thy_locs is not None:
        prefix_locs.append(('THEORY', thy_locs))
    if ts_locs is not None:
        prefix_locs.append(('TRANSITION STATE', ts_locs))
    if cnf_locs is not None:
        prefix_locs.append(('TAU', tau_locs))
    if cnf_locs is not None:
        prefix_locs.append(('CONFORMER', cnf_locs))
    if zma_locs is not None:
        prefix_locs.append(('ZMATRIX', zma_locs))
    if scn_locs is not None:
        prefix_locs.append(('SCAN', scn_locs))
    if cscn_locs is not None:
        prefix_locs.append(('CSCAN', cscn_locs))

    _fs = autofile.fs.manager(prefix, prefix_locs, end)

    return _fs


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
