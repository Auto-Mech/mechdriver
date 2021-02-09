""" New autofile build interface
"""

import os
import autofile
from mechanalyzer.inf import rxn as rinfo
from mechanalyzer.inf import spc as sinfo


def root_locs(spc_dct_i, saddle=False):
    """ Set the root locatores for the species and TS
    """

    if not saddle:
        spc_info = sinfo.from_dct(spc_dct_i)
        rxn_info = None
    else:
        spc_info = None
        rxn_info = rinfo.sort(spc_dct_i['rxn_info'])

    return {'spc_locs': spc_info, 'rxn_locs': rxn_info}


def build_fs(run_prefix, save_prefix, end,
             rxn_locs=None, spc_locs=None,
             thy_locs=None, ts_locs=None,
             cnf_locs=None, zma_locs=None,
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
                cnf_locs=cnf_locs, zma_locs=zma_locs,
                scn_locs=scn_locs, cscn_locs=cscn_locs)
        )

    return tuple(_fs)


def _build_fs(prefix, end,
              rxn_locs=None, spc_locs=None,
              thy_locs=None, ts_locs=None,
              cnf_locs=None, zma_locs=None,
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
        prefix_locs.append(('CONFORMER', cnf_locs))
    if zma_locs is not None:
        prefix_locs.append(('ZMATRIX', zma_locs))
    if scn_locs is not None:
        prefix_locs.append(('SCAN', scn_locs))
    if cscn_locs is not None:
        prefix_locs.append(('CSCAN', cscn_locs))

    _fs = autofile.fs.manager(prefix, prefix_locs, end)

    return _fs


if __name__ == '__main__':
    RUN_PREFIX = os.path.join(os.getcwd(), 'run')
    SAVE_PREFIX = os.path.join(os.getcwd(), 'save')

    RXN_LOCS = []
    SPC_LOCS = ['InChI=1S/H2O/h1H2', 0, 1]
    THY_LOCS = ['gaussian09', 'wb97xd', 'cc-pvtz', 'R']
    CNF_LOCS = [autofile.schema.generate_new_conformer_id()]
    ZMA_LOCS = [0]
    SCN_LOCS = [['D2'], [2.15]]

    CNF_FS = build_fs(
        RUN_PREFIX, SAVE_PREFIX, 'CONFORMER',
        spc_locs=SPC_LOCS,
        thy_locs=THY_LOCS[1:])

    CNF_RUN_PATH = CNF_FS[0][-1].path(CNF_LOCS)
    CNF_SAVE_PATH = CNF_FS[1][-1].path(CNF_LOCS)
    CNF_FS[0][-1].create(CNF_LOCS)
    CNF_FS[1][-1].create(CNF_LOCS)
    print('cnf run', CNF_RUN_PATH)
    print('cnf save', CNF_SAVE_PATH)

    ZMA_FS = build_fs(
        CNF_RUN_PATH, CNF_SAVE_PATH, 'ZMATRIX')

    ZMA_RUN_PATH = ZMA_FS[0][-1].path(ZMA_LOCS)
    ZMA_SAVE_PATH = ZMA_FS[1][-1].path(ZMA_LOCS)
    ZMA_FS[0][-1].create(ZMA_LOCS)
    ZMA_FS[1][-1].create(ZMA_LOCS)
    print('zma run', ZMA_RUN_PATH)
    print('zma save', ZMA_SAVE_PATH)

    SCN_FS = build_fs(
        CNF_RUN_PATH, CNF_SAVE_PATH, 'SCAN',
        zma_locs=ZMA_LOCS)

    SCN_RUN_PATH = SCN_FS[0][-1].path(SCN_LOCS)
    SCN_SAVE_PATH = SCN_FS[1][-1].path(SCN_LOCS)
    SCN_FS[0][-1].create(SCN_LOCS)
    SCN_FS[1][-1].create(SCN_LOCS)
    print('scn run', SCN_RUN_PATH)
    print('scn save', SCN_SAVE_PATH)
