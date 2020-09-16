"""
Executes the automation part of 1DMin
"""

from routines.trans._routines import lj
from routines.trans._routines import build


def run(tsk, spc_queue, spc_name, bath_name,
        spc_dct, pf_levels,
        thy_dct, etrans_keyword_dct,
        run_prefix, save_prefix):
    """ Run the task
    """

    if tsk == 'onedmin':
        lj.onedmin(spc_queue, spc_name, bath_name,
                   spc_dct, pf_levels,
                   thy_dct, etrans_keyword_dct,
                   run_prefix, save_prefix)
    elif tsk == 'write_transport':
        build.collate_properties(spc_queue, spc_name, bath_name,
                                 spc_dct, thy_dct, etrans_keyword_dct,
                                 save_prefix)
