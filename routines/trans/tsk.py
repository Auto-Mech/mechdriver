"""
Executes the automation part of 1DMin
"""

from routines.trans.lj import onedmin
from routines.trans._write import collate_properties


def run(tsk, spc_queue, spc_name, bath_name,
        spc_dct, pf_levels,
        thy_dct, etrans_keyword_dct,
        run_prefix, save_prefix):
    """ Run the task
    """

    if tsk == 'onedmin':
        onedmin(spc_queue, spc_name, bath_name,
                spc_dct, pf_levels,
                thy_dct, etrans_keyword_dct,
                run_prefix, save_prefix)
    elif tsk == 'write_transport':
        collate_properties(spc_queue, spc_name, bath_name,
                           spc_dct, pf_levels,
                           thy_dct, etrans_keyword_dct,
                           run_prefix, save_prefix)
