"""
Executes the automation part of 1DMin
"""

import os
import sys
from routines.trans.lj import onedmin


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
