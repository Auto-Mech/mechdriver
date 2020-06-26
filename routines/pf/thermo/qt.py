"""
Collates information for, and writes MESS files for rate calculations
"""

import importlib
import mess_io
from routines.pf.messf import blocks
from routines.pf.messf import models


BLOCK_MODULE = importlib.import_module('routines.pf.messf.blocks')


# Input string writer
def make_messpf_str(globkey_str, spc_str):
    """ Combine various MESS strings together to combined MESSPF
    """
    return '\n'.join([globkey_str, spc_str])


def get_pf_header(temps):
    """ prepare partition function header string
    """

    global_pf_str = mess_io.writer.global_pf(
        temperatures=temps,
        rel_temp_inc=0.001,
        atom_dist_min=0.6
    )

    return global_pf_str


def make_spc_mess_str(spc_dct_i, spc_name,
                      chn_pf_models, chn_pf_levels,
                      run_prefix, save_prefix):
    """ Write the MESS input file strings
    """

    # Read the filesystem for the information
    inf_dct = models.read_spc_data(
        spc_dct_i, spc_name,
        chn_pf_models, chn_pf_levels,
        run_prefix, save_prefix)

    # Write the mol block
    mess_writer = getattr(BLOCK_MODULE, inf_dct['writer'])
    mess_block = mess_writer(inf_dct)

    # Write the mess string
    spc_str = mess_io.writer.species(
        spc_label='TMP',
        spc_data=mess_block,
        zero_energy=inf_dct['zpe_chnlvl']
    )

    # dat string stuff
    dat_str_dct = {}

    return spc_str, dat_str_dct
