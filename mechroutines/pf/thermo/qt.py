"""
Collates information for, and writes MESS files for rate calculations
"""

import importlib
import mess_io
from mechroutines.pf.models import build

BLOCK_MODULE = importlib.import_module('mechroutines.pf.models.blocks')


# Input string writer
def make_messpf_str(temps, spc_dct, spc_name,
                    chn_pf_models, chn_pf_levels,
                    run_prefix, save_prefix):
    """ Combine various MESS strings together to combined MESSPF
    """

    # Write the header strings for the MESS input file
    globkey_str = make_pf_header(temps)

    # Write the molecules species strings
    spc_str, _ = make_spc_mess_str(
        spc_dct, spc_name,
        chn_pf_models, chn_pf_levels,
        run_prefix, save_prefix)

    # Combine the strings together
    mess_inp_str = mess_io.writer.messpf_inp_str(globkey_str, spc_str)

    return mess_inp_str


def make_pf_header(temps):
    """ prepare partition function header string
    """

    global_pf_str = mess_io.writer.global_pf(
        temperatures=temps,
        rel_temp_inc=0.001,
        atom_dist_min=0.6
    )

    return global_pf_str


def make_spc_mess_str(spc_dct, spc_name,
                      chn_pf_models, chn_pf_levels,
                      run_prefix, save_prefix):
    """ Write the MESS input file strings
    """
    # Read the filesystem for the information
    inf_dct, _ = build.read_spc_data(
        spc_dct, spc_name,
        chn_pf_models, chn_pf_levels,
        run_prefix, save_prefix, {}, calc_chn_ene=False)

    # Write the mol block
    mess_writer = getattr(BLOCK_MODULE, inf_dct['writer'])
    mess_block, dat_str_dct = mess_writer(inf_dct)
    if inf_dct['writer'] == 'tau_block':
        zero_energy = None
    else:
        zero_energy = inf_dct['zpe_chnlvl']

    # Write the mess string
    spc_str = mess_io.writer.species(
        spc_label=spc_name,
        spc_data=mess_block,
        zero_ene=zero_energy
    )

    return spc_str, dat_str_dct
