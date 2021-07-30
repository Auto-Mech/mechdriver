"""
Collates information for, and writes MESS files for rate calculations
"""

import mess_io
from mechroutines.models import build, blocks


def make_messpf_str(temps, spc_dct, spc_name, spc_locs,
                    pes_mod_dct_i, spc_mod_dct_i,
                    run_prefix, save_prefix):
    """ Reads and processes all information in the save filesys for
        a given species that are required for MESSPF calculations,
        as specified by the model dictionaries built from user input.

        :param temps: temperatures to calculate partition functions (in K)
        :type temps: tuple(float)
        :param spc_dct:
        :type spc_dct:
        :param spc_name: mechanism name of species to write MESSPF input for
        :type spc_name: str
        :param pes_mod_dct_i: keyword dict of specific kin model
        :type pes_mod_dct_i: dict[]
        :param spc_mod_dct_i: keyword dict of specific species model
        :type spc_mod_dct_i: dict[]
        :param run_prefix: root-path to the run-filesystem
        :type run_prefix: str
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
        :rtype: str
    """

    # Read the filesystem for the information
    inf_dct, _ = build.read_spc_data(
        spc_dct, spc_name,
        pes_mod_dct_i, spc_mod_dct_i,
        run_prefix, save_prefix, {}, calc_chn_ene=False,
        spc_locs=spc_locs)

    # Write the header string for the MESS input file
    globkey_str = mess_io.writer.global_pf_input(
        temperatures=temps,
        rel_temp_inc=0.001,
        atom_dist_min=0.6
    )

    # Write the species data block for the MESS input file
    mess_writer = getattr(blocks, inf_dct['writer'])
    mess_block, dat_str_dct = mess_writer(inf_dct)
    if inf_dct['writer'] == 'tau_block':
        zero_energy = None
    else:
        zero_energy = inf_dct['zpe_chnlvl']

    spc_str = mess_io.writer.species(
        spc_label=spc_name,
        spc_data=mess_block,
        zero_ene=zero_energy
    )

    # Combine the strings together to create full MESS input file string
    mess_inp_str = mess_io.writer.messpf_inp_str(globkey_str, spc_str)

    return mess_inp_str, dat_str_dct
