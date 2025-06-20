""" Wrapper functions to extract energy distributions from mess output
"""

from mess_io import reader

def ped_info(ped_inp_str, ped_ped_str, ped_ke_out_str):
    """ file strings
        wellreac: str if you want PEDs of well->bimol 
    """
    # Read: ped.inp:

    ped_spc, _ = reader.ped.ped_names(ped_inp_str)  # can supply
    energy_dct, _, _, _ = reader.pes(ped_inp_str)

    # Read ped.out file for product energy distributions
    ped_dct = reader.ped.get_ped(
        ped_ped_str, energy_dct, sp_labels='auto')

    # Read ke_ped.out file for energy density of each fragment
    dos_df = reader.rates.dos_rovib(ped_ke_out_str, sp_labels='auto')
    return ped_spc, ped_dct, dos_df, energy_dct


def hot_info(hot_inp_str, hot_log_str):
    """
    Extract required info for hotenergies
    """
    spc_blocks_hoten = reader.get_species(hot_inp_str)
    hot_frag_dct = reader.dct_species_fragments(spc_blocks_hoten)
    hot_spc_en = reader.hoten.get_hot_species(hot_inp_str)

    hoten_dct = reader.hoten.extract_hot_branching(
        hot_log_str, hot_spc_en, list(spc_blocks_hoten.keys()), sp_labels='auto')

    fne_bf = reader.hoten.extract_fne(hot_log_str, sp_labels='auto')

    return hot_frag_dct, hot_spc_en, hoten_dct, fne_bf

