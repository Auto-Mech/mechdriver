"""
utility functions
"""

import automol


def ini_elec_levels(spc_dct_i, spc_info):
    """ get initial elec levels
    """
    if 'elec_levels' in spc_dct_i:
        elec_levels = spc_dct_i['elec_levels']
    else:
        elec_levels = [[0., spc_info[2]]]

    return elec_levels


def atom_mass(spc_dct_i):
    """ write the atom string
    """
    geo = automol.chi.geometry(spc_dct_i['inchi'])
    return automol.geom.total_mass(geo)
