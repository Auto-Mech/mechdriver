"""
utility functions
"""

import automol
import autofile
from mechlib import structure


# Functions used by build to set pieces of information
def ini_elec_levels(spc_dct_i, spc_info):
    """ get initial elec levels
    """
    if 'elec_levels' in spc_dct_i:
        elec_levels = spc_dct_i['elec_levels']
    else:
        elec_levels = [[0., spc_info[2]]]

    return elec_levels


def combine_elec_levels(elec_levels_i, elec_levels_j):
    """ Put two elec levels together for two species
    """

    # Combine the energy levels
    init_elec_levels = []
    for _, elec_level_i in enumerate(elec_levels_i):
        for _, elec_level_j in enumerate(elec_levels_j):
            init_elec_levels.append(
                [elec_level_i[0]+elec_level_j[0],
                 elec_level_i[1]*elec_level_j[1]])

    # See if any levels repeat and thus need to be added together
    elec_levels = []
    for level in init_elec_levels:
        # Put level in in final list
        if level not in elec_levels:
            elec_levels.append(level)
        # Add the level to the one in the list
        else:
            idx = elec_levels.index(level)
            elec_levels[idx][1] += level[1]

    return elec_levels


def atom_mass(spc_dct_i):
    """ write the atom string
    """
    geo = automol.inchi.geometry(spc_dct_i['inchi'])
    return automol.geom.total_mass(geo)


def get_stoich(geom_i, geom_j):
    """ get the overall combined stoichiometry
    """

    form_i = automol.geom.formula(geom_i)
    form_j = automol.geom.formula(geom_j)
    form = automol.formula.join(form_i, form_j)
    stoich = ''
    for key, val in form.items():
        stoich += key + str(val)

    return stoich
