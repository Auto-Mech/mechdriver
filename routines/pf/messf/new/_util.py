"""
utility functions
"""

from qcelemental import periodictable as ptab
import automol


def ini_elec_levels(spc_dct, spc_info):
    """ get initial elec levels
    """
    if 'elec_levs' in spc_dct:
        elec_levels = spc_dct['elec_levs']
    else:
        elec_levels = [[0., spc_info[2]]]

    return elec_levels


def combine_elec_levels(spc_dct_i, spc_dct_j):
    """ Put two elec levels together for two species
    """

    if 'elec_levs' in spc_dct_i:
        elec_levels_i = spc_dct_i['elec_levs']
    else:
        elec_levels_i = [[0., spc_dct_i['mul']]]
    if 'elec_levs' in spc_dct_j:
        elec_levels_j = spc_dct_j['elec_levs']
    else:
        elec_levels_j = [[0., spc_dct_j['mul']]]

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


def get_bnd_keys(spc_dct, saddle):
    """ get bond broken and formed keys for a transition state
    """
    if saddle:
        frm_bnd_key = spc_dct['frm_bnd_key']
        brk_bnd_key = spc_dct['brk_bnd_key']
    else:
        frm_bnd_key = []
        brk_bnd_key = []

    return frm_bnd_key, brk_bnd_key


def is_atom(spc_dct_i):
    """ Check if species is an atom
    """
    geo = automol.inchi.geom(spc_dct_i['ich'])
    return automol.geom.is_atom(geo)


def atom_mass(spc_dct_i):
    """ write the atom string
    """
    geo = automol.inchi.geom(spc_dct_i['ich'])
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
