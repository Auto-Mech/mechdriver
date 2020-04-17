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


def is_atom(har_min_cnf_locs, har_cnf_save_fs):
    """ Check if species is an atom
    """
    if har_min_cnf_locs is not None:
        har_geo = har_cnf_save_fs[-1].file.geometry.read(har_min_cnf_locs)
        # print('This is an atom')
    return automol.geom.is_atom(har_geo)


def atom_mass(har_min_cnf_locs, har_cnf_save_fs):
    """ write the atom string
    """
    har_geo = har_cnf_save_fs[-1].file.geometry.read(har_min_cnf_locs)
    return ptab.to_mass(har_geo[0][0])


def get_stoich(harm_min_cnf_locs_i, harm_min_cnf_locs_j,
               harm_cnf_save_fs_i, harm_cnf_save_fs_j):
    """ get the overall combined stoichiometry
    """
    if harm_min_cnf_locs_i is not None:
        harm_geo_i = harm_cnf_save_fs_i[-1].file.geometry.read(
            harm_min_cnf_locs_i)
        if harm_min_cnf_locs_j is not None:
            harm_geo_j = harm_cnf_save_fs_j[-1].file.geometry.read(
                harm_min_cnf_locs_j)

    form_i = automol.geom.formula(harm_geo_i)
    form_j = automol.geom.formula(harm_geo_j)
    form = automol.formula.join(form_i, form_j)
    stoich = ''
    for key, val in form.items():
        stoich += key + str(val)

    return stoich
