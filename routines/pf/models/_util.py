"""
utility functions
"""

import automol
import autofile
from lib import structure


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


def set_dist_names(spc_dct_i, saddle):
    """ Set various things needed for TSs
    """

    dist_names = []
    if saddle:
        dist_names = []
        mig = 'migration' in spc_dct_i['class']
        elm = 'elimination' in spc_dct_i['class']
        if mig or elm:
            dist_names.append(spc_dct_i['dist_info'][0])
            dist_names.append(spc_dct_i['dist_info'][3])

    return dist_names


def get_bnd_keys(cnf_fs, cnf_locs, saddle, zma_locs=(0,)):
    """ get bond broken and formed keys for a transition state
    """
    if not saddle:
        frm_bnd_keys = []
        brk_bnd_keys = []
    else:
        frm_bnd_keys, brk_bnd_keys = structure.ts.rxn_bnd_keys(
            cnf_fs, cnf_locs, zma_locs=zma_locs)

    return frm_bnd_keys, brk_bnd_keys


def get_bnd_keys2(ts_path, saddle, zma_locs=(0,)):
    """ get bond broken and formed keys for a transition state
    """
    if not saddle:
        frm_bnd_keys = []
        brk_bnd_keys = []
    else:
        frm_bnd_keys, brk_bnd_keys = structure.ts.rxn_bnd_keys2(
            ts_path, zma_locs=zma_locs)

    return frm_bnd_keys, brk_bnd_keys



def get_rxn_coord_name(ts_path, frm_bnd_keys, sadpt=False, zma_locs=(0,)):
    """ get the name of the reaction coordinate
    """

    if sadpt:
        zma_fs = autofile.fs.zmatrix(ts_path)
        zma = zma_fs[-1].file.zmatrix.read(zma_locs)
        frm_name = automol.zmatrix.bond_key_from_idxs(
            zma, frm_bnd_keys)
    else:
        zma_fs = autofile.fs.zmatrix(ts_path)
        vma = zma_fs[-1].file.vmatrix.read(zma_locs)
        frm_name = automol.vmatrix.bond_key_from_idxs(
            vma, frm_bnd_keys)

    return frm_name


def set_ts_bnd(spc_dct_i, saddle):
    """ get ts bnd
    """
    if saddle:
        ts_bnd = spc_dct_i['ts_bnd']
    else:
        ts_bnd = None

    return ts_bnd


def set_rxn_class(spc_dct_i, saddle):
    """ get bond broken and formed keys for a transition state
    """
    if saddle:
        rxn_class = spc_dct_i['class']
    else:
        rxn_class = None

    return rxn_class


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
