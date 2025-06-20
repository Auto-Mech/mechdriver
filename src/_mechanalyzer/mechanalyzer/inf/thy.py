"""
  Handles data objects
"""

import elstruct
from mechanalyzer import par


THY_PROPS = [
    par.THY.PROGRAM,
    par.THY.METHOD,
    par.THY.BASIS,
    par.THY.ORB_RESTRICT
]


# Constructors
def from_data(program, method, basis, orb_restrict):
    """ thy info data structure from necessary data
    """

    # assert to put in a data for a real species
    return (program, method, basis, orb_restrict)


def from_dct(dct):
    """ thy info data structure from
    """

    assert set(THY_PROPS) <= set(dct), (
        f'Properties {" ".join(THY_PROPS)} not in dict'
    )

    inf_obj = tuple()
    for prop in THY_PROPS:
        inf_obj += (dct[prop],)

    return inf_obj


def modify_orb_label(thy_info, spc_info):
    """ Convert the orbital restriction label to one for a species
        using the multiplicity.
    """

    # Grab first three elements of thy_info object
    mod_thy_info = thy_info[:3]

    # Grab mult from spc info
    mult = spc_info[2]

    # Modify the label denoting the orbital restriction
    orb_label = value(thy_info, par.THY.ORB_RESTRICT)
    mod_thy_info += (
        elstruct.util.set_orbital_restriction_label(orb_label, mult),)

    return mod_thy_info


# Getters
def value(inf_obj, val):
    """ obtain a value
    """

    assert val in THY_PROPS, (
        f'Desired value {val} not in thy info object'
    )

    return inf_obj[THY_PROPS.index(val)]


def combine(mod_inf_obj1, mod_inf_obj2):
    """ Combine two modified thy info objects.

        :param mod_inf_obj1:
    """

    prog1, method1, basis1, orb_lbl1 = mod_inf_obj1
    prog2, method2, basis2, orb_lbl2 = mod_inf_obj2

    assert prog1 == prog2, ('program from info objs must be the same')
    assert method1 == method2, ('method from info objs must be the same')
    assert basis1 == basis2, ('basis from info objs must be the same')

    if (orb_lbl1, orb_lbl2) == ('R', 'R'):
        comb_orb_lbl = 'R'
    elif (orb_lbl1, orb_lbl2) in (('R', 'U'), ('U', 'R')):
        comb_orb_lbl = 'U'

    return (prog1, method1, basis1, comb_orb_lbl)


# I/O
def string(geo_obj, sp_obj=None):
    """ Convert objects to a string representative common thy representation
        in literature:

        sp_method/sp_basis//geo_method/geo_basis
    """
    
    thy_str = f'{geo_obj[3]}{geo_obj[1]}/{geo_obj[2]}'
    if sp_obj is not None:
        thy_str = f'{sp_obj[3]}{sp_obj[1]}/{sp_obj[2]}//' + thy_str

    return thy_str
