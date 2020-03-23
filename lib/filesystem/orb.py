""" new functions for filesystem stuff
"""


def mod_orb_restrict(spc_info, thy_info):
    """ append to the theory level the orb restricted stuff
    """
    orb_restr = set_orbital_restriction_label(spc_info, thy_info)
    thy_level = thy_info[0:3]
    thy_level.append(orb_restr)

    return thy_level


def set_orbital_restriction_label(spc_info, thy_level):
    """ orbital restriction logical
    """
    mul = spc_info[2]
    if thy_level[3] == 'RR':
        orb_restr = 'R'
    elif thy_level[3] == 'UU':
        orb_restr = 'U'
    elif thy_level[3] == 'RU':
        if mul == 1:
            orb_restr = 'R'
        else:
            orb_restr = 'U'
    return orb_restr


def mod_orbital_restrict(spc_info, thy_info):
    """ append to the theory level the orb restricted stuff
    """
    orb_restr = orbital_restriction(
        spc_info, thy_info)
    thy_level = thy_info[0:3]
    thy_level.append(orb_restr)
    return thy_level
