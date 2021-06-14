"""
Random functions that are needed in drivers and routines
"""


def make_rxn_str(rlst, prepend=''):
    """ convert list to string
    """
    return prepend + '+'.join(rlst)
