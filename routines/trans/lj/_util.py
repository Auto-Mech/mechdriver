""" util
"""

def trans_chg_mult(tgt_info, bath_info):
    """ Determine what the combined charge and mult is for a run
    """

    chg = tgt_info[1] + bath_info[1]
    mult = max([tgt_info[2], bath_info[2]])

    return chg, mult
