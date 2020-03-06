"""
Calculate an effective alpha parameter using certain rules
"""

import automol


def jasper_method():
    """ Calculate the alpha param using the method in Ahren's paper
    """
    pass


def calc_n_eff(geo,
               c_pp_ps_ss=1, c_pt_st=(2/3), c_pq_sq=(1/3),    
               c_tt_tq_qq=0, c_co_oo=(1/3), c_ss_ring=(1/2)):
    """ Calculate an effective N parameter using the given parametrization
    """

    # Count the rotors
    nlst = rotor_counts(geo)
    n_pp, n_ps, n_pt, n_pq, n_ss, n_st, n_sq, n_tt, n_tq, n_qq, n_co, n_oo = nlst

    # Count the rings
    n_ss_ring, n_rings = ring_counts(geo)

    # Use the rotor counts and the coefficients to calculate Neff
    n_eff = 1 + (
        c_pp_ps_ss * (n_pp + n_ps + n_ss) +
        c_pt_st * (n_pt + n_st) +
        c_pq_sq * (n_pq + n_sq) +
        c_tt_tq_qq * (n_tt + n_tq + n_qq) +
        c_co_oo * (n_co + n_oo) +
        c_ss_ring * n_ss_ring - n_rings
    )

    return n_eff


def rotor_counts(geo):
    """ Count up various types of rotors 
    """

    # Initialize the rotor counts
    n_pp, n_ps, n_pt, n_pq = 0, 0, 0, 0
    n_ss, n_st, n_sq = 0, 0, 0
    n_tt, n_tq = 0, 0
    n_qq = 0
    n_co, n_oo = 0, 0

    # Compile counts into a tuple
    n_lst = (
        n_pp, n_ps, n_pt, n_pq,
        n_ss, n_st, n_sq,
        n_tt, n_tq,
        n_qq,
        n_co, n_oo
    )

    return n_lst


def ring_counts(geo):
    """ Count up various types of rotors 
    """

    return n_ss_ring, n_rings
