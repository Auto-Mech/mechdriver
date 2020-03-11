"""
Calculate an effective alpha parameter using certain rules
"""

import automol


def jasper_method(geo, bath_model=None):
    """ Calculate the alpha param using the method in Ahren's paper
    """

    # Conver the geo to a graph
    gra = autmol.geom.graph(geo)

    # Determine the model for th target and the bath gas
    determine_model(gra)

    # Calculate Zalpa(Neff)
    n_eff = calc_n_eff(gra, automol.geom.symbols(geo))
    z_alpha_n_eff = calc_z_alpha(tgt_model, bath_model, temp, n_eff):

    # Normalize the Zalpha term with Z(N)
    n_heavy = automol.geom.heavy_count(geo)
    z_n_heavy = read_collision_freq(tgt_model, bath_model, n_heavy)

    return z_alpha_n_eff / z_n_heavy


def calc_n_eff(gra,
               c_pp_ps_ss=1, c_pt_st=(2/3), c_pq_sq=(1/3),    
               c_tt_tq_qq=0, c_co_oo=(1/3), c_ss_ring=(1/2)):
    """ Calculate an effective N parameter using the given parametrization
    """
    # Convert the geometry to a graph

    # Count the rotors
    nlst = rotor_counts(gra)
    n_pp, n_ps, n_pt, n_pq, n_ss, n_st, n_sq, n_tt, n_tq, n_qq, n_co, n_oo = nlst

    # Count the rings
    n_ss_ring, n_rings = ring_counts(gra)

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


def rotor_counts(geo, symbs):
    """ Count up various types of rotors 
    """

    # Initialize the rotor counts
    n_pp, n_ps, n_pt, n_pq = 0, 0, 0, 0
    n_ss, n_st, n_sq = 0, 0, 0
    n_tt, n_tq = 0, 0
    n_qq = 0
    n_co, n_oo = 0, 0

    # Loop over the bonds and count the number of atoms
    neighbors = automol.grpah.atom_neighbor_keys(gra)
    for bnd in automol.graph.bond_keys(gra):
        key1, key2 = bnd
        symb1, symb2 = symbols[key1], symbols[key2]
        if (symb1 == 'C' and symb2 == 'O') or (symb1 == 'O' and symb2 == 'C'):
            n_co += 1
        elif symb1 == 'O' and symb2 == 'O':
            n_oo += 1
        elif symb1 == 'C' and symb2 == 'C':
            # Figure out which neighbors are carbons and count the number
            atom1_neighbors = neighbors[key1]
            numc1 = 0
            for neighbor1 in atom1_neightbors:
                if symbols[neighbor1] == 'C':
                    numc1 += 1
            atom2_neighbors = neighbors[key2]
            numc2 = 0
            for neighbor2 in atom2_neightbors:
                if symbols[neighbor2] == 'C':
                    numc2 += 1
            # Determine appropriate term to increment
            if numc1 == 1 and numc2 == 1:
                npp += 1 
            elif (numc1 == 1 and numc2 == 2) or (numc1 == 2 and numc2 == 1):
                nps += 1
            elif (numc1 == 1 and numc2 == 3) or (numc1 == 3 and numc2 == 1):
                npt += 1
            elif (numc1 == 1 and numc2 == 4) or (numc1 == 4 and numc2 == 1):
                npq += 1
            elif numc1 == 2 and numc2 == 2:
                nss += 1 
            elif (numc1 == 2 and numc2 == 3) or (numc1 == 3 and numc2 == 2):
                nst += 1
            elif (numc1 == 2 and numc2 == 4) or (numc1 == 4 and numc2 == 2):
                nsq += 1
            elif numc1 == 3 and numc2 == 3:
                ntt += 1 
            elif (numc1 == 3 and numc2 == 4) or (numc1 == 4 and numc2 == 3):
                ntq += 1
            elif numc1 == 4 and numc2 == 4:
                nqq += 1 

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


def moiety_idx(geo):
    """ Determine the chemical moiety
    """
    pass


def calc_sig_eps(tgt_model, bath_model, n_cnt): 
    """ Returns in angstrom and cm-1
    """
    assert tgt_model in ('hydrocarbon', 'alcohol', 'peroxide')
    assert bath_model in ('Ar', 'He', 'H2', 'N2')

    if tgt_model == 'hydrocarbon':
        if bath_model == 'Ar':
            sigma = 3.40 * n_cnt**(0.18)
            epsilon = 113.0 * n_cnt**(0.31)
        elif bath_model == 'He':
            sigma = 3.33 * n_cnt**(0.17)
            epsilon = 21.3 * n_cnt**(0.31)
        elif bath_model == 'H2':
            sigma = 3.15 * n_cnt**(0.18)
            epsilon = 75.0 * n_cnt**(0.30)
        elif bath_model == 'N2':
            sigma = 3.68 * n_cnt**(0.16)
            epsilon = 100.0 * n_cnt**(0.25)
    elif tgt_model == 'alcohol':
        if bath_model == 'Ar':
            sigma = 3.05 * n_cnt**(0.20)
            epsilon = 150.0 * n_cnt**(0.29)
        if bath_model == 'He':
            sigma = 2.90 * n_cnt**(0.21)
            epsilon = 22.0 * n_cnt**(0.28)
    elif tgt_model == 'peroxide':
        if bath_model == 'Ar':
            sigma = 3.05 * n_cnt**(0.20)
            epsilon = 110.0 * n_cnt**(0.39)
        elif bath_model == 'He':
            sigma = 2.90 * n_cnt**(0.21)
            epsilon = 10.0 * n_cnt**(0.75)

    return sigma, epsilon


def calc_z_alpha(tgt_model, bath_model, temp, n_eff):
    """ Calculate the [Z*alpha](N_eff)
    """
    assert tgt_model in ('alkane', 'alcohol', 'peroxide')
    assert bath_model in ('Ar', 'He', 'H2', 'N2')
    assert float(300) <= temp <= float(2000) 

    # Set the fitting coefficients based on the model chosen
    fit_coeff_dct 
    if tgt_model == 'alkane':
        if bath_model == 'Ar':
            cff = {
                300: [-2.984767, 52.534594, -2.948084, 0.054986],
                1000: [31.934663, 125.075803, -10.208643, 0.273835],
                2000: [110.533659, 228.021958, -19.131921, 0.457981]}
        elif bath_model == 'He':
        elif bath_model == 'H2':
        elif bath_model == 'N2':
    elif tgt_model == 'alcohol':
        if bath_model == 'Ar':
        elif bath_model == 'He':
    elif tgt_model == 'peroxide':
        if bath_model == 'Ar':
        elif bath_model == 'He':
            sigma = 2.90 * n_cnt**(0.21)
            epsilon = 10.0 * n_cnt**(0.75)

    # Calculate the three alpha terms
    alpha_300 = (cff[300][0] * n_eff**(3) + cff[300][1] * n_eff**(2) + 
                 cff[300][2] * n_eff**(1) + cff[300][3])
    alpha_1000 = (cff[1000][0] * n_eff**(3) + cff[1000][1] * n_eff**(2) + 
                  cff[1000][2] * n_eff**(1) + cff[1000][3])
    alpha_2000 = (cff[2000][0] * n_eff**(3) + cff[2000][1] * n_eff**(2) + 
                 cff[2000][2] * n_eff**(1) + cff[2000][3])

    # Perform a linear interpolation

    return sigma, epsilon



