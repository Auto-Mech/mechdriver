""" Format utilities
"""

import automol.geom


CKIN_TRANS_HEADER_STR = """! THEORETICAL TRANSPORT PROPERTIES
!
! (1) Shape, index denotes atom (0), linear molec. (1), nonlinear molec. (2);
! (2) Epsilon, the Lennard-Jones well depth, in K;
! (3) Sigma, the Lennard-Jones collision diameter, in Angstrom;
! (4) Mu, total dipole moment, in Debye;
! (5) Alpha, mean static polarizability, in Angstrom^3; and
! (6) Z_rot, rotational relaxation collision number at 298 K."""


def max_rxn_length(rxn_param_dct):
    """ Gets the maximum length of the formatted reaction names in a mechanism
    """

    max_len = 0
    for rxn, params in rxn_param_dct.items():
        # Check if the pdep flag should be True or False
        forms = params.get_existing_forms()
        pdep = bool('troe' in forms or 'lind' in forms)
        rxn_name = format_rxn_name(rxn, pdep=pdep)
        max_len = max(len(rxn_name), max_len)

    return max_len

def max_spc_length(mech_spc_dct):
    """ Gets the maximum length of the species names in a mechanism

        (Will also work fine on any dict with species as the keys)
    """

    max_len = 0
    for spc in mech_spc_dct.keys():
        max_len = max(len(spc), max_len)

    return max_len


def name_column_length(names):
    """ Set the width of the name column
    """

    maxlen = 0
    for name in names:
        maxlen = max(maxlen, len(name))
    maxlen = max(maxlen, 9)
    names_len = str(maxlen + 3)

    return names_len


def format_rxn_name(rxn, pdep=False, rxnkeys = []):
    """ Receives a rxn and creates an appropriate string
        to be written in a Chemkin mech. Adds third body if applicable

        :param rxn: reaction names and third body
        :type rxn: tuple ((rct1, rct2), (prd1, prd2), (third_bod1,))
        :param rxnkeys: full set of reaction keys
        :type rxnkeys: list of rxns
        :return rxn_name: formatted reaction name for writing in the mech
        :rtype: str
    """
    # order names for consistent search
    rxnkeys_sorted = []
    for rxnk in rxnkeys:
        rxnkeys_sorted.append(
            tuple([tuple(sorted(rxnk[0])), tuple(sorted(rxnk[1])), rxnk[2]])
        )

    rcts = rxn[0]
    prds = rxn[1]
    thrbdy = rxn[2][0]

    # Convert to list if only one species
    if isinstance(rcts, str):
        rcts = [rcts]
    if isinstance(prds, str):
        prds = [prds]

    # Write the strings
    rct_str = ' + '.join(rcts)
    prd_str = ' + '.join(prds)

    # Add the +M or (+M) text if it is applicable
    if thrbdy is not None:
        rct_str += f' {thrbdy}'  # buffer space
        prd_str += f' {thrbdy}'
    elif pdep:  # if no tbody but pdep, add (+M)
        rct_str += ' (+M)'
        prd_str += ' (+M)'

    if len(prds) < 3:
        # check if backward reaction is also present
        rxnk_sorted = tuple([tuple(sorted(rxn[1])), tuple(sorted(rxn[0])), rxn[2]])
        if rxnk_sorted in rxnkeys_sorted:
            join_sign = ' => '  # if backward reaction found: write as irreversible
        else:
            join_sign = ' = '
    else:
        join_sign = ' => '

    rxn_name = f'{rct_str}{join_sign}{prd_str}'

    return rxn_name


def format_shape_idx(geo):
    """ Determine the shape index that signifies the overall
        molecular structure.

        :param geo: molecular geometry
        :type geo: automol.geom objects
        :rtype: int
    """

    if automol.geom.is_atom(geo):
        shape_idx = 0
    else:
        if automol.geom.is_linear(geo):
            shape_idx = 1
        else:
            shape_idx = 2

    return shape_idx
