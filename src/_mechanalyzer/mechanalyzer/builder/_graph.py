""" Handle the graph representation of a PES
"""

import networkx


def pes_graphs_dct(pes_dct):
    """ Convert PES dcts into graph
    """
    return {pes_inf: _graph(rxns) for pes_inf, rxns in pes_dct.items()}


def _graph(sub_pes_rxns):
    """ Convert part of a sub_PES dct into a graph
    """

    # Build the connections lst for edge and reagents lst for vertices
    conn_lst, rgts_lst = (), ()
    for rxn in sub_pes_rxns:
        _, (reacs, prods) = rxn

        # Flip order of reagents if reverse already in reagent list
        if reacs[::-1] in rgts_lst:
            reacs = reacs[::-1]
        if prods[::-1] in rgts_lst:
            prods = prods[::-1]

        # Now add to reagent list, if needed
        if reacs not in rgts_lst:
            rgts_lst += (reacs,)
        if prods not in rgts_lst:
            rgts_lst += (prods,)

        conn_lst += ((reacs, prods),)

    # Now loop over lists to join into strings
    rgts_lst = tuple('+'.join(rgts) for rgts in rgts_lst)
    conn_lst = tuple(('+'.join(reacs), '+'.join(prods))
                     for (reacs, prods) in conn_lst)

    nxg = networkx.Graph()
    nxg.add_nodes_from(rgts_lst)
    nxg.add_edges_from(conn_lst)

    return nxg
