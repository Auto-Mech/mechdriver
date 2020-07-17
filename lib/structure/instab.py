"""
 Library to deal unstable species
"""

# init: just check geom
# conf: check ratio of confs
# hr:   check geom along scan,

import automol


def get_gras(conn_geo, disconn_geo):
    """ get graphs
    """
    # maybe take zmas
    conn_gra = automol.geom.graph(conn_geo)
    disconn_gras = automol.graph.connected_components(
        automol.geom.graph(disconn_geo))

    return conn_gra, disconn_gras


# Write the instability files if needed
def write_instab(conn_gra, disconn_gras, save_fs, locs):
    """ write the instability files
    """

    # Obtain the transformation and reactants graph
    frm_bnd_keys, rcts_gra = _instab_info(conn_gra, disconn_gras)
    tra = (frozenset({frm_bnd_keys}),
           frozenset({frozenset({})}))

    # Write the files into the filesystem
    save_fs[-1].file.transformation.write(tra, locs)
    save_fs[-1].file.reactant_graph.write(rcts_gra, locs)


def _instab_info(conn_gra, disconn_gras):
    """ Obtain instability info
    """

    # Get the zma for the connected graph
    prd_zmas = [automol.geom.zmatrix(automol.graph.geometry(conn_gra))]

    # Get the zmas used for the identification
    rct_zmas = [automol.geom.zmatrix(automol.graph.geometry(gra))
                for gra in disconn_gras]

    # Get the keys
    ret = automol.zmatrix.ts.addition(rct_zmas, prd_zmas, ())
    _, _, frm_bnd_keys, _, rcts_gra = ret

    return frm_bnd_keys, rcts_gra


# DISCONN_GEO_STR = """
# C       -2.2947657629     -0.0356863088      0.0075513595
# O       -0.8974168549      0.0470354136      0.0255844602
# H       -2.6737732983     -0.7042322833      0.8113409886
# H       -2.6737771655     -0.3089516561     -1.0015842288
# O        0.3279198002     -2.2232984128     -0.2399706150
# H        1.3153969176     -2.1475164906     -0.2234378009"""
#
# CONN_GEO_STR = """
# C       -1.9332015404     -0.2451995079     -0.0114253692
# O       -0.5408704259     -0.3907355653     -0.0092590998
# H       -2.4063957726     -0.8260210019      0.8106129504
# H       -2.3645483729     -0.4721480489     -1.0108650702
# O       -0.3184096603     -1.6734963150     -0.2204733954
# H        0.6670094082     -1.7650492990     -0.1791058522"""
#
# DISCONN_GEO = automol.geom.from_string(DISCONN_GEO_STR)
# CONN_GEO = automol.geom.from_string(CONN_GEO_STR)

# if geo_disconnected(geo):
#     conn_gra, disconn_gras = get_gras(conn_geo, disconn_geo)
#     write_instab(conn_gra, disconn_gras, save_fs, locs):
#     _instab_info(CONN_GRA, DISCONN_GRA)
