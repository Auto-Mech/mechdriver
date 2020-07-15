"""
 Library to deal unstable species
"""

# init: just check geom
# conf: check ratio of confs
# hr:   check geom along scan, 


# Instability checkers
def identify_instability(scn_fs, scn_locs)
    """ Examine a geometry to see if it has disconnected into two species
    """
    
    # ret = es_runner.read_job(job=job, run_fs=run_fs)
    # if ret is not None:

    # Read the geometry and convert it to a graph anf componentns
    geo = scn_fs[-1].file.geometry.read(scn_locs)
    gra = automol.geom.graph(geo))
    gra_comp = automol.graph.connected_components(gra)
  
    return bool(len(gra_comp) == 1)


# Write the instability files if needed
def write_instab(gra, save_fs, locs):
    """ write the instability files
    """
    
    frm_bnd_keys, rcts_gra = _instab_info(gra)
    save_fs[-1].file.instab.write(instab, locs)
    save_fs[-1].file.instab_graph.write(rcts_gra, locs)


def _instab_info(disconn_gra):
    """ Obtain instability info
    """
   
    # Get the zmas used for the identification
    prd_zmas = []
    for gra in disconn_gra:
        prd_zmas.append(
            automol.geom.zmatrix(automol.graph.geometry(gra))
        )

    # Get the keys    
    ret = automol.zmatrix.ts.addition(rct_zmas, prd_zmas, rct_tors_names)
    _, _, frm_bnd_keys, _, rcts_gra = ret


    return frm_bnd_keys, rcts_gra
