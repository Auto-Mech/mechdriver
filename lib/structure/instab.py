"""
 Library to deal unstable species
"""

# init: just check geom
# conf: check ratio of confs
# hr:   check geom along scan,

import automol
import autofile
from lib import filesys
from automol.zmatrix._bimol_ts import addition


# Write the instability files
def write_instab(conn_zma, disconn_zma, save_fs, locs):
    """ write the instability files
    """

    # Obtain the transformation and reactants graph
    frm_bnd_keys, rcts_gra = _instab_info(conn_zma, disconn_zma)
    tra = (frozenset({frm_bnd_keys}),
           frozenset({frozenset({})}))

    # Write the files into the filesystem
    save_fs[-1].file.transformation.write(tra, locs)
    save_fs[-1].file.reactant_graph.write(rcts_gra, locs)


def _instab_info(conn_zma, disconn_zma):
    """ Obtain instability info
    """

    # Get the zma for the connected graph
    prd_zmas = [conn_zma]

    # Get the zmas used for the identification
    rct_zmas = _disconnected_zmas(disconn_zma)

    # print('prd')
    # for zma in prd_zmas:
    #     print(automol.zmatrix.string(zma))
    #     print()
    # print('\nrct')
    # for zma in rct_zmas:
    #     print(automol.zmatrix.string(zma))
    #     print()

    # Get the keys
    ret = addition(rct_zmas, prd_zmas, ())
    _, _, frm_bnd_keys, _, rcts_gra = ret

    return frm_bnd_keys, rcts_gra


def _disconnected_zmas(disconn_zma):
    """ get graphs
    """

    # Convert to disconnected component graph
    disconn_geo = automol.zmatrix.geometry(disconn_zma)
    disconn_gras = automol.graph.connected_components(
        automol.geom.graph(disconn_geo))

    print('disconn_gras', disconn_gras)

    # Get the zmas
    disconn_zmas = [automol.geom.zmatrix(automol.graph.geometry(gra))
                    for gra in disconn_gras]

    print('disconn_zmas', disconn_gras)

    return disconn_zmas


# Unstable check
def check_unstable_species(spc_dct, spc_name,
                           thy_info, ini_thy_info, save_prefix):
    """ see if a species and unstable and handle task management
    """

    if 'ts' not in spc_name:

        # Build filesystem
        spc_info = filesys.inf.get_spc_info(spc_dct[spc_name])
        _ = filesys.inf.modify_orb_restrict(spc_info, thy_info)
        mod_ini_thy_info = filesys.inf.modify_orb_restrict(
            spc_info, ini_thy_info)
        ini_thy_save_fs, _ = filesys.build.spc_thy_fs_from_root(
            save_prefix, spc_info, mod_ini_thy_info)

        # Check if the instability files exist
        thy_locs = mod_ini_thy_info[1:4]
        if (ini_thy_save_fs[-1].file.transformation.exists(thy_locs) and
                ini_thy_save_fs[-1].file.reactant_graph.exists(thy_locs)):
            stable = False
            thy_path = ini_thy_save_fs[-1].path(thy_locs)
            print('\nFound instability files for species {}'.format(spc_name),
                  'at path:\n{}'.format(thy_path))
        else:
            print('no inst', ini_thy_save_fs[-1].path(thy_locs))
            stable = True

    else:
        stable = True

    return stable


# Handle reaction lst
def break_all_unstable(rxn_lst, spc_dct, spc_model_dct, thy_dct, save_prefix):
    """ Loop over the reaction list and break up the unstable species
    """

    new_rxn_lst = []
    for rxn in rxn_lst:

        # Initialize dct
        new_rxn = {}

        # Get theory
        spc_model = rxn['model'][1]
        geo_model = spc_model_dct[spc_model]['es']['geo']
        ini_thy_info = filesys.inf.get_es_info(geo_model, thy_dct)

        # Asses the reactants for unstable species
        new_rxn['reacs'] = []
        for rct in rxn['reacs']:
            rct_stable = check_unstable_species(
                spc_dct, rct, ini_thy_info, ini_thy_info, save_prefix)
            if rct_stable:
                new_rct = rct
                new_rxn['reacs'].append(new_rct)
            else:
                new_rct = split_species(spc_dct, rct,
                                        ini_thy_info, save_prefix)
                new_rxn['reacs'].extend(new_rct)

        # Assess the products for unstable species
        new_rxn['prods'] = []
        for prd in rxn['prods']:
            prd_stable = check_unstable_species(
                spc_dct, prd, ini_thy_info, ini_thy_info, save_prefix)
            if prd_stable:
                new_prd = prd
                new_rxn['prods'].append(new_prd)
            else:
                new_prd = split_species(spc_dct, prd,
                                        ini_thy_info, save_prefix)
                new_rxn['prods'].extend(new_prd)

        # Add a print statement showing reaction being rewritten

        # Build rxn dct
        new_rxn.update(
            {'model': rxn['model'],
             'chn_idx': rxn['chn_idx'],
             'species': new_rxn['reacs']+new_rxn['prods']})

        # Flip the reaction if the reactants are unstable?

        # Append to list
        new_rxn_lst.append(new_rxn)

    return new_rxn_lst


def split_species(spc_dct, spc_name, thy_info, save_prefix,
                  zma_locs=(0,)):
    """  split up the unstable species
    """

    # Get filesys
    spc_info = filesys.inf.get_spc_info(spc_dct[spc_name])
    mod_thy_info = filesys.inf.modify_orb_restrict(spc_info, thy_info)
    thy_save_fs, thy_save_path = filesys.build.spc_thy_fs_from_root(
        save_prefix, spc_info, mod_thy_info)
    cnf_save_fs, cnf_save_locs = filesys.build.cnf_fs_from_thy(
        thy_save_path, mod_thy_info, cnf='min', saddle=False)
    cnf_save_paths = filesys.build.cnf_paths_from_locs(
        cnf_save_fs, cnf_save_locs)

    zma_save_fs = autofile.fs.manager(cnf_save_paths[0], 'ZMATRIX')

    # Read the zma for the unstable species
    instab_zma = zma_save_fs[-1].file.zmatrix.read(zma_locs)

    # Read the instability transformation information from the filesystem
    tra = thy_save_fs[-1].file.transformation.read(mod_thy_info[1:4])
    frm_bnd_key, brk_bnd_key = tra
    # rcts_gra = save_fs[-1].file.reactant_graph.write(locs)

    # Obtain the inchi strings for the species it breaks in to
    constituent_ichs = automol.zmatrix.ts.zmatrix_reactant_inchis(
        instab_zma, frm_bnd_key, brk_bnd_key)

    # Obtain the product names from the species dct
    prd_names = []
    for ich in constituent_ichs:
        for name, spc_dct_i in spc_dct.items():
            if ich == spc_dct_i.get('inchi'):
                prd_names.append(name)

    return prd_names
