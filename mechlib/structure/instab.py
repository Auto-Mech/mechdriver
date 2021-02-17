"""
 Library to deal unstable species
"""

# init: just check geom
# conf: check ratio of confs
# hr:   check geom along scan,

import automol
import autofile
import elstruct
from mechanalyzer.inf import thy as tinfo
# from mechanalyzer.inf import rxn as rinfo
from mechanalyzer.inf import spc as sinfo
from mechlib import filesys


# Write the instability files
def write_instab(conn_zma, disconn_zma,
                 instab_save_fs, thy_locs,
                 opt_ret,
                 zma_locs=(0,),
                 save_cnf=False):
    """ write the instability files
    """

    # Get a connected geometry
    conn_geo = automol.zmat.geometry(conn_zma)

    if opt_ret:

        # Obtain inf obj and inp str to write in filesys
        inf_obj, inp_str, out_str = opt_ret

        # Set and print the save path information
        save_path = instab_save_fs[-1].path()
        print(" - Saving...")
        print(" - Save path: {}".format(save_path))

        # Save the geometry information
        instab_fs = autofile.fs.instab(save_path)
        instab_fs[-1].create()
        instab_fs[-1].file.geometry_info.write(inf_obj)
        instab_fs[-1].file.geometry_input.write(inp_str)
        instab_fs[-1].file.geometry.write(conn_geo)
        instab_path = instab_fs[-1].path()

        # Grab the zma and instability transformation
        conn_zma, frm_bnd_keys, rcts_gra = _instab_info(conn_zma, disconn_zma)
        tra = (frozenset({frm_bnd_keys}),
               frozenset({frozenset({})}))

        # Save zma information seperately, if required
        zma_save_fs = autofile.fs.zmatrix(instab_path)
        zma_save_fs[-1].create(zma_locs)
        zma_save_fs[-1].file.geometry_info.write(inf_obj, zma_locs)
        zma_save_fs[-1].file.geometry_input.write(inp_str, zma_locs)
        zma_save_fs[-1].file.zmatrix.write(conn_zma, zma_locs)

        # Write the files into the filesystem
        zma_save_fs[-1].file.transformation.write(tra, zma_locs)
        zma_save_fs[-1].file.reactant_graph.write(rcts_gra, zma_locs)

        # Saving the energy to an SP filesys
        print(" - Saving energy...")
        prog = inf_obj.prog
        method = inf_obj.method
        ene = elstruct.reader.energy(prog, method, out_str)
        sp_save_fs = autofile.fs.single_point(save_path)
        sp_save_fs[-1].create(thy_locs)
        sp_save_fs[-1].file.input.write(inp_str, thy_locs)
        sp_save_fs[-1].file.info.write(inf_obj, thy_locs)
        sp_save_fs[-1].file.energy.write(ene, thy_locs)

        if save_cnf:
            # Save the geometry information
            cnf_fs = autofile.fs.conformer(save_path)
            cnf_locs = [autofile.schema.generate_new_conformer_id()]
            cnf_fs[-1].create(cnf_locs)
            cnf_fs[-1].file.geometry_info.write(inf_obj, cnf_locs)
            cnf_fs[-1].file.geometry_input.write(inp_str, cnf_locs)
            cnf_fs[-1].file.geometry.write(conn_geo, cnf_locs)
            cnf_path = cnf_fs[-1].path(cnf_locs)

            # Save zma information seperately, if required
            zma_save_fs = autofile.fs.zmatrix(cnf_path)
            zma_save_fs[-1].create(zma_locs)
            zma_save_fs[-1].file.geometry_info.write(inf_obj, zma_locs)
            zma_save_fs[-1].file.geometry_input.write(inp_str, zma_locs)
            zma_save_fs[-1].file.zmatrix.write(conn_zma, zma_locs)

            # Saving the energy to an SP filesys
            print(" - Saving energy...")
            prog = inf_obj.prog
            method = inf_obj.method
            ene = elstruct.reader.energy(prog, method, out_str)
            sp_save_fs = autofile.fs.single_point(cnf_path)
            sp_save_fs[-1].create(thy_locs)
            sp_save_fs[-1].file.input.write(inp_str, thy_locs)
            sp_save_fs[-1].file.info.write(inf_obj, thy_locs)
            sp_save_fs[-1].file.energy.write(ene, thy_locs)


# Write the instability files
def write_instab2(conn_zma, disconn_zmas,
                  instab_save_fs, thy_locs,
                  zma_locs=(0,),
                  save_cnf=False):
    """ write the instability files
    """

    # Get a connected geometry
    conn_geo = automol.zmat.geometry(conn_zma)

    # Set and print the save path information
    save_path = instab_save_fs[-1].path()
    print(" - Saving...")
    print(" - Save path: {}".format(save_path))

    # Save the geometry information
    instab_fs = autofile.fs.instab(save_path)
    instab_fs[-1].create()
    instab_fs[-1].file.geometry.write(conn_geo)
    instab_path = instab_fs[-1].path()

    # Grab the zma and instability transformation
    conn_zma, brk_bnd_keys, rcts_gra = _instab_info(conn_zma, disconn_zmas)
    tra = (frozenset({frozenset({})}),
           frozenset({brk_bnd_keys}))

    # Save zma information seperately, if required
    zma_save_fs = autofile.fs.zmatrix(instab_path)
    zma_save_fs[-1].create(zma_locs)
    zma_save_fs[-1].file.zmatrix.write(conn_zma, zma_locs)

    # Write the files into the filesystem
    zma_save_fs[-1].file.transformation.write(tra, zma_locs)
    zma_save_fs[-1].file.reactant_graph.write(rcts_gra, zma_locs)

    if save_cnf:

        # Save the geometry information
        cnf_fs = autofile.fs.conformer(save_path)
        cnf_locs = [autofile.schema.generate_new_conformer_id()]
        cnf_fs[-1].create(cnf_locs)
        cnf_fs[-1].file.geometry.write(conn_geo, cnf_locs)
        cnf_path = cnf_fs[-1].path(cnf_locs)

        # Save zma information seperately, if required
        zma_save_fs = autofile.fs.zmatrix(cnf_path)
        zma_save_fs[-1].create(zma_locs)
        zma_save_fs[-1].file.zmatrix.write(conn_zma, zma_locs)


# Unstable check
def check_unstable_species(tsk, spc_dct, spc_name,
                           thy_info, save_prefix):
    """ see if a species and unstable and handle task management
    """

    if 'ts' not in spc_name and tsk != 'init_geom':

        print('\nChecking filesystem if species {}'.format(spc_name),
              'is unstable...')

        # Build filesystem
        spc_info = sinfo.from_dct(spc_dct[spc_name])
        mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
        thy_save_fs, _ = filesys.build.spc_thy_fs_from_root(
            save_prefix, spc_info, mod_thy_info)

        thy_locs = mod_thy_info[1:4]
        thy_path = thy_save_fs[-1].path(thy_locs)

        instab_fs = autofile.fs.instab(thy_path)
        if instab_fs[-1].exists():

            instab_path = instab_fs[-1].path()
            zma_fs = autofile.fs.zmatrix(instab_path)
            zma_locs = (0,)

            # Check if the instability files exist
            if (zma_fs[-1].file.transformation.exists(zma_locs) and
                    zma_fs[-1].file.reactant_graph.exists(zma_locs)):
                stable = False
                print('- Found files denoting species instability at path')
                print('    {}'.format(zma_fs[-1].path(zma_locs)))
            else:
                stable = True
                print('- No files denoting instability were found at path')
                print('    {}'.format(zma_fs[-1].path(zma_locs)))

        else:
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
        ini_thy_info = tinfo.from_dct(geo_model)

        new_rxn['dummy'] = []

        # Asses the reactants for unstable species
        new_rxn['reacs'] = []
        for rct in rxn['reacs']:
            rct_stable = check_unstable_species(
                'rate', spc_dct, rct, ini_thy_info, save_prefix)
            if rct_stable:
                new_rct = rct
                new_rxn['reacs'].append(new_rct)
            else:
                print('\nSplitting species...')
                new_rct = split_species(spc_dct, rct,
                                        ini_thy_info, save_prefix)
                print('- New species: {}'.format(' '.join(new_rct)))
                new_rxn['reacs'].extend(new_rct)
                new_rxn['dummy'].append('reacs')

        # Assess the products for unstable species
        new_rxn['prods'] = []
        for prd in rxn['prods']:
            prd_stable = check_unstable_species(
                'rate', spc_dct, prd, ini_thy_info, save_prefix)
            if prd_stable:
                new_prd = prd
                new_rxn['prods'].append(new_prd)
            else:
                print('- Splitting species...')
                new_prd = split_species(spc_dct, prd,
                                        ini_thy_info, save_prefix)
                print('- New species: {}'.format(' '.join(new_prd)))
                new_rxn['prods'].extend(new_prd)
                new_rxn['dummy'].append('prods')

        if len(rxn['reacs']) > len(new_rxn['reacs']):
            print('WARNING: LIKELY MISSING DATA FOR REACTANTS FOR SPLIT')
        if len(rxn['prods']) > len(new_rxn['prods']):
            print('WARNING: LIKELY MISSING DATA FOR PRODUCTS FOR SPLIT')

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
    spc_info = sinfo.from_dct(spc_dct[spc_name])
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
    _, thy_save_path = filesys.build.spc_thy_fs_from_root(
        save_prefix, spc_info, mod_thy_info)

    instab_fs = autofile.fs.instab(thy_save_path)
    instab_path = instab_fs[-1].path()

    zma_save_fs = autofile.fs.zmatrix(instab_path)

    # Read the zma for the unstable species
    instab_zma = zma_save_fs[-1].file.zmatrix.read(zma_locs)

    # Read the instability transformation information from the filesystem
    tra = zma_save_fs[-1].file.transformation.read(zma_locs)
    frm_bnd_key, brk_bnd_key = tra
    # rcts_gra = save_fs[-1].file.reactant_graph.write(locs)

    # Obtain the inchi strings for the species it breaks in to
    constituent_ichs = automol.zmat.ts.zmatrix_product_inchis(
        instab_zma, frm_bnd_key, brk_bnd_key, stereo=True)
    print('constituent ichs', constituent_ichs)

    # Obtain the product names from the species dct
    prd_names = []
    prd_ichs = []
    for ich in constituent_ichs:
        print('constituent ichs:', ich, automol.inchi.smiles(ich))
        for name, spc_dct_i in spc_dct.items():
            if ich == spc_dct_i.get('inchi'):
                if ich not in prd_ichs:
                    prd_names.append(name)
                    prd_ichs.append(ich)

    return prd_names


def break_all_unstable2(rxn_lst, spc_dct, spc_model_dct, thy_dct, save_prefix):
    """ Loop over the reaction list and break up the unstable species
    """

    new_spc_queue = []
    for spc, spc_model in spc_queue:

        # Get theory
        geo_model = spc_model_dct[spc_model[1]]['es']['geo']
        ini_thy_info = filesys.inf.get_es_info(geo_model, thy_dct)

        # Asses the reactants for unstable species
        spc_stable = check_unstable_species(
            'thermo', spc_dct, spc, ini_thy_info, save_prefix)
        if spc_stable:
            new_spc_queue.append((spc, spc_model))
        else:
            print('\nSplitting species...')
            new_spcs = split_species(spc_dct, spc,
                                     ini_thy_info, save_prefix)
            for new_spc in new_spcs:
                new_spc_queue.append((new_spc, spc_model))

    # Remove redundant species
    new_spc_queue = list(set(new_spc_queue))

    return new_spc_queue
