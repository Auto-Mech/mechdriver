"""
 Library to deal unstable species
"""

import automol
import autofile
import elstruct
from mechanalyzer.inf import thy as tinfo
from mechlib import filesys


# Handle reaction lst
def split_unstable(rxn_lst, spc_dct, spc_model_dct, thy_dct, save_prefix):
    """ Loop over the reaction list and break up the unstable species
    """

    new_rxn_lst = []
    for rxn in rxn_lst:

        # Initialize dct
        new_rxn = {}

        # Get theory
        spc_model = rxn['model'][1]
        geo_model = spc_model_dct[spc_model]['es']['geo']
        ini_thy_info = tinfo.from_dct(thy_dct[geo_model])

        new_rxn['dummy'] = []

        # Asses the reactants for unstable species
        new_rxn['reacs'] = []
        for rct in rxn['reacs']:
            split_rct = _split_species(spc_dct, rct,
                                       ini_thy_info, save_prefix)
            if not split_rct:
                new_rct = rct
                new_rxn['reacs'].append(new_rct)
            else:
                print('\nSplitting species...')
                new_rct = split_rct
                print('- New species: {}'.format(' '.join(new_rct)))
                new_rxn['reacs'].extend(new_rct)
                new_rxn['dummy'].append('reacs')

        # Assess the products for unstable species
        new_rxn['prods'] = []
        for prd in rxn['prods']:
            split_prd = _split_species(spc_dct, prd,
                                       ini_thy_info, save_prefix)
            if not split_prd:
                new_prd = prd
                new_rxn['prods'].append(new_prd)
            else:
                print('- Splitting species...')
                new_prd = split_prd
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


def _split_dct(rxn_lst, spc_dct, thy_info, save_prefix,
               zma_locs=(0,)):
    """ Build a dictionry which maps the names of species into splits
        would like to build to just go over spc dct (good for mech pre-process)
        could do 
    """

    split_map = {}
    for rxn in rxn_lst:

        # Get theory
        spc_model = rxn['model'][1]
        geo_model = spc_model_dct[spc_model]['es']['geo']
        ini_thy_info = tinfo.from_dct(geo_model)
        
        for spc in rxn['species']:
            if spc not in split_map:
                split_names = _split_species(
                    spc_dct, spc, ini_thy_info, save_prefix)
                if split_names:
                    print('- Splitting species...')
                    print('- New species: {}'.format(' '.join(split_names)))
                    split_map[spc] = split_names
                else:
                    split_map[spc] = spc

    return split_map


def _split_species(spc_dct, spc_name, thy_info, save_prefix,
                   zma_locs=(0,)):
    """  split up the unstable species
    """

    # Initialize an empty list
    prd_names = []

    tra = filesys.read.instability_transformation(
        spc_dct, spc_name, thy_info, save_prefix, zma_locs=zma_locs)
    if tra is not None:
        zrxn, _ = tra
        prd_gras = automol.reac.product_graphs(zrxn)
        constituent_ichs = tuple(automol.graph.inchi(gra, stereo=True) 
                                 for gra in prd_gras)

        for ich in constituent_ichs:
            for name, spc_dct_i in spc_dct.items():
                if ich == spc_dct_i.get('inchi'):
                    prd_names.append(name)
                    break
        prd_names = list(set(prd_names))

    return prd_names
