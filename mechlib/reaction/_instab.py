"""
 Library to deal unstable species
"""

import automol
from mechlib import filesys
from mechlib.amech_io import printer as ioprinter


# Handle reaction lst
def split_unstable_full(pes_rlst, spc_rlst, spc_dct,
                        spc_model_dct_i, save_prefix):
    """ Loop over the pes reaction list and break up the unstable species
    """

    # Get split names from full PES run lst
    _split_rxn_names = ()
    if pes_rlst is not None:
        for _, rxn_lst in pes_rlst.items():
            split_rxn_lst, _ = split_unstable_pes(
                rxn_lst, spc_dct, spc_model_dct_i, save_prefix)
            for split_rxn in split_rxn_lst:
                _, (new_rcts, new_prds) = split_rxn
                _split_rxn_names += new_rcts
                _split_rxn_names += new_prds

    # Get split names from spc
    _split_spc_names = ()
    if spc_rlst is not None:
        _split_spc_names = split_unstable_spc(
             spc_rlst, spc_dct, spc_model_dct_i, save_prefix)
        _split_spc_names = tuple(_split_spc_names.values())[0]

    # Combine both and remove duplicates
    _split_names = _split_rxn_names + _split_spc_names
    split_names = tuple(i for n, i in enumerate(_split_names)
                        if i not in _split_names[:n])

    return {('SPC', 0, 0): split_names}


def split_unstable_pes(rxn_lst, spc_dct, spc_model_dct_i, save_prefix):
    """ Build a new list of reactions for a given PES where all of the
        reactant and product species of the channels are assessed for
        instability and broken up.

        :param rxn_lst:
        :type rxn_lst:
        :param spc_dct: species information
           dict[spc_name: spc_information]
        :param spc_mod_dct_i:
        :type spc_mod_dct_i: dict[]
    """

    # Get theory
    thy_info = spc_model_dct_i['vib']['geolvl'][1][1]

    # Loop over the reactions and split
    new_rxn_lst = ()
    unstable_chnl_idxs = ()
    for rxn in rxn_lst:

        # Unpack the reaction
        chnl_idx, (rcts, prds) = rxn

        # Build the mapping dictionary for the rxn species
        rxn_names = rcts + prds
        split_map = _split_mapping(spc_dct, thy_info, save_prefix,
                                   spc_names=rxn_names, zma_locs=(0,))

        # Assess and split the reactants and products for unstable species
        new_rcts = ()
        for rct in rcts:
            new_rcts += split_map[rct]

        new_prds = ()
        for prd in prds:
            new_prds += split_map[prd]

        # Flip the reaction if the reactants are unstable?
        # Append chnl idx if either reactants or products found to be unstable
        if (len(rcts) < len(new_rcts)) or (len(prds) < len(new_prds)):
            unstable_chnl_idxs += (chnl_idx,)

        # Append to list
        new_rxn = ((chnl_idx, (new_rcts, new_prds)),)
        new_rxn_lst += new_rxn

        # Check if the split species are in the spc dct
        if len(rcts) > len(new_rcts):
            print('WARNING: REACTANTS FROM SPLIT MISSING FROM SPC DCT')
        if len(prds) > len(new_prds):
            print('WARNING: PRODUCTS FROM SPLIT MISSING FROM SPC DCT')

    return new_rxn_lst, unstable_chnl_idxs


def split_unstable_spc(spc_rlst, spc_dct, spc_model_dct_i, save_prefix):
    """ Build a new list of species where each species of the input list
        has been assessed for instability and broken up.

        :param spc_rlst:
        :type spc_rlst:
        :param spc_dct: species information
           dict[spc_name: spc_information]
        :param spc_mod_dct_i:
        :type spc_mod_dct_i: dict[]
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
    """

    # Get theory
    thy_info = spc_model_dct_i['vib']['geolvl'][1][1]

    # Split each species
    _split_spc_names = ()
    for spc in list(spc_rlst.values())[0]:
        split_names = _split_species(
            spc_dct, spc, thy_info, save_prefix, zma_locs=(0,))
        if split_names:
            _split_spc_names += split_names
        else:
            _split_spc_names += (spc,)
    split_spc_names = tuple(i for n, i in enumerate(_split_spc_names)
                            if i not in _split_spc_names[:n])

    return {('SPC', 0, 0): split_spc_names}


def _split_mapping(spc_dct, thy_info, save_prefix,
                   spc_names=None, zma_locs=(0,)):
    """ Build a dictionary that describes how species decomposes into
        smaller species via some radical stability. Dictionary maps the
        species name to the names of the decomposition products.

        If no names provided, a mapping will be generated for all species
        in the species dictionary.

        :param spc_dct: species information
           dict[spc_name: spc_information]
        :param spc_names: mechanism names of species to assess
        :type spc_names: str
        :param thy_info: ???
        :type thy_info: ???
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
        :param zma_locs: locs for zma filesys (put in spc dct)
        :type zma_locs:
    """

    if spc_names is None:
        spc_names = tuple(name for name in spc_dct.keys() if 'ts_' not in name)

    split_map = {}
    for spc_name in spc_names:
        split_names = _split_species(
            spc_dct, spc_name, thy_info, save_prefix, zma_locs=zma_locs)
        if split_names:
            split_map[spc_name] = split_names
        else:
            split_map[spc_name] = (spc_name,)

    return split_map


def _split_species(spc_dct, spc_name, thy_info, save_prefix,
                   zma_locs=(0,)):
    """ Assess if a given species has an instability transformation
        file located in the save filesystem within a Z-Matrix layer:
        SPC/THY/CONFS/Z/ which are specified by the provided info.

        If file is found, use the contained information to break-up
        the species into products of the instability transformation. If
        no file found, return species.

        :param spc_dct: species information
           dict[spc_name: spc_information]
        :param spc_name: mechanism name of species to assess
        :type spc_name: str
        :param thy_info: ???
        :type thy_info: ???
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
        :param zma_locs: locs for zma filesys (put in spc dct)
        :type zma_locs:
    """

    # Initialize an empty list
    split_names = ()

    # Attempt to read the graph of the instability trans
    # Get the product graphs and inchis
    tra, path = filesys.read.instability_transformation(
        spc_dct, spc_name, thy_info, save_prefix, zma_locs=zma_locs)

    if tra is not None:
        ioprinter.info_message('\nFound instability files at path:')
        ioprinter.info_message(f'  {path}')

        zrxn, _ = tra
        prd_gras = automol.reac.product_graphs(zrxn)
        constituent_ichs = tuple(automol.graph.inchi(gra, stereo=True)
                                 for gra in prd_gras)

        _split_names = ()
        for ich in constituent_ichs:
            for name, spc_dct_i in spc_dct.items():
                # Try to match inchis with stereo included in checks
                if ich == spc_dct_i.get('inchi'):
                    _split_names += (name,)
                    break
                # Remove stereo since we used to not store this data
                # ich_noste1 = automol.inchi.standard_form(ich, stereo=False)
                # ich_noste2 = automol.inchi.standard_form(
                #     spc_dct_i.get('inchi'), stereo=False)
                # if ich_noste1 == ich_noste2:
                #     _split_names += (name,)
                #     break
        split_names = tuple(i for n, i in enumerate(_split_names)
                            if i not in _split_names[:n])

        ioprinter.info_message(f'- Splitting species {spc_name}'
                               f'into {split_names}')
        if len(split_names) < 2:
            ioprinter.warning_message(
                'Could not match all following InChI strings '
                '(corresponding to instability products) to\n'
                'to ones currently defined in the species.csv file:')
            for ich in constituent_ichs:
                ioprinter.info_message(f'  - {ich}')

    return split_names
