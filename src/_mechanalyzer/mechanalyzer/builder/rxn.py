""" Generate sets of reactions from reactants to construct mechanisms
"""

import itertools
import automol
from autorun import execute_function_in_parallel
from mechanalyzer.builder._update import update_spc_dct_from_reactions
from mechanalyzer.builder._update import update_rxn_dct
from mechanalyzer.builder._stereo import _add_third
from mechanalyzer.builder._stereo import _ste_rxn_lsts
from mechanalyzer.builder._stereo import _remove_enantiomer_reactions
from mechanalyzer.builder._stereo import _stereo_results
from chemkin_io.writer._util import format_rxn_name


# MAIN CALLABLE FUNCTIONS FOR GENERATING REACTION LISTS
def build_mechanism(mech_spc_dct, mech_rxn_dct, rxn_series, stereo=False, nprocs=1):
    """ Use the lst of reactions to build objects to describe mechanism
    """
    def _print_species_info(rcts_lst):
        # Use SMILES to print info message for what reactions generating
        mult_array = ['', 'singlet', 'doublet', 'triplet', 'quartet']
        for names in rcts_lst:
            _rct_smis = tuple(map(
                lambda name: automol.chi.smiles(mech_spc_dct[name]['inchi']), names))
            _rct_mults = tuple(map(
                lambda name: mult_array[int(mech_spc_dct[name]['mult'])], names))
            print(f'{names}: {_rct_smis}, {_rct_mults}')
        print('')
        if allowed_prds_lst:
            allowed_prd_smis = tuple(map(
                lambda name: automol.chi.smiles(mech_spc_dct[name]['inchi']), allowed_prds_lst))
            allow_str = ' '.join(allowed_prd_smis)
            print(f'\nEach reaction must produce {allow_str}')
            print('')

    # Initialize new_spc (needed for sequential steps
    # new_spc_lst = None
    print('---------------------------------------------------------\n')
    # Loop over the reaction series consisting of reactants and reaction type
    for sidx, series in enumerate(rxn_series):
        rct1_set, rct2_set, allowed_prds_lst, rxn_typs = series

        print(f'Generating Reactions for Series {sidx+1}')

        # Generate the reactions from the reactants
        ini_rxns_infos = ()
        for rtyp in rxn_typs:

            ReacInfo = automol.reac.reaction_info_from_string(rtyp)
            # Determine reactants to generate reactions for
            rcts_lst = _determine_reactants(
                mech_spc_dct, rct1_set, rct2_set, ReacInfo)
            # allowed_prds_lst = _determine_allowed_products(
            #     mech_spc_dct, allowed_prds)
            _print_species_info(rcts_lst)

            # Generate Reactions
            # for ichs in rct_ichs:
            #    ini_rxns += generate_reactions(ichs, allowed_prd_ichs, rtyp)
            ini_rxns_infos += execute_function_in_parallel(
                generate_reactions, rcts_lst, (mech_spc_dct, allowed_prds_lst, ReacInfo),
                nprocs=nprocs)
        # Add stereo
        if stereo:
            rxns = ()
            for _rxn_info in ini_rxns_infos:

                _rcts_ichs = (inf[0] for inf in _rxn_info[0])
                _prds_ichs = (inf[0] for inf in _rxn_info[1])
                log1 = ('\nExpanding Stereo for Reaction: '
                        f'{format_rxn_name((_rcts_ichs, _prds_ichs, _rxn_info[2],))}\n')
                _tmp_rxn = (_rcts_ichs, _prds_ichs)
                thrdbdy = _rxn_info[2]
                ste_rxns_lst, log2 = _ste_rxn_lsts(_tmp_rxn)
                ste_rxns_lst, rem_rxns, log3 = _remove_enantiomer_reactions(
                    ste_rxns_lst, reacs_stereo_inchi=_rxn[0])
                ste_infos_lst = ()
                # TODO
                # for ste_rxn in ste_rxns_lst:
                #     ste_infos_lst += (
                #         (rct_inf[0], ,
                # rxns += _add_third(ste_rxns_lst, thrdbdy)
                # log4 = _stereo_results(_rxn, ste_rxns_lst, rem_rxns)
                # print(log1 + log2 + log3 + log4)
        else:
            rxns_infos = ini_rxns_infos

        print('rxn_infos')
        print(rxns_infos)
        # Update the mechanism objects with unique spc and rxns
        mech_spc_dct = update_spc_dct_from_reactions(
            rxns_infos, mech_spc_dct)
        print(list(mech_spc_dct.keys()))
        mech_rxn_dct = update_rxn_dct(
            rxns_infos, mech_rxn_dct, mech_spc_dct)
        print(list(mech_rxn_dct.keys()))

        print('\n---------------------------------------------------------\n')

    return mech_spc_dct, mech_rxn_dct


def generate_reactions(
        mech_spc_dct, allowed_prds_lst, ReacInfo, rcts_lst, output_queue=None):
    """ For a given reactants
    """
    rxns_info = ()
    for rct_names in rcts_lst:
        print('rct names', rct_names)
        rct_ichs = [mech_spc_dct[name]['inchi'] for name in rct_names]
        # rct_mults = [mech_spc_dct[name]['mult'] for name in rct_names]
        _rct_smis = tuple(map(automol.chi.smiles, rct_ichs))
        log = f'Generating Reactions for {_rct_smis}...'

        # Generate the products with the desired reactant and reaction type
        rct_gras = _rct_gras(rct_ichs)
        prds_info_lst = _generated_prds_info_lst(rct_gras, ReacInfo)
        # prds_info_lst = _generated_prds_info_lst(rct_gras, rct_mults, ReacInfo)

        if allowed_prds_lst is not None:
            allowed_prds_info_lst = ((
                mech_spc_dct[name]['inchi'],
                mech_spc_dct[name]['charge'],
                mech_spc_dct[name]['mult'],) for name in allowed_prds_lst)
        else:
            allowed_prds_info_lst = ()
        rcts_info = tuple((
            mech_spc_dct[name]['inchi'],
            mech_spc_dct[name]['charge'],
            mech_spc_dct[name]['mult'],) for name in rct_names)

        # Build list of generated reactions including reactants and products
        for pidx, prds_info in enumerate(prds_info_lst):
            # Don't add if requested products not found
            if allowed_prds_lst:
                # if not all(prd in prds for prd in allowed_prd_ichs):
                if not all(prd in allowed_prds_info_lst for prd in prds_info):
                    continue
            # Don't add if self reaction was generated; weak check: ignores stereo
            if set(rct_ichs) == set((prd[0] for prd in prds_info)):
                continue
            # If continue not hit, save reaction to list and print to stdout
            _prd_smis = tuple(map(automol.chi.smiles, (prd[0] for prd in prds_info)))
            log += f'\nFound Product(s) {pidx+1}: {_prd_smis}'

            rxns_info += ((rcts_info, prds_info, (None,)),)

        if not prds_info_lst:
            log += '\nNO Product(s) Found'
        print(log)
    output_queue.put(rxns_info)


# Functions to set lists for mech building step
def _determine_reactants(spc_dct, rct1_set, rct2_set, ReacInfo):
    """ Determine a list of reactants for one or more reactions
        these could be unimolecular or bimolecular.
    """

    def _gen_set(spc_dct, spc_set):
        """ Make the list of species that will serve as
            reactant 1 or reactant 2
        """

        # Check if string id of class of species (e.g., 'all', 'radicals'), or
        # simply list of species name given
        if isinstance(spc_set, str):

            # Build ini list from input or spc dct
            spc_names, spc_info_lst = (), ()
            for name, dct in spc_dct.items():
                if (
                        (spc_set == 'radicals' and dct['mult'] % 2 == 0) or
                        spc_set != 'radicals'):
                    spc_names += (name,)
        else:
            spc_names = spc_set

        return spc_names

    rct1_names = _gen_set(spc_dct, rct1_set)
    if automol.ReactionClass.is_bimolecular(ReacInfo.reaction_class()):
        rct2_names = _gen_set(spc_dct, rct2_set)
        rcts_names = tuple(itertools.product(rct1_names, rct2_names))
    else:
        rcts_names = tuple((name,) for name in rct1_names)

    return rcts_names


def _determine_allowed_products(spc_dct, allowed_prds):
    """ Takes a list of species mechanism names for the products
        that must be produced in a reaction for a generated reaction
        to be saved and generates the corresponding InChI strings.
    """
    if allowed_prds is not None:
        spc_info_lst = tuple(
            (spc_dct[name]['inchi'], spc_dct[name]['mult'], spc_dct[name]['charge'])
            for name in allowed_prds)
    else:
        spc_info_lst = ()

    return spc_info_lst


# # Helper functions
# def _radicals(ich_lst, name_lst):
#     """ Determine the radicals
#     """
#     rad_ichs, rad_names = (), ()
#     for ich, name in zip(ich_lst, name_lst):
#         if automol.graph.is_radical_species(automol.chi.graph(ich)):
#             rad_ichs += (ich,)
#             rad_names += (name,)
#
#     return rad_ichs, rad_names


def _rct_gras(rct_ichs):
    """ Get reactant graphs from smiles
    """

    rct_gras = list(map(automol.amchi.graph, rct_ichs))
    rct_gras = list(map(automol.graph.without_stereo, rct_gras))
    rct_gras = list(map(automol.graph.explicit, rct_gras))
    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)

    return rct_gras


def _generated_prds_info_lst(rct_gras, ReacInfo, check=True):
    """ Check the products
    """

    def _product_multiplicities(gras, ReacInfo):
        return tuple(
            automol.graph.maximum_spin_multiplicity(gra)
            if  ReacInfo.is_high_spin() else
            automol.graph.possible_spin_multiplicities(gra)[0]
            for gra in gras)

    def _product_charges(gras):
        #TODO
        return tuple(0 for gra in gras)

    # Enumerate all possible reactions checking they are correct
    rxns = automol.reac.enumerate_reactions(
        rct_gras, rxn_type=ReacInfo.reaction_class())

    # Get InChI strings from the enumerated reaction objects.
    # Check for validity if requested
    prd_info_lst = ()
    for rxn in rxns:
        prd_gras_ = automol.reac.product_graphs(rxn, stereo=False)
        prd_ichs_ = tuple(map(automol.graph.chi, prd_gras_))
        prd_mults_ = _product_multiplicities(prd_gras_, ReacInfo)
        prd_chgs_ = _product_charges(prd_gras_)
        prd_info_lst += (tuple(zip(prd_ichs_, prd_chgs_, prd_mults_)),)

        if check:
            rct_gras_ = automol.reac.reactant_graphs(rxn, stereo=False)
            rxns_ = automol.reac.find(rct_gras_, prd_gras_)
            try:
                assert rct_gras_ == rct_gras
                assert any(automol.reac.class_(r) == ReacInfo.reaction_class() for r in rxns_)
            except AssertionError:
                print('WARNING: issue with reaction')

    # Remove duplicates (some reactions could have Reaction Objects, mult TS)
    # prd_ichs = automol.util.remove_duplicates_with_order(prd_ichs)

    return prd_info_lst
