""" functions for 
- calling rxn classification from automol
- broad rxn classification
"""

import sys
import numpy
import pandas as pd
import automol
from mechanalyzer.parser._util import get_mult
from mechanalyzer.calculator import formulas

# list of formulas for products/reactants identified for the sublcasses
FMLS_SET = numpy.array(['H1', 'O1', 'H1O1', 'O2', 'H1O2', 'C1H3'])

def classify_unimol(rcts, prds, spc_dct):
    """ Classifies unimolecular reaction from reactants and products names:
        - A=B: isomerization
        - A=C+H/O/OH/O2/HO2/CH3: addition
        - A=C+D: decomposition
        - A=C+D+E.. : decomposition(lumped)

        For recombination (depends on multiplicity of the products)
        it would be nice to distinguish type of bond that being broken

        :param rcts: reactant names
        :type rcts: tuple
        :param prds: product names
        :type prds: tuple
        :param spc_dct:
        :type spc_dct: dict[]
    """

    # extract formula dictionary
    fml_df = formulas.extract_fml_df(spc_dct)
    mult_rcts = get_mult(rcts, spc_dct)
    mult_prds = get_mult(prds, spc_dct)

    if len(prds) == 1:
        rxn_class_broad = 'Isomerization'
    elif len(prds) == 2:
        rxn_class_broad = 'Decomposition'
        # derive products multiplicity and formulas
        if mult_rcts == 1 and mult_prds >= 4:
            rxn_class_broad = 'Bond fission'
        elif mult_rcts > 1 and mult_prds == 2:
            rxn_class_broad = 'Beta-scission'

        prds_fmls = numpy.array(
            [fml_df['fml'][prds[0]], fml_df['fml'][prds[1]]])

        # check product composition
        # give priority to the second product (should be the lightest)
        if any(prds_fmls[1] == FMLS_SET):
            rxn_class_broad += f' +{prds[1]}'
        elif any(prds_fmls[0] == FMLS_SET):
            rxn_class_broad += f' +{prds[0]}'

    elif len(prds) > 2:
        rxn_class_broad = 'Decomposition(lumped)'

    return rxn_class_broad

def classify_bimol(rcts, prds, spc_dct):
    """ Classifies bimolecular reactions from reactants and products names
        with 1 product:

        - A+H/O/OH/O2/HO2/CH3 = B
        - A+R=B+RH: Habstraction-R (subclass indicates the abstractor)
        - A+R=A+R: isomerization-bim (isomerization aided by a radical.
          ex. CH2+H=CH2(S)+H, C6H6+H=FULV+H)
        - A+B=C+D: addition-decomposition - branch/prop/term
        - A+B=C+D+E..: addition-decomposition(lumped) - branch/prop/term

        For recombination (depends on the multiplicity of the reactants)
        with 2 products.

        :param rcts: reactant names
        :type rcts: tuple
        :param prds: product names
        :type prds: tuple
        :param spc_dct:
        :type spc_dct: dict[]
    """

    # extract formula dictionary
    fml_df = formulas.extract_fml_df(spc_dct)

    # extracts reactants and products multiplicity
    mult_rcts = get_mult(rcts, spc_dct)
    mult_prds = get_mult(prds, spc_dct)

    if mult_rcts < 4:
        rxn_class_broad = 'Addition'
    elif mult_rcts >= 4:
        rxn_class_broad = 'Recombination'

    if len(prds) == 1:

        rcts_fmls = numpy.array(
            [fml_df['fml'][rcts[0]], fml_df['fml'][rcts[1]]])

        # check reactant composition
        # give priority to the second reactant (should be the lightest)
        if any(rcts_fmls[1] == FMLS_SET):
            rxn_class_broad += f' {rcts[1]}'
        elif any(rcts_fmls[0] == FMLS_SET):
            rxn_class_broad += f' {rcts[0]}'

    elif len(prds) == 2:

        rxn_class_broad += '-decomposition'

        # H abstraction
        flag_habs = classify_habs(rcts, prds, fml_df, spc_dct)
        # Bimolecular isomerization
        flag_isom_bim = classify_isom_bim(rcts, prds, fml_df)

        if flag_habs == 1:
            rxn_class_broad = 'H abstraction'
        elif flag_isom_bim == 1:
            rxn_class_broad = 'Bimol Isomerization'

    elif len(prds) > 2:
        rxn_class_broad += '-decomposition(lumped)'

    # add subclass related to branching/propagation/termination
    # for addition-decomposition rxns
    if 'decomposition' in rxn_class_broad:
        # check if branching/propagation/termination
        rxn_class_broad += bran_prop_term(mult_rcts, mult_prds)

    return rxn_class_broad

def classify_isom_bim(rcts, prds, fml_df):
    """ Check if an A+B=C+D reaction is a bimolecular isomerization
        of the kind A+R=B+R by checking

        if the reactants and the products both match
        - 1 couple of rct/prd must be the exact same species
        - the other couple of rct/prd must match the stoichiometry

        :param rcts: reactant names
        :type rcts: tuple
        :param prds: product names
        :type prds: tuple
        :param fml_df:
        :type fml_df:
        :rtype: bool
    """
    rcts = numpy.array(rcts)
    prds = numpy.array(prds)
    rcts_fmls = numpy.array(
        [fml_df['fml'][rcts[0]], fml_df['fml'][rcts[1]]])
    prds_fmls = numpy.array(
        [fml_df['fml'][prds[0]], fml_df['fml'][prds[1]]])

    flag_species = (any(rcts[0] == prds) or any(rcts[1] == prds))
    flag_stoich = (any(rcts_fmls[0] == prds_fmls)
                   and any(rcts_fmls[1] == prds_fmls))

    return flag_species and flag_stoich

def classify_habs(rcts, prds, fml_df, spc_dct):
    """ Check if an A+B=C+D reaction is an hydrogen abstraction
        based on the stoichiometries and multiplicities of
        reactants and products.

        :rtype: bool
    """

    if not isinstance(rcts, tuple) or not isinstance(prds, tuple):
        print('error: reactants and products are not tuples')
        sys.exit()
    elif len(rcts) != 2 or len(prds) != 2:
        print('error: not A+B=C+D reaction')
        sys.exit()

    mult_rct = numpy.array([get_mult(rcts[0], spc_dct),
                            get_mult(rcts[1], spc_dct)])
    mult_prd = numpy.array([get_mult(prds[0], spc_dct),
                            get_mult(prds[1], spc_dct)])

    if (any(mult_rct == 1) and any(mult_rct > 1) and
            any(mult_prd == 1) and any(mult_prd > 1)):

        species_rct = rcts[numpy.where(mult_rct == 1)[0][0]]
        species_prd = prds[numpy.where(mult_prd > 1)[0][0]]
        stoich_add = [0, -1, 0]
        flag_try_prd2 = False

    # Habs with O2 and O: multiplicities are 1*3=2*2
    elif (any(mult_rct == 1) and any(mult_rct == 3) and
          all(mult_prd == 2)):

        species_rct = rcts[numpy.where(mult_rct == 3)[0][0]]
        species_prd = prds[0]
        stoich_add = [0, +1, 0]
        # try the second product in case the first is not right
        flag_try_prd2 = True

    try:
        flag_habs = set_flag_habs(species_rct, species_prd, fml_df, stoich_add)

        if flag_try_prd2 and not flag_habs:
            # try the second product
            species_prd = prds[1]
            flag_habs = set_flag_habs(
                species_rct, species_prd, fml_df, stoich_add)
    except NameError:
        flag_habs = False

    return flag_habs

def set_flag_habs(spc_rct, spc_prd, fml_df, stoich_add):
    """ Uses the formulae of the reactant and product to derive
        the product target and checks if the rct+stoich_add corresponds
        to the product one.

        :param spc_rct:
        :type spc_rct:
        :param spc_prd:
        :type spc_prd:
        :param fml_df:
        :type fml_df:
        :param stoich_add:
        :type: stoich_add:
    """

    stoich_rct = fml_df.loc[spc_rct][['nC', 'nH', 'nO']].values
    stoich_prd = fml_df.loc[spc_prd][['nC', 'nH', 'nO']].values
    stoich_prd_target = stoich_rct+stoich_add

    return all(stoich_prd == stoich_prd_target)

def bran_prop_term(rct_muls, prd_muls):
    """ Checks if a reaction can be further classified as a
        propagation, termination, or branching reaction using the
        reaction multiplicities.

        :param rct_muls: reactant multiplicities
        :type rct_muls: tuple(int)
        :param prd_muls: product multiplicities
        :type prd_muls: tuple(int)
    """

    if rct_muls == prd_muls:
        add = ' - propagation'
    elif rct_muls > prd_muls:
        add = ' - termination'
    elif rct_muls < prd_muls:
        add = ' - branching'
    else:
        add = ''

    return add

# FUNCTIONS FOR RXN GRAPH CLASSIFICATION #


def classify_graph(spc_dct, rct_names, prd_names):
    """ calls the graph classifier for a given reaction

    :param spc_dct: species dictionary
    :param rct_names: reactant names (r1, r2, )
    :param prd_names: product names (p1, p2, )

    :returns: reaction class (first of the possible identified classes)
    :rtype: str
    """

    # ID reaction
    rct_fmls = tuple(spc_dct[rct]['fml'] for rct in rct_names)
    prd_fmls = tuple(spc_dct[prd]['fml'] for prd in prd_names)

    rct_ichs = tuple(spc_dct[spc]['inchi'] for spc in rct_names)
    prd_ichs = tuple(spc_dct[spc]['inchi'] for spc in prd_names)

    if automol.form.reac.is_valid_reaction(rct_fmls, prd_fmls):
        try:
            rxn_objs = automol.reac.from_chis(
                rct_ichs, prd_ichs, stereo = False)
            rxn_classes = tuple(automol.reac.class_(obj) for obj in rxn_objs)
        except AssertionError:
            rxn_classes = ('AssertionError', )
        except TypeError:
            print('geoms of rxn classifier fail for rxn: '
                    f'{rct_ichs} = {prd_ichs}')
            rxn_classes = ('TypeError', )

        if rxn_classes:
            # save only the first possible reaction type
            # rclass = rxn_classes[0]
            rclass = '/'.join(set(rxn_classes))
        else:
            rclass = 'unclassified'

    else:
        rclass = 'unclassified - Wrong Stoichiometry'

    return rclass


def classify_ws(subpes_df, elem_reac_df, species_subpes, rxn):
    """ classifies well skipping channels of a given subpes
        WARNING: STILL UNDER CONSTRUCTION - SOME TEMPORARY FEATURES

    :param subpes_df: dataframe with subpes info
    :param elem_reac_df: dataframe with elementary reaction channels of subpes
    :param species_subpes: list of subpes species
    :param rxn: string with rxn belonging to the subpes

    :returns: reaction class of the WS channel considered
    :rtype: str
    """
    # derive unimolecular species list
    mult_species_subpes = pd.Series(
        list(map(len, species_subpes)), index=species_subpes)
    unimol_species = mult_species_subpes[mult_species_subpes == 1].index

    rct_names = subpes_df['rct_names_lst_ord'][rxn]
    prd_names = subpes_df['prd_names_lst_ord'][rxn]
    # reactants: if bimolecular, find the label of the elementary reaction
    # going to unimolecular species; if unimol, label is 'isom'
    # isolate A+B->C and C->A+B connections
    rxn_types = elem_reac_df[rct_names][unimol_species]
    rxn_types = rxn_types[rxn_types != 'unclassified']
    rxn_types_1 = rxn_types[rxn_types != '']

    rxn_types = elem_reac_df[prd_names][unimol_species]
    rxn_types = rxn_types[rxn_types != 'unclassified']
    rxn_types_2 = rxn_types[rxn_types != '']

    try:
        # TEMPORARY: SHOULD RECONSTRUCT FULL PATH FROM REACTANTS TO PRODUCTS
        rxn_type_1 = rxn_types_1.iloc[0] # rxn_types_1[0]
        # TEMPORARY: SHOULD RECONSTRUCT FULL PATH FROM REACTANTS TO PRODUCTS
        rxn_type_2 = rxn_types_2.iloc[0] # rxn_types_2[0]

        # WRITE THE REACTION TYPE STRING
        rxn_type_ws = rxn_type_1 + '-' + rxn_type_2 + ' (WS)'
        return rxn_type_ws

    except IndexError:
        return None

