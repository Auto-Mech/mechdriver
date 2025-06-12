""" Handles names for mechanisms
"""

import copy
import ioformat
import automol.chi
import automol.graph
from automol.graph import FunctionalGroup
import automol.form


# Name remaping function
def remap_mechanism_names(mech_spc_dct, rxn_param_dct, map_dct):
    """ Change all of the names in a spc_dct and rxn_param_dct
        according to the provided dictionary
    """
    def remap_single_rxn(rxn, map_dct):
        """ Remaps names for a single reaction
        """
        rcts, prds, thrd = rxn
        re_rcts = remap_rcts_or_prds(rcts, map_dct)
        re_prds = remap_rcts_or_prds(prds, map_dct)

        return (re_rcts, re_prds, thrd)

    def remap_rcts_or_prds(spcs, map_dct):
        """ Remaps names for a single set of species (rcts or prds)
        """
        re_spcs = []
        for spc in spcs:
            # Remap species name if in map_dct
            if map_dct.get(spc) is not None:
                re_spcs.append(map_dct.get(spc))
            # If species not in map dct, don't remap
            else:
                re_spcs.append(spc)
        re_spcs = tuple(re_spcs)

        return re_spcs

    # Fill the map dict if any names are missing
    for name in mech_spc_dct:
        if name not in map_dct:
            map_dct[name] = name

    # Alter the names of the dictionary
    re_mech_spc_dct = {}
    for name, dct in mech_spc_dct.items():
        re_mech_spc_dct[map_dct[name]] = dct

    # Alter the names of the reactions
    re_rxn_param_dct = {}
    for rxn, params in rxn_param_dct.items():
        #        rcts, prds, thrd = rxn
        #        re_rxn = (
        #            tuple(map_dct[rct] for rct in rcts),
        #            tuple(map_dct[prd] for prd in prds),
        #            thrd
        #        )
        re_rxn = remap_single_rxn(rxn, map_dct)
        re_rxn_param_dct[re_rxn] = params

    return re_mech_spc_dct, re_rxn_param_dct


# FUNCTIONAL NAME MAPPING
# Build various dictionaries
DEFAULT_FGRP_RENAME_RULE_DCT = {
    # Three Functional Groups
    'ADHY': (FunctionalGroup.ALKENE,
             FunctionalGroup.HYDROPEROXY,
             FunctionalGroup.HYDROPEROXY),
    'AAOH': (FunctionalGroup.ALCOHOL,
             FunctionalGroup.ALKENE,
             FunctionalGroup.ALKOXY),
    'AAK': (FunctionalGroup.ALCOHOL,
            FunctionalGroup.ALKENE,
            FunctionalGroup.KETONE),
    'EDHY': (FunctionalGroup.CYCLIC_ETHER,
             FunctionalGroup.HYDROPEROXY,
             FunctionalGroup.HYDROPEROXY),
    'PDHY': (FunctionalGroup.HYDROPEROXY,
             FunctionalGroup.HYDROPEROXY,
             FunctionalGroup.PEROXY),
    'THY': (FunctionalGroup.HYDROPEROXY,
            FunctionalGroup.HYDROPEROXY,
            FunctionalGroup.HYDROPEROXY),
    # Two Functional Groups
    'ADOH': (FunctionalGroup.ALCOHOL, FunctionalGroup.ALDEHYDE),
    'ALOH': (FunctionalGroup.ALCOHOL, FunctionalGroup.ALKENE),
    'KEOH': (FunctionalGroup.ALCOHOL, FunctionalGroup.KETONE),
    'DALD': (FunctionalGroup.ALDEHYDE, FunctionalGroup.ALDEHYDE),
    'ALAD': (FunctionalGroup.ALDEHYDE, FunctionalGroup.ALKENE),
    'AKOX': (FunctionalGroup.ALDEHYDE, FunctionalGroup.ALKOXY),
    'AHP': (FunctionalGroup.ALDEHYDE, FunctionalGroup.HYDROPEROXY),
    'KALD': (FunctionalGroup.ALDEHYDE, FunctionalGroup.KETONE),
    'ANHY': (FunctionalGroup.ALKENE, FunctionalGroup.HYDROPEROXY),
    'ALOX': (FunctionalGroup.ALKENE, FunctionalGroup.ALKOXY),
    'ALKE': (FunctionalGroup.ALKENE, FunctionalGroup.KETONE),
    'KKOX': (FunctionalGroup.ALKOXY, FunctionalGroup.KETONE),
    'CEHY': (FunctionalGroup.CYCLIC_ETHER, FunctionalGroup.HYDROPEROXY),
    'DHY': (FunctionalGroup.HYDROPEROXY, FunctionalGroup.HYDROPEROXY),
    'KHP': (FunctionalGroup.HYDROPEROXY, FunctionalGroup.KETONE),
    'OOHY': (FunctionalGroup.HYDROPEROXY, FunctionalGroup.PEROXY),
    'DKET': (FunctionalGroup.KETONE, FunctionalGroup.KETONE),
    # One Functional Group
    'ALD': (FunctionalGroup.ALDEHYDE,),
    'ALK': (FunctionalGroup.ALKENE,),
    'CETH': (FunctionalGroup.CYCLIC_ETHER,),
    'QOOH': (FunctionalGroup.HYDROPEROXY,),
    'KET': (FunctionalGroup.KETONE,),
    'RO2': (FunctionalGroup.PEROXY,),
}

# Functional groups to ignore when considering the naming scheme
IGNORE_FGRPS = ('methyl', 'alkane', 'alkoxy_oc',
                'allene', 'propyne', 'allyl')

# Dictionary to remap names to more common ones
NAME_EXCEPTION_DCT = {
    'HO': 'OH'
}


def functional_group_name_dct(mech_spc_dct, rename_rule_dct=None,
                              force_rename=False):
    """ Build a dictionary to map the names of a mechanism according
        to its functional groups and number of carbon atoms.

        The function detects which functional groups are in a species
        and then uses the list of these groups to assign a mechanism
        name according to a rule set bu the user in a supplied dictionary.

        Example:
            OOQOOH is name assigned to species with -OOH and -OO groups
            RO2 is name assigned to species with JUST -OO groups

        The hierarchy of names is determined by the order of the supplied
        dictionary.

        So the final name should correspond to
            <NCarbon>-<FUNC-GRP-RULE-NAME>

        In cases where there are less that 1 carbons or ID'd functional groups,
        no naming rule established

        :param rename_rule_dct: assigns string name to set of func. groups
        :type rename_rule_dct: dct[frozenset(grp1, grp2): str]
    """
    fgrp_map_dct = {}
    for name, dct in mech_spc_dct.items():
        fgrp_map_dct[name] = functional_group_name(
            dct['inchi'], name=name, rename_rule_dct=rename_rule_dct,
            force_rename=force_rename)
    return fgrp_map_dct


def functional_group_name(ich, name='', rename_rule_dct=None,
                          enant_label=True, force_rename=False):
    """ Assign the functional group name

        :param enant_label: Include the enantiomer label?
        :type enant_label: bool
    """

    def _conn_string(ich):
        """ Get the connectivity string
        """
        # OLD SCHEME
        # conn_string = automol.chi.connectivity(
        #     ich, parse_connection_layer=True, parse_h_layer=True)
        # return ioformat.hash_string(
        #   conn_string, 3, remove_char_lst=('-', '_'))

        # NEW SCHEME
        c_conn_str = automol.chi.connectivity(
            ich, parse_connection_layer=True, parse_h_layer=False)
        h_conn_str = automol.chi.connectivity(
            ich, parse_connection_layer=False, parse_h_layer=True)
        chash = ioformat.hash_string(c_conn_str, 3, remove_char_lst=('-', '_'))
        hhash = ioformat.hash_string(h_conn_str, 3, remove_char_lst=('-', '_'))
        lbl = chash + hhash
        return lbl

    def _fgrp_name_string(fgrp_cnt_dct, rule_dct):
        """ Determines what the new name for a species should
            be according the functional groups it has and the
            rule set by the user.

            grps and dct must be sets
        """

        # Build a list of all present func groups, if multiple found
        # ignore certain functional groups
        # put multiple iterations in that list
        # then sort it
        fgrps = []
        for fgrp, count in fgrp_cnt_dct.items():
            if fgrp not in IGNORE_FGRPS:
                for _ in range(count):
                    fgrps.append(fgrp)
        fgrps = tuple(sorted(fgrps))

        # Figure out which renaming rule to use
        if fgrps:
            _name = None
            for fgrp_name, fgrp_lst in rule_dct.items():
                if fgrp_lst == fgrps:
                    _name = fgrp_name
                    break
            _name = _name if _name is not None else ''
        else:
            _name = ''

        return _name

    # Rebuild the rename rule dct so that the functional group list is sorted
    # This is required to do the matching later
    if rename_rule_dct is None:
        rename_rule_dct = copy.deepcopy(DEFAULT_FGRP_RENAME_RULE_DCT)
    rename_rule_dct = {fgrp_name: tuple(sorted(list(fgrp_lst)))
                       for fgrp_name, fgrp_lst in rename_rule_dct.items()}
    print(rename_rule_dct)
    # Get the ich, geom, and gra and other info used for getting name
    gra = automol.chi.graph(ich)
    fml = automol.graph.formula(gra)

    # Get the number of atoms and functional groups
    hvy_atm_cnt = automol.graph.atom_count(gra, heavy_only=True)
    fgrp_cnt_dct = automol.graph.functional_group_count_dct(gra)
    print(fgrp_cnt_dct)

    if name and not force_rename:
        re_name = name
    else:
        # OLD SCHEME
        # # Get the labels
        # conn_lbl = _conn_string(ich)
        # ste_lbl = stereo_name_suffix(ich)
        # fgrp_lbl = _fgrp_name_string(fgrp_cnt_dct, rename_rule_dct)

        # # Build the names string
        # re_name = f'C{c_cnt}'
        # re_name += f'-{conn_lbl}'
        # if ste_lbl:
        #     re_name += f'{ste_lbl}'
        # re_name += f'-{fgrp_lbl}'

        # NEW SCHEME
        fml_lbl_full = automol.form.string(fml, hyd=True)
        fml_lbl_short = automol.form.string(fml, hyd=False)
        conn_lbl = _conn_string(ich)
        ste_lbl = stereo_name_suffix(ich, enant_label=enant_label)
        fgrp_lbl = _fgrp_name_string(fgrp_cnt_dct, rename_rule_dct)
        # Build the names string
        if hvy_atm_cnt == 1:
            re_name = f'{fml_lbl_full}'
        else:
            re_name = f'{fml_lbl_full}{fgrp_lbl}-{conn_lbl}{ste_lbl}'

            if len(re_name) > 16:
                re_name = f'{fml_lbl_full}{conn_lbl}{ste_lbl}'

            if len(re_name) > 16:
                re_name = f'{fml_lbl_short}{conn_lbl}{ste_lbl}'

            if len(re_name) > 16:
                print("WARNING! Name longer than 16 characters: {re_name}")

    # Put in name exception remapping
    re_name = NAME_EXCEPTION_DCT.get(re_name, re_name)

    return re_name


# FORMULA NAME MAPPING
def formula_name_dct(spc_dct):
    """ Build a name mapping dictionary that maps according to
        the stoichiometry of the species
    """

    # Generate the bookkeeping dictionaries to assign names
    fml_cnt_dct = formula_count_dct(spc_dct)
    _ich_name_dct = ich_name_dct(spc_dct)

    # Generate unique name from stoichiometry if species not present
    map_dct = {}
    for name, dct in spc_dct.items():
        _ich = dct['ich']
        if _ich not in _ich_name_dct:
            re_name, fml_cnt_dct = formula_name(_ich, fml_cnt_dct, spc_dct)
            map_dct[name] = re_name

    return map_dct


def formula_count_dct(spc_dct):
    """ get the spc_fml dct

        Loop over the spc dct, obtain the formula
        and use it to build a count of each formula
        in the dictionary
    """

    fml_count_dct = {}
    for dct in spc_dct.values():
        _, fml_count_dct = formula_name(
            dct['inchi'], fml_count_dct, spc_dct)

    return fml_count_dct


def formula_name(ich, fml_cnt_dct, spc_dct):
    """ Generate a unique name for a species with the
        given formula that corresponds to

        formula(N) where formula is the string N
        is the Nth iteration of said formula in the
        dictionary.

        Also, updates the overall formula dictionary
        which contains a count of how many times the
        formula appears in some mechanism.
    """

    fml_str = automol.chi.formula_layer(ich)

    if fml_str in fml_cnt_dct:

        fml_cnt_dct[fml_str] += 1

        fml_cnt = fml_cnt_dct[fml_str]
        name = fml_str + f'({fml_cnt})'

        name += f'({fml_cnt})'

        # If the number is in the dictionary, increase by one
        # happen when you have A(5), A(6), A(9)...A(7) miss throws off count
        while name in spc_dct:
            fml_cnt += 1
            name = fml_str + f'({fml_cnt})'
    else:
        fml_cnt_dct[fml_str] = 1
        name = fml_str

    return name, fml_cnt_dct


# Other helper functions
def rxn_ich_to_name(rxn, spc_dct):
    """ Set a reaction described by inchis to one where it is
        described by mechanism names
    """

    has_inf = False
    if rxn:
        if rxn[0]:
            if not isinstance(rxn[0], str):
                has_inf = True
    _ich_name_dct = ich_name_dct(spc_dct, incl_mult=has_inf, incl_chg=has_inf)
    return (
        tuple(_ich_name_dct[rgt] for rgt in rxn[0]),
        tuple(_ich_name_dct[rgt] for rgt in rxn[1]),
        tuple(rxn[2])
    )


def rxn_name_str(rxn, newline=False):
    """ get a reaction name string
    """
    if newline:
        rstr = ' =\n       '.join((' + '.join(rxn[0]), ' + '.join(rxn[1])))
    else:
        rstr = ' = '.join((' + '.join(rxn[0]), ' + '.join(rxn[1])))
    return rstr


def ich_name_dct(spc_dct, incl_mult=False, incl_chg=False):
    """ get dct[ich] = name
    """
    def set_key(dct):
        key = dct['inchi']
        if incl_chg:
            key = (key, dct['charge'],)
            if incl_mult:
                key += (dct['mult'],)
        elif incl_mult:
            key = (key, dct['mult'],)
        return key
    return {set_key(dct): name for name, dct in spc_dct.items()}


def stereo_name_suffix(ich, enant_label=True):
    """ Parse the stereo from the InChI and write a string describing the
        stereochemistry that is present.

        :param enant_label: Include the enantiomer label?
        :type enant_label: bool
    """

    ste_str = ''

    # Read the stereo chemistry from the InChI string
    bnd_ste_par_dct = automol.chi.bond_stereo_parities(ich)
    atm_ste_par_dct = automol.chi.atom_stereo_parities(ich)

    if bnd_ste_par_dct:
        ste_str += ''.join('E' if p else 'Z'
                           for _, p in sorted(bnd_ste_par_dct.items()))

    if len(atm_ste_par_dct) > 1:
        ste_str += ''.join('B' if p else 'A'
                           for _, p in sorted(atm_ste_par_dct.items()))

    if automol.chi.is_chiral(ich):
        if not enant_label:
            ste_str += '*'
        else:
            if automol.chi.is_inverted_enantiomer(ich):
                ste_str += '1'
            else:
                ste_str += '0'

    return ste_str


def remove_stereo_name_suffix(name):
    """ Removes a stereo suffix from a name

        Assumes name has the format:
        <CCNT><FGRP>-<CLYR><HLYR><STE>
    """

    if '-' in name:
        part1, part2 = name.split('-')
        # remove chars past the 6th char since clyr and hlyr are 3-char hash
        if len(part2) > 6:
            part2 = part2[:6]
        _name = '-'.join((part1, part2))
    else:
        _name = name

    return _name
