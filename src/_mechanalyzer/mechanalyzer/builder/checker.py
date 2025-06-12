""" Check mechanisms for various flaws/inconsistencies
"""

from collections import Counter as counter
from chemkin_io.writer._util import format_rxn_name


def run_all_checks(rxn_param_dct, rxn_ktp_dct, k_thresholds,
                   rxn_num_threshold):
    """ Run all mechanism checks and output a string describing the results.

        :param rxn_param_dct: rate constant parameters for a mechanism
        :type rxn_param_dct: dct
            {rxn1: (param_tuple1, param_tuple2, ...), rxn2: ...}
        :param rxn_ktp_dct: rate constant values for a mechanism
        :type rxn_ktp_dct: dct {rxn1: ktp_dct1, rxn2: ...}
        :param k_thresholds: rate constant thresholds
            for uni-, bi-, and ter-molecular reactions
        :type k_thresholds: list [float, float, float]
        :param rxn_num_threshold: # of reactions at and below
            which a species is considered "lone"
        :type rxn_num_threshold: int
        :return total_str: description of all the checks performed
        :rtype: str
    """

    def separator():
        """ Write a line of '+' to separate sections
        """
        return '\n' + '+' * 100 + '\n'

    total_str = separator()

    # Large rate constants
    large_rxn_ktp_dcts = get_large_kts(rxn_ktp_dct, k_thresholds)
    total_str += write_large_kts(large_rxn_ktp_dcts, k_thresholds)
    total_str += separator()

    # Negative rate constants
    negative_rxn_ktp_dct = get_negative_kts(rxn_ktp_dct)
    total_str += write_negative_kts(negative_rxn_ktp_dct)
    total_str += separator()

    # Currently doesn't work due to new RxnParams class
    # Duplicate (more than 2) reactions
    #duplicate_rxns = get_duplicates(rxn_param_dct)
    #total_str += write_duplicates(duplicate_rxns)
    #total_str += separator()

    # Mismatching rate expressions
    mismatched_rxns = get_mismatches(rxn_param_dct)
    total_str += write_mismatches(mismatched_rxns)
    total_str += separator()

    # Lone species
    lone_spcs = get_lone_spcs(rxn_param_dct, rxn_num_threshold)
    total_str += write_lone_spcs(lone_spcs, rxn_num_threshold)
    total_str += separator()

    # Sources and sinks
    source_spcs, sink_spcs = get_sources_and_sinks(rxn_param_dct)
    total_str += write_sources_and_sinks(source_spcs, sink_spcs)
    total_str += separator()

    return total_str


def get_sources_and_sinks(rxn_param_dct):
    """ Get species that only appear as reactants (sources) or
        only appear as products (sinks). Output sources and sinks and
        the reaction keys for all reactions in which each species
        participates.

        :param rxn_param_dct: rate constant parameters for a mechanism
        :type rxn_param_dct: dct
            {rxn1: (param_tuple1, param_tuple2, ...), rxn2: ...}
        :return source_species: species that only appear
            as reactants and associated reactions
        :rtype: dct {spc1: [rxn1, rxn2, ...], spc2: ...}
        :return sink_species: species that only appear
            as products and associated reactions
        :rtype: dct {spc1: [rxn1, rxn2, ...], spc2: ...}
    """

    # Get lists of source and sink species
    # Sorted to make the writer test succeed
    all_rcts, all_prds = _get_rcts_prds(rxn_param_dct)
    sources = sorted(list(set(all_rcts) - set(all_prds)))
    sinks = sorted(list(set(all_prds) - set(all_rcts)))

    # Store the reaction name(s) for each source or sink species
    source_spcs = {spc: [] for spc in sources}
    sink_spcs = {spc: [] for spc in sinks}
    for spc in sources:
        for rxn in rxn_param_dct.keys():
            if spc in rxn[0]:
                source_spcs[spc].append(rxn)
    for spc in sinks:
        for rxn in rxn_param_dct.keys():
            if spc in rxn[1]:
                sink_spcs[spc].append(rxn)

    return source_spcs, sink_spcs


def get_large_kts(rxn_ktp_dct, thresholds):
    """ Get the reactions and k(T,P) values for rate constants
        that exceed molecularity-specific thresholds

        :param rxn_ktp_dct: rate constant values for a mechanism
        :type rxn_ktp_dct: dct {rxn1: ktp_dct1, rxn2: ...}
        :param thresholds: rate constant thresholds for
            uni-, bi-, and ter-molecular reactions
        :type thresholds: list [float, float, float]
        :return large_rxn_ktp_dcts: list of rxn_ktp_dcts containing
            reactions that exceed the
            unimolecular, bimolecular, and termolecular thresholds
        :rtype: list(dict[], dict[], dict[])
    """

    unimolec_rxn_ktp_dct = {}
    bimolec_rxn_ktp_dct = {}
    termolec_rxn_ktp_dct = {}
    for rxn, ktp_dct in rxn_ktp_dct.items():
        molecularity = get_molecularity(rxn)
        unimolec_ktp_dct = {}
        bimolec_ktp_dct = {}
        termolec_ktp_dct = {}
        for pressure, (temps, kts) in ktp_dct.items():
            if max(kts) > thresholds[0] and molecularity == 1:
                unimolec_ktp_dct[pressure] = (temps, kts)
            if max(kts) > thresholds[1] and molecularity == 2:
                bimolec_ktp_dct[pressure] = (temps, kts)
            if max(kts) > thresholds[2] and molecularity == 3:
                termolec_ktp_dct[pressure] = (temps, kts)
        # Store the dictionaries if they contain anything
        if unimolec_ktp_dct != {}:
            unimolec_rxn_ktp_dct[rxn] = unimolec_ktp_dct
        if bimolec_ktp_dct != {}:
            bimolec_rxn_ktp_dct[rxn] = bimolec_ktp_dct
        if termolec_ktp_dct != {}:
            termolec_rxn_ktp_dct[rxn] = termolec_ktp_dct

    large_rxn_ktp_dcts = [unimolec_rxn_ktp_dct, bimolec_rxn_ktp_dct,
                          termolec_rxn_ktp_dct]

    return large_rxn_ktp_dcts


def get_negative_kts(rxn_ktp_dct):
    """ Finds any pressures within a mechanism that contain
        negative rate constant and outputs all rate constants at that pressure
        (including any neighboring positive values).

        :param rxn_ktp_dct: rate constant values for a mechanism
        :type rxn_ktp_dct: dct {rxn1: ktp_dct1, rxn2: ...}
        :return negative_rxn_ktp_dct: k(T)s at every pressure with
            one or more negative k(T)
        :rtype: dct {rxn1: ktp_dct1, rxn2: ...}
    """

    negative_rxn_ktp_dct = {}
    for rxn, ktp_dct in rxn_ktp_dct.items():
        negative_ktp_dct = {}
        for pressure, (temps, kts) in ktp_dct.items():
            if min(kts) < 0:
                negative_ktp_dct[pressure] = (temps, kts)
        # Store the ktp_dct if it contains anything
        if negative_ktp_dct != {}:
            negative_rxn_ktp_dct[rxn] = negative_ktp_dct

    return negative_rxn_ktp_dct


def get_lone_spcs(rxn_param_dct, threshold):
    """ Get species that are considered "lone" species based on only
        being included in a small number of reactions
        (the cutoff for which is set by threshold).

        Output each lone species and the reaction keys for
        all reactions in which each lone species participates.

        :param rxn_param_dct: rate constant parameters for a mechanism
        :type rxn_param_dct: dct
            {rxn1: (param_tuple1, param_tuple2, ...), rxn2: ...}
        :param threshold: number of reactions at and below which
            a species is considered "lone"
        :type threshold: int
        :return lone_spcs: dictionary containing
            each lone species and its reactions
        :rtype: dct {lone_spc1: [rxn1, rxn2, ...], lone_spc2: ...}

    """

    # Get the number of reactions that each species participates in
    all_rcts, all_prds = _get_rcts_prds(rxn_param_dct)
    spc_counts = dict(counter(all_rcts + all_prds))

    # Filter by the threshold
    lone_spcs = {}
    for spc, count in spc_counts.items():
        if count <= threshold:
            lone_spcs[spc] = []

    # Store the reaction name(s) for each lone species
    for spc, dct in lone_spcs.items():
        for rxn in rxn_param_dct.keys():
            if spc in rxn[0] or spc in rxn[1]:  # whether spc in prds or rcts
                dct.append(rxn)

    return lone_spcs


def get_duplicates(rxn_param_dct):
    """ Get reactions that have more than 2 rate expressions

        :param rxn_param_dct: rate constant parameters for a mechanism
        :type rxn_param_dct: dct
            {rxn1: (param_tuple1, param_tuple2, ...), rxn2: ...}
        :return duplicate_rxns: duplicate reactions and number of expressions
        :rtype: dct {rxn1: num_of_expressions1, rxn2: ...}
    """

    duplicate_rxns = {}
    for rxn, params in rxn_param_dct.items():
        
        if len(params) > 2:
            duplicate_rxns[rxn] = len(params)

    return duplicate_rxns


def get_mismatches(rxn_param_dct):
    """ Get reactions that have mismatching rate expressions
        (e.g., PLOG and Arrhenius).

        :param rxn_param_dct: rate constant parameters for a mechanism
        :type rxn_param_dct: dct
            {rxn1: (param_tuple1, param_tuple2, ...), rxn2: ...}
        :return mismatched_rxns: list of mismatched reactions with
            params and types of expressions
        :rtype: dct {rxn1: ((param_tuple1, param_tuple2, ...),
                            [type1, type2, ...], rxn2: ...}
    """

    def classify_params(param):
        """ Classify a rate based on the contents of the param
            
            Refactoring the RxnParams class should simplify this...
        """
        param_types = []
        if param.cheb is not None:
            param_types.append('Chebyshev')
        if param.plog is not None:
            param_types.append('PLOG')
        if param.troe is not None:
            param_types.append('Troe')
        if param.lind is not None:
            param_types.append('Lindemann')
        if param.arr is not None:
            param_types.append('Arrhenius')

        return param_types

    mismatched_rxns = {}
    for rxn, params in rxn_param_dct.items():
        param_types = classify_params(params)
        # If the reaction types are not all identical, save the rxn
        if len(set(param_types)) != 1:
            mismatched_rxns[rxn] = (params, param_types)

    return mismatched_rxns


def get_missing_spcs(rxn_param_dct, spc_dct):
    """ Compares a rxn_param_dct and a spc_dct to find missing species

        :param rxn_param_dct: rate constant parameters for a mechanism
        :type rxn_param_dct: dct
            {rxn1: (param_tuple1, param_tuple2, ...), rxn2: ...}
        :param spc_dct: info on species
        :type spc_dct: dct {spc1: info_dct1, spc2: ...}
        :return missing_from_csv: list of species missing from the spc_csv
        :rtype: list [spc1, spc2, ...]
        :return missing_from_mech: list of species missing from the mechanism
        :rtype: list [spc1, spc2, ...]
    """

    def strip_thrd_bod(thrd_bod, rxn):
        """ Strip the third body notation from a species name
        """

        if thrd_bod in ('(+M)', '+M'):
            stripped_thrd_bod = None
        elif thrd_bod[0] == '(':
            stripped_thrd_bod = thrd_bod[2:-1]
        elif thrd_bod[0] == '+':
            stripped_thrd_bod = thrd_bod[1:]
        else:
            stripped_thrd_bod = None
            print(f'The third body could not be read for the reaction {rxn}')

        return stripped_thrd_bod

    # Get the mechanism species
    mech_spcs = []
    for rxn in rxn_param_dct.keys():
        rcts, prds, thrd_bods = rxn
        for spc in rcts:
            mech_spcs.append(spc)
        for spc in prds:
            mech_spcs.append(spc)
        for thrd_bod in thrd_bods:
            if thrd_bod:  # if it's not None, try to read it
                stripped_thrd_bod = strip_thrd_bod(thrd_bod, rxn)
                if stripped_thrd_bod:  # will be None for '+M' or '(+M)'
                    mech_spcs.append(stripped_thrd_bod)

    mech_spcs = set(mech_spcs)

    # Get the spc_dct spcs
    csv_spcs = set(spc_dct.keys())

    # Get the differences between the two sets
    missing_from_csv = list(mech_spcs - csv_spcs)
    missing_from_mech = list(csv_spcs - mech_spcs)

    return missing_from_csv, missing_from_mech


def write_sources_and_sinks(source_spcs, sink_spcs):
    """ Write any species that occur as only sources or sinks,
        along with their associated reactions, to a string.

        :param source_spcs: species that only appear as
            reactants and associated reactions
        :type source_spcs: dct {spc1: [rxn1, rxn2, ...], spc2: ...}
        :param sink_spcs: species that only appear as
            products and associated reactions
        :type sink_spcs: dct {spc1: [rxn1, rxn2, ...], spc2: ...}
        :return source_sink_str: string with info on source and sink species
        :rtype: str
    """

    def write_dct(spc_dct, max_spc_len, buffer=7):
        """ Write either a source_spcs or sink_spcs dct to a string
        """
        new_str = 'Species' + ' ' * (max_spc_len - 7 + buffer) + 'Reactions\n'
        for spc, rxns in spc_dct.items():
            new_str += ('{0:<' + str(max_spc_len + buffer) + 's}').format(spc)
            for rxn_idx, rxn in enumerate(rxns):
                if rxn_idx != 0:  # don't add comma/space before first rxn
                    new_str += ', '
                new_str += format_rxn_name(rxn)
            new_str += '\n'
        new_str += '\n'

        return new_str

    source_sink_str = (
        '\nSOURCE AND SINK SPECIES\n\n' +
        'These species only appear as reactants:\n\n'
    )

    if source_spcs:
        longest_source_name = max(map(len, list(source_spcs.keys())))
        source_sink_str += write_dct(source_spcs, longest_source_name)
    else:  # if source_spcs is empty
        source_sink_str += 'No source species found\n\n'
    source_sink_str += 'These species only appear as products:\n\n'
    if sink_spcs:
        longest_sink_name = max(map(len, list(sink_spcs.keys())))
        source_sink_str += write_dct(sink_spcs, longest_sink_name)
    else:  # if sink_spcs is empty
        source_sink_str += 'No sink species found\n\n'
    source_sink_str += '\n'

    return source_sink_str


def write_large_kts(large_rxn_ktp_dcts, thresholds):
    """ Write the reactions and k(T,P) values for rate constants
        that exceed molecularity-specific thresholds.

        :param large_rxn_ktp_dcts: list of rxn_ktp_dcts containing
            reactions that exceed the
            unimolecular, bimolecular, and termolecular thresholds
        :type: lst [dict[], dict[], dict[]]
        :param thresholds: rate constant thresholds for
            uni-, bi-, and ter-molecular reactions
        :type thresholds: list [float, float, float]
        :return large_kts_str: string describing the large rate constants
        :rtype: str
    """

    unimolec_dct, bimolec_dct, termolec_dct = large_rxn_ktp_dcts

    # Header
    large_kts_str = '\nLARGE RATE CONSTANTS\n\n'
    large_kts_str += (
        f'Unimolecular threshold: {thresholds[0]:0.1E} s^-1\n' +
        f'Bimolecular threshold: {thresholds[1]:0.1E} cm^3 mol^-1 s^-1\n' +
        f'Termolecular threshold: {thresholds[2]:0.1E} cm^6 mol^-2 s^-1\n\n'
    )

    # Unimolecular
    if not unimolec_dct:
        large_kts_str += (
            'No unimolecular reactions exceed ' +
            f'{thresholds[0]:0.1E} s^-1\n\n')
    else:
        large_kts_str += (
            'Unimolecular rate constants that exceed ' +
            f'{thresholds[0]:0.1E} s^-1\n\n')
        large_kts_str += _write_rxn_ktp_dct(unimolec_dct)

    # Bimolecular
    if not bimolec_dct:
        large_kts_str += (
            'No bimolecular reactions exceed ' +
            f'{thresholds[1]:0.1E} cm^3 mol^-1 s^-1\n\n')
    else:
        large_kts_str += (
            'Bimolecular rate constants that exceed ' +
            f'{thresholds[1]:0.1E} cm^3 mol^-1 s^-1\n\n'
        )
        large_kts_str += _write_rxn_ktp_dct(bimolec_dct)

    # Termolecular
    if not termolec_dct:
        large_kts_str += (
            'No termolecular reactions exceed ' +
            f'{thresholds[2]:0.1E} cm^6 mol^-2 s^-1\n\n'
        )
    else:
        large_kts_str += (
            'Termolecular rate constants that exceed ' +
            f'{thresholds[2]:0.1E} cm^6 mol^-2 s^-1\n\n'
        )
        large_kts_str += _write_rxn_ktp_dct(bimolec_dct)

    return large_kts_str


def write_negative_kts(negative_rxn_ktp_dct):
    """ Write the rxn_ktp_dct containing negative rate constants to a string

    :param negative_rxn_ktp_dct: a dictionary containing
        any ktp_dcts with negative values
    :type negative_rxn_ktp_dct: dct {rxn1: ktp_dct1, rxn2: ...}
    :return negative_kts_str: string describing the negative_rxn_ktp_dct
    :rtype: str
    """

    negative_kts_str = '\nNEGATIVE RATE CONSTANTS\n\n'
    if negative_rxn_ktp_dct == {}:
        negative_kts_str += 'No reactions have negative rate constants\n\n\n'
    else:
        negative_kts_str += _write_rxn_ktp_dct(negative_rxn_ktp_dct)

    return negative_kts_str


def write_lone_spcs(lone_spcs, threshold):
    """ Write lone species and reactions in which they participate to a string.

        :param lone_spcs: dictionary containing
            each lone species and its reactions
        :type: dct {lone_spc1: [rxn1, rxn2, ...], lone_spc2: ...}
        :param threshold: number of reactions at and below which
            a species is considered "lone"
        :type threshold: int
        :return lone_spcs_str: string with lone species and their reactions
        :rtype: str
    """
    lone_spcs_str = (
        f'\nLONE SPECIES\n\nThese species appear in {threshold} ' +
        'or less reactions\n\n'
    )
    if lone_spcs:
        max_spc_len = max(map(len, list(lone_spcs.keys())))  # longest spc name
        buffer = 5
        lone_spcs_str += (
            'Species' + ' ' * (max_spc_len - 7 + buffer) + 'Reactions\n')
        for spc, rxns in lone_spcs.items():
            lone_spcs_str += (
                '{0:<' + str(max_spc_len + buffer) + 's}').format(spc)
            for rxn_idx, rxn in enumerate(rxns):
                if rxn_idx != 0:  # don't add comma/space before first rxn
                    lone_spcs_str += ', '
                lone_spcs_str += format_rxn_name(rxn)
            lone_spcs_str += '\n'
        lone_spcs_str += '\n\n'
    else:
        lone_spcs_str += 'No lone species found\n\n\n'

    return lone_spcs_str


def write_duplicates(duplicate_rxns):
    """ Write reactions with more than 2 duplicate expressions to a string

        :param duplicate_rxns: duplicate reactions and number of expressions
        :type: dct {rxn1: num_of_expressions1, rxn2: ...}
        :return dups_str: description of the duplicate reactions
        :rtype: str
    """

    dups_str = (
        '\nDUPLICATE REACTIONS\n\n' +
        'These reactions have more than 2 rate expressions:\n' +
        '(Number of rate expressions given in parentheses)\n\n'
    )

    if duplicate_rxns != {}:
        for rxn, num_dups in duplicate_rxns.items():
            rxn_name = format_rxn_name(rxn)
            dups_str += f'{rxn_name}     ({num_dups})\n'
    else:
        dups_str += 'No reactions with more than 2 expressions found\n'
    dups_str += '\n\n'

    return dups_str


def write_mismatches(mismatched_rxns):
    """ Write reactions with mismatching rate expressions to a string

        :param mismatched_rxns: list of mismatched reactions with
            params and types of expressions
        :type: dct {rxn1: ((param_tuple1, param_tuple2, ...),
                           [type1, type2, ...], rxn2: ...}
        :return mismatch_str: description of the mismatching reactions
        :rtype: str
    """

    mismatch_str = '\nMISMATCHED REACTIONS\n\n'

    if mismatched_rxns != {}:
        mismatch_str += (
            'The following reactions have mismatched rate expressions\n')
        for rxn, (_, param_types) in mismatched_rxns.items():
            rxn_name = format_rxn_name(rxn)
            mismatch_str += rxn_name + ': '
            for type_idx, rxn_type in enumerate(param_types):
                if type_idx != 0:
                    mismatch_str += ', '
                mismatch_str += rxn_type
            mismatch_str += '\n'
        mismatch_str += '\n\n'
    else:
        mismatch_str += (
            'No reactions with mismatching rate expressions found\n\n\n')

    return mismatch_str


def write_missing_spcs(missing_from_csv, missing_from_mech):
    """ Write a string saying what species are missing from the
        species.csv or mechanism.dat file.
    """

    missing_spc_str = '\nSPECIES MISSING FROM CSV OR MECHANISM\n\n'
    missing_spc_str += 'These species are missing from the spc.csv file:\n'
    for spc in missing_from_csv:
        missing_spc_str += spc + '\n'
    missing_spc_str += '\nThese species are missing from the mechanism file:\n'
    for spc in missing_from_mech:
        missing_spc_str += spc + '\n'
    missing_spc_str += '\n\n'

    return missing_spc_str


def _write_rxn_ktp_dct(rxn_ktp_dct):
    """ Write a rxn_ktp_dct in an easily readable string

    :param rxn_ktp_dct: dictionary containing k(T,P) values for reactions
    :type rxn_ktp_dct: dct {rxn1: ktp_dct1, rxn2: ...}
    :return output_str: string describing the rxn_ktp_dct
    :rtype: str
    """

    output_str = ''
    for rxn, ktp_dct in rxn_ktp_dct.items():
        output_str += format_rxn_name(rxn)
        for pressure, (temps, kts) in ktp_dct.items():
            output_str += f'\nPressure: {pressure} atm\n'
            output_str += '    Temperature (K)\n    '
            for temp in temps:
                output_str += (f'{temp:<12.1f}')
            output_str += '\n    Rate constant\n    '
            for rate in kts:
                output_str += (f'{rate:<12.3E}')
        output_str += '\n\n\n'

    return output_str


def _get_rcts_prds(rxn_param_dct):
    """ Gets lists of all reactants and products in a mechanism
        (both lists may contain duplicates)

    :param rxn_param_dct: rate constant parameters for a mechanism
    :return all_rcts: all reactants in the mechanism
    :rtype: list [rct1, rct2, ...]
    :return all_prds: all products in the mechanism
    :rtype: list [prd1, prd2, ...]
    """

    all_rcts = []
    all_prds = []
    for rxn in rxn_param_dct.keys():
        rcts, prds, _ = rxn
        for rct in rcts:
            all_rcts.append(rct)
        for prd in prds:
            all_prds.append(prd)

    return all_rcts, all_prds


def get_molecularity(rxn):
    """ Get the molecularity of a reaction

        :param rxn: tuple describing reactants, products, and third body
        :type rxn: tuple ((rct1, rct2, ...), (prd1, prd2, ...), (thirdbod,))
        :return molecularity: molecularity of reaction; 1=unimol, 2=bimol, etc.
        :rtype: int
    """

    rcts, _, third_bods = rxn
    num_rcts = len(rcts)
    third_bod = third_bods[0]
    if third_bod is not None:
        # Check if the third body counts toward molecularity
        if third_bod[0:2] == '(+':
            molecularity = num_rcts
        else:
            molecularity = num_rcts + 1
    # If no third body
    else:
        molecularity = num_rcts

    return molecularity
