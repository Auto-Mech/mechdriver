""" functions operating on the reactions block string
"""

import collections
import itertools
import numpy as np
import pyparsing as pp
from phydat import phycon
from autoreact.params import RxnParams
import autoparse.pattern as app
import autoparse.find as apf
from autoparse import cast as ap_cast
from ioformat import headlined_sections
from ioformat import remove_comment_lines


# gearing up to replace autoparse with pyparsing
PP_ARROW = pp.Combine(pp.Opt("<") + pp.Char("=") + pp.Opt(">"))
PP_THIRD_BODY = pp.Group("(" + pp.Word(pp.printables + "+", excludeChars=")") + ")")
PP_SPECIES_NAME = pp.Combine(
    pp.Word(pp.alphas, excludeChars="+=.<")  # Start with an alphabetical character
    #pp.Char(
    #    pp.printables, excludeChars=(pp.nums + "+=.<Ee")
    #)  # Exclude scientific notation- doesn't work with species starting with E/e
    + pp.ZeroOrMore(
        pp.Char(pp.printables, excludeChars="+=.<(")
        ^ "(" + ~pp.FollowedBy("+")
        ^ "<" + ~pp.FollowedBy("=")
    )
    )
# old, no accomodation for floats # PP_COEFF = pp.Word(pp.nums)
PP_COEFF = pp.Combine(
    pp.Optional(pp.Word(pp.nums) + ".")
    + pp.Word(pp.nums)  # Standard decimal with leading digit
    | pp.Word(pp.nums)  # Whole numbers
    | pp.Word(".") + pp.Word(pp.nums)  # Decimal without leading digit
)
PP_REAGENT = pp.Group(pp.Opt(PP_COEFF) + PP_SPECIES_NAME)
PP_SPECIES_LIST = pp.delimitedList(PP_REAGENT, delim="+")
PP_REACTION_EQUATION = (
    PP_SPECIES_LIST("reactants")
    + pp.Opt(PP_THIRD_BODY("thirdbody"))
    + PP_ARROW
    + PP_SPECIES_LIST("products")
    + pp.Opt(PP_THIRD_BODY("thirdbody"))
)


# Various strings needed to parse the data sections of the Reaction block
CHEMKIN_ARROW = (
    app.maybe(app.escape("<")) + app.escape("=") + app.maybe(app.escape(">"))
)
CHEMKIN_PLUS_EM = app.PLUS + "M"
CHEMKIN_PAREN_PLUS_EM = app.escape("(") + app.PLUS + "M" + app.escape(")")
CHEMKIN_PAREN_PLUS = app.escape("(") + app.PLUS
CHEMKIN_PAREN_CLOSE = app.escape(")")

SPECIES_NAME_PATTERN = (
    r"[^\s=+\-]"
    + app.zero_or_more(
        app.one_of_these(
            [
                app.LETTER,
                app.DIGIT,
                r"[#,()\-_]",
                app.escape("*"),
                app.escape("(+)"),
                app.escape("["),
                app.escape("]"),
            ]
        )
    )
    + app.zero_or_more(app.PLUS)
)
SPECIES_NAME_PATTERN_WCOEFFS = (
    r"\d*\.?\d*"  # Allow for optional numeric coefficients like 0.5, 1.5
    + r"[^\s=+\-]"  # Start of species name (non-whitespace, no plus, no equal signs)
    + app.zero_or_more(
        app.one_of_these(
            [
                app.LETTER,
                app.DIGIT,
                r"[#,()\-_]",
                app.escape("*"),
                app.escape("(+)"),
                app.escape("["),
                app.escape("]"),
            ]
        )
    )
    + app.zero_or_more(app.PLUS)
)

SPECIES_NAMES_PATTERN = app.series(
    app.padded(SPECIES_NAME_PATTERN), app.padded(app.PLUS)
) + app.maybe(app.padded(CHEMKIN_PAREN_PLUS_EM))
SPECIES_NAME_PATTERN_WCOEFFS = app.series(
    app.padded(SPECIES_NAME_PATTERN_WCOEFFS), app.padded(app.PLUS)
) + app.maybe(app.padded(CHEMKIN_PAREN_PLUS_EM))

# REACTION_PATTERN = (SPECIES_NAMES_PATTERN + app.padded(CHEMKIN_ARROW) +
#                    SPECIES_NAMES_PATTERN)
COEFF_PATTERN = app.NUMBER + app.LINESPACES + app.NUMBER + app.LINESPACES + app.NUMBER
COMMENTS_PATTERN = app.escape("!") + app.capturing(app.one_or_more(app.WILDCARD2))

BAD_STRS = ["inf", "INF", "nan"]


def get_rxn_param_dct(block_str, ea_units, a_units):
    """Parses all of the chemical equations and corresponding fitting
    parameters in the reactions block of the mechanism input file
    and subsequently pulls all of the species names and fitting
    parameters from the data string; this information is stored in a list.

    :param block_str: raw string for the entire reactions block
    :type block_str: str
    :param ea_units: units of activation energy
    :type ea_units: str
    :param a_units: units of rate constants; either 'moles' or 'molecules'
    :type a_units: str
    :return rxn_param_dct: dct {rxn1: params1, rxn2: ...}
    :rtype: dict
    """

    rxn_strs = get_rxn_strs(block_str)

    if rxn_strs is not None:
        # Loop over each reaction string, creating a RxnParams object for each
        rxns = []
        params_lst = []
        for rxn_str in rxn_strs:
            rxn = get_rxn_name(rxn_str)
            params = get_params(rxn_str, ea_units, a_units)
            rxns.append(rxn)
            params_lst.append(params)

        # Fix any duplicates
        rxns, params_lst = fix_duplicates(rxns, params_lst)
        # Zip into a dictionary
        rxn_param_dct = dict(zip(rxns, params_lst))

    else:
        rxn_param_dct = None

    return rxn_param_dct


def get_rxn_osclass_dct(block_str):
    """returns the reaction class according to the notation
    !#[REACTIONCLASS][species type][reaction type]
    function looks ugly but couldn't find a better way..
    :param block_str: raw string for the entire reactions block
    :type block_str: str
    :return block_dct: dct {rxn1: str, rxn2: ...}
    :rtype: dict
    """
    block_strs_lst = np.array(block_str.split("\n"))
    block_str_nocmt = apf.remove(COMMENTS_PATTERN, block_str)
    rxn_strs = get_rxn_strs(block_str_nocmt)

    rxn_osclass_dct = {}
    rxnclass_start = np.array(
        [i for i, line in enumerate(block_strs_lst) if "#[REACTIONCLASS]" in line]
    )
    rxnclassend = np.array(
        [i for i, line in enumerate(block_strs_lst) if "#[ENDREACTIONCLASS]" in line]
    )

    speciestype = [
        line.split("[")[2].split("]")[0] for line in block_strs_lst[rxnclass_start]
    ]
    reactiontype = [
        line.split("[")[3].split("]")[0] for line in block_strs_lst[rxnclass_start]
    ]

    for rxn_str in rxn_strs:
        rxn = get_rxn_name(rxn_str)
        rxn_str0 = rxn_str.split("\n")[0]
        idx = [i for i, line in enumerate(block_strs_lst) if rxn_str0 in line][0]
        idx_match = np.intersect1d(
            np.where(idx > rxnclass_start)[0], np.where(idx < rxnclassend)[0]
        )

        if len(idx_match) > 0:
            sptype, rxntype = [speciestype[idx_match[0]], reactiontype[idx_match[0]]]
        else:
            sptype, rxntype = ["", ""]
        rxn_osclass_dct[rxn] = {
            "speciestype": sptype,
            "reactiontype": rxntype,
        }

    return rxn_osclass_dct


def get_rxn_strs_dct(block_str):
    """returns the reaction strings
    function looks ugly but couldn't find a better way..
    :param block_str: raw string for the entire reactions block
    :type block_str: str
    :return block_dct: dct {rxn1: str, rxn2: ...}
    :rtype: dict
    """
    block_strs_lst = block_str.split("\n")
    block_str_nocmt = apf.remove(COMMENTS_PATTERN, block_str)
    rxn_strs = get_rxn_strs(block_str_nocmt)

    rxn_strs_dct = {}
    used_indices = [False] * len(block_strs_lst)  # check matched lines
    for rxn_str in rxn_strs:
        rxn = get_rxn_name(rxn_str)
        rxn_str_orig = []
        for line in rxn_str.split("\n"):
            # skip empty lines
            if line.strip() in ['', 'REACTIONS', 'END']:
                continue
            # skip lines 
            for idx, str_orig in enumerate(block_strs_lst):
                if line in str_orig and not used_indices[idx]:
                    rxn_str_orig.append(str_orig)  # add to original list
                    used_indices[idx] = True  # mark as used
                    break

        if rxn not in rxn_strs_dct.keys():
            rxn_strs_dct[rxn] = ["\n".join(rxn_str_orig)]
        else:
            rxn_strs_dct[rxn] = rxn_strs_dct[rxn] + ["\n".join(rxn_str_orig)]

    return rxn_strs_dct


def get_rxn_cmt_dct(block_str):
    """Parses all reactions inline comments if any
    keys of final dictionary are the same as rxn_param_dct

    :param block_str: raw string for the entire reactions block
    :type block_str: str
    :return cmts_dct: dct {rxn1: str, rxn2: ...}
    :rtype: dict
    """

    rxn_strs_dct = get_rxn_strs_dct(block_str)
    cmts_dct = {}
    # Loop over each reaction string, creating a rxn cmt dct for each
    for rxn, rxn_str in rxn_strs_dct.items():
        rxn_str = "\n".join(rxn_str)
        cmts = get_rxn_cmt(rxn_str)

        if rxn in cmts_dct.keys():
            cmts_dct[rxn]["inline"] = cmts_dct[rxn["inline"]] + cmts
        else:
            cmts_dct[rxn] = {"inline": cmts}

    return cmts_dct


def get_rxn_cmt(rxn_str):
    """get inline comments of a reaction string

    Args:
        rxn_str (str): raw reaction string
        cmt_str (str): associated comments
    """
    rxn_and_cmt = rxn_str.split("\n")[0].split("!")
    if len(rxn_and_cmt) > 1:
        cmt_str = "!" + "!".join(rxn_and_cmt[1:])
    else:
        cmt_str = ""
    return cmt_str


def get_pes_dct(block_str):
    """Parses all of the chemical equations
    and uses special comment line to parse them into PESs

    ! pes.subpes.channel  7.3.35

    :param block_str: raw string for the entire reactions block
    :type block_str: str
    :param ea_units: units of activation energy
    :type ea_units: str
    :param a_units: units of rate constants; either 'moles' or 'molecules'
    :type a_units: str
    :return rxn_param_dct: dct {rxn1: params1, rxn2: ...}
    :rtype: dict
    """

    rxn_strs = get_rxn_strs(block_str)
    rxns = []
    if rxn_strs is not None:
        pes_dct = {}
        for rxn_str in rxn_strs:
            rxn = get_rxn_name(rxn_str)
            if rxn in rxns:
                # channel already sorted, must be a dup entry
                continue
            idx_inf = get_pes_info(rxn_str)
            if idx_inf is not None:
                rxns.append(rxn)#add rxn found- 
                # dup rxns might not be classified and then throw an error
                pes_idx, subpes_idx, chnl_idx = idx_inf
                pes_inf = ("PES", pes_idx, subpes_idx)
                chnl_inf = (chnl_idx, rxn)
                if pes_inf in pes_dct:
                    pes_dct[pes_inf] += (chnl_inf,)
                else:
                    pes_dct[pes_inf] = (chnl_inf,)
            else:
                pes_dct = None
                break
    else:
        pes_dct = None

    return pes_dct


def get_rxn_name(rxn_str):
    """Parses a rxn_str to get the reaction key

    :param rxn_str: raw Chemkin string for a single reaction
    :type rxn_str: str
    :return rxn: tuple describing the reactants, products, and third body
    :rtype: tuple ((rct1, rct2, ...), (prd1, prd2, ...), (third_bod1, ...))
    """
    # split line; remove last 3 items (always rxn params) and re-join
    rxn_str = ' '.join(rxn_str.split('\n')[0].split()[:-3])
    parse_dct = PP_REACTION_EQUATION.parseString(rxn_str).asDict()
    rcts = list(
        itertools.chain(
            *(r if len(r) == 1 else int(r[0]) * [r[1]] for r in parse_dct["reactants"])
        )
    )

    if any(not r[0].isdigit() for r in parse_dct["products"] if len(r) > 1):
        prds = list(
            itertools.chain(
                *(r if len(r) == 1 else [r[0] + r[1]] for r in parse_dct["products"])
            )
        )
    else:
        prds = list(
            itertools.chain(
                *(
                    r if len(r) == 1 else int(r[0]) * [r[1]]
                    for r in parse_dct["products"]
                )
            )
        )

    if "thirdbody" in parse_dct:
        thd_bod = ("".join(parse_dct["thirdbody"]),)
    elif "M" in rcts or "M" in prds:
        assert "M" in rcts and "M" in prds, rxn_str
        rcts.remove("M")
        prds.remove("M")
        thd_bod = ("+M",)
    else:
        thd_bod = (None,)

    rxn = (tuple(rcts), tuple(prds), thd_bod)
    return rxn


def get_params(rxn_str, ea_units, a_units):
    """Extracts a RxnParams object from a str describing a rxn

    :param rxn_str: raw Chemkin string for a single reaction
    :type rxn_str: str
    :param ea_units: units of activation energies
    :type ea_units: string
    :param a_units: units of rate constants; either 'moles' or 'molecules'
    :type a_units: str
    :return params: object describing the rate parameters
    :rtype: autoreact.RxnParams object
    """

    rxn = get_rxn_name(rxn_str)
    param_tuple = (
        high_p(rxn_str, ea_units, a_units),
        low_p(rxn_str, ea_units, a_units),
        troe(rxn_str),
        cheb(rxn_str),
        plog(rxn_str, ea_units, a_units),
        colliders(rxn_str),
    )

    if param_tuple[3] is not None:  # Chebyshev
        cheb_dct = param_tuple[3]
        cheb_dct["one_atm_arr"] = param_tuple[0]  # might be None
        params = RxnParams(cheb_dct=cheb_dct)

    elif param_tuple[4] is not None:  # PLOG
        plog_dct = param_tuple[4]
        params = RxnParams(plog_dct=plog_dct)

    elif param_tuple[2] is not None:  # Troe
        assert (
            param_tuple[0] is not None
        ), f"For {rxn}, Troe params included, highP params absent"
        assert (
            param_tuple[1] is not None
        ), f"For {rxn}, Troe & highP params included, lowP params absent"
        assert rxn[2][0] is not None, (
            f"For {rxn}, Troe, highP, & lowP params included," + " third body absent"
        )
        troe_dct = {}
        troe_dct["highp_arr"] = param_tuple[0]
        troe_dct["lowp_arr"] = param_tuple[1]
        troe_dct["troe_params"] = param_tuple[2]
        troe_dct["collid"] = param_tuple[5]
        params = RxnParams(troe_dct=troe_dct)

    elif param_tuple[1] is not None:  # Lindemann
        assert (
            param_tuple[0] is not None
        ), f"For {rxn}, lowP params included, highP params absent"
        assert (
            rxn[2][0] is not None
        ), f"For {rxn}, highP & lowP params included, third body absent"
        lind_dct = {}
        lind_dct["highp_arr"] = param_tuple[0]
        lind_dct["lowp_arr"] = param_tuple[1]
        lind_dct["collid"] = param_tuple[5]
        params = RxnParams(lind_dct=lind_dct)

    else:  # simple Arrhenius
        assert param_tuple[0] is not None, f"For {rxn}, the highP params absent"
        arr_dct = {}
        arr_dct["arr_tuples"] = param_tuple[0]
        arr_dct["arr_collid"] = param_tuple[5]
        params = RxnParams(arr_dct=arr_dct)
    return params


def get_pes_info(rxn_str):
    """Get PES info"""

    ptt = app.capturing(
        app.INTEGER + app.escape(".") + app.INTEGER + app.escape(".") + app.INTEGER
    )

    pes_inf = None
    for line in rxn_str.splitlines():
        if "#" in line and "pes.subpes.channel" in line:
            cap = apf.first_capture(ptt, rxn_str)
            if cap is not None:
                pes_inf = tuple(int(x) - 1 for x in cap.strip().split("."))
            break

    return pes_inf


def get_rxn_strs(block_str, remove_bad_fits=False):
    """Parses all of the chemical equations and corresponding fitting
    parameters in the reactions block of the mechanism input file
    and stores them in a list.
    NB only works on rxn strs without comments
    
    :param block_str: raw string for the entire reactions block
    :type block_str: str
    :param remove_bad_fits: remove reactions with bad fits
    :type remove_bad_fits: bool
    :return rxn_strs: list of raw strings, one for each reaction
    :rtype: list(str)
    """

    if block_str.strip():
        rxn_strs = headlined_sections(
            string=block_str.strip(), headline_pattern=CHEMKIN_ARROW
        )
        if remove_bad_fits:
            rxn_strs = [
                dstr
                for dstr in rxn_strs
                if not any(string in dstr for string in BAD_STRS)
            ]
        # remove lines if they start with a comment !
        rxn_strs = [remove_comment_lines(dstr, delim_pattern=app.escape("!"))
                    for dstr in rxn_strs if dstr.strip()[0] != "!" and
                    'HTYPE' not in dstr] # ignore R+MOLEC types from CRECK model

    else:
        rxn_strs = None

    return rxn_strs


def fix_duplicates(rxns, params_lst):
    """This function finds any duplicates within the list of rxns. If any are
    found, combines the corresponding RxnParams objects

    :param rxns: all reaction keys
    :type rxns: list
    :param params: all reaction parameters
    :type params: list
    :return unique_rxns: unique reaction keys
    :rtype: list
    :return unique_params: combined reaction params to match unique_rxns
    :rtype: list
    """

    # Get the unique rxns and number of occurrences
    unique_rxns = list(set(rxns))
    dup_counts = collections.Counter(rxns)
    unique_params = []

    # Loop through each unique item, looking for duplicates
    for unique_rxn in unique_rxns:
        # Get the index of the first occurrence in the original list
        prev_match_idx = rxns.index(unique_rxn)
        # Get the params for the first occurence
        params = params_lst[prev_match_idx]
        # Get the number of occurrences of the rxn
        dup_count = dup_counts[unique_rxn]

        # Loop over each duplicate of the reaction; the 1 in range skips the
        # first loop and thus won't do any loops if dup_count == 1
        for _ in range(1, dup_count):
            # Get the current match index by looking just past the previous one
            match_idx = rxns.index(unique_rxn, prev_match_idx + 1)
            # Get the new params to be added to the existing ones
            new_params = params_lst[match_idx]
            # Combine the params with the existing ones
            params.combine_objects(new_params)
            # Update the search index
            prev_match_idx = match_idx

        # Append the combined params to the unique params
        unique_params.append(params)

    return unique_rxns, unique_params


def rct_names(rxn_str):
    """Parses the data string for a reaction in the reactions block
    for the line containing the chemical equation in order to
    read the names of the reactant species.

    :param rxn_str: raw Chemkin string for a single reaction
    :type rxn_str: str
    :return names: names of the reactants
    :rtype: tuple(str)
    """

    pattern = _first_line_pattern(
        rct_ptt=app.capturing(SPECIES_NAMES_PATTERN),
        prd_ptt=SPECIES_NAME_PATTERN_WCOEFFS,
        param_ptt=COEFF_PATTERN,
    )
    string = apf.first_capture(pattern, rxn_str)
    try:
        names = _split_reagent_string(string)
    except TypeError as exc:
        raise TypeError(
            f"Reaction line not formatted correctly:\n{rxn_str}\n"
            f"Check that there are three numbers after the "
            f"rxn equation"
        ) from exc

    return names


def prd_names(rxn_str):
    """Parses the data string for a reaction in the reactions block
    for the line containing the chemical equation in order to
    read the names of the product species.

    :param rxn_str: raw Chemkin string for a single reaction
    :type rxn_str: str
    :return names: names of the products
    :rtype: tuple(str)
    """

    pattern = _first_line_pattern(
        rct_ptt=SPECIES_NAMES_PATTERN,
        prd_ptt=app.capturing(SPECIES_NAMES_PATTERN),
        param_ptt=COEFF_PATTERN,
    )
    string = apf.first_capture(pattern, rxn_str)
    try:
        names = _split_reagent_string(string)
    except TypeError as exc:
        raise TypeError(
            f"Reaction line not formatted correctly:\n{rxn_str}\n"
            f"Check that there are three numbers after the "
            f"rxn equation"
        ) from exc

    return names


def third_body(rxn_str):
    """Parses the data string for a reaction in the reactions block
    for the line containing the chemical equation in order to
    read the names of the third body collider if present

    :param rxn_str: raw Chemkin string for a single reaction
    :type rxn_str: str
    :return trd_body: names of the colliders and corresponding fraction
    :rtype: tuple(str)
    """

    pattern = _first_line_pattern(
        rct_ptt=app.capturing(SPECIES_NAMES_PATTERN),
        prd_ptt=SPECIES_NAME_PATTERN_WCOEFFS,
        param_ptt=app.maybe(COEFF_PATTERN),
    )
    rgt_str = apf.first_capture(pattern, rxn_str)
    rgt_str = apf.remove(app.LINESPACES, rgt_str)
    rgt_split_paren = apf.split(CHEMKIN_PAREN_PLUS, rgt_str)
    rgt_split_plus = apf.split(app.PLUS, rgt_str)

    if len(rgt_split_paren) > 1:
        trd_body = "(+" + apf.split(CHEMKIN_PAREN_CLOSE, rgt_split_paren[1])[0] + ")"

    elif "M" in rgt_split_plus:
        trd_body = "+M"

    else:
        trd_body = None

    trd_body = (trd_body,)

    return trd_body


def high_p(rxn_str, ea_units, a_units):
    """Parses the data string for a reaction in the reactions block
    for the line containing the chemical equation in order to
    read the fitting parameters that are on the same line.

    :param rxn_str: raw Chemkin string for a single reaction
    :type rxn_str: str
    :param ea_units: units of activation energy
    :type ea_units: str
    :param a_units: units of rate constants; either 'moles' or 'molecules'
    :type a_units: str
    :return params: Arrhenius fitting parameters for high-P rates
    :rtype: list(list(float))
    """

    pattern = _first_line_pattern(
        rct_ptt=SPECIES_NAMES_PATTERN,
        prd_ptt=SPECIES_NAME_PATTERN_WCOEFFS,
        param_ptt=app.capturing(COEFF_PATTERN),
    )
    string_lst = apf.all_captures(pattern, rxn_str)
    if string_lst:
        fake_params = []
        for string in string_lst:
            fake_params.append(list(ap_cast(string.split())))
            params = fake_params[0]

        # Convert the units of Ea and A
        ea_conv_factor = get_ea_conv_factor(ea_units)
        a_conv_factor = get_a_conv_factor(rxn_str, a_units)
        params[2] = params[2] * ea_conv_factor
        params[0] = params[0] * a_conv_factor
        params = [params]  # convert to list inside a list
    else:
        params = None

    return params


def low_p(rxn_str, ea_units, a_units):
    """Parses the data string for a reaction in the reactions block
    for a line containing the low-pressure fitting parameters,
    then reads the parameters from this line.

    :param rxn_str: raw Chemkin string for a single reaction
    :type rxn_str: str
    :param ea_units: units of activation energies
    :type ea_units: string
    :param a_units: units of rate constants; either 'moles' or 'molecules'
    :type a_units: str
    :return params: Arrhenius fitting parameters for low-P rates
    :rtype: list(list(float))
    """

    pattern = (
        "LOW"
        + app.zero_or_more(app.SPACE)
        + app.escape("/")
        + app.zero_or_more(app.SPACE)
        + app.capturing(app.NUMBER)
        + app.one_or_more(app.SPACE)
        + app.capturing(app.NUMBER)
        + app.one_or_more(app.SPACE)
        + app.capturing(app.NUMBER)
        + app.zero_or_more(app.SPACE)
        + app.escape("/")
    )
    cap1 = apf.first_capture(pattern, rxn_str)
    if cap1 is not None:
        params = [float(val) for val in cap1]

        # Convert the units of Ea and A
        ea_conv_factor = get_ea_conv_factor(ea_units)
        a_conv_factor = get_a_conv_factor(rxn_str, a_units)
        params[2] = params[2] * ea_conv_factor
        params[0] = params[0] * a_conv_factor
        params = [params]  # convert to list inside a list

    else:
        params = None

    return params


def troe(rxn_str):
    """Parses the data string for a reaction in the reactions block
    for a line containing the Troe fitting parameters,
    then reads the parameters from this line.

    Only gets the 4 Troe-specific parameters: alpha, T***, T*, and T**

    :param rxn_str: raw Chemkin string for a single reaction
    :type rxn_str: str
    :return params: Troe fitting parameters
    :rtype: list(float)
    """

    pattern = (
        "TROE"
        + app.zero_or_more(app.SPACE)
        + app.escape("/")
        + app.zero_or_more(app.SPACE)
        + app.capturing(app.NUMBER)
        + app.one_or_more(app.SPACE)
        + app.capturing(app.NUMBER)
        + app.one_or_more(app.SPACE)
        + app.capturing(app.NUMBER)
        + app.maybe(app.one_or_more(app.SPACE) + app.capturing(app.NUMBER))
        + app.zero_or_more(app.SPACE)
        + app.escape("/")
    )
    cap1 = apf.first_capture(pattern, rxn_str)

    if cap1 is not None:
        params = []
        for val in cap1:
            if val is not None:
                params.append(float(val))
            else:
                params.append(None)
    else:
        params = None

    return params


def cheb(rxn_str):
    """Parses the data string for a reaction in the reactions block
    for the lines containing the Chebyshevs fitting parameters,
    then reads the parameters from these lines.

    :param rxn_str: raw Chemkin string for a single reaction
    :type rxn_str: str
    :return params: Chebyshev fitting parameters
    :rtype: dict[param: value]
    """

    original_rxn_str = rxn_str
    rxn_str = apf.remove(COMMENTS_PATTERN, rxn_str)

    tcheb_pattern = (
        "TCHEB"
        + app.zero_or_more(app.SPACE)
        + app.escape("/")
        + app.zero_or_more(app.SPACE)
        + app.capturing(app.NUMBER)
        + app.one_or_more(app.SPACE)
        + app.capturing(app.NUMBER)
        + app.zero_or_more(app.SPACE)
        + app.escape("/")
    )
    pcheb_pattern = (
        "PCHEB"
        + app.zero_or_more(app.SPACE)
        + app.escape("/")
        + app.zero_or_more(app.SPACE)
        + app.capturing(app.NUMBER)
        + app.one_or_more(app.SPACE)
        + app.capturing(app.NUMBER)
        + app.zero_or_more(app.SPACE)
        + app.escape("/")
    )
    cheb_pattern = (
        app.not_preceded_by(app.one_of_these(["T", "P"]))
        + "CHEB"
        + app.zero_or_more(app.SPACE)
        + app.escape("/")
        + app.capturing(app.one_or_more(app.WILDCARD2))
        + app.escape("/")
    )

    cheb_params_raw = apf.all_captures(cheb_pattern, rxn_str)

    if cheb_params_raw:
        params = {}
        # Get temp and pressure limits or use Chemkin defaults if non-existent
        cheb_temps = apf.first_capture(tcheb_pattern, rxn_str)
        cheb_pressures = apf.first_capture(pcheb_pattern, rxn_str)
        if cheb_temps is None:
            cheb_temps = ("300.00", "2500.00")
            print(
                "No Chebyshev temperature limits specified"
                + " for the below reaction."
                + f" Assuming 300 and 2500 K. \n \n {original_rxn_str}\n"
            )
        if cheb_pressures is None:
            cheb_pressures = ("0.001", "100.00")
            print(
                "No Chebyshev pressure limits specified"
                + " for the below reaction."
                + f" Assuming 0.001 and 100 atm. \n \n {original_rxn_str}\n"
            )

        # Get all the numbers from the CHEB parameters
        cheb_params = []
        for cheb_line in cheb_params_raw:
            cheb_params.extend(cheb_line.split())

        # Get alpha dimensions N and M, which are the first two CHEB entries
        cheb_n = int(float(cheb_params[0]))
        cheb_m = int(float(cheb_params[1]))

        # Start on third value (after N and M) and get all polynomial coeffs
        coeffs = []
        for idx, coeff in enumerate(cheb_params[2:]):
            # extra coefficients are allowed but ignored
            if idx + 1 > (cheb_n * cheb_m):
                break
            coeffs.append(coeff)
        assert len(coeffs) == (cheb_n * cheb_m), (
            f"For the below reaction, there should be {cheb_n*cheb_m}"
            + " Chebyshev polynomial"
            + f" coefficients, but there are only {len(coeffs)}."
            + f" \n \n {original_rxn_str}\n"
        )
        alpha = np.array(list(map(float, coeffs)))

        params["tlim"] = tuple(float(val) for val in cheb_temps)
        params["plim"] = tuple(float(val) for val in cheb_pressures)
        params["alpha"] = alpha.reshape([cheb_n, cheb_m])

    else:
        params = None

    return params


def plog(rxn_str, ea_units, a_units):
    """Parses the data string for a reaction in the reactions block
    for the lines containing the PLOG fitting parameters,
    then reads the parameters from these lines.

    :param rxn_str: raw Chemkin string for a single reaction
    :type rxn_str: str
    :param ea_units: units of activation energies
    :type ea_units: string
    :param a_units: units of rate constants; either 'moles' or 'molecules'
    :type a_units: str
    :return params: PLOG fitting parameters
    :rtype: dict[pressure: params]
    """

    pattern = (
        "PLOG"
        + app.zero_or_more(app.SPACE)
        + app.escape("/")
        + app.zero_or_more(app.SPACE)
        + app.capturing(app.NUMBER)
        + app.one_or_more(app.SPACE)
        + app.capturing(app.NUMBER)
        + app.one_or_more(app.SPACE)
        + app.capturing(app.NUMBER)
        + app.one_or_more(app.SPACE)
        + app.capturing(app.NUMBER)
        + app.zero_or_more(app.SPACE)
        + app.escape("/")
    )
    params_lst = apf.all_captures(pattern, rxn_str)

    # Build dictionary of parameters, indexed by parameter
    if params_lst:

        # Get the Ea and A conversion factors
        ea_conv_factor = get_ea_conv_factor(ea_units)
        a_conv_factor = get_a_conv_factor(rxn_str, a_units)
        params = {}
        for param in params_lst:
            pressure = float(param[0])
            vals = list(map(float, param[1:]))
            vals[2] = vals[2] * ea_conv_factor
            vals[0] = vals[0] * a_conv_factor
            if pressure not in params:
                params[pressure] = [vals]
            else:
                params[pressure].append(vals)  # add duplicate expressions

    else:
        params = None

    return params


def colliders(rxn_str):
    """Parses the data string for a reaction in the reactions block
    for the line containing the names of several bath gases and
    their corresponding collider efficiencies

    :param rxn_str: raw Chemkin string for a single reaction
    :type rxn_str: str
    :return params: collider efficiencies for each bath gas
    :rtype: dict {spc1: eff1, spc2: ...}
    """

    bad_strings = ("DUP", "LOW", "TROE", "CHEB", "PLOG", CHEMKIN_ARROW)
    species_char = app.one_of_these(
        [app.LETTER, app.DIGIT, app.escape("("), app.escape(")"), app.UNDERSCORE]
    )
    species_name = app.one_or_more(species_char)

    # Loop over the lines and search for string with collider facts
    if "LOW" in rxn_str or "TROE" in rxn_str or "M=" in rxn_str or "M =" in rxn_str:
        params = {}
        for line in rxn_str.splitlines():
            if not any(apf.has_match(string, line) for string in bad_strings):
                factor_pattern = (
                    app.capturing(species_name)
                    + app.zero_or_more(app.SPACE)
                    + app.escape("/")
                    + app.zero_or_more(app.SPACE)
                    + app.capturing(app.NUMBER)
                    + app.zero_or_more(app.SPACE)
                    + app.escape("/")
                    + app.zero_or_more(app.SPACE)
                )
                baths = apf.all_captures(factor_pattern, line)
                if baths:
                    for bath in baths:
                        params[bath[0]] = float(bath[1])
        # If nothing was put into the dictionary, set it to None
        if not params:
            params = None
    else:
        params = None

    return params


def _first_line_pattern(rct_ptt, prd_ptt, param_ptt):
    """Defines the pattern for the first line in a reaction data
    string that contains the chemical equation and high-pressure
    fitting parameters for the reaction.

    :param rct_ptt: string pattern for the reactant species
    :type rct_ptt: str
    :param prd_ptt: string pattern for the product species
    :type prd_ptt: str
    :param param_ptt: string pattern for high-pressure parameters
    :type param_ptt: str
    :rtype: str
    """
    return rct_ptt + app.padded(CHEMKIN_ARROW) + prd_ptt + app.LINESPACES + param_ptt


def _split_reagent_string(rgt_str):
    """Parses out the names of all the species given in a string with
    the chemical equation within the reactions block.

    :param rgt_str: string with the reaction chemical equation
    :type rgt_str: str
    :return rgts: names of the species in the reaction
    :type rgts: tuple(str)
    """

    def _interpret_reagent_count(rgt_cnt_str):
        """Count the species in a string containing one side
        of a chemical equation.

        :param rgt_cnt_str: string of one side of chemcial equation
        :type rgt_cnt_str: str
        :return: rgts: names of species from string
        :rtype: tuple(str)
        """
        _pattern = (
            app.STRING_START
            + app.capturing(app.maybe(app.DIGIT))
            + app.capturing(app.one_or_more(app.NONSPACE))
        )
        cnt, rgt = apf.first_capture(_pattern, rgt_cnt_str)
        cnt = int(cnt) if cnt else 1
        rgts = (rgt,) * cnt
        return rgts

    rgt_str = apf.remove(app.LINESPACES, rgt_str)
    rgt_str = apf.remove(CHEMKIN_PAREN_PLUS_EM, rgt_str)
    rgt_split = apf.split(CHEMKIN_PAREN_PLUS, rgt_str)

    if len(rgt_split) > 1:
        rgt_str = rgt_split[0]

    pattern = app.PLUS + app.not_followed_by(app.PLUS)
    rgt_cnt_strs = apf.split(pattern, rgt_str)

    rgts = tuple(itertools.chain(*map(_interpret_reagent_count, rgt_cnt_strs)))
    # remove '+M':
    if "M" in rgts:
        rgts_lst = list(rgts)
        rgts_lst.remove("M")
        rgts = tuple(rgts_lst)

    return rgts


def get_ea_conv_factor(ea_units):
    """Get the factor for converting Ea to the desired units of kcal/mole

    :param ea_units: units of activation energies
    :type ea_units: string
    :return ea_conv_factor: conversion factor from ea_units to cal/mol
    :rtype: float
    """

    if ea_units == "cal/mole":
        ea_conv_factor = 1
    elif ea_units == "kcal/mole":
        ea_conv_factor = phycon.KCAL2CAL
    elif ea_units == "joules/mole":
        ea_conv_factor = phycon.J2CAL
    elif ea_units == "kjoules/mole":
        ea_conv_factor = phycon.KJ2CAL
    elif ea_units == "kelvins":
        ea_conv_factor = phycon.KEL2CAL
    else:
        raise NotImplementedError(
            f"Invalid ea_units: {ea_units}."
            + " Options: 'kcal/mole', 'cal/mole', 'joules/mole',"
            + " 'kjoules/mole', 'kelvins'"
        )

    return ea_conv_factor


def get_a_conv_factor(rxn_str, a_units):
    """Get the factor for converting A to the desired basis of moles

    :param rxn_str: data string for species in reaction block
    :type rxn_str: str
    :param a_units: units of rate constants; either 'moles' or 'molecules'
    :type a_units: str
    :return a_conv_factor: conversion factor from a_units to moles
    :rtype: float
    """

    # Get the molecularity
    rcts = rct_names(rxn_str)
    if not isinstance(rcts, tuple):  # convert to list to avoid mistake
        rcts = [rcts]
    molecularity = len(rcts)

    # Find out whether there is a third body
    trd_body = third_body(rxn_str)[0]
    # if 3rd body has '(', no effect on units
    if trd_body is not None and "(+" not in trd_body:
        molecularity += 1

    if a_units == "moles":
        a_conv_factor = 1
    elif a_units == "molecules":
        a_conv_factor = phycon.NAVO ** (molecularity - 1)
    else:
        raise NotImplementedError(
            f"Invalid a_units: {a_units}. Options: 'moles' or 'molecules'"
        )

    return a_conv_factor


# Not used here, but called elsewhere
def ratek_fit_info(rxn_str):
    """Read the information describing features of the fits to the
    rate constants
    """

    # Read the temperatures and the Errors from the lines
    pressure_ptt = (
        "Pressure:" + app.SPACES + app.capturing(app.one_of_these([app.FLOAT, "High"]))
    )
    trange_ptt = (
        "Temps: "
        + app.SPACES
        + app.capturing(app.INTEGER)
        + "-"
        + app.capturing(app.INTEGER)
        + app.SPACES
        + "K"
    )
    mean_ptt = (
        "MeanAbsErr:" + app.SPACES + app.capturing(app.FLOAT) + app.escape("%") + ","
    )
    max_ptt = "MaxErr:" + app.SPACES + app.capturing(app.FLOAT) + app.escape("%")
    pressure_caps = apf.all_captures(pressure_ptt, rxn_str)
    trange_caps = apf.all_captures(trange_ptt, rxn_str)
    mean_caps = apf.all_captures(mean_ptt, rxn_str)
    max_caps = apf.all_captures(max_ptt, rxn_str)

    pressures = []
    for pressure in pressure_caps:
        if pressure != "High":
            pressures.append(float(pressure))
        elif pressure == "High":
            pressures.append(pressure)
    trange_vals = []
    for cap in trange_caps:
        temp1, temp2 = cap
        trange_vals.append([int(temp1), int(temp2)])
    if mean_caps is not None:
        mean_vals = [float(val) for val in mean_caps]
    else:
        mean_vals = []
    if max_caps is not None:
        max_vals = [float(val) for val in max_caps]
    else:
        max_vals = []

    # Build the inf_dct
    inf_dct = {}
    for idx, pressure in enumerate(pressures):
        inf_dct[pressure] = {"temps": trange_vals[idx]}
        if mean_vals:
            inf_dct[pressure].update({"mean_err": mean_vals[idx]})
        if max_vals:
            inf_dct[pressure].update({"max_err": max_vals[idx]})

    return inf_dct


if __name__ == "__main__":
    RXN_STR = "I3C6OOH2-1O2 = I3C6Q12-6"
    RXN_STR = "C7H121OOH6-7 = C7H121-6 + HO2"
    RXN_STR = "C7H121-6 + HO2 = C7H121OOH6-7"
    RXN_STR = "C7H121-6 + HO2 (+ M) = C7H121OOH6-7 (+M)"
    RXN_STR = "2OH = H2O2"

    print(get_rxn_name("I3C6OOH2-1O2 = I3C6Q12-6"))
    print(get_rxn_name("C7H121OOH6-7 = C7H121-6 + HO2"))
    print(get_rxn_name("C7H121-6 + HO2 = C7H121OOH6-7"))
    print(get_rxn_name("C7H121-6 + HO2 (+ M) = C7H121OOH6-7 (+M)"))
    print(get_rxn_name("C7H121-6 + HO2 (+ M) = C7H121OOH6-7 ( + M )"))
    print(get_rxn_name("C7H121-6 + HO2 + M = C7H121OOH6-7 +M"))
    print(get_rxn_name("2OH = H2O2"))
    print(get_rxn_name("IC6D2 + O = IC6D1-3 + OH"))
    print(get_rxn_name("H2+M<=>H+H+M"))
    print(get_rxn_name("H2(+M)<=>H+H(+M)"))
    print(get_rxn_name("C(O)CC(+M)<=>H+H(+M)"))
    print(get_rxn_name("C3H2(S)+M<=>C3H2+M"))

    # print(PP_SPECIES_NAME.parseString('I3C6OOH2-1O2 = I3C6Q12-6').asList())
    # print(PP_SPECIES_NAME.parseString('C(O)CC(+M)<=>H+H(+M)').asList())
