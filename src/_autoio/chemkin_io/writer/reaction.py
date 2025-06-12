"""
Writes Chemkin-formatted strings containing the rate parameters
"""

from chemkin_io.writer import _util as util

BUFFER = 10  # buffer between longest reaction name and Arr params
INLINE_BUFFER = 2  # buffer between Arr params and inline comment
INDENT = 2  # indent for things like TROE, PLOG


def replace_rxn_in_cki(ckin_str, rxn_strs_dct, rxn_param_dct_new, rxn_cmts_dct_new={}):
    """Replace the reactions of ckin_str with parameters identified in rxn_strs_dct
        with the new parameters assigned in rxn_param_dct_new

    Args:
        ckin_str (str): original chemkin string to edit
        rxn_strs_dct (dict {rxn: strings}): strings of each reaction in the original file
        rxn_param_dct_new (dict {rxn: params}): reaction parameter dictionary
        rxn_cmts_dct_new (dict {rxn: str}, optional): comments dictionary
    """
    # Get the length of the longest reaction name
    max_len = util.max_rxn_length(rxn_param_dct_new)
    n_of_rxns_replaced = len(rxn_param_dct_new)

    for rxn, params in rxn_param_dct_new.items():
        # ckin string of single new reaction: rxn keys must be from the whole dct
        cmts_dct = rxn_cmts_dct_new.get(rxn)
        ckin_str_rxn = single_rxn(
            rxn, params, cmts_dct=cmts_dct, max_len=max_len, rxnkeys=rxn_strs_dct.keys()
        )
        if rxn in rxn_strs_dct.keys():
            # get the position of the first occurrence of the reaction
            rxn_strs = rxn_strs_dct[rxn]
            ckin_str_lst = ckin_str.split("\n")
            first_rxn_line = rxn_strs[0].split("\n")[0]
            idx = ckin_str_lst.index(first_rxn_line)
            connector = (
                "=>" if "=>" in first_rxn_line and "<=>" not in first_rxn_line else "="
            )
            # delete strings from initial file
            for rxn_str in rxn_strs:
                ckin_str = ckin_str.replace(rxn_str, "")

            # convert connector to => if necessary and print warning
            if connector == "=>" and "=>" not in ckin_str_rxn:
                ckin_str_rxn = ckin_str_rxn.replace("=", "=>")
                print(
                    "*Warning: replaced = connector with => in reaction {} \n \
                        according to the connector in the original mechanism. \n \
                        check !!!".format(
                        first_rxn_line
                    )
                )
            # add new strings
            ckin_str_lst = ckin_str.split("\n")
            ckin_str_lst.insert(idx, ckin_str_rxn)
            # reconstruct string
            ckin_str = "\n".join(ckin_str_lst)
        else:  # add the reaction if not present
            print(
                "Using replace function, but rxn {} \
                was not present before- check".format(
                    rxn
                )
            )
            ckin_str += "\n" + ckin_str_rxn

    # replace additional '\n'
    ckin_str = ckin_str.replace("\n" * n_of_rxns_replaced, "\n")

    return ckin_str


def write_rxn_param_dct(rxn_param_dct, rxn_cmts_dct=None, sortrxns=False):
    """Write all reactions in a rxn_param_dct to a Chemkin string

    :param rxn_param_dct: fitting parameters for all reactions
    :type rxn_param_dct: dict {rxn: params}
    :param rxn_cmts_dct: comments for all reactions
    :type rxn_cmts_dct: dict {rxn: cmts_dct}
    :param sortrxns: reorder dict keys alphabetically
                    (otherwise, order may change every time if other operations have been done before)
    :type sortrxns: bool
    """
    if sortrxns:
        rxn_param_dct = dict(sorted(rxn_param_dct.items()))

    rxn_cmts_dct = rxn_cmts_dct or {}  # sets to empty dict if None

    # Get the length of the longest reaction name
    max_len = util.max_rxn_length(rxn_param_dct)

    # Write each reaction
    ckin_str = ""
    for rxn, params in rxn_param_dct.items():
        cmts_dct = rxn_cmts_dct.get(rxn)
        sing_rxn_str = single_rxn(
            rxn,
            params,
            cmts_dct=cmts_dct,
            max_len=max_len,
            rxnkeys=rxn_param_dct.keys(),
        )
        ckin_str += sing_rxn_str

    return ckin_str


def single_rxn(rxn, params, cmts_dct=None, max_len=45, rxnkeys=[]):
    """Writes a reaction string from a RxnParams object

    :param rxn: tuple describing reactants, products, and third body
    :type rxn: tuple ((rct1, rct2, ...), (prd1, prd2, ...), (thirdbod,))
    :param params: parameters for a reaction
    :type params: autoreact.RxnParams object
    :param cmts_dct: comments for a single reaction
    :type cmt_dct: dict {'header': cmt, 'inline': cmt, 'footer': cmt}
    :param max_len: length of the longest reaction name in the mechanism
    :type max_len: int
    :param rxnkeys: all keys of the reaction dictionary to be fitted (needed to assing "=" or "=>" sign in header)
    :type rxnkeys: list(tuple)
    :return ckin_str: Chemkin-formatted string describing the reaction
    :rtype: str
    """

    # Get the Chemkin string for the chemical reaction
    # rxn = util.format_rxn_name(rxn, fullset_rxns = rxnkeys)
    # ^ AVC: this doesn't work, because the arr, plog, etc. functions format their own reactions
    # (Also, TROE needs to add +Ms if missing)

    # Get the functional forms to write (usually will only be one)
    forms = params.get_existing_forms()
    # Get the comment strings
    if cmts_dct is None:
        header_cmt = ""
        inline_cmt = ""
        footer_cmt = ""
    else:
        header_cmt = cmts_dct["header"] if "header" in cmts_dct.keys() else ""
        inline_cmt = cmts_dct["inline"] if "inline" in cmts_dct.keys() else ""
        footer_cmt = cmts_dct["footer"] if "footer" in cmts_dct.keys() else ""
    # Get information on unusual duplicates (i.e., dup other than Arrhenius)
    dups, dup_counts = params.check_for_dups()  # dups is a Boolean

    # Loop over each functional form and write each one (usually only one)
    ckin_str = header_cmt
    for form in forms:
        if form == "arr":
            ckin_str += arr(
                rxn,
                params.arr,
                colliders=params.arr_collid,
                max_len=max_len,
                inline_cmt=inline_cmt,
                rxnkeys=rxnkeys,
            )
        elif form == "plog":
            ckin_str += plog(
                rxn,
                params.plog,
                max_len=max_len,
                inline_cmt=inline_cmt,
                rxnkeys=rxnkeys,
            )
        elif form == "cheb":
            ckin_str += cheb(
                rxn,
                params.cheb["alpha"],
                params.cheb["tlim"],
                params.cheb["plim"],
                params.cheb["one_atm_arr"],
                max_len=max_len,
                inline_cmt=inline_cmt,
                rxnkeys=rxnkeys,
            )
        elif form == "troe":
            ckin_str += troe(
                rxn,
                params.troe["highp_arr"],
                params.troe["lowp_arr"],
                params.troe["troe_params"],
                colliders=params.troe["collid"],
                max_len=max_len,
                inline_cmt=inline_cmt,
                rxnkeys=rxnkeys,
            )
        elif form == "lind":
            ckin_str += lind(
                rxn,
                params.lind["highp_arr"],
                params.lind["lowp_arr"],
                colliders=params.lind["collid"],
                max_len=max_len,
                inline_cmt=inline_cmt,
                rxnkeys=rxnkeys,
            )
        if dups:
            ckin_str += "  DUP\n"

    ckin_str += footer_cmt

    # Handle any duplicates that may exist (this is not usual)
    dup_ckin_str = handle_duplicates(rxn, params, dup_counts, max_len, rxnkeys=rxnkeys)
    ckin_str += dup_ckin_str
    ckin_str += "\n"

    return ckin_str


def handle_duplicates(rxn, params, dup_counts, max_len, rxnkeys=[]):
    """Writes any unusual duplicate cases. These only occur in strange cases when
    duplicates of a functional form (e.g., two PLOGs) are described for
    the same reaction. These should not really occur, but are nonetheless
    handled here.

    :param rxn: Chemkin-formatted reaction name
    :type rxn: str
    :param params: parameters for a reaction
    :type params: autoreact.RxnParams object
    :param dup_counts: the number of duplicates for each form
    :type dup_counts: dct {'plog': num_dups, 'cheb': ...} ('arr' excluded)
    :param max_len: length of the longest reaction name in the mechanism
    :type max_len: int
    :return ckin_str: Chemkin-formatted string describing the reaction
    :rtype: str
    """

    ckin_str = ""
    for form, dup_count in dup_counts.items():
        # Note that Arrhenius is absent because this is allowed to have dups
        # and is thus not an unusual case
        if form == "plog":
            for dup_idx in range(dup_count):
                ckin_str += plog(
                    rxn, params.plog_dups[dup_idx], max_len=max_len, rxnkeys=rxnkeys
                )
                ckin_str += "  DUP\n"
        elif form == "cheb":
            for dup_idx in range(dup_count):
                ckin_str += cheb(
                    rxn,
                    params.cheb_dups[dup_idx]["alpha"],
                    params.cheb_dups[dup_idx]["tlim"],
                    params.cheb_dups[dup_idx]["plim"],
                    params.cheb_dups[dup_idx]["one_atm_arr"],
                    max_len=max_len,
                    rxnkeys=rxnkeys,
                )
                ckin_str += "  DUP\n"
        elif form == "troe":
            for dup_idx in range(dup_count):
                ckin_str += troe(
                    rxn,
                    params.troe_dups[dup_idx]["highp_arr"],
                    params.troe_dups[dup_idx]["lowp_arr"],
                    params.troe_dups[dup_idx]["troe_params"],
                    colliders=params.troe_dups[dup_idx]["collid"],
                    max_len=max_len,
                    rxnkeys=rxnkeys,
                )
                ckin_str += "  DUP\n"
        elif form == "lind":
            for dup_idx in range(dup_count):
                ckin_str += lind(
                    rxn,
                    params.lind_dups[dup_idx]["highp_arr"],
                    params.lind_dups[dup_idx]["lowp_arr"],
                    colliders=params.lind_dups[dup_idx]["collid"],
                    max_len=max_len,
                    rxnkeys=rxnkeys,
                )
                ckin_str += "  DUP\n"

    return ckin_str


def arr(rxn, arr_tuples, colliders=None, max_len=45, inline_cmt=None, rxnkeys=[]):
    """Writes a reaction in the Arrhenius form

    :param rxn: tuple describing reactants, products, and third body
    :type rxn: tuple ((rct1, rct2, ...), (prd1, prd2, ...), (thirdbod,))
    :param arr_tuples: Arrhenius high-P (i.e., high-P) parameters
    :type arr_tuples: tuple of tuples ((A1, n1, Ea1), (A2, n2, Ea2), ...)
    :param max_len: length of the longest reaction name in the mechanism
    :type max_len: int
    :return arr_str: Chemkin reaction string with Arrhenius parameters
    :rtype: str
    """

    reaction = util.format_rxn_name(rxn, rxnkeys=rxnkeys)

    # Write the Arrhenius parameter string
    arr_str = ""
    for idx, arr_tuple in enumerate(arr_tuples):
        if idx != 0:
            inline_cmt = None  # this prevents duplicate inline comments
        arr_str += _highp_str(
            reaction, arr_tuple, max_len=max_len, inline_cmt=inline_cmt
        )
        if len(arr_tuples) > 1:
            arr_str += "  DUP\n"

    # Write the collider efficiencies string
    if colliders is not None:
        arr_str += _format_collider_string(colliders)

    return arr_str


def plog(rxn, plog_dct, max_len=45, inline_cmt=None, rxnkeys=[]):
    """Writes a reaction in the PLOG form

    :param rxn: tuple describing reactants, products, and third body
    :type rxn: tuple ((rct1, rct2, ...), (prd1, prd2, ...), (thirdbod,))
    :param plog_dct: Arrhenius fitting parameters at all pressures
    :type plog_dct: dict{pressure: [Arrhenius params]}
    :param max_len: length of the longest reaction name in the mechanism
    :type max_len: int
    :return plog_str: Chemkin reaction string with PLOG parameters
    :rtype: str
    """

    def _pressure_str(pressure, params, max_len=45):
        """Write a line in a PLOG string

        :param pressure: pressure at which to write the string
        :type pressure: float
        :param params: Arrhenius parameters at the specified pressure
        :type params: list of floats
        :param max_len: length of longest reaction name in the mechanism
        :type max_len: int
        :return single_str: Chemkin reaction string with PLOG parameters
            at a single pressure
        :rtype: str
        """

        plog_buffer = max_len + BUFFER - 12
        if plog_buffer < 0:
            plog_buffer = "0"
        else:
            plog_buffer = str(plog_buffer)

        [a_par, n_par, ea_par] = params

        single_str = (
            "{0:<" + plog_buffer + "s}{1:<12.3E}{2:<10.3E}{3:>9.3f}{4:>9.0f} /\n"
        ).format("  PLOG /", pressure, a_par, n_par, ea_par)

        return single_str

    # Obtain a list of the pressures and sort from low to high pressure
    unsorted_pressures = plog_dct.keys()
    pressures = sorted(unsorted_pressures)

    # Write the header for the reaction
    reaction = util.format_rxn_name(rxn, rxnkeys=rxnkeys)
    if 1 in pressures:
        header_params = plog_dct[1][0]  # 1 atm fit if available
    else:
        header_params = [1.0, 0.0, 0.0]  # otherwise, fake params
    plog_str = _highp_str(
        reaction, header_params, max_len=max_len, inline_cmt=inline_cmt
    )

    # Loop over each pressure
    for pressure in pressures:
        plog_params = plog_dct[pressure]
        for param_set in plog_params:
            assert len(param_set) % 3 == 0, (
                f"Arr params should be a multiple of 3, is {len(param_set)}"
                + f" for {reaction}"
            )

            # Loop over however many Arrhenius sets there are,
            # writing a PLOG line for each
            num_arr_sets = int(len(param_set) / 3)
            for idx in range(num_arr_sets):
                current_param_set = param_set[3 * idx : 3 * (idx + 1)]
                plog_str += _pressure_str(pressure, current_param_set, max_len=max_len)

    return plog_str


def cheb(
    rxn, alpha, tlim, plim, one_atm_arr=None, max_len=45, inline_cmt=None, rxnkeys=[]
):
    """Writes a reaction in the Chebyshev form

    :param rxn: tuple describing reactants, products, and third body
    :type rxn: tuple ((rct1, rct2, ...), (prd1, prd2, ...), (thirdbod,))
    :param alpha: Chebyshev coefficient matrix
    :type alpha: numpy.ndarray
    :param tlim: Chebyshev temperature limits
    :type tlim: tuple (tmin, tmax)
    :param plim: Chebyshev pressure limits
    :type plim: tuple (pmin, pmax)
    :param one_atm_arr: Arrhenius parameters at 1 atm
    :type one_atm_arr: tuple of tuples ((A1, n1, Ea1), (A2, n2, Ea2), ...)
    :param max_len: length of the longest reaction name in the mechanism
    :type max_len: int
    :return cheb_str: Chemkin reaction string with Chebyshev parameters
    :rtype: str
    """

    # Write reaction header
    reaction = util.format_rxn_name(rxn, rxnkeys=rxnkeys)
    if one_atm_arr is None:
        one_atm_arr = [[1, 0, 0]]
    cheb_str = _highp_str(
        reaction, one_atm_arr[0], max_len=max_len, inline_cmt=inline_cmt
    )

    # Write the temperature and pressure ranges
    cheb_str += _misc_troe_cheb("TCHEB", tlim, newline=True, val="float")
    cheb_str += _misc_troe_cheb("PCHEB", plim, newline=True, val="float")

    # Write the dimensions of the alpha matrix
    nrows = len(alpha)
    ncols = len(alpha[0])
    cheb_str += _misc_troe_cheb("CHEB", (nrows, ncols), newline=True, val="int")

    # Write the parameters from the alpha matrix
    for row in alpha:
        cheb_str += _misc_troe_cheb("CHEB", row, newline=True, val="exp")

    return cheb_str


def troe(
    rxn,
    high_params,
    low_params,
    troe_params,
    colliders=None,
    max_len=45,
    inline_cmt=None,
    rxnkeys=[],
):
    """Writes a reaction in the Troe form

    :param rxn: tuple describing reactants, products, and third body
    :type rxn: tuple ((rct1, rct2, ...), (prd1, prd2, ...), (thirdbod,))
    :param high_params: Arrhenius high-P parameters
    :type high_params: list of floats
    :param low_params: Arrhenius low-P parameters
    :type low_params: list of floats
    :param troe_params: Troe parameters: alpha, T***, T*, and T**
        (T** is optional)
    :type troe_params: list of floats
    :param colliders: collision enhancement factors for bath gases
    :type colliders: list((str, float))
    :param max_len: length of the longest reaction name in the mechanism
    :type max_len: int
    :return troe_str: Chemkin reaction string with Troe parameters
    :rtype: str
    """

    reaction = util.format_rxn_name(
        rxn, pdep=True, rxnkeys=rxnkeys
    )  # pdep: add (+M) if absent

    assert (
        len(high_params[0]) == 3
    ), f"{len(high_params)} highP params for {reaction}, should be 3"
    assert (
        len(low_params[0]) == 3
    ), f"{len(low_params)} lowP params for {reaction}, should be 3"
    assert len(troe_params) in (
        3,
        4,
    ), f"{len(troe_params)} Troe params for {reaction}, should be 3 or 4"

    troe_str = _highp_str(
        reaction, high_params[0], max_len=max_len, inline_cmt=inline_cmt
    )
    troe_str += _lowp_str(low_params[0], max_len=max_len)
    troe_str += _misc_troe_cheb("TROE", troe_params, newline=True, val="exp")

    # Write the collider efficiencies string
    if colliders:
        troe_str += _format_collider_string(colliders)

    return troe_str


def lind(
    rxn,
    high_params,
    low_params,
    colliders=None,
    max_len=45,
    inline_cmt=None,
    rxnkeys=[],
):
    """Writes a reaction in the Lindemann form

    :param rxn: tuple describing reactants, products, and third body
    :type rxn: tuple ((rct1, rct2, ...), (prd1, prd2, ...), (thirdbod,))
    :param high_params: Arrhenius high-P parameters
    :type high_params: list of floats
    :param low_params: Arrhenius low-P parameters
    :type low_params: list of floats
    :param colliders: names and collision enhancement factors for bathgases
    :type colliders: list((str, float))
    :param max_len: length of the longest reaction name in the mechanism
    :type max_len: int
    :return lind_str: Chemkin reaction string with Lindemann parameters
    :rtype: str
    """

    reaction = util.format_rxn_name(
        rxn, pdep=True, rxnkeys=rxnkeys
    )  # pdep: add (+M) if absent
    lind_str = _highp_str(
        reaction, high_params[0], max_len=max_len, inline_cmt=inline_cmt
    )
    lind_str += _lowp_str(low_params[0], max_len=max_len)

    # Write the collider efficiencies string
    if colliders:
        lind_str += _format_collider_string(colliders)

    return lind_str


def _highp_str(reaction, arr_tuple, max_len=45, inline_cmt=None):
    """Write a single set of high-pressure Arrhenius parameters

    :param reaction: Chemkin-formatted reaction name
    :type reaction: str
    :param arr_tuple: Arrhenius high-P parameters
    :type arr_tuple: tuple (or list) of floats (A, n, Ea)
    :param max_len: length of the longest reaction name in the mechanism
    :type max_len: int
    :param inline_cmt: optional inline comment
    :type inline_cmt: str
    :return highp_str: Chemkin reaction string with high-P Arrhenius params
    :rtype: str
    """

    highp_buffer = str(max_len + BUFFER)
    [a_par, n_par, ea_par] = arr_tuple
    highp_str = ("{0:<" + highp_buffer + "s}{1:<10.3E}{2:>9.3f}{3:>9.0f}").format(
        reaction, a_par, n_par, ea_par
    )

    if inline_cmt:
        highp_str += " " * INLINE_BUFFER + inline_cmt

    highp_str += "\n"

    return highp_str


def _lowp_str(arr_tuple, max_len=45):
    """Write a single set of low-pressure Arrhenius parameters (i.e., w/'LOW')

    :param reaction: Chemkin-formatted reaction name
    :type reaction: str
    :param arr_tuple: Arrhenius low-P parameters
    :type arr_tuple: tuple (or list) of floats (A, n, Ea)
    :param max_len: length of the longest reaction name in the mechanism
    :type max_len: int
    :return lowp_str: Chemkin reaction string with low-P Arrhenius params
    :rtype: str
    """

    lowp_buffer = str(max_len + BUFFER)
    [a_par, n_par, ea_par] = arr_tuple
    lowp_str = ("{0:<" + lowp_buffer + "s}{1:<10.3E}{2:>9.3f}{3:>9.0f}  /").format(
        "  LOW  /", a_par, n_par, ea_par
    )

    lowp_str += "\n"

    return lowp_str


def _misc_troe_cheb(header, params, newline=False, val="exp"):
    """Write a string containing slash-delimited params used for
    Troe and Chebyshev functional forms

    :param header: the header str to be written before the slashes
    :type header: str
    :param params: params to be written
    :type params: list or numpy array or tuple of floats
    :param newline: whether or not to add a new line after the params
    :type newline: Bool
    :param val: indicates how the parameters are to be formatted;
        should be 'exp' or 'float' or 'int'
    :type val: str
    :return params_str: a string containing the formatted params
    :rtype: str
    """

    if val == "exp":
        val_str = "{0:12.3E}"
    elif val == "float":
        val_str = "{0:12.2f}"
    elif val == "int":
        val_str = "{0:12.0f}"

    params_str = f"  {header.upper():5s}/"
    for param in params:
        if param:  # skip any None entries
            params_str += "".join(val_str.format(param))
    params_str += " /"
    if newline:
        params_str += "\n"

    return params_str


def _format_collider_string(colliders):
    """Write the string for the bath gas collider and their efficiencies
    for the Lindemann and Troe functional expressions:

    :param colliders: the {collider: efficiency} dct
    :type colliders: dct {str: float}
    :return: collider_str: Chemkin string with colliders and efficiencies
    :rtype: str
    """

    collider_str = "  "  # two blank spaces
    for collider, efficiency in colliders.items():
        collider_str += f"{collider:s}/{efficiency:4.3f}/   "
    collider_str += "\n"

    return collider_str
