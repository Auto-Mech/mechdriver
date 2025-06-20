""" functions operating on mechanism strings
"""

import autoparse.pattern as app
import autoparse.find as apf
from ioformat import remove_comment_lines
from ioformat import remove_whitespace_from_string
from chemkin_io.parser.reaction import get_rxn_param_dct


def species_block(mech_str, remove_comments=True):
    """Parses the species block out of the mechanism input file.

    :param mech_str: string of mechanism input file
    :type mech_str: str
    :param remove_comments: elect to remove comment liness from string
    :type remove_comments: bool
    :return block_str: string containing species block
    :rtype: string
    """

    block_str = _block(
        string=_clean_up(mech_str, remove_comments=remove_comments),
        start_pattern=app.one_of_these(["SPECIES", "SPEC"]),
        end_pattern="\nEND",
    )

    _check_if_empty(block_str, "SPECIES")

    return block_str


def reactions(mech_str):
    """Parses all of the chemical equations and corresponding fitting from the
    mechanism file.

    :param mech_str: string of mechanism input file
    :type mech_str: str
    :return rxn_param_dct: dct {rxn1: params1, rxn2: ...}
    :rtype: dict
    """
    ea_units, a_units = reaction_units(mech_str)
    block_str = reaction_block(mech_str, remove_comments=True)
    rxn_param_dct = get_rxn_param_dct(block_str, ea_units, a_units)

    return rxn_param_dct


def reaction_block(mech_str, remove_comments=True):
    """Parses the reaction block out of the mechanism input file.

    :param mech_str: string of mechanism input file
    :type mech_str: str
    :param remove_comments: elect to remove comment liness from string
    :type remove_comments: bool
    :return block_str: string containing reaction block
    :rtype: string
    """

    block_str = _block(
        string=_clean_up(mech_str, remove_comments=remove_comments),
        start_pattern=app.one_of_these(["REACTIONS", "REAC"]),
        end_pattern="\nEND",
        # at least ensure that "end" is on a new line,
        # otherwise "end" might be written anywhere
    )

    _check_if_empty(block_str, "REACTIONS")

    return block_str


def thermo_block(mech_str, remove_comments=True):
    """Parses the thermo block out of the mechanism input file.

    :param mech_str: string of mechanism input file
    :type mech_str: str
    :param remove_comments: elect to remove comment liness from string
    :type remove_comments: bool
    :return block_str: string containing thermo block
    :rtype: string
    """

    block_str = _block(
        string=_clean_up(mech_str, remove_comments=remove_comments),
        start_pattern=app.one_of_these(
            ["THERMO ALL", "THERM ALL", "THER ALL", "THERMO", "THERM", "THER"]
        ),
        end_pattern="\nEND",
    )

    _check_if_empty(block_str, "THERMO")

    return block_str


def element_block(mech_str, remove_comments=True):
    """Parses the species block out of the mechanism input file.

    :param mech_str: string of mechanism input file
    :type mech_str: str
    :param remove_comments: elect to remove comment liness from string
    :type remove_comments: bool
    :return block_str: string containing species block
    :rtype: string
    """

    block_str = _block(
        string=_clean_up(mech_str, remove_comments=remove_comments),
        start_pattern=app.one_of_these(["ELEMENTS"]),
        end_pattern="\nEND",
    )

    _check_if_empty(block_str, "ELEMENTS")

    return block_str


def _block(string, start_pattern, end_pattern):
    """return a block delimited by start and end patterns"""
    contents_pattern = app.capturing(app.one_or_more(app.WILDCARD, greedy=False))
    pattern = start_pattern + contents_pattern + end_pattern
    contents = apf.first_capture(pattern, string)

    return contents


def reaction_units(mech_str):
    """Parses from the mechanism input file, the units of the
    pre-exponential (A) and activation enery (Ea) fitting parameter.
    :param mech_str: string of mechanism input file
    :type mech_str: str
    :return units: units for fitiing parameters (Ea unit, A unit)
    :rtype: list(float)
    """

    def _reaction_units(mech_str, start_pattern, units_pattern):
        """Helper function used to parse the units at the head of
        the reaction block of the mechanism file string.
        :param mech_str: string of mechanism input file
        :type mech_str: str
        :param start_pattern: start pattern at line at head
        :type start_pattern: str
        :param units_pattern: patterns for various unit strings
        :type units_pattern: str
        :return units: units for fitiing parameters (A unit, Ea unit)
        :rtype: list
        """

        rxn_line_pattern = start_pattern + app.capturing(app.LINE_FILL)
        units_str = apf.first_capture(rxn_line_pattern, mech_str)
        units_lst = apf.all_captures(units_pattern, units_str)

        ckin_ea_units = [
            "CAL/MOLE",
            "KCAL/MOLE",
            "JOULES/MOLE",
            "KJOULES/MOLE",
            "KELVINS",
        ]
        ckin_a_units = ["MOLES", "MOLECULES"]

        if units_lst:
            if any(unit in ckin_ea_units for unit in units_lst):
                for unit in ckin_ea_units:
                    if unit in units_lst:
                        ea_unit = unit.lower()
            else:
                ea_unit = "cal/mole"
            if any(unit in ckin_a_units for unit in units_lst):
                for unit in ckin_a_units:
                    if unit in units_lst:
                        a_unit = unit.lower()
            else:
                a_unit = "moles"
            units = (ea_unit, a_unit)
        else:
            units = ("cal/mole", "moles")

        return units

    units = _reaction_units(
        mech_str=_clean_up(mech_str),
        start_pattern=app.one_of_these(["REACTIONS", "REAC"]),
        units_pattern=app.one_or_more(app.one_of_these([app.LETTER, app.escape("/")])),
    )

    return units


# Clean up the ChemKin mechanism strings
def _clean_up(mech_str, remove_comments=True):
    """Cleans up mechanism input string by converting specific comment
    lines that are used later and removes other comments and
    whitespace from mech string.

    :param mech_str: string of mechanism input file
    :param mech_str: str
    :return mech_str: string with altered comment lines
    :rtype: string
    """
    mech_str = _convert_comment_lines(mech_str)
    if remove_comments:
        mech_str = remove_comment_lines(mech_str, delim_pattern=app.escape("!"))
    mech_str = remove_whitespace_from_string(mech_str)

    return mech_str


def _convert_comment_lines(mech_str):
    """alter based on above...
    Reads a string for the mechanism input file and alters certain
    comment lines, by removing the comment symbols. This is so they
    are not removed later by functions which remove comments from string.

    :param mech_str: string of mechanism input file
    :type mech_str: str
    :return mech_str: string with altered comment lines
    :rtype: string
    """

    # Set the lines in the file (in_lines) and their replacements (out_lines)
    inlines = [
        app.escape("!") + app.SPACES + app.escape("Pressure:"),
    ]
    outlines = [
        app.escape("Pressure:"),
    ]

    # Loop over lines and make all the replacements in the mech string
    for inline, outline in zip(inlines, outlines):
        mech_str = apf.replace(inline, outline, mech_str, case=True)

    return mech_str


def _check_if_empty(block_str, block_type):
    """Checks if a block str is None and prints a message if so

    :param block_str: string containing some kind of Chemkin block
    :type block_str: str
    :param block_type: the type of block ('species', 'reaction', etc.)
    :type block_type: str
    """

    if block_str is None:
        print(f"No Chemkin {block_type} block detected. Check start and end.")
