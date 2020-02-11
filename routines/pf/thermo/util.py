"""
utility functions that MAY not exist yet
"""

import re


def get_atom_counts_dict(formula):
    """ get the atom types and numbers out of
    """

    # Search for alpha-integer pairs
    search_str = r"([A-Z][a-z]?)(\d+)?"

    # Obtain a dictionary for the number associated with atom symbol
    atom_counts_dict = {k: int(v) if v else 1
                        for k, v in re.findall(search_str, formula)}

    return atom_counts_dict


def parse_line16(string):
    """
    Parse string containing exponental numbers of length 16 chars.
    Note: Numbers may not have a space in between.
    """

    assert len(string) % 16 == 0, 'Given string should have 16n chararacters'

    # Replace the exponent D with E
    string2 = string.replace('D', 'E')

    # Build a list of values from the string
    nchunks = len(string2) // 16
    vals = [0.0 for i in range(nchunks)]
    for i in range(nchunks):
        vals[i] = float(string2[i*16: (i+1)*16])

    return vals
