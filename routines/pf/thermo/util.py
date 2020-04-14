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
