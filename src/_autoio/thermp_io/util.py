"""
  Utility function for ThermP writing
"""

import re


def get_atom_counts_dict(formula):
    """ Obtain the atom symbols and counts from a forumla.

        :param formula: stoichiometry of species
        :type formula: str
        :return atom_counts_dct: dct containing the number of each atom
        :rtype: dict[str: int]
    """

    # Search for alpha-integer pairs
    search_str = r"([A-Z][a-z]?)(\d+)?"

    # Obtain a dictionary for the number associated with atom symbol
    atom_counts_dct = {k: int(v) if v else 1
                       for k, v in re.findall(search_str, formula)}

    return atom_counts_dct
