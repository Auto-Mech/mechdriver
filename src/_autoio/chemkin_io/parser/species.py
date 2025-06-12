""" functions operating on the species block string
"""

import autoparse.find as apf


def names(block_str, exclude_names=()):
    """ Parses the names of species from the species block
        of the mechanism input file.

        :param block_str: string for species blockA
        :type: block_str: str
        :param exclude_names: names of species to ignore during parsing
        :type exclude_names: list(str)
        :return spc_names: names of species that were not excluded
        :rtype: tuple
    """

    if block_str is not None:
        spc_names = apf.split_words(block_str)
        spc_names = tuple(filter(lambda x: x not in exclude_names, spc_names))
    else:
        spc_names = None

    return spc_names
