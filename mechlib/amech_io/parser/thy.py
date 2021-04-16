""" theory
"""

import ioformat
from mechlib.amech_io.parser._keywrd import check_val_dictionary2
from mechlib.amech_io.parser._keywrd import THY_VAL_DCT


def theory_dictionary(thy_str):
    """ Parse the theory.dat file
    """

    # Format input to keyword-value theory dcts
    thy_blocks = ioformat.ptt.named_end_blocks(
        thy_str, 'level', footer='level')
    thy_dct = ioformat.ptt.keyword_dcts_from_blocks(thy_blocks)

    # Assess if the theory.dat input is valid
    for dct in thy_dct.values():
        check_val_dictionary2(dct, THY_VAL_DCT, 'thy.dat')

    return thy_dct
