""" theory
"""

import automol
import ioformat
from mechlib.amech_io.parser._keywrd import defaults_from_val_dct
from mechlib.amech_io.parser._keywrd import check_dct1


# DCTS
# rquired, Type, allowed, default,
# maybe set defaults using the qchem params script?
THY_REQ = ('program', 'method', 'basis', 'orb_res')
THY_VAL_DCT = {
    'program': (str, (), None),
    'method': (str, (), None),
    'basis': (str, (), None),
    'orb_res': (str, ('RR', 'UU', 'RU'), None),
    'ncycles': (int, (), None),
    'mem': (float, (), None),
    'nprocs': (int, (), None),
    'econv': (float, (), None),
    'gconv': (float, (), None)
}


def theory_dictionary(thy_str):
    """ Parse the theory.dat file
    """

    # Format input to keyword-value theory dcts
    thy_blocks = ioformat.ptt.named_end_blocks(
        thy_str, 'level', footer='level')
    thy_dct = ioformat.ptt.keyword_dcts_from_blocks(thy_blocks)

    # Add defaults to each thy dictionary
    full_thy_dct = {}
    for lvl, dct in thy_dct.items():
        full_thy_dct[lvl] = automol.util.dict_.right_update(
            defaults_from_val_dct(THY_VAL_DCT), dct)

    # Check each dictionary
    for lvl, dct in full_thy_dct.items():
        check_dct1(dct, THY_VAL_DCT, THY_REQ, 'Thy-{}'.format(lvl))

    return full_thy_dct
