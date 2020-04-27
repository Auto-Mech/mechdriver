"""
  Handle parsing of the class.dat file used to set reaction types
"""

import os
import sys
import chemkin_io
from lib.amech_io.reader import ptt
from lib.amech_io.cleaner import remove_whitespace


CLA_INP = 'inp/class.csv'


def set_class_with_dct(cla_dct, reacs, prods):
    """ set the class using the class dictionary
    """
    rxn = (reacs, prods)
    rxn_rev = (prods, reacs)
    if rxn in cla_dct:
        inp_class = cla_dct[rxn]
        flip_rxn = False
    elif rxn_rev in cla_dct:
        inp_class = cla_dct[rxn_rev]
        flip_rxn = True
    else:
        inp_class = None
        flip_rxn = False

    return inp_class, flip_rxn


def parse_rxn_class_file(job_path):
    """ Read the class dictionary
    """

    if os.path.exists(os.path.join(job_path, CLA_INP)):
        print('  class.dat found. Reading contents...')
        cla_str = ptt.read_inp_str(job_path, CLA_INP)
        cla_dct = _build_cla_dct(cla_str)
    else:
        print('  No class.dat found.')
        cla_dct = {}

    return cla_dct


def _build_cla_dct(cla_str):
    """ read file
    """
    cla_dct = {}
    cla_str = remove_whitespace(cla_str)
    for line in cla_str.splitlines():
        try:
            [rxn_line, rclass] = line.split('||')
            reacs = chemkin_io.parser.reacion.reactant_names(rxn_line)
            prods = chemkin_io.parser.reacion.product_names(rxn_line)
            cla_dct[(reacs, prods)] = rclass
        except:
            print('*ERROR: Error in formatting line')
            print(line)
            sys.exit()

    return cla_dct
