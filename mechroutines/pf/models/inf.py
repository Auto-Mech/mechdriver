"""
Random functions that are needed in drivers and routines
"""

from mechroutines.pf.models import typ
from mechlib.amech_io import parser
from mechlib.amech_io import printer as ioprinter
from mechanalyzer.inf import rxn as rinfo


def make_rxn_str(rlst, prepend=''):
    """ convert list to string
    """
    return prepend + '+'.join(rlst)
