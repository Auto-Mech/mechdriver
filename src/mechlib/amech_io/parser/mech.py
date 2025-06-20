""" Parses the `mechanism.dat` input file for MechDriver that contains
    all of the reactions of mechanism. The format of this file corresponds
    to some user-specified format.
"""

from mechanalyzer.parser import pes
from mechanalyzer.parser.mech import parse_mechanism
from mechanalyzer.builder.sorter import sorting

# default values used for the basic PES-SUBPES sorting
SORT_STR = ['pes', 'subpes', 0]
ISOLATE_SPECIES = ()


def pes_dictionary(mech_str, mech_type, spc_dct, printlog=True):
    """ Calls mechanalyzer to do the dictionary
    """

    pes_dct = pes.pes_dictionary(
        mech_str, spc_dct)

    if pes_dct is None and mech_str is not None:
        print('No # pes.subpes.channel comment type found in mech: resorting ...')
        print('reminder: ! pes.subpes.channel format does not work!! # mandatory')
        rxn_param_dct = parse_mechanism(mech_str, mech_type)
        if rxn_param_dct is not None:
            srt_mch = sorting(
                rxn_param_dct, spc_dct, SORT_STR, ISOLATE_SPECIES)
            pes_dct = srt_mch.return_pes_dct()

    if pes_dct is not None and printlog:
        pes.print_pes_channels(pes_dct)

    return pes_dct
