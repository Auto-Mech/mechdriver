""" Parses the `mechanism.dat` input file for MechDriver that contains
    all of the reactions of mechanism. The format of this file corresponds
    to some user-specified format.
"""

from mechanalyzer.parser import pes
from mechanalyzer.parser.mech import parse_mechanism
from mechanalyzer.builder import sorter


def pes_dictionary(mech_str, mech_type, spc_dct):
    """ Constructs the Potential-Energy-Surface dictionary for all of the
        channels of the user input utilizing the sorter functionality
        from mechanalyzer. Currently, we sort just via PES and then SUB-PES.

        Format: {(formula, pes idx, sub pes idx): ((chnl_idx, rcts, prds),)

        Also, currently prints the PES channels.

        :param mech_str: mechanism.dat input file string
        :type mech_str: str
        :param mech_type:
        :type mech_type: str
        :param spc_dct:
        :type spc_dct: dict[str: ____]
        :rtype: dict[tuple(str, int, int)] = tuple(int, tuple(str))
    """

    # Initialize values used for the basic PES-SUBPES sorting
    sort_str = ['pes', 'subpes', 0]
    isolate_species = ()

    # Build and print the full sorted PES dict
    _, mech_info, _ = parse_mechanism(mech_str, mech_type, spc_dct)
    srt_mch = sorter.sorting(mech_info, spc_dct, sort_str, isolate_species)
    pes_dct = srt_mch.return_pes_dct()

    pes.print_pes_channels(pes_dct)

    return pes_dct
