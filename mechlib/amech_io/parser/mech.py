""" Parses the `mechanism.dat` input file for MechDriver that contains
    all of the reactions of mechanism. The format of this file corresponds
    to some user-specified format.
"""

import automol.inchi
import automol.formula
from mechanalyzer.parser import pes
from mechanalyzer.parser import parse_pes_dct
from mechanalyzer.parser.mech import parse_mechanism
from mechanalyzer.builder import sorter


def pes_dictionary(mech_str, mech_type, spc_dct):
    """ Constructs the Potential-Energy-Surface dictionary for all of the
        channels of the user input utilizing the sorter functionality
        from mechanalyzer. Currently, we sort just via PES and then SUB-PES.

        Also, currently prints the PES channels.

        :param mech_str: mechanism.dat input file string
        :type mech_str: str
        :param mech_type: format of the mechanism string
        :type mech_type: str
        :param spc_dct: info for mechanism species and transition states
        :type spc_dct: mechdriver.spc_dct object
        :rtype: mechanalyer.pes.pes_dct object
    """

    def _fix_formula(pes_dct, spc_dct):
        """ have to fix fml since it is not parsed from file
        """
        new_pes_dct = {}
        for (_, pes_idx, subpes_idx), chnls in pes_dct.items():
            rcts = chnls[0][1][0]
            rct_ichs = tuple(spc_dct[rct]['inchi'] for rct in rcts)
            rct_ich = automol.inchi.join(rct_ichs)
            fml = automol.formula.string(automol.inchi.formula(rct_ich))
            new_pes_dct[(fml, pes_idx, subpes_idx)] = chnls
        return new_pes_dct

    # Initialize values used for the basic PES-SUBPES sorting
    sort_str = ['pes', 'subpes', 0]
    isolate_species = ()

    # Build and print the full sorted PES dict
    pes_dct = None
    if mech_str is not None:

        # Try and just read a PES dictionary from a file
        pes_dct = parse_pes_dct(mech_str)

        # If that fails then try and read a file normally and sort it
        if pes_dct is None:
            _, mech_info, _ = parse_mechanism(mech_str, mech_type, spc_dct)
            if mech_info is not None:
                srt_mch = sorter.sorting(
                    mech_info, spc_dct, sort_str, isolate_species)
                pes_dct = srt_mch.return_pes_dct()
        else:
            pes_dct = _fix_formula(pes_dct, spc_dct)
            print('Building PES dictionary from input file specification')

    if pes_dct is not None:
        pes.print_pes_channels(pes_dct)

    return pes_dct
