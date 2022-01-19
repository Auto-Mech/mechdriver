""" Parses the `mechanism.dat` input file for MechDriver that contains
    all of the reactions of mechanism. The format of this file corresponds
    to some user-specified format.
"""

import mechanalyzer.parser.pes


def pes_dictionary(mech_str, mech_type, spc_dct, printlog=True):
    """ Calls mechanalyzer to do the dictionary
    """
    return mechanalyzer.parser.pes.pes_dictionary(
        mech_str, mech_type, spc_dct, printlog=printlog)
