"""
Sort by selected criteria written in this script
"""

import os
import sys

import chemkin_io.writer
import numpy as np
from ioformat import pathtools

from mechanalyzer.builder import sorter
from mechanalyzer.parser import ckin_ as ckin_parser
from mechanalyzer.parser import mech as mparser
from mechanalyzer.parser import spc as sparser

# Set useful global variables
CWD = os.getcwd()


def main(
    mech: str = "mechanism.dat",
    spc: str = "species.csv",
    therm: str = "therm.dat",
    sortopts: str = "sort.dat",
    outmech: str = "outmech.dat",
    outspc: str = "outspc.csv",
    outgroups: str = "pes_groups.dat",
):
    """Sort the reactions in a mechanism

    :param mech: Input mechanism file name, defaults to "mechanism.dat"
    :param spc: Input species file name, defaults to "species.csv"
    :param therm: Input thermo file name, defaults to "therm.dat"
    :param sort: Input sort file name, defaults to "sort.dat"
    :param outmech: Output mechanism file name, defaults to "outmech.dat"
    :param outspc: Output species file name, defaults to "outspc.csv"
    :param outgroups: Output PES groups file name, defaults to "pes_groups.dat"
    """

    # Read the input files
    spc_str = pathtools.read_file(CWD, spc, remove_comments="!")
    mech_str = pathtools.read_file(CWD, mech, remove_comments="!")
    sort_str = pathtools.read_file(CWD, sortopts, remove_comments="#")

    # Check if the input strings exist
    if any(string is None for string in (spc_str, mech_str, sort_str)):
        print("ERROR: Input file(s) species.csv, mechanism.dat, sort.dat missing")
        sys.exit()

    # read thermo and filter pes groups
    if os.path.exists(therm):
        therm_str = pathtools.read_file(CWD, therm, remove_comments="!")
        spc_therm_dct = ckin_parser.parse_spc_therm_dct(
            therm_str, np.arange(300, 2010, 10)
        )
    else:
        print("Warning: No therm file found. Running sort without it...")
        spc_therm_dct = None

    # Build sorted mechanism files
    isolate_spc, sort_lst, prompt_filter_dct = mparser.parse_sort(sort_str)
    param_dct_sort, mech_spc_dct, cmts_dct, pes_groups, rxns_filter = (
        sorter.sorted_mech(
            spc_str,
            mech_str,
            isolate_spc,
            sort_lst,
            spc_therm_dct=spc_therm_dct,
            dct_flt_grps=prompt_filter_dct,
        )
    )
    rxn_cmts_dct = chemkin_io.writer.comments.get_rxn_cmts_dct(rxn_sort_dct=cmts_dct)

    # Write the output files (need to make general at some point)
    headers = sparser.csv_headers(mech_spc_dct)
    sortd_csv_str = sparser.csv_string(mech_spc_dct, headers)
    sortd_mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
        mech_spc_dct=mech_spc_dct,
        rxn_param_dct=param_dct_sort,
        rxn_cmts_dct=rxn_cmts_dct,
    )
    # replace all '!' with '#' in mech string
    sortd_mech_str = sortd_mech_str.replace("! pes", "# pes")
    pathtools.write_file(sortd_csv_str, CWD, outspc)
    pathtools.write_file(sortd_mech_str, CWD, outmech)

    if pes_groups:
        pes_groups_str = chemkin_io.writer.pesgroups.write_pes_groups(pes_groups)
        pathtools.write_file(pes_groups_str, CWD, outgroups)

    try:
        np.savetxt(
            CWD + "/rxns_prompt_dh.out",
            rxns_filter[2:, :],
            delimiter="\t\t",
            header="\t".join(rxns_filter[1, :]),
            fmt=list(rxns_filter[0, :]),
        )
    except TypeError:
        print("no filtering selected, prompt list not generated")
    # pathtools.write_file(rxns_filter, CWD, 'rxns_prompt_dh.out')
