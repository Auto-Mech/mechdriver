""" Stereoexpansion script
"""

import os
import time

import chemkin_io
import ioformat

import mechanalyzer
from mechanalyzer.builder import sorter


def main():
    """Main function for CLI"""

    # Initialize the start time for script execution
    t0 = time.time()

    # Set useful global variables
    oscwd = os.getcwd()

    # READ AND PARSE INPUT
    print("\n---- Parsing the species and mechanism files ---\n")

    # Read the build input file and set up info
    bld_str = ioformat.pathtools.read_file(
        oscwd, "build.dat", remove_comments="#", remove_whitespace=True
    )
    file_dct, _ = mechanalyzer.parser.build_input_file(bld_str)
    # Read input species and mechanism files into dictionary
    out_loc = oscwd
    out_spc = file_dct["out_spc"]
    out_mech = file_dct["out_mech"]
    inp_spc_str = ioformat.pathtools.read_file(
        oscwd, file_dct["inp_spc"], remove_comments="!", remove_whitespace=True
    )
    inp_mech_str = ioformat.pathtools.read_file(
        oscwd, file_dct["inp_mech"], remove_comments="!", remove_whitespace=True
    )
    sort_str = ioformat.pathtools.read_file(
        oscwd, file_dct["sort"], remove_comments="#", remove_whitespace=True
    )
    enant = file_dct["enant"] if "enant" in file_dct else False
    rename = file_dct["rename"] if "rename" in file_dct else False
    enant_label = file_dct["enant_label"] if "enant_label" in file_dct else True

    # old main(): carry out all the mechanism things you wanna do
    # Build the initial dictionaries
    mech_spc_dct = mechanalyzer.parser.spc.build_spc_dct(inp_spc_str, "csv")
    rxn_param_dct = mechanalyzer.parser.mech.parse_mechanism(inp_mech_str, "chemkin")
    isolate_spc, sort_lst, _ = mechanalyzer.parser.mech.parse_sort(sort_str)

    # ADD STEREO TO SPC AND REACTIONS
    print(
        "\n---- Adding stereochemistry to InChIs of mechanism"
        " species where needed ---\n"
    )
    mech_spc_dct = mechanalyzer.parser.spc.stereochemical_spc_dct(
        mech_spc_dct, nprocs="auto", all_stereo=False, enant=enant
    )
    print("mechspc dct2", mech_spc_dct)
    print("Mechanism species with stereo added")
    for name, dct in mech_spc_dct.items():
        print(f'Name: {name:<25s} InChI: {dct["inchi"]}')

    # Expand the reactions in the mechanism to include stereochemical variants
    print(
        "\n---- Expanding the list of mechanism reactions to include all"
        " valid, stereoselective permutations ---\n"
    )
    full_rxn_lst = mechanalyzer.builder.expand_mech_stereo(
        rxn_param_dct, mech_spc_dct, nprocs="auto", enant=enant
    )
    spc_orig_name_dct = {}

    print("turning reaction list into mechanism dictionary")
    ste_mech_spc_dct, ste_rxn_dct = {}, {}
    ste_mech_spc_dct = mechanalyzer.builder.update_spc_dct_from_reactions(
        full_rxn_lst,
        ste_mech_spc_dct,
        rename=rename,
        enant_label=enant_label,
        spc_orig_name_dct=spc_orig_name_dct,
    )
    ste_rxn_dct = mechanalyzer.builder.update_rxn_dct(
        full_rxn_lst, ste_rxn_dct, ste_mech_spc_dct
    )
    ste_mech_spc_dct = mechanalyzer.builder.remove_spc_not_in_reactions(
        ste_rxn_dct, ste_mech_spc_dct
    )
    print("Writing expanded stereomechanism")
    write_mechanism(
        ste_mech_spc_dct,
        ste_rxn_dct,
        out_loc,
        out_spc,
        out_mech,
        sort_lst.copy(),
        isolate_spc,
        list(file_dct["headers"]),
        "_full",
    )
    print("Writing original mechanism")
    write_mechanism(
        mech_spc_dct,
        rxn_param_dct,
        out_loc,
        out_spc,
        out_mech,
        sort_lst.copy(),
        isolate_spc,
        list(file_dct["headers"]),
        "_orig",
    )

    # Compute script run time and print to screen
    tf = time.time()
    print("\n\nScript executed successfully.")
    print(f"Time to complete: {tf-t0:.2f}")
    print("Exiting...")


def write_mechanism(
    ste_mech_spc_dct,
    ste_rxn_dct,
    out_loc,
    out_spc,
    out_mech,
    sort_lst,
    isolate_spc,
    header_lst,
    suffix="",
):
    """write mechanism"""
    # Write the new (sorted) species dictionary to a string
    headers = ("smiles", "inchi", "inchikey", "mult", "charge")
    ste_mech_spc_dct_sort = mechanalyzer.parser.spc.reorder_by_atomcount(
        ste_mech_spc_dct
    )
    csv_str = mechanalyzer.parser.spc.csv_string(ste_mech_spc_dct_sort, headers)

    # Write initial string to call the sorter
    mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
        elem_tuple=None,
        mech_spc_dct=ste_mech_spc_dct_sort,
        spc_nasa7_dct=None,
        rxn_param_dct=ste_rxn_dct,
        rxn_cmts_dct=None,
    )

    # Use strings to generate ordered objects
    param_dct_sort, ste_mech_spc_dct_sort, cmts_dct, _, _ = sorter.sorted_mech(
        csv_str, mech_str, isolate_spc, sort_lst
    )
    rxn_cmts_dct = chemkin_io.writer.comments.get_rxn_cmts_dct(rxn_sort_dct=cmts_dct)

    # Write the dictionaries to ordered strings
    csv_str = mechanalyzer.parser.spc.csv_string(ste_mech_spc_dct_sort, header_lst)
    mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
        elem_tuple=(),
        mech_spc_dct=ste_mech_spc_dct_sort,
        spc_nasa7_dct=None,
        rxn_param_dct=param_dct_sort,
        rxn_cmts_dct=rxn_cmts_dct,
    )

    print("Writing files...")
    print(f"Writing {out_spc + suffix}")
    ioformat.pathtools.write_file(csv_str, out_loc, out_spc + suffix)
    print(f"Writing {out_mech + suffix}")
    ioformat.pathtools.write_file(mech_str, out_loc, out_mech + suffix)
