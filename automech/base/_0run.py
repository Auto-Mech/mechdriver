""" Main AutoMech execution script
"""

import autofile
from drivers import esdriver, ktpdriver, procdriver, thermodriver, transdriver
from mechlib.amech_io import parser as ioparser
from mechlib.amech_io import printer as ioprinter

# import argparse
from mechlib.filesys import prefix_fs


def run(path: str = ".", safemode_off: bool = False):
    """Central Execution script to launch a MechDriver process which will
    parse all of the user-supplied input files in a specified directory, then
    launches all of the requested electronic structure, transport,
    thermochemistry and kinetics calculations via their associated
    sub-drivers.
    """
    if safemode_off:
        autofile.turn_off_safemode()
        ioprinter.info_message("Running with safemode turned OFF...")

    # Print the header message and host name
    ioprinter.program_header("amech")
    ioprinter.random_cute_animal()
    ioprinter.host_name()

    # Parse all of the input
    ioprinter.program_header("inp")

    ioprinter.info_message("\nReading files provided in the inp directory...")
    input = ioparser.read_amech_input(path)

    ioprinter.info_message("\nParsing input files for runtime parameters...")
    thy_dct = ioparser.thy.theory_dictionary(input["thy"])
    kmod_dct, smod_dct = ioparser.models.models_dictionary(input["mod"], thy_dct)
    inp_key_dct = ioparser.run.input_dictionary(input["run"])
    pes_idx_dct, spc_idx_dct = ioparser.run.chem_idxs(input["run"])
    tsk_lst_dct = ioparser.run.tasks(input["run"], thy_dct)
    spc_dct, glob_dct = ioparser.spc.species_dictionary(
        input["spc"], input["dat"], input["geo"], input["act"], inp_key_dct, "csv"
    )
    pes_dct = ioparser.mech.pes_dictionary(input["mech"], "chemkin", spc_dct)
    pes_rlst, spc_rlst = ioparser.rlst.run_lst(
        pes_dct, spc_dct, pes_idx_dct, spc_idx_dct
    )
    # Do a check
    ioprinter.info_message("\nFinal check if all required input provided...")
    ioparser.run.check_inputs(tsk_lst_dct, pes_dct, kmod_dct, smod_dct)

    # Build the Run-Save Filesystem Directories
    prefix_fs(inp_key_dct["run_prefix"], inp_key_dct["save_prefix"])

    # Run Drivers Requested by User
    es_tsks = tsk_lst_dct.get("es")
    if es_tsks is not None:
        ioprinter.program_header("es")
        esdriver.run(
            pes_rlst,
            spc_rlst,
            es_tsks,
            spc_dct,
            glob_dct,
            thy_dct,
            inp_key_dct["run_prefix"],
            inp_key_dct["save_prefix"],
            print_debug=inp_key_dct["print_debug"],
        )
        ioprinter.program_exit("es")

    therm_tsks = tsk_lst_dct.get("thermo")
    if therm_tsks is not None:
        ioprinter.program_header("thermo")
        thermodriver.run(
            pes_rlst,
            spc_rlst,
            therm_tsks,
            kmod_dct,
            smod_dct,
            spc_dct,
            thy_dct,
            inp_key_dct["run_prefix"],
            inp_key_dct["save_prefix"],
            path,
        )
        ioprinter.program_exit("thermo")

    trans_tsks = tsk_lst_dct.get("trans")
    if trans_tsks is not None:
        ioprinter.program_header("trans")
        if pes_dct:
            transdriver.run(
                pes_rlst,
                spc_rlst,
                trans_tsks,
                kmod_dct,
                smod_dct,
                spc_dct,
                thy_dct,
                inp_key_dct["run_prefix"],
                inp_key_dct["save_prefix"],
                path,
            )
        ioprinter.program_exit("trans")

    ktp_tsks = tsk_lst_dct.get("ktp")
    if ktp_tsks is not None:
        ioprinter.program_header("ktp")
        ktpdriver.run(
            pes_rlst,
            input["pesgrp"],
            ktp_tsks,
            spc_dct,
            glob_dct,
            thy_dct,
            kmod_dct,
            smod_dct,
            inp_key_dct["run_prefix"],
            inp_key_dct["save_prefix"],
            path,
        )
        ioprinter.program_exit("ktp")

    proc_tsks = tsk_lst_dct.get("proc")
    if proc_tsks is not None:
        ioprinter.program_header("proc")
        procdriver.run(
            pes_rlst,
            spc_rlst,
            proc_tsks,
            spc_dct,
            thy_dct,
            kmod_dct,
            smod_dct,
            inp_key_dct["run_prefix"],
            inp_key_dct["save_prefix"],
            path,
        )
        ioprinter.program_exit("proc")

    # Check if any drivers were requested to be run
    if all(
        tsks is None for tsks in (es_tsks, therm_tsks, trans_tsks, ktp_tsks, proc_tsks)
    ):
        ioprinter.warning_message(
            "User did not provide (uncommented) driver tasks lists in run.dat"
        )

    # Exit Program
    ioprinter.obj("vspace")
    ioprinter.program_exit("amech")
