""" Central Execution script to launch a MechDriver process which will
    parse all of the user-supplied input files in a specified directory, then
    launches all of the requested electronic structure, transport,
    thermochemistry and kinetics calculations via their associated
    sub-drivers.
"""

import sys
from mechlib.filesys import prefix_fs
from mechlib.amech_io import parser as ioparser
from mechlib.amech_io import printer as ioprinter
from drivers import esdriver, thermodriver, ktpdriver, transdriver, procdriver


# Set runtime options based on user input
JOB_PATH = sys.argv[1]  # Add a check to see if [1] exists; path exits

# Print the header message and host name (probably combine into one function)
ioprinter.program_header('amech')
ioprinter.random_cute_animal()
ioprinter.host_name()

# Parse all of the input
ioprinter.program_header('inp')

INP_STRS = ioparser.read_amech_input(JOB_PATH)

THY_DCT = ioparser.thy.theory_dictionary(INP_STRS['thy'])
KMOD_DCT, SMOD_DCT = ioparser.models.models_dictionary(
    INP_STRS['mod'], THY_DCT)
INP_KEY_DCT = ioparser.run.input_dictionary(INP_STRS['run'])
PES_IDX_DCT = ioparser.run.pes_idxs(INP_STRS['run'])
SPC_IDX_DCT = ioparser.run.spc_idxs(INP_STRS['run'])
TSK_LST_DCT = ioparser.run.tasks(INP_STRS['run'], THY_DCT, KMOD_DCT, SMOD_DCT)
SPC_DCT, GLOB_DCT = ioparser.spc.species_dictionary(
    INP_STRS['spc'], INP_STRS['dat'], INP_STRS['geo'], 'csv')
PES_DCT = ioparser.mech.pes_dictionary(
    INP_STRS['mech'], 'chemkin', SPC_DCT)

PES_RLST, SPC_RLST = ioparser.rlst.run_lst(
    PES_DCT, SPC_DCT, PES_IDX_DCT, SPC_IDX_DCT)

# Build the Run-Save Filesystem Directories
prefix_fs(INP_KEY_DCT['run_prefix'], INP_KEY_DCT['save_prefix'])

# Run Drivers Requested by User
ES_TSKS = TSK_LST_DCT.get('es')
if ES_TSKS is not None:
    ioprinter.program_header('es')
    esdriver.run(
        PES_RLST, SPC_RLST,
        ES_TSKS,
        SPC_DCT, GLOB_DCT, THY_DCT,
        INP_KEY_DCT['run_prefix'], INP_KEY_DCT['save_prefix']
    )
    ioprinter.program_exit('es')

THERM_TSKS = TSK_LST_DCT.get('thermo')
if THERM_TSKS is not None:
    ioprinter.program_header('thermo')
    thermodriver.run(
        SPC_RLST,
        THERM_TSKS,
        KMOD_DCT, SMOD_DCT,
        SPC_DCT,
        INP_KEY_DCT['run_prefix'], INP_KEY_DCT['save_prefix']
    )
    ioprinter.program_exit('thermo')

TRANS_TSKS = TSK_LST_DCT.get('trans')
if TRANS_TSKS is not None:
    ioprinter.program_header('trans')
    if PES_DCT:
        transdriver.run(
            SPC_RLST,
            SPC_DCT,
            THY_DCT,
            TRANS_TSKS,
            INP_KEY_DCT
        )
    ioprinter.program_exit('trans')

KTP_TSKS = TSK_LST_DCT.get('ktp')
if KTP_TSKS is not None:
    ioprinter.program_header('ktp')
    ktpdriver.run(
        PES_RLST,
        KTP_TSKS,
        SPC_DCT, GLOB_DCT,
        KMOD_DCT, SMOD_DCT, THY_DCT,
        INP_KEY_DCT['run_prefix'], INP_KEY_DCT['save_prefix']
    )
    ioprinter.program_exit('ktp')

PROC_TSKS = TSK_LST_DCT.get('proc')
if PROC_TSKS is not None:
    ioprinter.program_header('proc')
    PES_IDX = None
    procdriver.run(
        PES_RLST, SPC_RLST,
        SPC_DCT,
        PROC_TSKS,
        THY_DCT,
        INP_KEY_DCT,
        SMOD_DCT
    )
    ioprinter.program_exit('proc')

# Exit Program
ioprinter.obj('vspace')
ioprinter.program_exit('amech')
