"""
   Main Driver to parse and sort the mechanism input files and
   launch the desired drivers
"""

import sys
from mechlib.amech_io.parser import read_run, read_spc, read_thy, read_mech, read_model
from mechlib.amech_io import printer as ioprinter


# Set runtime options based on user input
JOB_PATH = sys.argv[1]

# Print the header message and host name (probably combine into one function)
ioprinter.program_header('amech')
ioprinter.random_cute_animal()
ioprinter.host_name()

# Parse all of the input
ioprinter.program_header('inp')
a = read_run(JOB_PATH) 
b = read_spc(JOB_PATH) 
c = read_thy(JOB_PATH) 
d = read_model(JOB_PATH) 
e = read_mech(JOB_PATH) 


# THY_DCT = parser.theory.build_thy_dct(JOB_PATH)
# PES_MODEL_DCT, SPC_MODEL_DCT = parser.model.read_models_sections(JOB_PATH)
# SPC_DCT = parser.species.build_spc_dct(JOB_PATH, 'csv')
# RUN_PES_DCT = parser.mechanism.build_pes_dct(
# SUB_SCRIPT_DCT = build_sub_script_dct(JOB_PATH)


