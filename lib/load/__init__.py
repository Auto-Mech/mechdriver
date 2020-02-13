""" Libraries of functions that parse the moldriver input files
"""

from lib.load import run
from lib.load import model
from lib.load import theory
from lib.load import species
from lib.load import mechanism
from lib.load import tsks
from lib.load import ptt
from lib.load import keywords


__all__ = [
    'run',
    'model',
    'theory',
    'species',
    'mechanism',
    'tsks',
    'ptt',
    'keywords'
]

# def parse_automech_input():
#     """ Parse all of the AutoMechanic input files
#     """
#
#     # Parse the Parse the run file
#     print('\n\nrun.dat file contents')
#     run_inp_file_str = run.read_run_file()
#
#     # Get the input section from the run,dat file
#     run_inp_block_str = run.input_block(run_inp_file_str)
#     run_inp_dct = run.build_run_input_keyword_dct(run_inp_block_str)
#     print(run_inp_dct)
#     ids = run_inp_dct['ids']
#     mech_typ = run_inp_dct['mech']
#     spc_typ = run_inp_dct['spc']
#
#     # Read the obj section from the run.dat file
#     run_inp_obj_str = run.object_block(run_inp_file_str)
#     run_obj_lst = run.objects_lst(run_inp_obj_str)
#     print(run_obj_lst)
#
#     # Read the proc sections from the run.dat file
#     proc_dct = run.read_proc_sections(run_inp_file_str)
#     print(proc_dct)
#
#     # Read in the contents of the species file and get dct
#     print('\n\nspecies.dat file contents')
#     spc_str = species.read_species_file()
#     spc_dct = species.build_spc_dct(spc_str, spc_typ)
#     print(spc_dct)
#
#     # Read in the contents of the mechanism file
#     mech_str = mechanism.read_mech_file()
#     pes_dct = species.build_spc_dct(mech_str, mech_typ)
#     print(pes_dct)
#
#     # Read the contents of the model.dat file
#     print('\n\nmodel.dat file contents')
#     model_inp_file_str = model.read_model_file()
#     model_dct = model.read_models_sections(model_inp_file_str)
#     print(model_dct)
#
#     # Read contents of theory.dat file and build dct of els methods
#     print('\n\ntheory.dat file contents')
#     thy_inp_file_str = theory.read_theory_file()
#     thy_dct = theory.read_theory_sections(thy_inp_file_str)
#     print(thy_dct)
#
#     return spc_dct, pes_dct, model_dct, thy_dct, run_obj_lst, proc_dct
