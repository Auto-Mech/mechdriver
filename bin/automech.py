"""
   Main Driver to parse and sort the mechanism input files and
   launch the desired drivers
"""

import sys
from mechlib.filesys import prefix_fs
from mechlib.amech_io import parser as ioparser
from mechlib.amech_io import printer as ioprinter


# Set runtime options based on user input
JOB_PATH = sys.argv[1]

# Print the header message and host name (probably combine into one function)
ioprinter.program_header('amech')
ioprinter.random_cute_animal()
ioprinter.host_name()

# Parse all of the input
ioprinter.program_header('inp')

INP_STRS = ioparser.read_amech_input(JOB_PATH)

THY_DCT = ioparser.thy.theory_dictionary(INP_STRS['thy'])
KMOD_DCT, SMOD_DCT = ioparser.models.models_dictionary(INP_STRS['mod'])
INP_KEY_DCT = ioparser.run.input_dictionary(INP_STRS['run'])
PES_IDX_DCT = ioparser.run.pes_idxs(INP_STRS['run'])
SPC_IDX_DCT = ioparser.run.spc_idxs(INP_STRS['run'])
TSK_LST_DCT = ioparser.run.tasks(INP_STRS['run'], THY_DCT, KMOD_DCT, SMOD_DCT)
SPC_DCT = ioparser.spc.species_dictionary(
    INP_STRS['spc'], INP_STRS['dat'], INP_STRS['geo'])
PES_DCT = ioparser.mech.pes_dictionary(INP_STRS['mech'], 'chemkin', SPC_DCT, PES_IDX_DCT)

print('kin_mod\n', KMOD_DCT)
print('spc_mod\n', SMOD_DCT)
print('thy dct\n', THY_DCT)
print('spc_dct\n', SPC_DCT)
print('inp_dct\n', INP_KEY_DCT)
print('pes_dct\n', PES_IDX_DCT)
print('spc_dct\n', SPC_IDX_DCT)

# Initialize the filesystem
prefix_fs(INP_KEY_DCT['run_prefix'], INP_KEY_DCT['save_prefix'])

print('FINISH')
import sys
sys.exit()

# Run Drivers
ES_TSKS = TSK_LST_DCT.get('es')
if ES_TSKS is not None:

    ioprinter.program_header('es')
    if PES_IDX_DCT:
        for (formula, pes_idx, sub_pes_idx), rxn_lst in RUN_PES_DCT.items():
            ioprinter.pes(pes_idx, formula, sub_pes_idx)
            for rxn in rxn_lst:
                ioprinter.channel(rxn['chn_idx'], rxn['reacs'], rxn['prods'])
            esdriver.run(
                pes_idx,
                rxn_lst,
                SPC_DCT,
                ES_TSK_LST,
                THY_DCT,
                RUN_INP_DCT
            )
    else:
        PES_IDX = 0
        esdriver.run(
            PES_IDX,
            RUN_SPC_LST_DCT,
            SPC_DCT,
            ES_TSK_LST,
            THY_DCT,
            RUN_INP_DCT
        )
    ioprinter.program_exit('es')

# ThermoDriver
THERM_TSKS = TSK_LST_DCT.get('thermo')
if THERM_TSKS is not None:

    ioprinter.program_header('thermo')

    # Call ThermoDriver for spc in PES
    for spc in RUN_SPC_LST_DCT:
        ioprinter.info_message(
            'Calculating Thermochem for species: {}'.format(spc),
            newline=1)
    thermodriver.run(
        SPC_DCT,
        PES_MODEL_DCT, SPC_MODEL_DCT,
        THY_DCT,
        RUN_SPC_LST_DCT,
        RUN_INP_DCT,
    )

    ioprinter.program_exit('thermo')

# TransportDriver
TRANS_TSKS = TSK_LST_DCT.get('trans')
if TRANS_TSKS is not None:

    ioprinter.program_header('trans')

    # Call ThermoDriver for spc in PES
    if RUN_OBJ_DCT['pes']:
        for spc in RUN_SPC_LST_DCT:
            ioprinter.message(
                'Calculating Transport for species: {}'.format(spc), newline=1)
            transdriver.run(
                SPC_DCT,
                THY_DCT,
                RUN_SPC_LST_DCT,
                TRANS_TSK_LST,
                RUN_INP_DCT
            )

# kTPDriver
KTP_TSKS = TSK_LST_DCT.get('ktp')
if KTP_TSKS is not None:

    ioprinter.program_header('ktp')

    if RUN_OBJ_DCT['pes']:
        for (formula, pes_idx, sub_pes_idx), rxn_lst in RUN_PES_DCT.items():

            # Print PES form and SUB PES Channels
            ioprinter.pes(pes_idx, formula, sub_pes_idx)
            for chn_idx, rxn in enumerate(rxn_lst):
                ioprinter.channel(
                    rxn['chn_idx'],
                    '+'.join(rxn['reacs']),
                    '+'.join(rxn['prods']))

            ktpdriver.run(
                formula, pes_idx, sub_pes_idx,
                SPC_DCT,
                THY_DCT,
                rxn_lst,
                PES_MODEL_DCT, SPC_MODEL_DCT,
                RUN_INP_DCT,
            )
    else:
        ioprinter.error_message("Can't run kTPDriver without a PES being specified")

    ioprinter.program_exit('ktp')

# Proc Driver
PROC_TSKS = TSK_LST_DCT.get('proc')
if PROC_TSKS is not None:

    ioprinter.program_header('proc')

    PES_IDX = 0
    printdriver.run(
        PES_IDX,
        RUN_SPC_LST_DCT,
        SPC_DCT,
        PRINT_TSK_LST,
        THY_DCT,
        RUN_INP_DCT,
        SPC_MODEL_DCT
    )

    ioprinter.program_exit('proc')

# Exit Program
ioprinter.obj('vspace')
ioprinter.program_exit('amech')
