""" reaction list test
"""
import os
import pandas
from qcelemental import constants as qcc
import chemkin_io
import automol
import moldr
import scripts

ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')
WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

# 0. choose which mechanism to run

# MECHANISM_NAME = 'ch4+nh2'  # options: syngas, natgas, heptane
MECHANISM_NAME = 'test'  # options: syngas, natgas, heptane
#MECHANISM_NAME = 'vhp'  # options: syngas, natgas, heptane
# MECHANISM_NAME = 'onereac'  # options: syngas, natgas, heptane
# MECHANISM_NAME = 'estoktp/add30'  # options: syngas, natgas
# MECHANISM_NAME = 'estoktp/habs65'  # options: syngas, natgas

# 1. script control parameters

# a. Strings to launch executable
# script_strings for electronic structure are obtained from run_qchem_par since
# they vary with method

PROJROT_SCRIPT_STR = ("#!/usr/bin/env bash\n"
                      "RPHt.exe")
PF_SCRIPT_STR = ("#!/usr/bin/env bash\n"
                 "messpf pf.inp build.out >> stdout.log &> stderr.log")
RATE_SCRIPT_STR = ("#!/usr/bin/env bash\n"
                   "mess mess.inp build.out >> stdout.log &> stderr.log")
VARECOFMCFLUX_SCRIPT_STR = {"#!/usr/bin/env bash\n"
                     "/home/ygeorgi/build/rotd/multi ")
MCFLUX_SCRIPT_STR = {"#!/usr/bin/env bash\n"
                     "/home/ygeorgi/build/rotd/mc_flux ")
CONV_MULTI_SCRIPT_STR = {"#!/usr/bin/env bash\n"
                     "/home/ygeorgi/build/rotd/mc_flux ")
TST_CHECK_SCRIPT_STR = {"#!/usr/bin/env bash\n"
                     "/home/ygeorgi/build/rotd/tst_check ")
MOLPRO_PATH_STR = ('/home/sjklipp/bin/molpro')
#NASA_SCRIPT_STR = ("#!/usr/bin/env bash\n"
#                   "cp ../PF/build.out pf.dat\n"
#                   "cp /tcghome/sjklipp/PACC/nasa/new.groups .\n"
#                   "python /tcghome/sjklipp/PACC/nasa/makepoly.py"
#                   " >> stdout.log &> stderr.log")

# b. Electronic structure parameters; code, method, basis, convergence control

# use reference to determine starting geometries, which are used to define
# z-matrices also used for reference geometries in HL calculations

# code is designed to run
# (i) an arbitrary number of low level optimizations
# (ii) an arbitrary number of high level energies at geometries corresponding to a single reference method
# the code presumes that the optimizations for the reference method already exist
# the low level optimizations are used to determine the temperature dependence of the partition functions
# the high level energies are used to determine the energy part of the 0 K Heat of formation
# the latter requires the specification of some reference species and reference Heats of formation
# set up standard prog, method, basis for opts
# RU => RHF for singlet, UHF for multiplets
# RR => RHF for singlet, RHF for multiplets
# UU => UHF for singlet, UHF for multiplets

RUN_OPT_LEVELS = []
#RUN_OPT_LEVELS.append(['molpro', 'casscf', 'cc-pVDZ', 'RR'])
RUN_OPT_LEVELS.append(['g09', 'wb97xd', '6-31g*', 'RU'])
# RUN_OPT_LEVELS.append(['g09', 'wb97xd', 'cc-pVTZ', 'RU'])
# RUN_OPT_LEVELS.append(['g09', 'b2plypd3', 'cc-pVTZ', 'RU'])
# RUN_OPT_LEVELS.append(['g09', 'wb97xd', 'cc-pVTZ', 'RU'])
# RUN_OPT_LEVELS.append(['g09', 'm062x', 'cc-pVTZ', 'RU'])
# RUN_OPT_LEVELS.append(['g09', 'b3lyp', '6-31g*', 'RU'])

# set up a set of standard hl methods
THEORY_REF_HIGH_LEVEL = ['', 'wb97xd', '6-31g*', 'RU']
RUN_HIGH_LEVELS = []
#RUN_HIGH_LEVELS.append(['molpro', 'mp2', 'cc-pVTZ'])
#RUN_HIGH_LEVELS.append(['molpro', 'CCSD(T)', 'cc-pVDZ', 'RR'])
#RUN_HIGH_LEVELS.append(['molpro', 'CCSD(T)', 'cc-pVTZ', 'RR'])
#RUN_HIGH_LEVELS.append(['molpro', 'CCSD(T)', 'cc-pVQZ', 'RR'])
#RUN_HIGH_LEVELS.append(['molpro', 'CCSD(T)', 'cc-pV5Z', 'RR'])
#RUN_HIGH_LEVELS.append(['molpro', 'CCSD(T)-F12', 'cc-pVDZ-F12', 'RR'])
#RUN_HIGH_LEVELS.append(['molpro', 'CCSD(T)-F12', 'cc-pVTZ-F12', 'RR'])
#RUN_HIGH_LEVELS.append(['molpro', 'CCSD(T)-F12', 'cc-pVQZ-F12', 'RR'])
#RUN_HIGH_LEVELS.append(['molpro', 'CCSD(T)', 'aug-cc-pVDZ', 'RR'])
#RUN_HIGH_LEVELS.append(['molpro', 'CCSD(T)', 'aug-cc-pVTZ', 'RR'])
#RUN_HIGH_LEVELS.append(['molpro', 'CCSD(T)', 'aug-cc-pVQZ', 'RR'])
#RUN_HIGH_LEVELS.append(['molpro', 'CCSD(T)', 'aug-cc-pV5Z', 'RR'])

PF_LEVELS = []
#PF_LEVELS.append([['', 'wb97xd', '6-31g*', 'RU'], ['', 'wb97xd', '6-31g*', 'RU'], ['', 'wb97xd', '6-31g*', 'RU']])
# PF_LEVELS contains the elec. struc. levels for the harmonic, torsional, and anharmonic analyses

# c. What type of electronic structure calculations to run
RUN_SPC_QCHEM = True
RUN_TS_QCHEM = True
RUN_TS_KICKS_QCHEM = False
RUN_VDW_QCHEM = False

RUN_INI_GEOM = True
RUN_REMOVE_IMAG = True

KICKOFF_SADDLE = True

RUN_CONF_SAMP = True
RUN_CONF_OPT = True
RUN_MIN_GRAD = False
RUN_MIN_HESS = True
RUN_MIN_VPT2 = False
RUN_CONF_GRAD = False
RUN_CONF_HESS = False
RUN_CONF_VPT2 = False
RUN_CONF_SCAN = True
RUN_CONF_SCAN_GRAD = False
RUN_CONF_SCAN_HESS = False
RUN_TAU_SAMP = False
RUN_TAU_GRAD = False
RUN_TAU_HESS = False
RUN_HL_MIN_ENE = True

RUN_TS_CONF_SAMP = True
RUN_TS_CONF_OPT = False
RUN_TS_MIN_GRAD = False
RUN_TS_MIN_HESS = False
RUN_TS_MIN_VPT2 = False
RUN_TS_CONF_GRAD = False
RUN_TS_CONF_HESS = False
RUN_TS_CONF_VPT2 = False
RUN_TS_CONF_SCAN = True
RUN_TS_CONF_SCAN_GRAD = False
RUN_TS_CONF_SCAN_HESS = False
RUN_TS_TAU_SAMP = False
RUN_TS_TAU_GRAD = False
RUN_TS_TAU_HESS = False
RUN_TS_HL_MIN_ENE = True

RUN_VDW_CONF_SAMP = False
RUN_VDW_CONF_OPT = False
RUN_VDW_MIN_GRAD = False
RUN_VDW_MIN_HESS = False
RUN_VDW_MIN_VPT2 = False
RUN_VDW_CONF_GRAD = False
RUN_VDW_CONF_HESS = False
RUN_VDW_CONF_VPT2 = False
RUN_VDW_CONF_SCAN = False
RUN_VDW_CONF_SCAN_GRAD = False
RUN_VDW_CONF_SCAN_HESS = False
RUN_VDW_TAU_SAMP = False
RUN_VDW_TAU_GRAD = False
RUN_VDW_TAU_HESS = False
RUN_VDW_HL_MIN_ENE = True

# d. Parameters for number of torsional samplings
NSAMP_CONF = 5
NSAMP_CONF_EXPR = True
NSAMP_CONF_A = 3
NSAMP_CONF_B = 1
NSAMP_CONF_C = 4
NSAMP_CONF_D = 100
NSAMP_CONF_PAR = [NSAMP_CONF_EXPR, NSAMP_CONF_A, NSAMP_CONF_B, NSAMP_CONF_C,
                  NSAMP_CONF_D, NSAMP_CONF]

NSAMP_TAU = 10
NSAMP_TAU_EXPR = False
NSAMP_TAU_A = 3
NSAMP_TAU_B = 1
NSAMP_TAU_C = 3
NSAMP_TAU_D = 15
NSAMP_TAU_PAR = [NSAMP_TAU_EXPR, NSAMP_TAU_A, NSAMP_TAU_B, NSAMP_TAU_C,
                 NSAMP_TAU_D, NSAMP_TAU]

NSAMP_TS_CONF = 5
NSAMP_TS_CONF_EXPR = True
NSAMP_TS_CONF_A = 3
NSAMP_TS_CONF_B = 1
NSAMP_TS_CONF_C = 3
NSAMP_TS_CONF_D = 100
NSAMP_TS_CONF_PAR = [NSAMP_TS_CONF_EXPR, NSAMP_TS_CONF_A, NSAMP_TS_CONF_B, NSAMP_TS_CONF_C,
                     NSAMP_TS_CONF_D, NSAMP_TS_CONF]

NSAMP_TS_TAU = 5
NSAMP_TS_TAU_EXPR = True
NSAMP_TS_TAU_A = 3
NSAMP_TS_TAU_B = 1
NSAMP_TS_TAU_C = 3
NSAMP_TS_TAU_D = 100
NSAMP_TS_TAU_PAR = [NSAMP_TS_TAU_EXPR, NSAMP_TS_TAU_A, NSAMP_TS_TAU_B, NSAMP_TS_TAU_C,
                    NSAMP_TS_TAU_D, NSAMP_TS_TAU]

NSAMP_VDW = 10
NSAMP_VDW_EXPR = False
NSAMP_VDW_A = 3
NSAMP_VDW_B = 1
NSAMP_VDW_C = 3
NSAMP_VDW_D = 15
NSAMP_VDW_PAR = [NSAMP_VDW_EXPR, NSAMP_VDW_A, NSAMP_VDW_B, NSAMP_VDW_C,
                 NSAMP_VDW_D, NSAMP_VDW]

# e. What to run for thermochemical kinetics
RUN_SPC_THERMO = True
SPC_MODELS = [['RIGID', 'HARM']]
#SPC_MODELS = [['1DHR', 'HARM']]
#SPC_MODELS = [['RIGID', 'HARM'], ['1DHR', 'HARM']]
# The first component specifies the torsional model - TORS_MODEL.
# It can take 'RIGID', '1DHR', or 'TAU'
# The second component specifies the vibrational model - VIB_MODEL.
# It can take 'HARM', or 'VPT2' values.
RUN_REACTION_RATES = True
RUN_VDW_RCT_RATES = False
RUN_VDW_PRD_RATES = False

# f. Partition function parameters
TAU_PF_WRITE = True

# Defaults and dictionaries
SCAN_INCREMENT = 30. * qcc.conversion_factor('degree', 'radian')
KICKOFF_SIZE = 0.1
KICKOFF_BACKWARD = False
RESTRICT_OPEN_SHELL = False
OVERWRITE = False
# Temperatue and pressure
TEMP_STEP = 100.
NTEMPS = 30
TEMPS = [300., 500., 750., 1000., 1250., 1500., 1750., 2000.]
PRESS = [0.1, 1., 10., 100.]
# Collisional parameters
EXP_FACTOR = 150.0
EXP_POWER = 50.0
EXP_CUTOFF = 80.0
EPS1 = 100.0
EPS2 = 200.0
SIG1 = 10.0
SIG2 = 20.0
MASS1 = 15.0
MASS2 = 25.0
ETSFR_PAR = [EXP_FACTOR, EXP_POWER, EXP_CUTOFF, EPS1, EPS2, SIG1, SIG2, MASS1, MASS2]

ELC_DEG_DCT = {
    ('InChI=1S/B', 2): [[0., 2], [16., 4]],
    ('InChI=1S/C', 3): [[0., 1], [16.4, 3], [43.5, 5]],
    ('InChI=1S/N', 2): [[0., 6], [8., 4]],
    ('InChI=1S/O', 3): [[0., 5], [158.5, 3], [226.5, 1]],
    ('InChI=1S/F', 2): [[0., 4], [404.1, 2]],
    ('InChI=1S/Cl', 2): [[0., 4], [883.4, 2]],
    ('InChI=1S/Br', 2): [[0., 4], [685.2, 2]],
    ('InChI=1S/HO/h1H', 2): [[0., 2], [138.9, 2]],
    ('InChI=1S/NO/c1-2', 2): [[0., 2], [123.1, 2]],
    ('InChI=1S/O2/c1-2', 1): [[0., 2]]
}

ELC_SIG_LST = {'InChI=1S/CN/c1-2', 'InChI=1S/C2H/c1-2/h1H'}
#SMILES_TEST = '[C]#C'
#for smiles in SMILES_TEST:
#    ICH_TEST = (automol.convert.smiles.inchi(SMILES_TEST))
#    print(SMILES_TEST, ICH_TEST)

# 2. create run and save directories
RUN_PREFIX = '/lcrc/project/PACC/run'
if not os.path.exists(RUN_PREFIX):
    os.mkdir(RUN_PREFIX)

SAVE_PREFIX = '/lcrc/project/PACC/save'
if not os.path.exists(SAVE_PREFIX):
    os.mkdir(SAVE_PREFIX)

# 3. read in data from the mechanism directory
DATA_PATH = os.path.dirname(os.path.realpath(__file__))
MECH_PATH = os.path.join(DATA_PATH, 'data', MECHANISM_NAME)
MECH_STR = open(os.path.join(MECH_PATH, 'mechanism.txt')).read()
SPC_TAB = pandas.read_csv(os.path.join(MECH_PATH, 'smiles.csv'))

# 4. process species data from the mechanism file
# Also add in basis set species

SPC_TAB['charge'] = 0
SMI_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['smiles']))
SMI_DCT['REF_H2'] = '[H][H]'
SMI_DCT['REF_CH4'] = 'C'
SMI_DCT['REF_H2O'] = 'O'
SMI_DCT['REF_NH3'] = 'N'
CHG_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['charge']))
CHG_DCT['REF_H2'] = 0
CHG_DCT['REF_CH4'] = 0
CHG_DCT['REF_H2O'] = 0
CHG_DCT['REF_NH3'] = 0
MUL_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['mult']))
MUL_DCT['REF_H2'] = 1
MUL_DCT['REF_CH4'] = 1
MUL_DCT['REF_H2O'] = 1
MUL_DCT['REF_NH3'] = 1
SPC_BLK_STR = chemkin_io.species_block(MECH_STR)
SPC_NAMES = chemkin_io.species.names(SPC_BLK_STR)

# You need one species reference for each element in the set of references and species
SPC_REF_NAMES = ('REF_H2', 'REF_CH4', 'REF_H2O', 'REF_NH3')
SPC_REF_ICH = []
for ref_name in SPC_REF_NAMES:
    smi = SMI_DCT[ref_name]
    SPC_REF_ICH.append(automol.smiles.inchi(smi))
SPC_NAMES += SPC_REF_NAMES
print('SPC_NAMES')
print(SPC_NAMES)
SPC_INFO = {}
for name in SPC_NAMES:
    smi = SMI_DCT[name]
    ich = automol.smiles.inchi(smi)
    chg = CHG_DCT[name]
    mul = MUL_DCT[name]
    SPC_INFO[name] = [ich, chg, mul]

SPC_STR = {}

RXN_BLOCK_STR = chemkin_io.reaction_block(MECH_STR)
RXN_STRS = chemkin_io.reaction.data_strings(RXN_BLOCK_STR)
RCT_NAMES_LST = list(
    map(chemkin_io.reaction.DataString.reactant_names, RXN_STRS))
PRD_NAMES_LST = list(
    map(chemkin_io.reaction.DataString.product_names, RXN_STRS))

GEOM_PATH = os.path.join(DATA_PATH, 'data', 'geoms')
print(GEOM_PATH)
GEOM_DCT = moldr.util.geometry_dictionary(GEOM_PATH)

# take starting geometry from saved directory if possible, otherwise get it
# from inchi via rdkit
SPC_QCHEM_RUN_FLAGS = [
    RUN_INI_GEOM, RUN_REMOVE_IMAG, RUN_CONF_SAMP, RUN_MIN_GRAD, RUN_MIN_HESS,
    RUN_MIN_VPT2, RUN_CONF_SCAN, RUN_CONF_GRAD, RUN_CONF_HESS, RUN_TAU_SAMP,
    RUN_TAU_GRAD, RUN_TAU_HESS, RUN_HL_MIN_ENE]

TS_QCHEM_RUN_FLAGS = [
    RUN_TS_CONF_SAMP, RUN_TS_MIN_GRAD, RUN_TS_MIN_HESS,
    RUN_TS_MIN_VPT2, RUN_TS_CONF_SCAN, RUN_TS_CONF_GRAD, RUN_TS_CONF_HESS, RUN_TS_TAU_SAMP,
    RUN_TS_TAU_GRAD, RUN_TS_TAU_HESS, RUN_TS_HL_MIN_ENE, RUN_TS_KICKS_QCHEM]

VDW_QCHEM_RUN_FLAGS = [
    RUN_VDW_CONF_SAMP, RUN_VDW_MIN_GRAD, RUN_VDW_MIN_HESS,
    RUN_VDW_MIN_VPT2, RUN_VDW_CONF_SCAN, RUN_VDW_CONF_GRAD, RUN_VDW_CONF_HESS, RUN_VDW_TAU_SAMP,
    RUN_VDW_TAU_GRAD, RUN_VDW_TAU_HESS, RUN_VDW_HL_MIN_ENE]

SPC_NSAMP_PARS = [NSAMP_CONF_PAR, NSAMP_TAU_PAR]
TS_NSAMP_PARS = [NSAMP_TS_CONF_PAR, NSAMP_TS_TAU_PAR]
VDW_NSAMP_PARS = [NSAMP_VDW_PAR]
#NSAMP_PARS = [NSAMP_CONF_PAR, NSAMP_TAU_PAR, NSAMP_VDW_PAR, NSAMP_TS_CONF_PAR, NSAMP_TS_TAU_PAR]
KICKOFF_PARS = [KICKOFF_BACKWARD, KICKOFF_SIZE]

#    nsamp_vdw_par = nsamp_pars[2]
#    nsamp_ts_conf_par = nsamp_pars[3]
#    nsamp_ts_taupar = nsamp_pars[4]

if RUN_SPC_QCHEM:
    scripts.es.species_qchem(
        spc_names=SPC_NAMES,
        spc_info=SPC_INFO,
        run_opt_levels=RUN_OPT_LEVELS,
        ref_high_level=THEORY_REF_HIGH_LEVEL,
        run_high_levels=RUN_HIGH_LEVELS,
        geom_dct=GEOM_DCT,
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,
        qchem_flags=SPC_QCHEM_RUN_FLAGS,
        nsamp_pars=SPC_NSAMP_PARS,
        scan_increment=SCAN_INCREMENT,
        kickoff_pars=KICKOFF_PARS,
        overwrite=OVERWRITE,
        )

if RUN_TS_QCHEM:
    scripts.es.ts_qchem(
        rct_names_lst=RCT_NAMES_LST,
        prd_names_lst=PRD_NAMES_LST,
        smi_dct=SMI_DCT,
        chg_dct=CHG_DCT,
        mul_dct=MUL_DCT,
        run_opt_levels=RUN_OPT_LEVELS,
        ref_high_level=THEORY_REF_HIGH_LEVEL,
        run_high_levels=RUN_HIGH_LEVELS,
        geom_dct=GEOM_DCT,
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,
        qchem_flags=TS_QCHEM_RUN_FLAGS,
        nsamp_pars=TS_NSAMP_PARS,
        scan_increment=SCAN_INCREMENT,
        kickoff_pars=KICKOFF_PARS,
        overwrite=OVERWRITE,
        )

if RUN_VDW_QCHEM:
    scripts.es.vdw_qchem(
        rct_names_lst=RCT_NAMES_LST,
        prd_names_lst=PRD_NAMES_LST,
        smi_dct=SMI_DCT,
        chg_dct=CHG_DCT,
        mul_dct=MUL_DCT,
        run_opt_levels=RUN_OPT_LEVELS,
        ref_high_level=THEORY_REF_HIGH_LEVEL,
        run_high_levels=RUN_HIGH_LEVELS,
        geom_dct=GEOM_DCT,
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,
        qchem_flags=VDW_QCHEM_RUN_FLAGS,
        nsamp_pars=VDW_NSAMP_PARS,
        scan_increment=SCAN_INCREMENT,
        kickoff_pars=KICKOFF_PARS,
        overwrite=OVERWRITE,
        )

TEMP_PAR = [TEMP_STEP, NTEMPS]

print('RUN_SPC_THERMO test:')
print(RUN_SPC_THERMO)
if RUN_SPC_THERMO:
    scripts.ktp.species_thermo(
        spc_names=SPC_NAMES,
        spc_info=SPC_INFO,
        spc_ref_names=SPC_REF_NAMES,
        elc_deg_dct=ELC_DEG_DCT,
        temp_par=TEMP_PAR,
        ref_high_level=THEORY_REF_HIGH_LEVEL,
        run_high_levels=RUN_HIGH_LEVELS,
        spc_models=SPC_MODELS,
        pf_levels=PF_LEVELS,
        save_prefix=SAVE_PREFIX,
        projrot_script_str=PROJROT_SCRIPT_STR,
        pf_script_str=PF_SCRIPT_STR,
        )

if RUN_REACTION_RATES:
    scripts.ktp.reaction_rates(
        rct_names_lst=RCT_NAMES_LST,
        prd_names_lst=PRD_NAMES_LST,
        smi_dct=SMI_DCT,
        chg_dct=CHG_DCT,
        mul_dct=MUL_DCT,
        elc_deg_dct=ELC_DEG_DCT,
        temperatures=TEMPS,
        pressures=PRESS,
        etsfr_par=ETSFR_PAR,
        run_opt_levels=RUN_OPT_LEVELS,
        ref_high_level=THEORY_REF_HIGH_LEVEL,
        run_high_levels=RUN_HIGH_LEVELS,
        spc_models=SPC_MODELS,
        pf_levels=PF_LEVELS,
        save_prefix=SAVE_PREFIX,
        projrot_script_str=PROJROT_SCRIPT_STR,
        rate_script_str=RATE_SCRIPT_STR,
        )
