""" reaction list test
"""
import os
import sys
import numpy
import pandas
from qcelemental import constants as qcc
import chemkin_io
import projrot_io
import thermo
import automol
import elstruct
import autofile
import moldr
import mess_io.writer

ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')
WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

# 0. choose which mechanism to run

# MECHANISM_NAME = 'ch4+nh2'  # options: syngas, natgas, heptane
MECHANISM_NAME = 'test'  # options: syngas, natgas, heptane
#MECHANISM_NAME = 'onereac'  # options: syngas, natgas, heptane
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
#NASA_SCRIPT_STR = ("#!/usr/bin/env bash\n"
#                   "cp ../PF/build.out pf.dat\n"
#                   "cp /tcghome/sjklipp/PACC/nasa/new.groups .\n"
#                   "python /tcghome/sjklipp/PACC/nasa/makepoly.py"
#                   " >> stdout.log &> stderr.log")

# b. Electronic structure parameters; code, method, basis, convergence control

# use reference to determine starting geometries, which are used to define
# z-matrices also used for reference geometries in HL calculations
#LL_METHOD_REF = 'wb97xd'
#LL_BASIS_REF = '6-31g*'
#LL_ORB_RESTR_REF = False
THEORY_REF_LOW_LEVEL = ['', 'wb97d', '6-31g*', 'RU']

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
RUN_OPT_LEVELS.append(['g09', 'wb97xd', '6-31g*', 'RU'])
# RUN_OPT_LEVELS.append(['g09', 'wb97xd', 'cc-pVTZ', 'RU'])
# RUN_OPT_LEVELS.append(['g09', 'b2plypd3', 'cc-pVTZ', 'RU'])
# RUN_OPT_LEVELS.append(['g09', 'wb97xd', 'cc-pVTZ', 'RU'])
# RUN_OPT_LEVELS.append(['g09', 'm062x', 'cc-pVTZ', 'RU'])
# RUN_OPT_LEVELS.append(['g09', 'b3lyp', '6-31g*', 'RU'])

# set up a set of standard hl methods
#HL_METHOD_REF = 'wb97xd'
#HL_BASIS_REF = '6-31g*'
#HL_ORB_RESTR_REF = 'RU'
THEORY_REF_HIGH_LEVEL = ['', 'wb97xd', '6-31g*', 'RU']
RUN_HIGH_LEVELS = []
#RUN_HIGH_LEVELS.append(['molpro', 'mp2', 'cc-pVTZ'])
#RUN_HIGH_LEVELS.append(['molpro', 'CCSD(T)', 'cc-pVDZ', 'RR'])
RUN_HIGH_LEVELS.append(['molpro', 'CCSD(T)', 'cc-pVTZ', 'RR'])
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
PF_LEVELS.append([['', 'wb97xd', '6-31g*', 'RU'], ['', 'wb97xd', '6-31g*', 'RU'], ['', 'wb97xd', '6-31g*', 'RU']])
# PF_LEVELS contains the elec. struc. levels for the harmonic, torsional, and anharmonic analyses
#HAR_LEVELS = []
#HAR_LEVELS.append(['wb97xd','6-31g*'])
#TORS = []
#TORS.append(['wb97xd','6-31g*'])
#ANH_LEVELS = []
#ANH_LEVELS.append(['wb97xd','6-31g*'])

# c. What type of electronic structure calculations to run
RUN_SPECIES_QCHEM = True
RUN_REACTIONS_QCHEM = True
RUN_TS_KICKS_QCHEM = True
RUN_VDW_QCHEM = False

RUN_INI_GEOM = True
RUN_REMOVE_IMAG = False

KICKOFF_SADDLE = True

RUN_CONF_SAMP = True
RUN_CONF_OPT = True
RUN_MIN_GRAD = True
RUN_MIN_HESS = True
RUN_MIN_VPT2 = False
RUN_CONF_GRAD = True
RUN_CONF_HESS = True
RUN_CONF_VPT2 = False

RUN_CONF_SCAN = False
RUN_CONF_SCAN_GRAD = False
RUN_CONF_SCAN_HESS = False

RUN_TAU_SAMP = False
RUN_TAU_GRAD = False
RUN_TAU_HESS = False

RUN_TS_CONF_SAMP = False
RUN_TS_CONF_OPT = False
RUN_TS_MIN_GRAD = False
RUN_TS_MIN_HESS = False
RUN_TS_MIN_VPT2 = False
RUN_TS_CONF_SCAN = False
RUN_TS_CONF_SCAN_GRAD = False
RUN_TS_CONF_SCAN_HESS = False
RUN_TS_TAU_SAMP = False

RUN_HL_MIN_ENE = True

# setting these to true turns on corresponding run for min, conf, conf_scan,
# and tau
RUN_GRAD = False
if RUN_GRAD:
    RUN_MIN_GRAD = True
    RUN_CONF_GRAD = True
    RUN_CONF_SCAN_GRAD = True
    RUN_TAU_GRAD = True

RUN_HESS = False
if RUN_HESS:
    RUN_MIN_HESS = True
    RUN_CONF_HESS = True
    RUN_CONF_SCAN_HESS = True
    RUN_TAU_HESS = True
RUN_GRAD_PF = False
RUN_HESS_PF = False

# d. Parameters for number of torsional samplings
NSAMP_CONF = 5
NSAMP_CONF_EXPR = True
NSAMP_CONF_A = 3
NSAMP_CONF_B = 1
NSAMP_CONF_C = 3
NSAMP_CONF_D = 100
NSAMP_CONF_PAR = [NSAMP_CONF_EXPR, NSAMP_CONF_A, NSAMP_CONF_B, NSAMP_CONF_C,
                  NSAMP_CONF_D, NSAMP_CONF]

# NSAMP_TAU = 100
NSAMP_TAU = 10
NSAMP_TAU_EXPR = False
NSAMP_TAU_A = 3
NSAMP_TAU_B = 1
NSAMP_TAU_C = 3
NSAMP_TAU_D = 15
NSAMP_TAU_PAR = [NSAMP_TAU_EXPR, NSAMP_TAU_A, NSAMP_TAU_B, NSAMP_TAU_C,
                 NSAMP_TAU_D, NSAMP_TAU]

NSAMP_VDW = 10
NSAMP_VDW_EXPR = False
NSAMP_VDW_A = 3
NSAMP_VDW_B = 1
NSAMP_VDW_C = 3
NSAMP_VDW_D = 15
NSAMP_VDW_PAR = [NSAMP_VDW_EXPR, NSAMP_VDW_A, NSAMP_VDW_B, NSAMP_VDW_C,
                 NSAMP_VDW_D, NSAMP_VDW]

# e. What to run for thermochemical kinetics
RUN_SPECIES_THERMO = False
#SPECIES_MODELS = [['1DHR', 'HARM']]
SPECIES_MODELS = [['RIGID', 'HARM']]
#SPECIES_MODELS = [['RIGID', 'HARM'], ['1DHR', 'HARM']]
#SPECIES_MODELS = [['1DHR', 'HARM']]
# The first component specifies the torsional model - TORS_MODEL.
# It can take 'RIGID', '1DHR', or 'TAU'
# The second component specifies the vibrational model - VIB_MODEL.
# It can take 'HARM', or 'VPT2' values.
RUN_REACTIONS_RATES = True
RUN_VDW_RCT_RATES = False
RUN_VDW_PRD_RATES = False

# f. Partition function parameters
TAU_PF_WRITE = True

# Defaults
SCAN_INCREMENT = 30. * qcc.conversion_factor('degree', 'radian')
KICKOFF_SIZE = 0.1
KICKOFF_BACKWARD = False
RESTRICT_OPEN_SHELL = False
OVERWRITE = False
# Temperatue and pressure
TEMPS = [300., 500., 750., 1000., 1250., 1500., 1750., 2000.]
TEMP_STEP = 100.
NTEMPS = 30
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
species_info = {}
for name in SPC_NAMES:
    smi = SMI_DCT[name]
    ich = automol.smiles.inchi(smi)
    chg = CHG_DCT[name]
    mul = MUL_DCT[name]
    species_info[name] = [ich, chg, mul]

species_str = {}

GEOM_PATH = os.path.join(DATA_PATH, 'data', 'geoms')
print(GEOM_PATH)
GEOM_DCT = moldr.util.geometry_dictionary(GEOM_PATH)

# take starting geometry from saved directory if possible, otherwise get it
# from inchi via rdkit

if RUN_SPECIES_QCHEM:
    for name in SPC_NAMES:
        # species
        print("Species: {}".format(name))
        ich = species_info[name][0]
        smi = automol.inchi.smiles(ich)
        print("smiles: {}".format(smi), "inchi: {}".format(ich))
        for opt_level_idx, _ in enumerate(RUN_OPT_LEVELS):
            # theory
            prog = RUN_OPT_LEVELS[opt_level_idx][0]
            SP_SCRIPT_STR, OPT_SCRIPT_STR, KWARGS, OPT_KWARGS = moldr.util.run_qchem_par(prog)

            orb_restr = moldr.util.orbital_restriction(
                species_info[name], RUN_OPT_LEVELS[opt_level_idx])
            thy_level = RUN_OPT_LEVELS[opt_level_idx][1:3]
            thy_level.append(orb_restr)

            # a. conformer sampling
            spc_run_fs = autofile.fs.species(RUN_PREFIX)
            spc_run_fs.leaf.create(species_info[name])
            spc_run_path = spc_run_fs.leaf.path(species_info[name])

            spc_save_fs = autofile.fs.species(SAVE_PREFIX)
            spc_save_fs.leaf.create(species_info[name])
            spc_save_path = spc_save_fs.leaf.path(species_info[name])

            thy_run_fs = autofile.fs.theory(spc_run_path)
            thy_run_fs.leaf.create(thy_level)
            thy_run_path = thy_run_fs.leaf.path(thy_level)

            thy_save_fs = autofile.fs.theory(spc_save_path)
            thy_save_fs.leaf.create(thy_level) 
            thy_save_path = thy_save_fs.leaf.path(thy_level)

            geo_init = moldr.util.reference_geometry(
                species_info=species_info[name],
                theory_level=RUN_OPT_LEVELS[opt_level_idx],
                prefix=SAVE_PREFIX,
                geom_dct=GEOM_DCT)

            # this uses theory run path - should start with a check in save path to see if initial geometry has already been saved
            # eventually theory data will be removed
            # also may need to remove hessian etc from saved geometry ...
            if RUN_INI_GEOM:
                geo = moldr.driver.run_initial_geometry_opt(
                    species_info=species_info[name],
                    theory_level=RUN_OPT_LEVELS[opt_level_idx],
                    run_prefix=spc_run_path,
                    save_prefix=spc_save_path,
                    script_str=OPT_SCRIPT_STR,
                    overwrite=OVERWRITE,
                    geo_init=geo_init,
                    **OPT_KWARGS,
                )

                if RUN_REMOVE_IMAG:
                    imag, geo, disp_xyzs = moldr.driver.run_check_imaginary(
                        species_info=species_info[name],
                        theory_level=RUN_OPT_LEVELS[opt_level_idx],
                        run_prefix=spc_run_path,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )
                    if imag:
                        moldr.driver.run_kickoff_saddle(
                            geo, disp_xyzs,
                            species_info=species_info[name],
                            theory_level=RUN_OPT_LEVELS[opt_level_idx],
                            run_path=thy_run_path,
                            script_str=OPT_SCRIPT_STR,
                            kickoff_backward=KICKOFF_BACKWARD,
                            kickoff_size=KICKOFF_SIZE,
                            opt_cart=False,
                            **OPT_KWARGS)
                        print('removing saddlepoint hessian')

                        run_fs = autofile.fs.run(thy_run_path)
                        run_fs.leaf.remove([elstruct.Job.HESSIAN])


                moldr.driver.save_initial_geometry(
                    species_info=species_info[name],
                    theory_level=RUN_OPT_LEVELS[opt_level_idx],
                    run_prefix=spc_run_path,
                    save_prefix=spc_save_path,
                )

            if RUN_CONF_SAMP:
                moldr.driver.conformer_sampling(
                    species_info=species_info[name],
                    theory_level=RUN_OPT_LEVELS[opt_level_idx],
                    run_prefix=spc_run_path,
                    save_prefix=spc_save_path,
                    script_str=OPT_SCRIPT_STR,
                    overwrite=OVERWRITE,
                    nsamp_par=NSAMP_CONF_PAR,
                    **OPT_KWARGS,
                )

                if RUN_MIN_GRAD:
                    moldr.driver.run_minimum_energy_gradient(
                        species_info=species_info[name],
                        theory_level=RUN_OPT_LEVELS[opt_level_idx],
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

                if RUN_MIN_HESS:
                    moldr.driver.run_minimum_energy_hessian(
                        species_info=species_info[name],
                        theory_level=RUN_OPT_LEVELS[opt_level_idx],
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

                if RUN_MIN_VPT2:
                    moldr.driver.run_minimum_energy_vpt2(
                        species_info=species_info[name],
                        theory_level=RUN_OPT_LEVELS[opt_level_idx],
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

                if RUN_CONF_SCAN:
                    moldr.driver.hindered_rotor_scans(
                        species_info=species_info[name],
                        theory_level=RUN_OPT_LEVELS[opt_level_idx],
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=OVERWRITE,
                        scan_increment=SCAN_INCREMENT,
                        **OPT_KWARGS,
                    )

                if RUN_CONF_GRAD:
                    moldr.driver.run_conformer_gradients(
                        species_info=species_info[name],
                        theory_level=RUN_OPT_LEVELS[opt_level_idx],
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

                if RUN_CONF_HESS:
                    moldr.driver.run_conformer_hessians(
                        species_info=species_info[name],
                        theory_level=RUN_OPT_LEVELS[opt_level_idx],
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

            if RUN_TAU_SAMP:
                moldr.driver.tau_sampling(
                    species_info=species_info[name],
                    theory_level=RUN_OPT_LEVELS[opt_level_idx],
                    run_prefix=spc_run_path,
                    save_prefix=spc_save_path,
                    script_str=OPT_SCRIPT_STR,
                    overwrite=OVERWRITE,
                    nsamp_par=NSAMP_TAU_PAR,
                    **OPT_KWARGS,
                )

                if RUN_TAU_GRAD:
                    moldr.driver.run_tau_gradients(
                        species_info=species_info[name],
                        theory_level=RUN_OPT_LEVELS[opt_level_idx],
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )
 
                if RUN_TAU_HESS:
                    moldr.driver.run_tau_hessians(
                        species_info=species_info[name],
                        theory_level=RUN_OPT_LEVELS[opt_level_idx],
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=OPT_SCRIPT_STR,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

            if TAU_PF_WRITE:
                moldr.driver.tau_pf_write(
                    name=name,
                    save_prefix=thy_save_path,
                    run_grad=RUN_GRAD_PF,
                    run_hess=RUN_HESS_PF,
                )

        for high_level_idx, _ in enumerate(RUN_HIGH_LEVELS):

            orb_restr = moldr.util.orbital_restriction(
                    species_info[name], THEORY_REF_HIGH_LEVEL)
            ref_level = THEORY_REF_HIGH_LEVEL[1:3]
            ref_level.append(orb_restr)
            print('ref level test:', ref_level)

            spc_run_fs = autofile.fs.species(RUN_PREFIX)
            spc_run_fs.leaf.create(species_info[name])
            spc_run_path = spc_run_fs.leaf.path(species_info[name])
            spc_save_fs = autofile.fs.species(SAVE_PREFIX)
            spc_save_fs.leaf.create(species_info[name])
            spc_save_path = spc_save_fs.leaf.path(species_info[name])

            ref_run_fs = autofile.fs.theory(spc_run_path)
            ref_run_fs.leaf.create(ref_level)
            ref_run_path = ref_run_fs.leaf.path(ref_level)
            ref_save_fs = autofile.fs.theory(spc_save_path)
            ref_save_fs.leaf.create(ref_level)
            ref_save_path = ref_save_fs.leaf.path(ref_level)

            min_cnf_locs = moldr.util.min_energy_conformer_locators(
                ref_save_path)
            cnf_run_fs = autofile.fs.conformer(ref_run_path)
            cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
            cnf_save_fs = autofile.fs.conformer(ref_save_path)
            cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
            min_cnf_geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)

            # evaluate the high level energy and save it

            prog = RUN_HIGH_LEVELS[high_level_idx][0]
            SP_SCRIPT_STR, OPT_SCRIPT_STR, KWARGS, OPT_KWARGS = (
                moldr.util.run_qchem_par(prog))
            if RUN_HL_MIN_ENE:
                moldr.driver.run_single_point_energy(
                    geo=min_cnf_geo,
                    species_info=species_info[name],
                    theory_level=RUN_HIGH_LEVELS[high_level_idx],
                    run_prefix=cnf_run_path,
                    save_prefix=cnf_save_path,
                    script_str=SP_SCRIPT_STR,
                    overwrite=OVERWRITE,
                    **KWARGS,
                )

if RUN_SPECIES_THERMO:
    # get the high level basis energies 
    ene_ref = []

    # generate and store the energies so that the reference energies can be used as needed
    ene_hl = {}
    print('Evaluating thermo for following species')
    print(SPC_NAMES)
    for name in SPC_NAMES:
        # set up species information
        smi = SMI_DCT[name]
        ich = automol.smiles.inchi(smi)
        print("smiles: {}".format(smi), "inchi: {}".format(ich))

        # get and store the high level energies
        print('name in high level energy routine')
        print(name)
        for high_level_idx, _ in enumerate(RUN_HIGH_LEVELS):
            min_ene = moldr.driver.get_high_level_energy(
                species_info=species_info[name],
                theory_low_level=THEORY_REF_HIGH_LEVEL,
                theory_high_level=RUN_HIGH_LEVELS[high_level_idx],
                save_prefix=SAVE_PREFIX,)
            ene_hl[(name, high_level_idx)] = min_ene

    spc_save_fs = autofile.fs.species(SAVE_PREFIX)
    for tors_model, vib_model in SPECIES_MODELS:
    # rot_model = RRHO or 1DHR or MDHR or TAU
        print('tors_model:',tors_model)
        print('vib_model:',vib_model)
        for har_level, tors_level, vpt2_level in PF_LEVELS:
            # get the zero-point energy for each species
            print('harmonic_level:',har_level)
            print('tors_level:',tors_level)
            print('vpt2_level:',vpt2_level)
            print('Calculating zpe')
            species_zpe = {}
            for name in SPC_NAMES:
                smi = SMI_DCT[name]
                ich = automol.smiles.inchi(smi)
                print("smiles: {}".format(smi), "inchi: {}".format(ich))
                spc_save_fs.leaf.create(species_info[name])
                spc_save_path = spc_save_fs.leaf.path(species_info[name])

                species_zpe[name] = moldr.driver.get_zero_point_energy(
                    species_info[name],
                    tors_model, vib_model,
                    har_level, tors_level, vpt2_level,
                    script_str=PF_SCRIPT_STR,
                    elec_levels=[[0., 1]], sym_factor=1.,
                    save_prefix=spc_save_path)
                print(name, species_zpe[name])

            pf_inp_str = {}
            print('finished zpe')
            # get the partition function for each species
            
            chemkin_poly_strs = ['' for i in range(len(RUN_HIGH_LEVELS))]
            for name in SPC_NAMES:
                # set up species information
                smi = SMI_DCT[name]
                ich = automol.smiles.inchi(smi)
                print("smiles: {}".format(smi), "inchi: {}".format(ich))
                spc_save_fs.leaf.create(species_info[name])
                spc_save_path = spc_save_fs.leaf.path(species_info[name])

                # to be generalized
                sym_factor = 1.
                elec_levels = [[0., mul]]
                if (ich, mul) in ELC_DEG_DCT:
                    elec_levels = ELC_DEG_DCT[(ich, mul)]

                # cycle through the low levels generating partition functions  for each
                species_str = moldr.driver.species_block(
                    species_info=species_info[name],
                    tors_model=tors_model,
                    vib_model=vib_model,
                    har_level=har_level,
                    tors_level=tors_level,
                    vpt2_level=vpt2_level,
                    script_str=PROJROT_SCRIPT_STR,
                    elec_levels=elec_levels,
                    sym_factor=sym_factor,
                    save_prefix=spc_save_path,
                    )

                # create a messpf input file
                temp_step = TEMP_STEP
                ntemps = NTEMPS
                global_pf_str = mess_io.writer.write_global_pf(
                    [], temp_step, ntemps, rel_temp_inc=0.001,
                    atom_dist_min=0.6)
                print(global_pf_str)
                species_head_str = 'Species ' + name
                print(species_head_str)
                pf_inp_str[name] = '\n'.join(
                    [global_pf_str, species_head_str,
                     species_str])
                print(species_str)

                orb_restr = moldr.util.orbital_restriction(
                        species_info[name], tors_level)
                tors_levelp = tors_level[1:3]
                tors_levelp.append(orb_restr)

                thy_save_fs = autofile.fs.theory(spc_save_path)
                thy_save_fs.leaf.create(tors_levelp)
                thy_save_path = thy_save_fs.leaf.path(tors_levelp)
                bld_locs = ['PF', 0]
                bld_save_fs = autofile.fs.build(thy_save_path)
                bld_save_fs.leaf.create(bld_locs)
                pf_path = bld_save_fs.leaf.path(bld_locs)
                print('Build Path for Partition Functions')
                print(pf_path)

                # run messpf
                with open(os.path.join(pf_path, 'pf.inp'), 'w') as pf_file:
                    pf_file.write(pf_inp_str[name])
                moldr.util.run_script(PF_SCRIPT_STR, pf_path)
                print('finished partition function')

                formula = thermo.util.inchi_formula(ich)
                print('\nformula:')
                print(formula)

                # Get atom count dictionary
                atom_dict = thermo.util.get_atom_counts_dict(formula)
                print('\natom dict:')
                print(atom_dict)

                # Get the list of the basis
                spc_bas = thermo.heatform.get_reduced_basis(SPC_REF_ICH, formula)
                print('\nbasis:')
                print(spc_bas)

                # Get the coefficients for the balanced heat-of-formation eqn
                coeff = thermo.heatform.calc_coefficients(spc_bas, atom_dict)
                print('\ncoeff:')
                print(coeff)

                # prepare NASA polynomials
                nasa_inp_str = ('nasa')
                bld_locs = ['NASA_POLY', 0]
                bld_save_fs.leaf.create(bld_locs)
                nasa_path = bld_save_fs.leaf.path(bld_locs)
                print('NASA build path')
                print(nasa_path)

                h_basis = []
                for hl_idx, _ in enumerate(RUN_HIGH_LEVELS):
                    for  name_ref in SPC_REF_NAMES:
                        smi = SMI_DCT[name_ref]
                        ich = automol.convert.smiles.inchi(smi)
                        if ich in spc_bas:
                            tmp = ene_hl[(name_ref, hl_idx)] + species_zpe[name_ref]/EH2KCAL
                            h_basis.append(tmp)
                    print('\ne_basis:')
                    print(h_basis)

                    # Get the 0 K heat of formation
                    species_ene = ene_hl[(name, high_level_idx)] + species_zpe[name]/EH2KCAL
                    h0form = thermo.heatform.calc_hform_0k(species_ene, h_basis, spc_bas, coeff, ref_set='ATcT')
                    print('h0form = ',h0form)

                    # need to change back to starting directory after running thermp and pac99 or rest of code is confused
                    starting_path = os.getcwd()
                    os.chdir(nasa_path)

                    # Write thermp input file
                    ENTHALPYT = 0.
                    BREAKT = 1000.
                    thermo.runner.write_thermp_input(
                        formula=formula,
                        deltaH=h0form,
                        enthalpyT=ENTHALPYT,
                        breakT=BREAKT,
                        thermp_file_name='thermp.dat')

                    # Run thermp
                    thermo.runner.run_thermp(
                        pf_path=pf_path,
                        thermp_path=nasa_path,
                        thermp_file_name='thermp.dat',
                        pf_file_name='pf.dat'
                        )

                    # Run pac99
                    print('formula test')
                    print(formula)
                    print(nasa_path)
                    thermo.runner.run_pac99(nasa_path, formula)

                    with open(os.path.join(nasa_path, 'thermp.out'), 'r') as thermp_outfile:
                        thermp_out_str = thermp_outfile.read()

                    with open(os.path.join(nasa_path, formula+'.c97'), 'r') as pac99_file:
                        pac99_str = pac99_file.read()

                    # Get the pac99 polynomial
                    pac99_poly_str = thermo.nasapoly.get_pac99_polynomial(pac99_str)
                    print('\nPAC99 Polynomial:')
                    print(pac99_poly_str)

                    # Convert the pac99 polynomial to chemkin polynomial
                    comment_str = '! tors model: {0}\n'.format(tors_model)
                    comment_str += '! vib model: {0}\n'.format(vib_model)
                    comment_str += '! har level: {0}\n'.format(har_level)
                    comment_str += '! tors level: {0}\n'.format(tors_level)
                    comment_str += '! vpt2 level: {0}\n'.format(vpt2_level)
                    comment_str += '! hl_method level: {0}\n'.format(RUN_HIGH_LEVELS[hl_idx][1])
                    comment_str += '! hl_basis level: {0}\n'.format(RUN_HIGH_LEVELS[hl_idx][2])

                    chemkin_poly_str = thermo.nasapoly.convert_pac_to_chemkin(name, atom_dict, comment_str, pac99_poly_str)
                    print('\nCHEMKIN Polynomial:')
                    print(chemkin_poly_str)
                    print(hl_idx)
                    chemkin_poly_strs[hl_idx] += chemkin_poly_str
                    
                    print('starting_path in thermo')
                    print(starting_path)
                    os.chdir(starting_path) 
                    with open(os.path.join(nasa_path, formula+'.ckin'), 'w') as nasa_file:
                        nasa_file.write(chemkin_poly_str)
                    with open(os.path.join(starting_path, formula+'.ckin'), 'w') as nasa_file:
                        nasa_file.write(chemkin_poly_str)

            for hl_idx, _ in enumerate(RUN_HIGH_LEVELS):
                with open(os.path.join(starting_path, formula+'.ckin'), 'w') as nasa_file:
                    nasa_file.write(chemkin_poly_strs[hl_idx])


# 5. process reaction data from the mechanism file
RXN_BLOCK_STR = chemkin_io.reaction_block(MECH_STR)
RXN_STRS = chemkin_io.reaction.data_strings(RXN_BLOCK_STR)
RCT_NAMES_LST = list(
    map(chemkin_io.reaction.DataString.reactant_names, RXN_STRS))
PRD_NAMES_LST = list(
    map(chemkin_io.reaction.DataString.product_names, RXN_STRS))

if RUN_REACTIONS_QCHEM:
    for rct_names, prd_names in zip(RCT_NAMES_LST, PRD_NAMES_LST):
        # print the CHEMKIN reaction name for reference
        rxn_name = '='.join(['+'.join(rct_names), '+'.join(prd_names)])
        print()
        print("Reaction: {}".format(rxn_name))

        # determine inchis, charges, and multiplicities

        rct_smis = list(map(SMI_DCT.__getitem__, rct_names))
        prd_smis = list(map(SMI_DCT.__getitem__, prd_names))
        rct_ichs = list(map(automol.smiles.inchi, rct_smis))
        prd_ichs = list(map(automol.smiles.inchi, prd_smis))
        rct_chgs = list(map(CHG_DCT.__getitem__, rct_names))
        prd_chgs = list(map(CHG_DCT.__getitem__, prd_names))
        rct_muls = list(map(MUL_DCT.__getitem__, rct_names))
        prd_muls = list(map(MUL_DCT.__getitem__, prd_names))

        # determine the transition state multiplicity
        ts_mul = automol.mult.ts.low(rct_muls, prd_muls)
        ts_chg = sum(rct_chgs)
        print('ts_chg test:',ts_chg)
        ts_info = ('', ts_chg, ts_mul)

        # theory
        for opt_level_idx, _ in enumerate(RUN_OPT_LEVELS):
#        for prog, method, basis in RUN_OPT_LEVELS:
            prog = RUN_OPT_LEVELS[opt_level_idx][0]
            SCRIPT_STR, OPT_SCRIPT_STR, KWARGS, OPT_KWARGS = (
                moldr.util.run_qchem_par(prog))

            ts_orb_restr = moldr.util.orbital_restriction(
                ts_info, RUN_OPT_LEVELS[opt_level_idx])
            thy_level = RUN_OPT_LEVELS[opt_level_idx][1:3]
            thy_level.append(orb_restr)

#            ts_orb_restr = moldr.util.orbital_restriction(
#                ts_mul, RESTRICT_OPEN_SHELL)

            # check direction of reaction
            rxn_ichs = [rct_ichs, prd_ichs]
            rxn_chgs = [rct_chgs, prd_chgs]
            rxn_muls = [rct_muls, prd_muls]
#            rxn_info = zip(rxn_ichs, rxn_chgs, rxn_muls)
#            print('rxn_info test:', rxn_info)
            rxn_exo = moldr.util.reaction_energy(
                    SAVE_PREFIX, rxn_ichs, rxn_chgs, rxn_muls, RUN_OPT_LEVELS[opt_level_idx])
#                    SAVE_PREFIX, rxn_ichs, rxn_chgs, rxn_muls, thy_level)
#                SAVE_PREFIX, rxn_ichs, rxn_chgs, rxn_muls, method, basis,
#                RESTRICT_OPEN_SHELL)
            print(rxn_exo)
            if rxn_exo > 0:
                rct_ichs, prd_ichs = prd_ichs, rct_ichs
                rct_chgs, prd_chgs = prd_chgs, rct_chgs
                rct_muls, prd_muls = prd_muls, rct_muls
                print('ts search will be performed in reverse direction')

            # obtain geometries from a hierachy of (i) data directory and (ii)
            # previous species calculation
            rct_geos = []
            for ich, chg, mul in zip(rct_ichs, rct_chgs, rct_muls):
                rct_info = [ich, chg, mul]
#                orb_restr = moldr.util.orbital_restriction(
#                    rct_info, RUN_OPT_LEVELS[opt_level_idx])
#                rct_level = RUN_OPT_LEVELS[opt_level_idx][1:3]
#                rct_level.append(orb_restr)
#                print('geo test in rxn_qchem:', rct_info, rct_level, SAVE_PREFIX, GEOM_DCT)
                geo = moldr.util.reference_geometry(
                    rct_info, RUN_OPT_LEVELS[opt_level_idx], SAVE_PREFIX,
                    GEOM_DCT)
                rct_geos.append(geo)

            prd_geos = []
            for ich, chg, mul in zip(prd_ichs, prd_chgs, prd_muls):
                prd_info = [ich, chg, mul]
#                orb_restr = moldr.util.orbital_restriction(
#                    prd_info, RUN_OPT_LEVELS[opt_level_idx])
#                prd_level = RUN_OPT_LEVELS[opt_level_idx][1:3]
#                prd_level.append(orb_restr)
                geo = moldr.util.reference_geometry(
                    prd_info, RUN_OPT_LEVELS[opt_level_idx], SAVE_PREFIX,
                    GEOM_DCT)
                prd_geos.append(geo)

            # determine the transition state z-matrix
            # replace this with save values if they are available
            rct_zmas = list(map(automol.geom.zmatrix, rct_geos))
            prd_zmas = list(map(automol.geom.zmatrix, prd_geos))

            typ = None

            # # (migrations are not yet implemented)
            # ret = automol.zmatrix.ts.hydrogen_migration(rct_zmas, prd_zmas)
            # if ret and typ is None:
            #     typ = 'hydrogen migration'

            ret = automol.zmatrix.ts.beta_scission(rct_zmas, prd_zmas)
            if ret and typ is None:
                typ = 'beta scission'
                ts_zma, dist_name, tors_names = ret

            ret = automol.zmatrix.ts.addition(rct_zmas, prd_zmas)
            if ret and typ is None:
                typ = 'addition'
                ts_zma, dist_name, tors_names = ret

            # fix this later
            # ret = automol.zmatrix.ts.hydrogen_abstraction(rct_zmas, prd_zmas,
            #                                               sigma=True)
            ret = automol.zmatrix.ts.hydrogen_abstraction(rct_zmas, prd_zmas,
                                                          sigma=False)
            if ret and typ is None:
                typ = 'hydrogen abstraction'
                ts_zma, dist_name, tors_names = ret

            if typ is None:
                print("Failed to classify reaction.")
            else:
                print("Type: {}".format(typ))

                # determine the grid
                dist_coo, = automol.zmatrix.coordinates(ts_zma)[dist_name]
                syms = automol.zmatrix.symbols(ts_zma)
                bnd_len_key = tuple(sorted(map(syms.__getitem__, dist_coo)))

                bnd_len_dct = {
                    ('C', 'C'): 1.54 * ANG2BOHR,
                    ('C', 'H'): 1.09 * ANG2BOHR,
                    ('H', 'H'): 0.74 * ANG2BOHR,
                    ('N', 'N'): 1.45 * ANG2BOHR,
                    ('O', 'O'): 1.48 * ANG2BOHR,
                    ('C', 'N'): 1.47 * ANG2BOHR,
                    ('C', 'O'): 1.43 * ANG2BOHR,
                    ('H', 'O'): 1.20 * ANG2BOHR,
                    ('H', 'N'): 0.99 * ANG2BOHR,
                }

                if typ in ('beta scission', 'addition'):
                    rmin = 1.4 * ANG2BOHR
                    rmin = 2.8 * ANG2BOHR
                    if bnd_len_key in bnd_len_dct:
                        bnd_len = bnd_len_dct[bnd_len_key]
                        rmin = bnd_len + 0.2 * ANG2BOHR
                        rmax = bnd_len + 1.6 * ANG2BOHR
                elif typ == 'hydrogen abstraction':
                    rmin = 0.7 * ANG2BOHR
                    rmax = 2.2 * ANG2BOHR
                    if bnd_len_key in bnd_len_dct:
                        bnd_len = bnd_len_dct[bnd_len_key]
                        rmin = bnd_len
                        rmax = bnd_len + 1.0 * ANG2BOHR

                npoints = 8
                grid = numpy.linspace(rmin, rmax, npoints)

                # construct the filesystem
                rxn_ichs = [rct_ichs, prd_ichs]
                rxn_chgs = [rct_chgs, prd_chgs]
                rxn_muls = [rct_muls, prd_muls]

                # set up the filesystem
                is_rev = autofile.system.reaction_is_reversed(
                    rxn_ichs, rxn_chgs, rxn_muls)
                rxn_ichs, rxn_chgs, rxn_muls = autofile.system.sort_together(
                    rxn_ichs, rxn_chgs, rxn_muls)
                print(" - The reaction direction is {}"
                      .format('backward' if is_rev else 'forward'))

                ts_chg = 0
                for rct_chg in rct_chgs:
                    ts_chg += rct_chg
                ts_info = ['', ts_chg, ts_mul]

                rxn_run_fs = autofile.fs.reaction(RUN_PREFIX)
                rxn_run_fs.leaf.create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
                rxn_run_path = rxn_run_fs.leaf.path(
                    [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

                rxn_ichs = tuple(map(tuple, rxn_ichs))
                rxn_chgs = tuple(map(tuple, rxn_chgs))
                rxn_muls = tuple(map(tuple, rxn_muls))
                print('rxn_save test0', rxn_ichs, rxn_chgs, rxn_muls, ts_mul)
                print(SAVE_PREFIX)
                rxn_save_fs = autofile.fs.reaction(SAVE_PREFIX)
                rxn_save_fs.leaf.create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
                rxn_save_path = rxn_save_fs.leaf.path(
                    [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

                orb_restr = moldr.util.orbital_restriction(
                    ts_info, RUN_OPT_LEVELS[opt_level_idx])
                ref_level = RUN_OPT_LEVELS[opt_level_idx][1:3]
                ref_level.append(orb_restr)
                print('ref level test:', ref_level)

                thy_run_fs = autofile.fs.theory(rxn_run_path)
                thy_run_fs.leaf.create(ref_level)
                thy_run_path = thy_run_fs.leaf.path(
                    ref_level)

                thy_save_fs = autofile.fs.theory(rxn_save_path)
                thy_save_fs.leaf.create(ref_level)
                thy_save_path = thy_save_fs.leaf.path(ref_level)

                print('entering run_scan:')

                moldr.driver.run_scan(
                    zma=ts_zma,
                    species_info=ts_info,
                    theory_level=RUN_OPT_LEVELS[opt_level_idx],
                    grid_dct={dist_name: grid},
                    run_prefix=thy_run_path,
                    save_prefix=thy_save_path,
                    script_str=SCRIPT_STR,
                    overwrite=OVERWRITE,
                    update_guess=False,
                    reverse_sweep=False,
                    **OPT_KWARGS
                )

                moldr.driver.save_scan(
                    run_prefix=thy_run_path,
                    save_prefix=thy_save_path,
                    coo_names=[dist_name],
                )

                scn_save_fs = autofile.fs.scan(thy_save_path)
                locs_lst = [
                    locs for locs in scn_save_fs.leaf.existing([[dist_name]])
                    if scn_save_fs.leaf.file.energy.exists(locs)]
                print(locs_lst)
                enes = [scn_save_fs.leaf.file.energy.read(locs)
                        for locs in locs_lst]
                max_locs = locs_lst[enes.index(max(enes))]
                max_ene = max(enes)
                max_zma = scn_save_fs.leaf.file.zmatrix.read(max_locs)

                print('optimizing ts')
                # find saddlepoint from maximum on the grid opt scan
                print('thy_run_path in ts_opt:', thy_run_path)
                ts_run_fs = autofile.fs.ts(thy_run_path)
                ts_run_fs.trunk.create()
                ts_run_path = ts_run_fs.trunk.path()
                print('ts_run_path:', ts_run_path)
#                ts_run_path = ts_run_fs.trunk.path(ref_level)
#                print('ts_run_path:', ts_run_path)

                ts_save_fs = autofile.fs.ts(thy_save_path)
                ts_save_fs.trunk.create()
                ts_save_path = ts_save_fs.trunk.path()
                print('ts_save_path:', ts_save_path)
#                ts_save_fs.leaf.create()

#                cnf_run_fs = autofile.fs.conformer(ts_run_path)
#                cnf_run_fs.trunk.create()
#                cnf_save_fs = autofile.fs.conformer(ts_save_path)
#                cnf_save_fs.trunk.create()

#                cnf_locs_lst = cnf_save_fs.leaf.existing()
#                print('cnf_locs_lst:', cnf_locs_lst)
#                for locs in cnf_locs_lst:
#                    cnf_run_path = cnf_run_fs.leaf.path(locs)
#                    cnf_save_path = cnf_run_fs.leaf.path(locs)
#                    print(locs)

#                print('cnf_run_path:', cnf_run_path)
#                print('cnf_save_path:', cnf_save_path)
                moldr.driver.run_job(
                    job='optimization',
                    script_str=SCRIPT_STR,
                    prefix=ts_run_path,
                    geom=max_zma,
                    species_info=ts_info,
                    theory_level=RUN_OPT_LEVELS[opt_level_idx],
                    saddle=True,
                    overwrite=OVERWRITE,
                    **OPT_KWARGS,
                )
                opt_ret = moldr.driver.read_job(
                    job='optimization',
                    prefix=ts_run_path,
                )
                if opt_ret is not None:
                    inf_obj, inp_str, out_str = opt_ret
                    prog = inf_obj.prog
                    method = inf_obj.method
                    ene = elstruct.reader.energy(prog, method, out_str)
                    geo = elstruct.reader.opt_geometry(prog, out_str)
                    zma = elstruct.reader.opt_zmatrix(prog, out_str)

                    print(" - Saving...")
                    print(" - Save path: {}".format(ts_save_path))

#                   thy_save_fs.leaf.file.geometry_info.write(inf_obj)
#                   thy_save_fs.leaf.file.geometry_input.write(inp_str)
                    ts_save_fs.trunk.file.energy.write(ene)
                    ts_save_fs.trunk.file.geometry.write(geo)
                    ts_save_fs.trunk.file.zmatrix.write(zma)

#                    moldr.driver.run_job(
#                        job='hessian',
#                        script_str=SCRIPT_STR,
#                        prefix=thy_run_path,
#                        geom=geo,
#                        species_info=ts_info,
#                        theory_level=RUN_OPT_LEVELS[opt_level_idx],
#                        overwrite=OVERWRITE,
#                        **KWARGS,
#                    )
#                    hess_ret = moldr.driver.read_job(
#                        job='hessian',
#                        prefix=thy_run_path,
#                    )
#                    if hess_ret is not None:
#                        inf_obj, inp_str, out_str = hess_ret
#                        prog = inf_obj.prog
#                        method = inf_obj.method
#                        hess = elstruct.reader.hessian(prog, out_str)
#                        freqs = elstruct.util.harmonic_frequencies(geo, hess)
#
#                        print(" - Saving hessian...")
#                        print(" - Save path: {}".format(thy_save_path))
#
#                        thy_save_fs.leaf.file.hessian_info.write(inf_obj)
#                        thy_save_fs.leaf.file.hessian_input.write(inp_str)
#                        thy_save_fs.leaf.file.hessian.write(hess)
#                        thy_save_fs.leaf.file.harmonic_frequencies.write(freqs)

            if RUN_TS_CONF_SAMP:
                moldr.driver.conformer_sampling(
                    species_info=ts_info,
                    theory_level=RUN_OPT_LEVELS[opt_level_idx],
                    run_prefix=ts_run_path,
                    save_prefix=ts_save_path,
                    script_str=SCRIPT_STR,
                    overwrite=OVERWRITE,
                    saddle=True,
                    nsamp_par=NSAMP_CONF_PAR,
                    **OPT_KWARGS,
                )

                if RUN_TS_MIN_GRAD:
                    moldr.driver.run_minimum_energy_gradient(
                        species_info=ts_info,
                        theory_level=RUN_OPT_LEVELS[opt_level_idx],
                        run_prefix=ts_run_path,
                        save_prefix=ts_save_path,
                        script_str=SCRIPT_STR,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

                if RUN_TS_MIN_HESS:
                    moldr.driver.run_minimum_energy_hessian(
                        species_info=ts_info,
                        theory_level=RUN_OPT_LEVELS[opt_level_idx],
                        run_prefix=ts_run_path,
                        save_prefix=ts_save_path,
                        script_str=SCRIPT_STR,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

                if RUN_TS_MIN_VPT2:
                    moldr.driver.run_minimum_energy_vpt2(
                        species_info=ts_info,
                        theory_level=RUN_OPT_LEVELS[opt_level_idx],
                        run_prefix=ts_run_path,
                        save_prefix=ts_save_path,
                        script_str=SCRIPT_STR,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

                if RUN_TS_TAU_SAMP:

                    moldr.driver.save_tau(
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                    )

                    zma = ts_save_fs.trunk.file.zmatrix.read()
                    tors_ranges = automol.zmatrix.torsional_sampling_ranges(
                        zma, tors_names)
                    tors_range_dct = dict(zip(tors_names, tors_ranges))

                    moldr.driver.run_tau(
                        zma=zma,
                        species_info=ts_info,
                        theory_level=RUN_OPT_LEVELS[opt_level_idx],
                        nsamp=nsamp,
                        tors_range_dct=tors_range_dct,
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=SCRIPT_STR,
                        overwrite=OVERWRITE,
                        **OPT_KWARGS,
                    )
#                        saddle=True,
# used to have saddle=True in call, but this is not used. Probably a bug.

                    moldr.driver.save_tau(
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                    )

                if RUN_TS_CONF_SCAN:
                    zma = ts_save_fs.trunk.file.zmatrix.read(thy_save_path)
                    val_dct = automol.zmatrix.values(zma)
                    tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
                        zma, tors_names, SCAN_INCREMENT)
                    tors_grids = [
                        numpy.linspace(*linspace) + val_dct[name]
                        for name, linspace in zip(tors_name, tors_linspaces)]
                    for tors_name, tors_grid in zip(tors_names, tors_grids):
                        moldr.driver.run_scan(
                            zma=zma,
                            species_info=ts_info,
                            theory_level=RUN_OPT_LEVELS[opt_level_idx],
                            grid_dct={tors_name: tors_grid},
                            run_prefix=ts_run_path,
                            save_prefix=ts_save_path,
                            script_str=SCRIPT_STR,
                            overwrite=OVERWRITE,
                            **OPT_KWARGS,
                        )
#                            saddle=True,
# ??

                        moldr.driver.save_scan(
                            run_prefix=ts_run_path,
                            save_prefix=ts_save_path,
                            coo_names=[tors_name],
                        )
# AVC: debugged to here
                    hind_rot_dct = {}
                    scn_run_fs = autofile.fs.scan(ts_run_path)
                    scn_save_fs = autofile.fs.scan(ts_save_path)

                    min_ene = ts_save_fs.leaf.file.energy.read()
                    for tors_name in tors_names:
                        enes = [scn_save_fs.leaf.file.energy.read(locs)
                                for locs in scn_save_fs.leaf.existing()]
                        enes = numpy.subtract(enes, min_ene)
                        hind_rot_dct[tors_name] = enes*EH2KCAL

                    print('ts hindered rotor potential')
                    print(hind_rot_dct)

                if RUN_TS_KICKS_QCHEM:
                    ret = moldr.driver.read_job(
                        job=elstruct.Job.HESSIAN, prefix=ts_run_path)
                    if ret:
                        inf_obj, _, out_str = ret
                        prog = inf_obj.prog
                        hess = elstruct.reader.hessian(prog, out_str)
                        freqs = elstruct.util.harmonic_frequencies(geo, hess, project=True)
                        norm_coos = elstruct.util.normal_coordinates(geo, hess, project=True)
                        assert freqs[0] < -100

                        print('Kicking off from saddle in forward direction')
                        im_norm_coo = numpy.array(norm_coos)[:, 0]
                        disp_xyzs = numpy.reshape(im_norm_coo, (-1, 3))
                        dir_afs = autofile.fs.direction()
                        fwd_run_path = dir_afs.direction.dir.path(thy_run_path, [True])
                        dir_afs.direction.dir.create(thy_run_path, [True])
                        fwd_save_path = dir_afs.direction.dir.path(thy_save_path, [True])
                        dir_afs.direction.dir.create(thy_save_path, [True])
                        print(automol.geom.string(geo))
                        print(disp_xyzs)
                        moldr.driver.run_kickoff_saddle(
                            geo, disp_xyzs, chg, mul, method, basis, orb_restr,
                            fwd_run_path, SCRIPT_STR, prog, OVERWRITE,
                            kickoff_size=KICKOFF_SIZE, kickoff_backward=False,
                            opt_cart=True, **OPT_KWARGS)
                        print('Saving product of kick off from saddle in forward direction')
                        ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, prefix=fwd_run_path)
                        if ret:
                            inf_obj, inp_str, out_str = ret
                            prog = inf_obj.prog
                            method = inf_obj.method
                            ene = elstruct.reader.energy(prog, method, out_str)
                            geo = elstruct.reader.opt_geometry(prog, out_str)
                            fwd_save_path = dir_afs.direction.dir.path(thy_save_path, [True])
                            print('save path', fwd_save_path)
                            dir_afs.direction.file.geometry_info.write(inf_obj, thy_save_path, [True])
                            dir_afs.direction.file.geometry_input.write(inp_str, thy_save_path, [True])
                            dir_afs.direction.file.geometry.write(geo, thy_save_path, [True])
                            dir_afs.direction.file.energy.write(ene, thy_save_path, [True])

                        print('Kicking off from saddle in backward direction')
                        bwd_run_path = dir_afs.direction.dir.path(thy_run_path, [False])
                        dir_afs.direction.dir.create(thy_run_path, [False])
                        bwd_save_path = dir_afs.direction.dir.path(thy_save_path, [False])
                        dir_afs.direction.dir.create(thy_save_path, [False])
                        moldr.driver.run_kickoff_saddle(
                            geo, disp_xyzs, chg, mul, method, basis,
                            orb_restr, bwd_run_path, SCRIPT_STR, prog,
                            OVERWRITE, kickoff_size=KICKOFF_SIZE,
                            kickoff_backward=True, **OPT_KWARGS)
                        print('Saving product of kick off from saddle in backward direction')
                        ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, prefix=bwd_run_path)
                        if ret:
                            inf_obj, inp_str, out_str = ret
                            prog = inf_obj.prog
                            method = inf_obj.method
                            ene = elstruct.reader.energy(prog, method, out_str)
                            geo = elstruct.reader.opt_geometry(prog, out_str)
                            bwd_save_path = dir_afs.direction.dir.path(thy_save_path, [False])
                            print('save path', bwd_save_path)
                            dir_afs.direction.file.geometry_info.write(inf_obj, thy_save_path, [False])
                            dir_afs.direction.file.geometry_input.write(inp_str, thy_save_path, [False])
                            dir_afs.direction.file.geometry.write(geo, thy_save_path, [False])
                            dir_afs.direction.file.energy.write(ene, thy_save_path, [False])

if RUN_VDW_QCHEM:
    if NSAMP_VDW_EXPR:
        nsamp = min(NSAMP_VDW_A + NSAMP_VDW_B * NSAMP_VDW_C**ntaudof, NSAMP_VDW_D)
    else:
        nsamp = NSAMP_VDW

    VDW_NAMES_LST = []
    for rct_names, prd_names in zip(RCT_NAMES_LST, PRD_NAMES_LST):
        rct_muls = list(map(MUL_DCT.__getitem__, rct_names))
        prd_muls = list(map(MUL_DCT.__getitem__, prd_names))
        ts_mul = automol.mult.ts.low(rct_muls, prd_muls)
        if len(rct_names) == 2:
            if sorted(rct_names) not in VDW_NAMES_LST:
                VDW_NAMES_LST.append([sorted(rct_names), ts_mul])
        if len(prd_names) == 2:
            if sorted(prd_names) not in VDW_NAMES_LST:
                VDW_NAMES_LST.append([sorted(prd_names), ts_mul])

    for names, ts_mul in VDW_NAMES_LST:
        smis = list(map(SMI_DCT.__getitem__, names))
        ichs = list(map(automol.smiles.inchi, smis))
        chgs = list(map(CHG_DCT.__getitem__, names))
        muls = list(map(MUL_DCT.__getitem__, names))

        method = METHOD
        basis = BASIS
        geos = []
        for ich, chg, mul in zip(ichs, chgs, muls):
            orb_restr = moldr.util.orbital_restriction(mul, RESTRICT_OPEN_SHELL)
            geo = moldr.util.reference_geometry(
                ich, chg, mul, method, basis, orb_restr, SAVE_PREFIX, GEOM_DCT)
            geos.append(geo)
           
        geo1, geo2 = geos
        geo1 = automol.geom.mass_centered(geo1)
        geo2 = automol.geom.mass_centered(geo2)
        for idx in range(nsamp):
            print('Optimizing vdw geometry {}/{}'.format(idx+1,nsamp))
            angs1 = numpy.multiply(
                numpy.random.rand(3), [1*numpy.pi, 2*numpy.pi, 2*numpy.pi])
            angs2 = numpy.multiply(
                numpy.random.rand(3), [1*numpy.pi, 2*numpy.pi, 2*numpy.pi])
            angs12 = numpy.multiply(
                numpy.random.rand(2), [1*numpy.pi, 2*numpy.pi])
            geo1 = automol.geom.euler_rotated(geo1, *angs1)
            geo2 = automol.geom.euler_rotated(geo2, *angs2)
            dist_cutoff = 3.*qcc.conversion_factor('angstrom', 'bohr')

            geo = automol.geom.join(geo1, geo2, dist_cutoff, *angs12)
            print("Species: {}".format('+'.join(names)))
            print('vdw starting geometry')
            print(automol.geom.xyz_string(geo))

    # set up the filesystem
            ich = automol.inchi.recalculate(automol.inchi.join(ichs))
            chg = sum(chgs)
            mul = ts_mul
            orb_restr = moldr.util.orbital_restriction(mul, RESTRICT_OPEN_SHELL)
            spc_run_path = moldr.util.species_path(ich, chg, mul, RUN_PREFIX)
            spc_save_path = moldr.util.species_path(ich, chg, mul, SAVE_PREFIX)
            thy_run_path = moldr.util.theory_path(method, basis, orb_restr, spc_run_path)
            thy_save_path = moldr.util.theory_path(method, basis, orb_restr, spc_save_path)

    # generate reference geometry
    # generate the z-matrix and sampling ranges

            moldr.driver.run_job(
                job=elstruct.Job.OPTIMIZATION,
                geom=geo,
                species_info=ts_info,
                theory_level=RUN_OPT_LEVELS[opt_level_idx],
                prefix=thy_run_path,
                script_str=SCRIPT_STR,
                overwrite=OVERWRITE,
                **OPT_KWARGS,
            )

    # save info for the initial geometry (from inchi or from save directory)
            ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, prefix=thy_run_path)
            if ret:
                print('Saving reference geometry')
                print(" - Save path: {}".format(thy_save_path))

                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                geo = elstruct.reader.opt_geometry(prog, out_str)
                thy_afs = autofile.fs.theory()
                thy_afs.theory.file.geometry.write(geo, spc_save_path, [method, basis, orb_restr])
                ene = elstruct.reader.energy(prog, method, out_str)
                print('ene test in vdw')
                print(ene)
                thy_afs.theory.file.energy.write(ene, spc_save_path, [method, basis, orb_restr])


if RUN_REACTIONS_RATES:
    for rct_names, prd_names in zip(RCT_NAMES_LST, PRD_NAMES_LST):
        # print the CHEMKIN reaction name for reference
        rxn_name = '='.join(['+'.join(rct_names), '+'.join(prd_names)])
        print()
        print('Mess Input for')
        print("Reaction: {}".format(rxn_name))

        # determine inchis, charges, and multiplicities

        rct_smis = list(map(SMI_DCT.__getitem__, rct_names))
        prd_smis = list(map(SMI_DCT.__getitem__, prd_names))
        rct_ichs = list(map(automol.smiles.inchi, rct_smis))
        prd_ichs = list(map(automol.smiles.inchi, prd_smis))
        rct_chgs = list(map(CHG_DCT.__getitem__, rct_names))
        prd_chgs = list(map(CHG_DCT.__getitem__, prd_names))
        rct_muls = list(map(MUL_DCT.__getitem__, rct_names))
        prd_muls = list(map(MUL_DCT.__getitem__, prd_names))
        rxn_name = '='.join(['+'.join(rct_names), '+'.join(prd_names)])

        ts_mul = automol.mult.ts.low(rct_muls, prd_muls)
        ts_chg = sum(rct_chgs)
        print('ts_chg test:',ts_chg)
        ts_info = ('', ts_chg, ts_mul)

# header section
        temperatures = TEMPS
        pressures = PRESS
        header_str = mess_io.writer.write_global_reaction(temperatures, pressures)
        print(header_str)

# energy transfer section
        exp_factor = EXP_FACTOR 
        exp_power = EXP_POWER 
        exp_cutoff = EXP_CUTOFF 
        eps1 = EPS1 
        eps2 = EPS2 
        sig1 = SIG1 
        sig2 = SIG2 
        mass1 = MASS1 
        mass2 = MASS2 
        energy_trans_str = mess_io.writer.write_energy_transfer(
            exp_factor, exp_power, exp_cutoff, eps1, eps2, sig1, sig2, mass1, mass2)
        print(energy_trans_str)

    for tors_model, vib_model in SPECIES_MODELS:
    # rot_model = RRHO or 1DHR or MDHR or TAU
        print('tors_model:',tors_model)
        print('vib_model:',vib_model)
        for har_level, tors_level, vpt2_level in PF_LEVELS:
            # get the zero-point energy for each species
            print('harmonic_level:',har_level)
            print('tors_level:',tors_level)
            print('vpt2_level:',vpt2_level)
        for opt_level_idx, _ in enumerate(RUN_OPT_LEVELS):
#        for prog, method, basis in RUN_OPT_LEVELS:
            prog = RUN_OPT_LEVELS[opt_level_idx][0]
            SCRIPT_STR, OPT_SCRIPT_STR, KWARGS, OPT_KWARGS = (
                moldr.util.run_qchem_par(prog))

            ts_orb_restr = moldr.util.orbital_restriction(
                ts_info, RUN_OPT_LEVELS[opt_level_idx])
            thy_level = RUN_OPT_LEVELS[opt_level_idx][1:3]
            thy_level.append(orb_restr)

            # check direction of reaction
            rxn_ichs = [rct_ichs, prd_ichs]
            rxn_chgs = [rct_chgs, prd_chgs]
            rxn_muls = [rct_muls, prd_muls]
            rxn_exo = moldr.util.reaction_energy(
                SAVE_PREFIX, rxn_ichs, rxn_chgs, rxn_muls, RUN_OPT_LEVELS[opt_level_idx])
            print(rxn_exo)
            if rxn_exo > 0:
                rct_ichs, prd_ichs = prd_ichs, rct_ichs
                rct_chgs, prd_chgs = prd_chgs, rct_chgs
                rct_muls, prd_muls = prd_muls, rct_muls
                print('ts search will be performed in reverse direction')

            # construct the filesystem
            rxn_ichs = [rct_ichs, prd_ichs]
            rxn_chgs = [rct_chgs, prd_chgs]
            rxn_muls = [rct_muls, prd_muls]

            # set up the filesystem
            is_rev = autofile.system.reaction_is_reversed(
                rxn_ichs, rxn_chgs, rxn_muls)
            rxn_ichs, rxn_chgs, rxn_muls = autofile.system.sort_together(
                rxn_ichs, rxn_chgs, rxn_muls)
            print(" - The reaction direction is {}"
                  .format('backward' if is_rev else 'forward'))

            ts_mul = automol.mult.ts.low(rct_muls, prd_muls)
            ts_chg = 0
            for rct_chg in rct_chgs:
                ts_chg += rct_chg
            ts_info = ['', ts_chg, ts_mul]

            rxn_ichs = tuple(map(tuple, rxn_ichs))
            rxn_chgs = tuple(map(tuple, rxn_chgs))
            rxn_muls = tuple(map(tuple, rxn_muls))
            print('rxn_save test', rxn_ichs, rxn_chgs, rxn_muls, ts_mul)
            print(SAVE_PREFIX)
            rxn_save_fs = autofile.fs.reaction(SAVE_PREFIX)
#            rxn_save_fs.leaf.create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
            rxn_save_path = rxn_save_fs.leaf.path(
                [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

            orb_restr = moldr.util.orbital_restriction(
                ts_info, RUN_OPT_LEVELS[opt_level_idx])
            ref_level = RUN_OPT_LEVELS[opt_level_idx][1:3]
            ref_level.append(orb_restr)
            print('ref level test:', ref_level)
            thy_save_fs = autofile.fs.theory(rxn_save_path)
            thy_save_fs.leaf.create(ref_level)
            thy_save_path = thy_save_fs.leaf.path(ref_level)

    # cycle over reactant and product species
    # check if unimolecular or bimolecular species
            print('Unimolecular or bimolecular test')
            indxw = 0
            indxp = 0
            well_str = ''
            bim_str = ''
            if not rct_ichs[1]:
                print(rct_ichs)
                indxw += 1
                species_info = (rct_ichs[0], rct_chgs[0], rct_muls[0])
                well_data = moldr.driver.species_block(
                    species_info=species_info,
                    tors_model=tors_model,
                    vib_model=vib_model,
                    har_level=har_level,
                    tors_level=tors_level,
                    vpt2_level=vpt2_level,
                    script_str=PROJROT_SCRIPT_STR,
                    elec_levels=elec_levels,
                    sym_factor=sym_factor,
                    save_prefix=spc_save_path,
                    )
                well_label = 'W'+str(indxw)
                reac_label = well_label
                zero_energy = 0.0
                well_str += mess_io.writer.write_well(well_label, well_data, zero_energy)
                well_data.replace()
    # write W_indxw 
            else:
                indxp += 1
                species_data = ['', '']
                species_label = ['', '']
                bimol_label = 'P'+str(indxp)
                reac_label = bimol_label
                for idx, (rct_ich, rct_chg, rct_mul) in enumerate(zip(rct_ichs, rct_chgs, rct_muls)):
                    species_label[idx] = automol.inchi.smiles(rct_ich)
                    species_info = (rct_ich, rct_chg, rct_mul)
                    species_data[idx] = moldr.driver.species_block(
                        species_info=species_info,
                        tors_model=tors_model,
                        vib_model=vib_model,
                        har_level=har_level,
                        tors_level=tors_level,
                        vpt2_level=vpt2_level,
                        script_str=PROJROT_SCRIPT_STR,
                        elec_levels=elec_levels,
                        sym_factor=sym_factor,
                        save_prefix=spc_save_path,
                        )
                ground_energy = 0.0
                bim_str += mess_io.writer.write_bimolecular(
                    bimol_label, species_label[0], species_data[0],
                    species_label[1], species_data[1], ground_energy)

    # write P_indxp
            if not prd_ichs[1]:
                print(prd_ichs)
                indxw += 1
                for prd_ich, prd_chg, prd_mul in zip(prd_ichs, prd_chgs, prd_muls):
                    species_info = (prd_ich, prd_chg, prd_mul)
                    well_data = moldr.driver.species_block(
                        species_info=species_info,
                        tors_model=tors_model,
                        vib_model=vib_model,
                        har_level=har_level,
                        tors_level=tors_level,
                        vpt2_level=vpt2_level,
                        script_str=PROJROT_SCRIPT_STR,
                        elec_levels=elec_levels,
                        sym_factor=sym_factor,
                        save_prefix=spc_save_path,
                        )
                well_label = 'W'+str(indxw)
                prod_label = well_label
                zero_energy = 0.0
                well_str += mess_io.writer.write_well(well_label, well_data, zero_energy)
    # write W_indxw 
            else:
                indxp += 1
                species_data = ['', '']
                species_label = ['', '']
                bimol_label = 'P'+str(indxp)
                prod_label = bimol_label
                for idx, (rct_ich, rct_chg, rct_mul) in enumerate(zip(rct_ichs, rct_chgs, rct_muls)):
                    species_info = (rct_ich, rct_chg, rct_mul)
                    species_data[idx] = moldr.driver.species_block(
                        species_info=species_info,
                        tors_model=tors_model,
                        vib_model=vib_model,
                        har_level=har_level,
                        tors_level=tors_level,
                        vpt2_level=vpt2_level,
                        script_str=PROJROT_SCRIPT_STR,
                        elec_levels=elec_levels,
                        sym_factor=sym_factor,
                        save_prefix=spc_save_path,
                        )
                    species_label[idx] = automol.inchi.smiles(rct_ich)
                ground_energy = 0.0
                bim_str += mess_io.writer.write_bimolecular(
                    bimol_label, species_label[0], species_data[0],
                    species_label[1], species_data[1], ground_energy)

            elec_levels = [[0., mul]]
            symfactor = 1.
#            prd_info=[prd_ichs[i], prd_chgs[i], prd_muls[i]]
            # is it OK to use product info??
#            spc_save_path = moldr.util.species_path(
#                prd_ichs[i], prd_chgs[i], prd_muls[i], SAVE_PREFIX)
#            thy_save_path = moldr.util.theory_path(method, basis, ts_orb_restr, spc_save_path)
            ts_data_str = moldr.driver.species_block(
                species_info=ts_info,
                tors_model=tors_model,
                vib_model=vib_model,
                har_level=har_level,
                tors_level=tors_level,
                vpt2_level=vpt2_level,
                script_str=PROJROT_SCRIPT_STR,
                elec_levels=elec_levels,
                sym_factor=sym_factor,
                save_prefix=thy_save_path,
                )
            ts_label = 'B1'

            ts_str = mess_io.writer.write_ts_sadpt(ts_label, reac_label, prod_label, ts_data_str, zero_energy)
            print(ts_str)

#                    ts_sadpt_writer
            core = mess_io.writer.write_core_rigidrotor(geo, symfactor)
            molecule_section_str1 = mess_io.writer.write_molecule(
                core, freqs, zpe, elec_levels, hind_rot='',
            )
            print(molecule_section_str1)
    #
            mess_inp_str = '/n'.join([header_str, energy_trans_str, well_str, bim_str, ts_str])

            bld_locs = ['MESS', 0]
            bld_save_fs = autofile.fs.build(thy_save_path)
            bld_save_fs.leaf.create(bld_locs)
            mess_path = bld_save_fs.leaf.path(bld_locs)
            print('Build Path for MESS rate files:')
            print(mess_path)
            with open(os.path.join(mess_path, 'mess.inp'), 'w') as mess_file:
                mess_file.write(mess_inp_str[name])
            moldr.util.run_script(RATE_SCRIPT_STR, mess_path)
                
#            rct_block_str = []
#            for i, _ in enumerate(rct_ichs):
#                elec_levels = [[0., mul]]
#                if (ich, mul) in ELC_DEG_DCT:
#                    elec_levels = ELC_DEG_DCT[(ich, mul)]
#                symfactor = 1.
#                rct_info=[rct_ichs[i], rct_chgs[i], rct_muls[i]]
#                orb_restr = moldr.util.orbital_restriction(rct_muls[i], RESTRICT_OPEN_SHELL)
#                spc_save_path = moldr.util.species_path(
#                    rct_ichs[i], rct_chgs[i], rct_muls[i], SAVE_PREFIX)
#                thy_save_path = moldr.util.theory_path(method, basis, orb_restr, spc_save_path)
#                rct_block_str.append(moldr.driver.species_block(
#                        species_info=rct_info,
#                        tors_model=tors_model,
#                        vib_model=vib_model,
#                        har_level=har_level,
#                        tors_level=tors_level,
#                        vpt2_level=vpt2_level,
#                        orb_restr=orb_restr,
#                        script_str=PROJROT_SCRIPT_STR,
#                        elec_levels=elec_levels,
#                        sym_factor=sym_factor,
#                        save_prefix=thy_save_path,
#                        ))
#            indx_well = 0
#            indx_bim = 0
#            if rct_ichs[1]:
#                print(rct_ichs)
#                indx_well += 1
#                reac_label = 'W{}'.format(indx_well)
#
## write W_indxw
#            else:
#                indx_bim += 1
#
#                prd_block_str = []
#                for i, _ in enumerate(prd_ichs):
#                    elec_levels = [[0., mul]]
#                    if (ich, mul) in ELC_DEG_DCT:
#                        elec_levels = ELC_DEG_DCT[(ich, mul)]
#                    symfactor = 1.
#                    prd_info=[prd_ichs[i], prd_chgs[i], prd_muls[i]]
#                    orb_restr = moldr.util.orbital_restriction(prd_muls[i], RESTRICT_OPEN_SHELL)
#                    spc_save_path = moldr.util.species_path(
#                            prd_ichs[i], prd_chgs[i], prd_muls[i], SAVE_PREFIX)
#                    thy_save_path = moldr.util.theory_path(method, basis, orb_restr, spc_save_path)
#                    prd_block_str.append(moldr.driver.species_block(
#                            species_info=prd_info,
#                            tors_model=tors_model,
#                            vib_model=vib_model,
#                            har_level=har_level,
#                            tors_level=tors_level,
#                            vpt2_level=vpt2_level,
#                            script_str=PROJROT_SCRIPT_STR,
#                            elec_levels=elec_levels,
#                            sym_factor=sym_factor,
#                            save_prefix=thy_save_path,
#                            ))


        # write P_indxp
        # cycle over vdw Species
        # cycle over transition states
        # for now - have only one TS

        #        if reaction_typ = addition:
        # reactants
        #            print(species_str(rct1))
        #            print(species_str(rct2))
        # well
        #            print(species_str(prod1))
        #        if reaction_typ = abstraction:
        # reactants
        #            print(species_str(rct1))
        #            print(species_str(rct2))
        # vdw
        #           for vdw_species in ...
        #                   print(species_str(vdwi))
        # products
        #            print(species_str(prod1))
        #            print(species_str(prod2))

        #        if reaction_typ = abstraction
        # reactants
        #            print(species_str(rct1))
        #            print(species_str(rct2))
        # vdw
        #            print(species_str(vdw1))
        #            print(species_str(vdw2))
        # products
        #            print(species_str(prod1))
        #            print(species_str(prod2))
        # ts

sys.exit()
