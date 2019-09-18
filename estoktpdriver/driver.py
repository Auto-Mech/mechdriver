""" reaction list test
"""
import os
from itertools import chain
import collections
import json
import pandas
from qcelemental import constants as qcc
import thermo
import chemkin_io
import automol
from automol import formula
import moldr
import thermodriver

ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')
WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

# 0. choose which mechanism to run

# MECHANISM_NAME = 'ch4+nh2'  # options: syngas, natgas, heptane
#MECHANISM_NAME = 'isooctane'  # options: syngas, natgas, heptane
#MECHANISM_NAME = 'natgas'  # options: syngas, natgas, heptane
#MECHANISM_NAME = 'butane'  # options: syngas, natgas, heptane
MECHANISM_NAME = 'syngas'  # options: syngas, natgas, heptane
#MECHANISM_NAME = 'vhp'  # options: syngas, natgas, heptane
# MECHANISM_NAME = 'onereac'  # options: syngas, natgas, heptane
# MECHANISM_NAME = 'estoktp/add30'  # options: syngas, natgas
# MECHANISM_NAME = 'estoktp/habs65'  # options: syngas, natgas

# 1. create run and save directories
RUN_PREFIX = '/lcrc/project/PACC/run'
if not os.path.exists(RUN_PREFIX):
    os.mkdir(RUN_PREFIX)

SAVE_PREFIX = '/lcrc/project/PACC/save'
if not os.path.exists(SAVE_PREFIX):
    os.mkdir(SAVE_PREFIX)

# 2. Prepare species and reaction lists and dictionaries

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

# read in data from the mechanism directory
DATA_PATH = os.path.dirname(os.path.realpath(__file__))
MECH_PATH = os.path.join(DATA_PATH, 'data', MECHANISM_NAME)
MECH_TYPE = 'CHEMKIN'
#MECH_TYPE = 'json'
#MECH_FILE = 'mech.json'
RAD_RAD_SORT = True
if MECH_TYPE == 'CHEMKIN':
    MECH_STR = open(os.path.join(MECH_PATH, 'mechanism.txt')).read()
    SPC_TAB = pandas.read_csv(os.path.join(MECH_PATH, 'smiles.csv'))
# 4. process species data from the mechanism file
# Also add in basis set species

    SPC_TAB['charge'] = 0
    #print('SPC_TAB:', SPC_TAB)
    #print('SPC_TAB TEST', SPC_TAB['name'])
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
    SPC_DCT = {}
    for name in SPC_NAMES:
        SPC_DCT[name] = {}
        smi = SMI_DCT[name]
        ich = automol.smiles.inchi(smi)
        chg = CHG_DCT[name]
        mul = MUL_DCT[name]
        SPC_INFO[name] = [ich, chg, mul]
        SPC_DCT[name]['ich'] = ich
        SPC_DCT[name]['chg'] = chg
        SPC_DCT[name]['mul'] = mul

    RXN_BLOCK_STR = chemkin_io.reaction_block(MECH_STR)
    RXN_STRS = chemkin_io.reaction.data_strings(RXN_BLOCK_STR)
    RCT_NAMES_LST = list(
        map(chemkin_io.reaction.DataString.reactant_names, RXN_STRS))
    PRD_NAMES_LST = list(
        map(chemkin_io.reaction.DataString.product_names, RXN_STRS))

    # Sort reactant and product name lists by formula to facilitate
    # multichannel, multiwell rate evaluations

    FORMULA_STR = ''
    RXN_NAME_LST = []
    FORMULA_STR_LST = []
    for rct_names, prd_names in zip(RCT_NAMES_LST, PRD_NAMES_LST):
        rxn_name = '='.join(['+'.join(rct_names), '+'.join(prd_names)])
        RXN_NAME_LST.append(rxn_name)
        rct_smis = list(map(SMI_DCT.__getitem__, rct_names))
        rct_ichs = list(map(automol.smiles.inchi, rct_smis))
        prd_smis = list(map(SMI_DCT.__getitem__, prd_names))
        prd_ichs = list(map(automol.smiles.inchi, prd_smis))
        formula_dict = ''
        for rct_ich in rct_ichs:
            formula_i = thermo.util.inchi_formula(rct_ich)
            formula_i_dict = thermo.util.get_atom_counts_dict(formula_i)
            formula_dict = automol.formula._formula.join(formula_dict, formula_i_dict)
        formula_dict = collections.OrderedDict(sorted(formula_dict.items()))
        FORMULA_STR = ''.join(map(str, chain.from_iterable(formula_dict.items())))
        FORMULA_STR_LST.append(FORMULA_STR)


    #print('formula string test list', formula_str_lst, 'rxn name list', rxn_name_lst)
    #for rct_names in RCT_NAMES_LST:
    #    print('first rct names test:', rct_names)
    RXN_INFO_LST = list(zip(FORMULA_STR_LST, RCT_NAMES_LST, PRD_NAMES_LST, RXN_NAME_LST))
    #for _, rct_names_lst, _, _ in RXN_INFO_LST:
    #    print('second rct names test:', rct_names_lst)
    #    for rct_names in rct_names_lst:
    #        print('rct names test:', rct_names)
    RXN_INFO_LST.sort()
    FORMULA_STR_LST, RCT_NAMES_LST, PRD_NAMES_LST, RXN_NAME_LST = zip(*RXN_INFO_LST)
    for rxn_name in RXN_NAME_LST:
        print("Reaction: {}".format(rxn_name))
    #for rct_names in RCT_NAMES_LST:
    #    print('rct names test:', rct_names)
    for formula in FORMULA_STR_LST:
        print("Reaction: {}".format(formula))

elif MECH_TYPE == 'json':
    with open(os.path.join(MECH_PATH, MECH_FILE)) as f:
        MECH_DATA_IN = json.load(f, object_pairs_hook=collections.OrderedDict)
        MECH_DATA = []
    for reaction in MECH_DATA_IN:
        if isinstance(reaction, dict):
            MECH_DATA = MECH_DATA_IN
            break
        else:
            for entry in MECH_DATA_IN[reaction]:
                MECH_DATA.append(entry)

#    print(mech_data)
#    SENS_SORT = sorted(mech_data, key = lambda i: i['Sensitivity'])
    FORMULA_STR = ''
#    RXN_NAME_LST = []
    FORMULA_STR_LST = []
#    print ('mech_data test', mech_data)
    RXN_NAME_LST = []
    RCT_ICHS_LST = []
    RCT_MULS_LST = []
    RCT_NAMES_LST = []
    PRD_ICHS_LST = []
    PRD_MULS_LST = []
    PRD_NAMES_LST = []
    RXN_SENS = []
    RXN_UNC = []
    RXN_VAL = []
    RXN_FAM = []
    for reaction in MECH_DATA:
        # set up reaction info
#        print('reaction test:', reaction)
        rct_ichs = []
        rct_muls = []
        rct_names = []
        prd_ichs = []
        prd_muls = []
        prd_names = []
#        print('reaction name json test:', reaction['name'])
        if 'Reactants' in reaction and 'Products' in reaction:
            for i, rct in enumerate(reaction['Reactants']):
                rct_ichs.append(rct['InChi'])
                rct_muls.append(rct['multiplicity'])
                rct_names.append(rct['name'])
            rad_rad_reac = True
            if len(rct_ichs) == 1:
                rad_rad_reac = False
            else:
                if min(rct_muls) == 1:
                    rad_rad_reac = False
            for i, prd in enumerate(reaction['Products']):
                prd_names.append(prd['name'])
                prd_ichs.append(prd['InChi'])
                prd_muls.append(prd['multiplicity'])
#                print('prd_muls test:', prd['name'], prd['multiplicity'])
            rad_rad_prod = True
            if len(prd_ichs) == 1:
                rad_rad_prod = False
            else:
                if min(prd_muls) == 1:
                    rad_rad_prod = False
            if RAD_RAD_SORT and not rad_rad_reac and not rad_rad_prod:
                continue
            RCT_ICHS_LST.append(rct_ichs)
            RCT_MULS_LST.append(rct_muls)
            RCT_NAMES_LST.append(rct_names)
            PRD_ICHS_LST.append(prd_ichs)
            PRD_MULS_LST.append(prd_muls)
            PRD_NAMES_LST.append(prd_names)
#        print(reaction['name'], RAD_RAD_SORT, rad_rad_reac, rad_rad_prod)
        RXN_NAME_LST.append(reaction['name'])
        if 'Sensitivity' in reaction:
            RXN_SENS.append(reaction['Sensitivity'])
        else:
            RXN_SENS.append('')
        if 'Uncertainty' in reaction:
            RXN_UNC.append(reaction['Uncertainty'])
        else:
            RXN_UNC.append('')
        if 'Value' in reaction:
            RXN_VAL.append(reaction['Value'])
        else:
            RXN_VAL.append('')
        if 'Family' in reaction:
            RXN_FAM.append(reaction['Family'])
        else:
            RXN_FAM.append('')

        formula = ''
#        formula_dict = ''
        for rct_ich in rct_ichs:
            formula_i = automol.inchi.formula(rct_ich)
#            formula_i = thermo.util.inchi_formula(rct_ich)
#            formula_i_dict = thermo.util.get_atom_counts_dict(formula_i)
#            print('formula test:', formula_i_test, formula_i, formula_i_dict)
            formula = automol.formula._formula.join(formula, formula_i)
        formula = collections.OrderedDict(sorted(formula.items()))
        FORMULA_STR = ''.join(map(str, chain.from_iterable(formula.items())))
#            formula_dict = automol.formula._formula.join(formula_dict, formula_i_dict)
#        formula_dict = collections.OrderedDict(sorted(formula_dict.items()))
#        FORMULA_STR = ''.join(map(str, chain.from_iterable(formula_dict.items())))
        FORMULA_STR_LST.append(FORMULA_STR)

    RXN_INFO_LST = list(zip(
        FORMULA_STR_LST, RCT_NAMES_LST, PRD_NAMES_LST,
        RXN_NAME_LST, RXN_SENS, RXN_UNC, RXN_VAL, RXN_FAM,
        RCT_ICHS_LST, RCT_MULS_LST, PRD_ICHS_LST,
        PRD_MULS_LST))
    #for _, rct_names_lst, _, _ in RXN_INFO_LST:
    #    print('second rct names test:', rct_names_lst)
    #    for rct_names in rct_names_lst:
    #        print('rct names test:', rct_names)
#    print('sens test:', RXN_SENS)
    RXN_INFO_LST = sorted(RXN_INFO_LST, key=lambda x: (x[0]))
    OLD_FORMULA = RXN_INFO_LST[0][0]
    SENS = RXN_INFO_LST[0][4]
    ORDERED_FORMULA = []
    ORDERED_SENS = []
    for entry in RXN_INFO_LST:
        formula = entry[0]
        if formula == OLD_FORMULA:
            SENS = max(SENS, entry[4])
        else:
            ORDERED_SENS.append(SENS)
            ORDERED_FORMULA.append(OLD_FORMULA)
            SENS = entry[4]
            OLD_FORMULA = formula
    ORDERED_SENS.append(SENS)
    ORDERED_FORMULA.append(OLD_FORMULA)
    SENS_DCT = {}
    for i, sens in enumerate(ORDERED_SENS):
        SENS_DCT[ORDERED_FORMULA[i]] = sens
    RXN_INFO_LST = sorted(RXN_INFO_LST, key=lambda x: (SENS_DCT[x[0]], x[4]), reverse=True)

#    SORTED_RXN_LST = list(zip(ordered_sens, ordered_formula))
#    SORTED_RXN_LST.sort()
#    ordered_sens, ordered_formula = zip(*SORTED_RXN_LST)
#    for ordered

#    RXN_INFO_LST = sorted(RXN_INFO_LST, key=lambda x: (x[0], x[4]))
#    RXN_INFO_LST.sort()
    FORMULA_STR_LST, RCT_NAMES_LST, PRD_NAMES_LST, RXN_NAME_LST, RXN_SENS, RXN_UNC, RXN_VAL, RXN_FAM, RCT_ICHS_LST, RCT_MULS_LST, PRD_ICHS_LST, PRD_MULS_LST = zip(*RXN_INFO_LST)
#    print('rxn_names test:', RXN_NAME_LST)
    for i, rxn_name in enumerate(RXN_NAME_LST):
        rct_smis = []
        for rct_ich in RCT_ICHS_LST[i]:
            rct_smis.append(automol.inchi.smiles(rct_ich))
        rct_smis = ' + '.join(rct_smis)
        prd_smis = []
        for prd_ich in PRD_ICHS_LST[i]:
            prd_smis.append(automol.inchi.smiles(prd_ich))
        prd_smis = ' + '.join(prd_smis)
        print("Reaction: {} {:.2f}".format(rxn_name, RXN_SENS[i]))
        print("          {}: {} <=> {}".format(FORMULA_STR_LST[i], rct_smis, prd_smis))

    # set up species info
    SPC_NAMES = []
    SPC_INFO = {}
    SMI_DCT = {}
    CHG_DCT = {}
    MUL_DCT = {}
    SPC_DCT = {}
    for i, spc_names_lst in enumerate(RCT_NAMES_LST):
        for j, spc_name in enumerate(spc_names_lst):
            chg = 0
            if spc_name not in SPC_NAMES:
                SPC_NAMES.append(spc_name)
                SPC_INFO[spc_name] = [RCT_ICHS_LST[i][j], chg, RCT_MULS_LST[i][j]]
                SMI_DCT[spc_name] = RCT_ICHS_LST[i][j]
                CHG_DCT[spc_name] = chg
                MUL_DCT[spc_name] = RCT_MULS_LST[i][j]
                SPC_DCT[spc_name]['chg'] = [chg]
                SPC_DCT[spc_name]['smi'] = RCT_ICHS_LST[i][j]
                SPC_DCT[spc_name]['mul'] = RCT_MULS_LST[i][j]
    for i, spc_names_lst in enumerate(PRD_NAMES_LST):
        for j, spc_name in enumerate(spc_names_lst):
            chg = 0
            if spc_name not in SPC_NAMES:
                SPC_NAMES.append(spc_name)
                SPC_INFO[spc_name] = [PRD_ICHS_LST[i][j], chg, PRD_MULS_LST[i][j]]
                SMI_DCT[spc_name] = PRD_ICHS_LST[i][j]
                CHG_DCT[spc_name] = chg
                MUL_DCT[spc_name] = PRD_MULS_LST[i][j]
                SPC_DCT[spc_name]['chg'] = [chg]
                SPC_DCT[spc_name]['smi'] = RCT_ICHS_LST[i][j]
                SPC_DCT[spc_name]['mul'] = RCT_MULS_LST[i][j]
    RXN_INFO_LST = list(zip(FORMULA_STR_LST, RCT_NAMES_LST, PRD_NAMES_LST, RXN_NAME_LST))

REF_MOLS = [automol.smiles.inchi('[H][H]'), automol.smiles.inchi('C'), automol.smiles.inchi('O')]
#os.sys.exit()

GEOM_PATH = os.path.join(DATA_PATH, 'data', 'geoms')
#print(GEOM_PATH)
GEOM_DCT = moldr.util.geometry_dictionary(GEOM_PATH)

# 2. script control parameters

# a. Strings to launch executable
# script_strings for electronic structure are obtained from run_qchem_par since
# they vary with method

PROJROT_SCRIPT_STR = ("#!/usr/bin/env bash\n"
                      "RPHt.exe")
PF_SCRIPT_STR = ("#!/usr/bin/env bash\n"
                 "messpf pf.inp build.out >> stdout.log &> stderr.log")
RATE_SCRIPT_STR = ("#!/usr/bin/env bash\n"
                   "mess mess.inp build.out >> stdout.log &> stderr.log")
VARECOF_SCRIPT_STR = ("#!/usr/bin/env bash\n"
                      "/home/ygeorgi/build/rotd/multi ")
MCFLUX_SCRIPT_STR = ("#!/usr/bin/env bash\n"
                     "/home/ygeorgi/build/rotd/mc_flux ")
CONV_MULTI_SCRIPT_STR = ("#!/usr/bin/env bash\n"
                         "/home/ygeorgi/build/rotd/mc_flux ")
TST_CHECK_SCRIPT_STR = ("#!/usr/bin/env bash\n"
                        "/home/ygeorgi/build/rotd/tst_check ")
MOLPRO_PATH_STR = ('/home/sjklipp/bin/molpro')
#NASA_SCRIPT_STR = ("#!/usr/bin/env bash\n"
#                   "cp ../PF/build.out pf.dat\n"
#                   "cp /tcghome/sjklipp/PACC/nasa/new.groups .\n"
#                   "python /tcghome/sjklipp/PACC/nasa/makepoly.py"
#                   " >> stdout.log &> stderr.log")

# b. Electronic structure parameters; code, method, basis, convergence control

mc_nsamp0 = [True, 3, 1, 3, 100]

ES_DCT = {
        'lvl_wbs': {
            'orb_res': 'RU', 'program': 'gaussian09', 'method': 'wb97xd', 'basis': '6-31g*',
            'mc_nsamp': mc_nsamp0
            },
        'lvl_wbm': {
            'orb_res': 'RU', 'program': 'gaussian09', 'method': 'wb97xd', 'basis': '6-31+g*',
            'mc_nsamp': mc_nsamp0
            },
        'lvl_wbt': {
            'orb_res': 'RU', 'program': 'gaussian09', 'method': 'wb97xd', 'basis': 'cc-pvtz',
            'mc_nsamp': mc_nsamp0
            },
        'lvl_b2t': {
            'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b2plypd3', 'basis': 'cc-pvtz',
            'mc_nsamp': mc_nsamp0
            },
        'lvl_b3s': {
            'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp', 'basis': '6-31g*',
            'mc_nsamp': mc_nsamp0
            },
        'lvl_b3t': {
            'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp', 'basis': '6-31g*',
            'mc_nsamp': mc_nsamp0
            },
        'cc_lvl_d': {
            'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)', 'basis': 'cc-pvdz',
            'mc_nsamp': mc_nsamp0
            },
        'cc_lvl_t': {
            'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)', 'basis': 'cc-pvtz',
            'mc_nsamp': mc_nsamp0
            },
        'cc_lvl_q': {
            'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)', 'basis': 'cc-pvqz',
            'mc_nsamp': mc_nsamp0
            },
        'cc_lvl_df': {
            'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)-f12',
            'basis': 'cc-pvdz-f12',
            'mc_nsamp': mc_nsamp0
            },
        'cc_lvl_tf': {
            'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)-f12',
            'basis': 'cc-pvtz-f12',
            'mc_nsamp': mc_nsamp0
            },
        'cc_lvl_qf': {
            'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)-f12',
            'basis': 'cc-pvqz-f12',
            'mc_nsamp': mc_nsamp0
            },
        }

OPT_LVL0 = 'lvl_wbs'
OPT_LVL1 = 'lvl_wbm'
OPT_LVL2 = 'lvl_wbt'
SCAN_LVL1 = OPT_LVL1
SP_LVL1 = 'cc_lvl_df'
SP_LVL2 = 'cc_lvl_tf'
SP_LVL3 = 'cc_lvl_qf'

TSK_INFO_LST = [
    ['conf_samp', OPT_LVL1, OPT_LVL0, False],
    ['conf_hess', OPT_LVL1, OPT_LVL1, False],
    ['hr_scan', SCAN_LVL1, OPT_LVL1, False],
    ['conf_energy', SP_LVL1, OPT_LVL1, False],
    ['conf_energy', SP_LVL2, OPT_LVL1, False],
    ]

HIGH_LEVEL_COEFF = [-0.5, 1.5]

OPT_ES = True
OPT_MESS = False
OPT_THERMO = False
OPT_ALLPF = False
OPTIONS = [OPT_ES, OPT_MESS, OPT_THERMO, OPT_ALLPF]

SPC_QUEUE = list(SPC_NAMES)

REF_MOLS='basic'
thermodriver.driver.run(
        TSK_INFO_LST, ES_DCT, SPC_DCT, SPC_QUEUE, REF_MOLS, RUN_PREFIX, SAVE_PREFIX, OPTIONS)

OPT_ES = False
OPT_MESS = True
OPT_THERMO = True
OPT_ALLPF = True
OPTIONS = [OPT_ES, OPT_MESS, OPT_THERMO, OPT_ALLPF]

thermodriver.driver.run(
        TSK_INFO_LST, ES_DCT, SPC_DCT, SPC_QUEUE, REF_MOLS, RUN_PREFIX, SAVE_PREFIX, OPTIONS)
# set up a combination of energies
# E_HL = sum_i E_HL(i) * Coeff(i)
#HIGH_LEVEL_COEFF = [1]
#HIGH_LEVEL_COEFF = None

# for now pf_levels and spc_models are automatically determined on basis of electronic structure levels
#PF_LEVELS = []
#PF_LEVELS.append([['', 'wb97xd', '6-31g*', 'RU'], ['', 'wb97xd', '6-31g*', 'RU'], ['', 'wb97xd', '6-31g*', 'RU']])
# PF_LEVELS contains the elec. struc. lvlels for the harmonic, torsional, and anharmonic analyses

# e. What to run for thermochemical kinetics
#RUN_SPC_THERMO = True
#SPC_MODELS = [['RIGID', 'HARM']]
#SPC_MODELS = [['1DHR', 'HARM']]
#SPC_MODELS = [['RIGID', 'HARM'], ['1DHR', 'HARM']]
# The first component specifies the torsional model - TORS_MODEL.
# It can take 'RIGID', '1DHR', or 'TAU' and eventually 'MDHR'
# The second component specifies the vibrational model - VIB_MODEL.
# It can take 'HARM', or 'VPT2' values.
#RUN_REACTION_RATES = True
#RUN_VDW_RCT_RATES = False
#RUN_VDW_PRD_RATES = False

# f. Partition function parameters
#TAU_PF_WRITE = True

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
EXP_POWER = 0.85
EXP_CUTOFF = 15.
EPS1 = 100.0
EPS2 = 200.0
SIG1 = 6.
SIG2 = 6.
MASS1 = 15.0
ETSFR_PAR = [EXP_FACTOR, EXP_POWER, EXP_CUTOFF, EPS1, EPS2, SIG1, SIG2, MASS1]
