""" reaction list test
"""
import os
import sys
from itertools import chain
import collections
import json
import numpy
from qcelemental import constants as qcc
import thermo
import chemkin_io
import automol
from automol import formula
import moldr
import thermodriver
import ktpdriver

ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')
WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')
#PES_NUMS= 'all'
#PES_NUMS = [2, 3]
PESNUMS= 'all'
CHANNELS= 'all'
#PESNUMS=  [16,17,18,19,20]
# 0. choose which mechanism to run

MECHANISM_NAME = sys.argv[1]
MECH_TYPE = sys.argv[2]
if len(sys.argv) > 3:
    PESNUMS = sys.argv[3]
    if len(sys.argv) > 4:
        CHANNELS = sys.argv[4]
RUN_THERMO = False
RUN_RATES = True
# 1. create run and save directories
RUN_PREFIX = '/lcrc/project/PACC/run'
if not os.path.exists(RUN_PREFIX):
    os.mkdir(RUN_PREFIX)

SAVE_PREFIX = '/lcrc/project/PACC/save'
if not os.path.exists(SAVE_PREFIX):
    os.mkdir(SAVE_PREFIX)

# 2. Prepare special species and reaction dictionaries

ELC_SIG_LST = {'InChI=1S/CN/c1-2', 'InChI=1S/C2H/c1-2/h1H'}

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

SYMM_DCT = {
    ('InChI=1S/HO/h1H', 2): 1.
}

DATA_PATH = '/home/sjklipp/PACC/mech_test'
GEOM_PATH = os.path.join(DATA_PATH, 'data', 'geoms')
GEOM_DCT = moldr.util.geometry_dictionary(GEOM_PATH)

# 3. Prepare species and reaction lists 
#SMILES_TEST = '[C]#C'
#for smiles in SMILES_TEST:
#    ICH_TEST = (automol.convert.smiles.inchi(SMILES_TEST))
#    print(SMILES_TEST, ICH_TEST)

# read in data from the mechanism directory
MECH_PATH = os.path.join(DATA_PATH, 'data', MECHANISM_NAME)
MECH_FILE = 'mech.json'
RAD_RAD_SORT = False

# 4. process species data from the mechanism file
# Also add in basis set species
SORT_RXNS = False

if MECH_TYPE == 'CHEMKIN':

    MECH_STR = open(os.path.join(MECH_PATH, 'mechanism.txt')).read()
    SPC_STR = open(os.path.join(MECH_PATH, 'species.csv')).read()

    SMI_DCT = chemkin_io.mechparser.mechanism.species_name_smiles_dct(SPC_STR)
    ICH_DCT = chemkin_io.mechparser.mechanism.species_name_inchi_dct(SPC_STR)
    MUL_DCT = chemkin_io.mechparser.mechanism.species_name_mult_dct(SPC_STR)
    CHG_DCT = chemkin_io.mechparser.mechanism.species_name_charge_dct(SPC_STR)
    SENS_DCT = chemkin_io.mechparser.mechanism.species_name_sens_dct(SPC_STR)

    CHECK_STEREO = False

    if CHECK_STEREO:
        SPC_STR = 'name,SMILES,InChI,mult,charge,sens \n'
        for name in ICH_DCT:
            ich = ICH_DCT[name]
            # print('ich test:', ich)
            smi = SMI_DCT[name]
            mul = MUL_DCT[name]
            chg = CHG_DCT[name]
            sens = SENS_DCT[name]
            if not automol.inchi.is_complete(ich):
                print('adding stereochemistry for {0}, {1}, {2}'.format(name, smi, ich))
                # note this returns a list of ich's with the different possible stereo values
                # for now just taking the first of these
                ich = automol.inchi.add_stereo(ich)
                ich = ich[0]
                print('new ich:', ich)
                ICH_DCT[name] = ich
                # print('new dct:', ICH_DCT[name])
            SPC_STR += '{0},\'{1}\',\'{2}\',{3},{4},{5} \n'.format(
                name, smi, ich, mul, chg, sens)

        with open(os.path.join(MECH_PATH, 'species_stereo.csv'), 'w') as stereo_csv_file:
            stereo_csv_file.write(SPC_STR)

        # ICH_DCT = chemkin_io.mechparser.mechanism.species_name_inchi_dct(SPC_STR)

    #SPC_TAB['charge'] = 0
    #print('SPC_TAB:', SPC_TAB)
    #print('SPC_TAB TEST', SPC_TAB['name'])
    #SMI_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['smiles']))
    #SMI_DCT['REF_H2'] = '[H][H]'
    #SMI_DCT['REF_CH4'] = 'C'
    #SMI_DCT['REF_H2O'] = 'O'
    #SMI_DCT['REF_NH3'] = 'N'
    #CHG_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['charge']))
    #CHG_DCT['REF_H2'] = 0
    #CHG_DCT['REF_CH4'] = 0
    #CHG_DCT['REF_H2O'] = 0
    #CHG_DCT['REF_NH3'] = 0
    #MUL_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['mult']))
    #MUL_DCT['REF_H2'] = 1
    #MUL_DCT['REF_CH4'] = 1
    #MUL_DCT['REF_H2O'] = 1
    #MUL_DCT['REF_NH3'] = 1
    #SPC_BLK_STR = chemkin_io.mechparser.mechanism.species_block(MECH_STR)
    #SPC_NAMES = chemkin_io.mechparser.species.names(SPC_BLK_STR)

    # You need one species reference for each element in the set of references and species
    #SPC_REF_NAMES = ('REF_H2', 'REF_CH4', 'REF_H2O', 'REF_NH3')
    #SPC_REF_ICH = []
    #for ref_name in SPC_REF_NAMES:
    #    smi = SMI_DCT[ref_name]
    #    SPC_REF_ICH.append(automol.smiles.inchi(smi))
    #SPC_NAMES += SPC_REF_NAMES

    #print('SPC_NAMES')
    #print(SPC_NAMES)
    # SPC_INFO = {}
    SPC_DCT = {}
    for name in MUL_DCT:
        SPC_DCT[name] = {}
        # SPC_INFO[name] = [ICH_DCT[name], CHG_DCT[name], MUL_DCT[name]] 
        SPC_DCT[name]['smi'] = SMI_DCT[name]
        SPC_DCT[name]['ich'] = ICH_DCT[name]
        SPC_DCT[name]['chg'] = CHG_DCT[name]
        SPC_DCT[name]['mul'] = MUL_DCT[name]

    RXN_BLOCK_STR = chemkin_io.mechparser.mechanism.reaction_block(MECH_STR)
    RXN_STRS = chemkin_io.mechparser.reaction.data_strings(RXN_BLOCK_STR)
    RCT_NAMES_LST = list(
        map(chemkin_io.mechparser.reaction.reactant_names, RXN_STRS))
    PRD_NAMES_LST = list(
        map(chemkin_io.mechparser.reaction.product_names, RXN_STRS))

    # Sort reactant and product name lists by formula to facilitate
    # multichannel, multiwell rate evaluations

    FORMULA_STR = ''
    RXN_NAME_LST = []
    FORMULA_STR_LST = []
    for rct_names, prd_names in zip(RCT_NAMES_LST, PRD_NAMES_LST):
        rxn_name = '='.join(['+'.join(rct_names), '+'.join(prd_names)])
        RXN_NAME_LST.append(rxn_name)
        rct_smis = list(map(SMI_DCT.__getitem__, rct_names))
        rct_ichs = list(map(ICH_DCT.__getitem__, rct_names))
        prd_smis = list(map(SMI_DCT.__getitem__, prd_names))
        prd_ichs = list(map(ICH_DCT.__getitem__, prd_names))
        formula_dct = ''
        for rct_ich in rct_ichs:
            formula_i_dct = automol.inchi.formula_dct(rct_ich)
            formula_dct = automol.formula._formula.join(formula_dct, formula_i_dct)
        FORMULA_STR = automol.formula._formula.string(formula_dct)
        FORMULA_STR_LST.append(FORMULA_STR)

    RXN_INFO_LST = list(zip(FORMULA_STR_LST, RCT_NAMES_LST, PRD_NAMES_LST, RXN_NAME_LST))
    if SORT_RXNS:
        RXN_INFO_LST.sort()
        FORMULA_STR_LST, RCT_NAMES_LST, PRD_NAMES_LST, RXN_NAME_LST = zip(*RXN_INFO_LST)
        for rxn_name in RXN_NAME_LST:
            print("Reaction: {}".format(rxn_name))
        for formula in FORMULA_STR_LST:
            print("Reaction: {}".format(formula))

elif MECH_TYPE == 'json':
    CHECK_STEREO = False
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

    # first convert the essential pieces of the json file to chemkin formatted data so
    # (i) can easily remove species that don't really exist
    # (ii) revise products of reactions for species that don't exist
    # (iii) do the stereochemistry generation only one

    FORMULA_STR = ''
    FORMULA_STR_LST = []
    RXN_NAME_LST = []
    RCT_NAMES_LST = []
    PRD_NAMES_LST = []
    RCT_SMIS_LST = []
    RCT_ICHS_LST = []
    RCT_MULS_LST = []
    PRD_SMIS_LST = []
    PRD_ICHS_LST = []
    PRD_MULS_LST = []
    PRD_NAMES_LST = []
    RXN_SENS = []
    RXN_UNC = []
    RXN_VAL = []
    RXN_FAM = []
    UNQ_RXN_LST = []
    FLL_RXN_LST = []
    idxp = 0
    for idx, reaction in enumerate(MECH_DATA):
        if 'Reactants' in reaction and 'Products' in reaction:
            print(idx, reaction['name'])
            if reaction['name'] in FLL_RXN_LST:
                print('duplicate reaction found:', reaction['name'], idx)
            else:
                UNQ_RXN_LST.append(reaction['name'])
            FLL_RXN_LST.append(reaction['name'])
    print('reaction duplicate test:', len(UNQ_RXN_LST), len(FLL_RXN_LST))

    for ridx, reaction in enumerate(MECH_DATA):
        # set up reaction info
        rct_smis = []
        rct_ichs = []
        rct_muls = []
        rct_names = []
        prd_smis = []
        prd_ichs = []
        prd_muls = []
        prd_names = []
        if 'Reactants' in reaction and 'Products' in reaction:
            for rct in reaction['Reactants']:
                rct_names.append(rct['name'])
                rct_smis.append(rct['SMILES'][0])
                # rct_ichs_smi = automol.smiles.inchi(rct['SMILES'][0])
                # rct_ichs_rmg = rct['InChI']
                #if rct_ichs_smi != rct_ichs_rmg:
                    # print('Warning: RMG inchi {} differs from conversion from
                    # smiles {}{}:'.format(
                        # rct_ichs_rmg, rct_ichs_smi, rct['name']))
                ich = rct['InChi']
                # print('{} , {} , {}'.format(rct['name'], rct['SMILES'][0], rct['multiplicity']))
                if CHECK_STEREO:
                    if not automol.inchi.is_complete(ich):
                        print('adding stereochemsiry for {}'.format(ich))
                        ich = automol.inchi.add_stereo(rct['InChi'])[0]
                rct_ichs.append(ich)
                #rct_ichs.append(rct_ichs_smi)
                rct_muls.append(rct['multiplicity'])
            rad_rad_reac = True
            if len(rct_ichs) == 1:
                rad_rad_reac = False
            else:
                if min(rct_muls) == 1:
                    rad_rad_reac = False
            for prd in reaction['Products']:
                prd_names.append(prd['name'])
                prd_smis.append(prd['SMILES'][0])
                # prd_ichs_smi = automol.smiles.inchi(prd['SMILES'][0])
                # prd_ichs_rmg = prd['InChi']
                # if prd_ichs_smi != prd_ichs_rmg:
                    # print('Warning: RMG inchi {} differs from conversion from
                    # smiles {}{}:'.format(
                        # prd_ichs_rmg, prd_ichs_smi, prd['name']))
                # prd_ichs.append(prd_ichs_smi)
                ich = prd['InChi']
                # print('{} , {} , {}'.format(prd['name'], prd['SMILES'][0], prd['multiplicity']))
                if CHECK_STEREO:
                    if not automol.inchi.is_complete(ich):
                        print('adding stereochemsiry for {}'.format(ich))
                        ich = automol.inchi.add_stereo(prd['InChi'])[0]
                prd_ichs.append(ich)
                prd_muls.append(prd['multiplicity'])
            rad_rad_prod = True
            if len(prd_ichs) == 1:
                rad_rad_prod = False
            else:
                if min(prd_muls) == 1:
                    rad_rad_prod = False
            if RAD_RAD_SORT and not rad_rad_reac and not rad_rad_prod:
                continue
            RCT_SMIS_LST.append(rct_smis)
            RCT_ICHS_LST.append(rct_ichs)
            RCT_MULS_LST.append(rct_muls)
            RCT_NAMES_LST.append(rct_names)
            PRD_SMIS_LST.append(prd_smis)
            PRD_ICHS_LST.append(prd_ichs)
            PRD_MULS_LST.append(prd_muls)
            PRD_NAMES_LST.append(prd_names)
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

        # formula = ''
        formula_dct = ''
        for rct_ich in rct_ichs:
            formula_i_dct = automol.inchi.formula_dct(rct_ich)
            formula_dct = automol.formula._formula.join(formula_dct, formula_i_dct)
        FORMULA_STR = automol.formula._formula.string(formula_dct)
            # formula_i = automol.inchi.formula(rct_ich)
            # formula = automol.formula._formula.join(formula, formula_i)
        # formula = collections.OrderedDict(sorted(formula.items()))
        # FORMULA_STR = ''.join(map(str, chain.from_iterable(formula.items())))
        FORMULA_STR_LST.append(FORMULA_STR)

    UNQ_ICH_LST = []
    UNQ_MUL_LST = []
    UNQ_SMI_LST = []
    UNQ_LAB_LST = []
    UNQ_LAB_IDX_LST = []
    csv_str = 'name,SMILES,mult'
    csv_str += '\n'
    spc_str = 'SPECIES'
    spc_str += '\n'
    for ichs, muls, smis in zip(RCT_ICHS_LST, RCT_MULS_LST, RCT_SMIS_LST):
        for ich, mul, smi in zip(ichs, muls, smis):
            unique = True
            for unq_ich, unq_mul in zip(UNQ_ICH_LST, UNQ_MUL_LST):
                if ich == unq_ich and mul == unq_mul:
                    unique = False
            if unique:
                UNQ_ICH_LST.append(ich)
                UNQ_MUL_LST.append(mul)
                UNQ_SMI_LST.append(smi)
                #formula = ''
                #formula_i = automol.inchi.formula(ich)
                #formula = automol.formula._formula.join(formula, formula_i)
                #formula = collections.OrderedDict(sorted(formula.items()))
                #lab = ''.join(map(str, chain.from_iterable(formula.items())))

                formula_dct = automol.inchi.formula_dct(ich)
                lab = automol.formula._formula.string(formula_dct)

                UNQ_LAB_LST.append(lab)
                lab_idx = -1
                for lab_i in UNQ_LAB_LST:
                    if lab == lab_i:
                        lab_idx += 1
                UNQ_LAB_IDX_LST.append(lab_idx)
                if lab_idx == 0:
                    label = lab 
                else:
                    label = lab + '(' + str(lab_idx) + ')'
                csv_str += ','.join([label, smi, str(mul)])
                csv_str += '\n'
                spc_str += label
                spc_str += '\n'
    for ichs, muls, smis in zip(PRD_ICHS_LST, PRD_MULS_LST, PRD_SMIS_LST):
        for ich, mul, smi in zip(ichs, muls, smis):
            unique = True
            for unq_ich, unq_mul in zip(UNQ_ICH_LST, UNQ_MUL_LST):
                if ich == unq_ich and mul == unq_mul:
                    unique = False
            if unique:
                UNQ_ICH_LST.append(ich)
                UNQ_MUL_LST.append(mul)
                UNQ_SMI_LST.append(smi)
                #formula = ''
                #formula_i = automol.inchi.formula(ich)
                #formula = automol.formula._formula.join(formula, formula_i)
                #formula = collections.OrderedDict(sorted(formula.items()))
                #lab = ''.join(map(str, chain.from_iterable(formula.items())))

                formula_dct = automol.inchi.formula_dct(ich)
                lab = automol.formula._formula.string(formula_dct)

                UNQ_LAB_LST.append(lab)
                lab_idx = -1
                for lab_i in UNQ_LAB_LST:
                    if lab == lab_i:
                        lab_idx += 1
                UNQ_LAB_IDX_LST.append(lab_idx)
                if lab_idx == 0:
                    label = lab 
                else:
                    label = lab + '(' + str(lab_idx) + ')'
                csv_str += ','.join([label, smi, str(mul)])
                csv_str += '\n'
                spc_str += label
                spc_str += '\n'

    spc_str += 'END'
    spc_str += '\n'
    spc_str += '\n'

    with open(os.path.join(MECH_PATH, 'smiles_sort.csv'), 'w') as sorted_csv_file:
        sorted_csv_file.write(csv_str)

    # UNQ_SMI_LST = []
    # rsmi = 0
    # psmi = 0
    # spc_str = 'SPECIES'
    # spc_str += '\n'
    # for smis in RCT_SMIS_LST:
        # for smi in smis:
            # rsmi += 1
            # print('rct_smi test:', smi)
            # if smi not in UNQ_SMI_LST:
                # UNQ_SMI_LST.append(smi)
                # print('unq_smi test:', smi)
                # spc_str += smi
                # spc_str += '\n'
    # for smis in PRD_SMIS_LST:
        # for smi in smis:
            # psmi += 1
            # print('prd_smi test:', smi)
            # if smi not in UNQ_SMI_LST:
                # UNQ_SMI_LST.append(smi)
                # print('unq_smi test:', smi)
                # spc_str += smi
                # spc_str += '\n'
    # spc_str += 'END'
    # spc_str += '\n'
        
    RXN_INFO_LST = list(zip(
        FORMULA_STR_LST, RCT_NAMES_LST, PRD_NAMES_LST, RXN_NAME_LST, RXN_SENS,
        RXN_UNC, RXN_VAL, RXN_FAM, RCT_SMIS_LST, RCT_ICHS_LST, RCT_MULS_LST,
        PRD_SMIS_LST, PRD_ICHS_LST, PRD_MULS_LST))
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

    FORMULA_STR_LST, RCT_NAMES_LST, PRD_NAMES_LST, RXN_NAME_LST, RXN_SENS, RXN_UNC, RXN_VAL, RXN_FAM, RCT_SMIS_LST, RCT_ICHS_LST, RCT_MULS_LST, PRD_SMIS_LST, PRD_ICHS_LST, PRD_MULS_LST = zip(*RXN_INFO_LST)

    RXN_NAMEP_LST = []
    rxn_namep_str = 'REACTIONS   KCAL/MOLE   MOLES'
    rxn_namep_str += '\n'
    for i, rxn_name in enumerate(RXN_NAME_LST):
        rxn_namep = []
        rct_labs = []
        for rct_smi, rct_ich, rct_mul in zip(RCT_SMIS_LST[i], RCT_ICHS_LST[i], RCT_MULS_LST[i]):
            for ich, mul, lab, lab_idx in zip(UNQ_ICH_LST, UNQ_MUL_LST, UNQ_LAB_LST, UNQ_LAB_IDX_LST):
                if rct_ich == ich and rct_mul == mul:
                    if lab_idx == 0:
                        rct_lab = lab
                    else:
                        rct_lab = lab + '(' + str(lab_idx) + ')'
                    break
            rct_labs.append(rct_lab)
        rct_label = '+'.join(rct_labs)
        prd_labs = []
        for prd_smi, prd_ich, prd_mul in zip(PRD_SMIS_LST[i], PRD_ICHS_LST[i], PRD_MULS_LST[i]):
            for ich, mul, lab, lab_idx in zip(UNQ_ICH_LST, UNQ_MUL_LST, UNQ_LAB_LST, UNQ_LAB_IDX_LST):
                if prd_ich == ich and prd_mul == mul:
                    if lab_idx == 0:
                        prd_lab = lab
                    else:
                        prd_lab = lab + '(' + str(lab_idx) + ')'
                    break
            prd_labs.append(prd_lab)
        prd_label = '+'.join(prd_labs)
        rate_str = str('  1.e10   1.0   10000.  ! Sens = ')
        rxn_namep = rct_label + ' <=> ' + prd_label + rate_str + str(RXN_SENS[i])
        rxn_namep_str += rxn_namep 
        rxn_namep_str += '\n'
        RXN_NAMEP_LST.append(rxn_namep)

    mech_str = spc_str + rxn_namep_str
    mech_str += 'END'
    mech_str += '\n'

    with open(os.path.join(MECH_PATH, 'mech_sort.txt'), 'w') as sorted_mech_file:
        sorted_mech_file.write(mech_str)

    # set up species info
    SPC_NAMES = []
    # SPC_INFO = {}
    CHG_DCT = {}
    MUL_DCT = {}
    SPC_DCT = {}
    for i, spc_names_lst in enumerate(RCT_NAMES_LST):
        for j, spc_name in enumerate(spc_names_lst):
            chg = 0
            if spc_name not in SPC_NAMES:
                SPC_NAMES.append(spc_name)
                # SPC_INFO[spc_name] = [RCT_ICHS_LST[i][j], chg, RCT_MULS_LST[i][j]]
                CHG_DCT[spc_name] = chg
                MUL_DCT[spc_name] = RCT_MULS_LST[i][j]
                SPC_DCT[spc_name] = {}
                SPC_DCT[spc_name]['chg'] = chg
                SPC_DCT[spc_name]['ich'] = RCT_ICHS_LST[i][j]
                SPC_DCT[spc_name]['mul'] = RCT_MULS_LST[i][j]
    for i, spc_names_lst in enumerate(PRD_NAMES_LST):
        for j, spc_name in enumerate(spc_names_lst):
            chg = 0
            if spc_name not in SPC_NAMES:
                SPC_NAMES.append(spc_name)
                # SPC_INFO[spc_name] = [PRD_ICHS_LST[i][j], chg, PRD_MULS_LST[i][j]]
                CHG_DCT[spc_name] = chg
                MUL_DCT[spc_name] = PRD_MULS_LST[i][j]
                SPC_DCT[spc_name] = {}
                SPC_DCT[spc_name]['chg'] = chg
                SPC_DCT[spc_name]['ich'] = PRD_ICHS_LST[i][j]
                SPC_DCT[spc_name]['mul'] = PRD_MULS_LST[i][j]
    RXN_INFO_LST = list(zip(FORMULA_STR_LST, RCT_NAMES_LST, PRD_NAMES_LST, RXN_NAME_LST))

PES_LST = {}
current_formula = ''
for fidx, formula in enumerate(FORMULA_STR_LST):
    if current_formula == formula:
        PES_LST[formula]['RCT_NAMES_LST'].append(RCT_NAMES_LST[fidx])
        PES_LST[formula]['PRD_NAMES_LST'].append(PRD_NAMES_LST[fidx])
        PES_LST[formula]['RXN_NAME_LST'].append(RXN_NAME_LST[fidx])
    else:
        current_formula = formula
        PES_LST[formula] = {}
        PES_LST[formula]['RCT_NAMES_LST'] = [RCT_NAMES_LST[fidx]]
        PES_LST[formula]['PRD_NAMES_LST'] = [PRD_NAMES_LST[fidx]]
        PES_LST[formula]['RXN_NAME_LST'] = [RXN_NAME_LST[fidx]]

for spc in SPC_DCT:
    if tuple([SPC_DCT[spc]['ich'], SPC_DCT[spc]['mul']]) in ELC_DEG_DCT:
        SPC_DCT[spc]['elec_levs'] = ELC_DEG_DCT[SPC_DCT[spc]['ich'], SPC_DCT[spc]['mul']]
    if tuple([SPC_DCT[spc]['ich'], SPC_DCT[spc]['mul']]) in SYMM_DCT:
        SPC_DCT[spc]['sym'] = SYMM_DCT[SPC_DCT[spc]['ich'], SPC_DCT[spc]['mul']]

for spc in SPC_DCT:
    ich = SPC_DCT[spc]['ich']
    if ich in GEOM_DCT:
        SPC_DCT[spc]['geo_obj'] = GEOM_DCT[ich]

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

MC_NSAMP0 = [True, 4, 1, 3, 100]

ES_DCT = {
        'lvl_wbs': {
            'orb_res': 'RU', 'program': 'gaussian09', 'method': 'wb97xd', 'basis': '6-31g*',
            'mc_nsamp': MC_NSAMP0
            },
        'lvl_wbm': {
            'orb_res': 'RU', 'program': 'gaussian09', 'method': 'wb97xd', 'basis': '6-31+g*',
            'mc_nsamp': MC_NSAMP0
            },
        'lvl_wbt': {
            'orb_res': 'RU', 'program': 'gaussian09', 'method': 'wb97xd', 'basis': 'cc-pvtz',
            'mc_nsamp': MC_NSAMP0
            },
        'lvl_b2d': {
            'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b2plypd3', 'basis': 'cc-pvdz',
            'mc_nsamp': MC_NSAMP0
            },
        'lvl_b2t': {
            'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b2plypd3', 'basis': 'cc-pvtz',
            'mc_nsamp': MC_NSAMP0
            },
        'lvl_b2q': {
            'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b2plypd3', 'basis': 'cc-pvqz',
            'mc_nsamp': MC_NSAMP0
            },
        'lvl_b3s': {
            'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp', 'basis': '6-31g*',
            'mc_nsamp': MC_NSAMP0
            },
        'lvl_b3t': {
            'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp', 'basis': '6-31g*',
            'mc_nsamp': MC_NSAMP0
            },
        'cc_lvl_d': {
            'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)', 'basis': 'cc-pvdz',
            'mc_nsamp': MC_NSAMP0
            },
        'cc_lvl_t': {
            'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)', 'basis': 'cc-pvtz',
            'mc_nsamp': MC_NSAMP0
            },
        'cc_lvl_q': {
            'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)', 'basis': 'cc-pvqz',
            'mc_nsamp': MC_NSAMP0
            },
        'cc_lvl_df': {
            'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)-f12',
            'basis': 'cc-pvdz-f12',
            'mc_nsamp': MC_NSAMP0
            },
        'cc_lvl_tf': {
            'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)-f12',
            'basis': 'cc-pvtz-f12',
            'mc_nsamp': MC_NSAMP0
            },
        'cc_lvl_qf': {
            'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)-f12',
            'basis': 'cc-pvqz-f12',
            'mc_nsamp': MC_NSAMP0
            },
        }

OPT_LVL0 = 'lvl_wbs'
OPT_LVL1 = 'lvl_wbm'
OPT_LVL2 = 'lvl_b2t'
SCAN_LVL1 = OPT_LVL0
SP_LVL1 = 'cc_lvl_df'
SP_LVL2 = 'cc_lvl_tf'
SP_LVL3 = 'cc_lvl_qf'

# The logic key in tsk_info_lst is for overwrite
OVERWRITE = False

if RUN_THERMO:
    TSK_INFO_LST = [
        ['find_geom', OPT_LVL0, OPT_LVL0, OVERWRITE],
        ['conf_samp', OPT_LVL1, OPT_LVL0, OVERWRITE],
        ['conf_hess', OPT_LVL1, OPT_LVL1, OVERWRITE],
        ['find_geom', OPT_LVL2, OPT_LVL1, OVERWRITE],
        ['conf_hess', OPT_LVL2, OPT_LVL2, OVERWRITE],
        ['hr_scan', SCAN_LVL1, OPT_LVL1, OVERWRITE],
        ['conf_energy', SP_LVL1, OPT_LVL2, OVERWRITE],
        ['conf_energy', SP_LVL2, OPT_LVL2, OVERWRITE],
        ['conf_energy', SP_LVL3, OPT_LVL2, OVERWRITE],
        ['sym_samp', OPT_LVL1, OPT_LVL1, OVERWRITE],
        # ['conf_vpt2', OPT_LVL1, OPT_LVL1, OVERWRITE],
        ]

    ENE_COEFF = [0.0, -0.693, 1.693]
    #ENE_COEFF = [-0.4, 1.4]

    OPT_ES = True
    OPT_MESS = False
    OPT_THERMO = False
    OPT_ALLPF = False
    OPTIONS = [OPT_ES, OPT_MESS, OPT_THERMO, OPT_ALLPF]

    SPC_QUEUE = list(SPC_NAMES)
    SPC_QUEUE = ['H2O2(11)']
    # SPC_QUEUE = ['CHOH']
    # SPC_QUEUE = ['CH2OH']

    REF_MOLS = 'basic'
    # REF_MOLS = [
        # automol.smiles.inchi('[H][H]'), automol.smiles.inchi('O')
        # automol.smiles.inchi('C'),
        # automol.smiles.inchi('O'), automol.smiles.inchi('N')]
    # REF_MOLS = [
        # automol.smiles.inchi('[H][H]'), automol.smiles.inchi('C'),
        # automol.smiles.inchi('O'), automol.smiles.inchi('N')]
    # REF_MOLS = [
        # automol.smiles.inchi('[H]'), automol.smiles.inchi('[C]'),
        # automol.smiles.inchi('[O]'), automol.smiles.inchi('[N]')]
    # REF_MOLS = [
        # automol.smiles.inchi('[H][H]'), automol.smiles.inchi('C'),
        # automol.smiles.inchi('O=O'), automol.smiles.inchi('N#N')]
        # ]
    thermodriver.driver.run(
        TSK_INFO_LST, ES_DCT, SPC_DCT, SPC_QUEUE, REF_MOLS, RUN_PREFIX,
        SAVE_PREFIX, options=OPTIONS)

    OPT_ES = False
    OPT_MESS = True
    OPT_THERMO = True
    OPT_ALLPF = True
    OPTIONS = [OPT_ES, OPT_MESS, OPT_THERMO, OPT_ALLPF]

    thermodriver.driver.run(
        TSK_INFO_LST, ES_DCT, SPC_DCT, SPC_QUEUE, REF_MOLS, RUN_PREFIX,
        SAVE_PREFIX, ENE_COEFF, options=OPTIONS)

if RUN_RATES:

    OPT_ES = True
    OPT_MESS = False
    OPT_THERMO = False
    OPT_ALLPF = False
    OPTIONS = [OPT_ES, OPT_MESS, OPT_THERMO, OPT_ALLPF]

    TSK_INFO_LST = [
        ['find_geom', OPT_LVL0, OPT_LVL0, OVERWRITE],
        ['conf_samp', OPT_LVL0, OPT_LVL0, OVERWRITE],
        ['conf_hess', OPT_LVL0, OPT_LVL0, OVERWRITE],
        ['hr_scan', SCAN_LVL1, OPT_LVL0, OVERWRITE],
        # ['sym_samp', OPT_LVL0, OPT_LVL0, OVERWRITE],
        ['find_ts', OPT_LVL0, OPT_LVL0, OVERWRITE],
        ['conf_samp', OPT_LVL0, OPT_LVL0, OVERWRITE],
        ['conf_hess', OPT_LVL0, OPT_LVL0, OVERWRITE],
        ['hr_scan', SCAN_LVL1, OPT_LVL0, OVERWRITE],
        ['conf_energy', SP_LVL1, OPT_LVL0, OVERWRITE],
        # ['conf_energy', SP_LVL1, OPT_LVL1, OVERWRITE],
        # ['conf_energy', SP_LVL2, OPT_LVL1, OVERWRITE],
        # ['conf_energy', SP_LVL3, OPT_LVL1, OVERWRITE],
        # ['sym_samp', OPT_LVL0, OPT_LVL0, OVERWRITE],
        # ['conf_vpt2', OPT_LVL1, OPT_LVL1, OVERWRITE],
        ]
    if isinstance(PESNUMS, str):
        if PESNUMS== 'all':
            PESNUMS = numpy.arange(len(PES_LST)+1)
        elif '-' in PESNUMS:
            start, end = PESNUMS.split('-')
            PESNUMS = numpy.arange(int(start), int(end)+1)
        elif '[' in PESNUMS:
            nums = PESNUMS.replace('[','').replace(']','').split(',')
            PESNUMS = [int(num) for num in nums]
    # loop over PESs
    for pes_idx, PES in enumerate(PES_LST, start=1):
        if pes_idx in PESNUMS:
            PES_RCT_NAMES_LST = PES_LST[PES]['RCT_NAMES_LST']
            PES_PRD_NAMES_LST = PES_LST[PES]['PRD_NAMES_LST']
            PES_RXN_NAME_LST = PES_LST[PES]['RXN_NAME_LST']
            if isinstance(CHANNELS, str):
                if CHANNELS == 'all':
                    print(len(PES_RXN_NAME_LST))
                    pes_chns = numpy.arange(len(PES_RXN_NAME_LST)+1)
                elif '-' in CHANNELS:
                    start, end = CHANNELS.split('-')
                    pes_chns = numpy.arange(int(start), int(end)+1)
                elif '[' in CHANNELS:
                    nums = CHANNELS.replace('[','').replace(']','').split(',')
                    pes_chns = [int(num) for num in nums]
            RCT_NAMES_LST = []
            PRD_NAMES_LST = []
            RXN_NAME_LST = []
            print('running ktp on PES{}: {} for the following channels...'.format(str(pes_idx), PES))
            for chn_idx, _ in enumerate(PES_RXN_NAME_LST):
                if chn_idx+1 in pes_chns:
                    RCT_NAMES_LST.append(PES_RCT_NAMES_LST[chn_idx])
                    PRD_NAMES_LST.append(PES_PRD_NAMES_LST[chn_idx])
                    RXN_NAME_LST.append(PES_RXN_NAME_LST[chn_idx])
                print('{} = {}'.format(
                    '+ '.join(PES_RCT_NAMES_LST[chn_idx]),
                    '+ '.join(PES_PRD_NAMES_LST[chn_idx])))
            # set up name lists and ts species dictionary for a given PES
           # print('RCT_NAMES_LST test:', RCT_NAMES_LST)
           # print('PRD_NAMES_LST test:', PRD_NAMES_LST)
           # print('RXT_NAME_LST test:', RXN_NAME_LST)
            RXN_LST = []
            for rxn, _ in enumerate(RCT_NAMES_LST):
                RXN_LST.append(
                    {'species': [], 'reacs': list(RCT_NAMES_LST[rxn]), 'prods':
                     list(PRD_NAMES_LST[rxn])})
            ts_idx = 0
            for idx, rxn in enumerate(RXN_LST):
                reacs = rxn['reacs']
                prods = rxn['prods']
                tsname = 'ts_{:g}'.format(ts_idx)
                SPC_DCT[tsname] = {}
                if reacs and prods:
                    SPC_DCT[tsname] = {'reacs': reacs, 'prods': prods}
                SPC_DCT[tsname]['ich'] = ''
                ts_chg = 0
                for rct in RCT_NAMES_LST[idx]:
                    ts_chg += SPC_DCT[rct]['chg']
                SPC_DCT[tsname]['chg'] = ts_chg
                ts_mul_low, ts_mul_high, rad_rad = moldr.util.ts_mul_from_reaction_muls(
                    RCT_NAMES_LST[idx], PRD_NAMES_LST[idx], SPC_DCT)
                SPC_DCT[tsname]['mul'] = ts_mul_low
                SPC_DCT[tsname]['rad_rad'] = rad_rad
                ts_idx += 1
                if ts_mul_low != ts_mul_high and rad_rad:
                    spc_dct = SPC_DCT[tsname].copy()
                    tsname = 'ts_{:g}'.format(ts_idx)
                    SPC_DCT[tsname] = spc_dct
                    SPC_DCT[tsname]['mul'] = ts_mul_high
                    ts_idx += 1
                    
            print('RUNNING WITH MESS')
            # run ktp for a given PES
            ktpdriver.driver.run(
                TSK_INFO_LST, ES_DCT, SPC_DCT, RCT_NAMES_LST, PRD_NAMES_LST,
                '/lcrc/project/PACC/run', '/lcrc/project/PACC/save', options=OPTIONS)

# f. Partition function parameters determined internally
# TORS_MODEL can take values: 'RIGID', '1DHR', or 'TAU' and eventually 'MDHR'
# VIB_MODEL can take values: 'HARM', or 'VPT2' values.

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
