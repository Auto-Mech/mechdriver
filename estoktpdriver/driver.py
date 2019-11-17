""" reaction list test
"""
import os
import sys
import collections
import json
import numpy
import chemkin_io
import automol
from automol import formula
import moldr
import thermodriver
import ktpdriver
from datalibs import phycon, eleclvl, symm
from submission import read_dat

# Set mechanism and type to be read based on user input
DATA_PATH = sys.argv[1]
MECHANISM_NAME = sys.argv[2]
MECH_TYPE = sys.argv[3]
MECH_PATH = os.path.join(DATA_PATH, 'data', MECHANISM_NAME)
MECH_FILE = 'mech.json'
# Set path to the mechanims data files
#DATA_PATH = '/home/sjklipp/PACC/mech_test'
#DATA_PATH = '/lcrc/project/CMRP/pacc/automech'

# Obtain DCT containing geometries to use as input
GEOM_PATH = os.path.join(DATA_PATH, 'data', 'geoms')
GEOM_DCT = moldr.util.geometry_dictionary(GEOM_PATH)

# Read parameters that dictate job running options
PARAMS = read_dat.params(os.path.join(MECH_PATH, 'params.dat'))

# Set further parameters for what reactions and PESs to be run
if len(sys.argv) > 4:
    PARAMS.PESNUMS = sys.argv[4]
    if len(sys.argv) > 5:
        PARAMS.CHANNELS = sys.argv[5]
print('PESNUMS and PARAMS.CHANNELS:', PARAMS.PESNUMS, PARAMS.CHANNELS)

# Create the AutoMech run and save directory file systems
if not os.path.exists(PARAMS.RUN_PREFIX):
    os.mkdir(PARAMS.RUN_PREFIX)
if not os.path.exists(PARAMS.SAVE_PREFIX):
    os.mkdir(PARAMS.SAVE_PREFIX)
# Run EStokTPDriver for ChemKin inputs
if MECH_TYPE == 'CHEMKIN':

    # Process species data from the mechanism file
    # Also add in basis set species
    MECH_STR = open(os.path.join(MECH_PATH, 'mechanism.txt')).read()
    SPC_STR = open(os.path.join(MECH_PATH, 'species.csv')).read()
    SMI_DCT = chemkin_io.parser.mechanism.spc_name_dct(SPC_STR, 'smiles')
    ICH_DCT = chemkin_io.parser.mechanism.spc_name_dct(SPC_STR, 'inchi')
    MUL_DCT = chemkin_io.parser.mechanism.spc_name_dct(SPC_STR, 'mult')
    CHG_DCT = chemkin_io.parser.mechanism.spc_name_dct(SPC_STR, 'charge')
    SENS_DCT = chemkin_io.parser.mechanism.spc_name_dct(SPC_STR, 'sens')

    print('CHECK_STEREO TEST:', PARAMS.CHECK_STEREO, type(PARAMS.CHECK_STEREO))
    if PARAMS.CHECK_STEREO:
        SPC_STR = 'name,SMILES,InChI,mult,charge,sens \n'
        for name in ICH_DCT:
            ich = ICH_DCT[name]
            smi = SMI_DCT[name]
            mul = MUL_DCT[name]
            chg = CHG_DCT[name]
            sens = SENS_DCT[name]
            if not automol.inchi.is_complete(ich):
                print('adding stereochemistry for {0}, {1}, {2}'.format(
                    name, smi, ich))
                # this returns a list of ichs w/ different possible stereo vals
                # for now just taking the first of these
                ich = automol.inchi.add_stereo(ich)
                print('new ich possibilities:', ich)
                ich = ich[-1]
                print('new ich:', ich)
                ICH_DCT[name] = ich
            SPC_STR += '{0},\'{1}\',\'{2}\',{3},{4},{5} \n'.format(
                name, smi, ich, mul, chg, sens)

        STEREO_PATH = os.path.join(MECH_PATH, 'species_stereo.csv')
        with open(STEREO_PATH, 'w') as stereo_csv_file:
            stereo_csv_file.write(SPC_STR)

    SPC_NAMES = []
    SPC_DCT = {}
    for name in MUL_DCT:
        SPC_DCT[name] = {}
        SPC_DCT[name]['smi'] = SMI_DCT[name]
        SPC_DCT[name]['ich'] = ICH_DCT[name]
        SPC_DCT[name]['chg'] = CHG_DCT[name]
        SPC_DCT[name]['mul'] = MUL_DCT[name]
        SPC_NAMES.append(name)

    RXN_BLOCK_STR = chemkin_io.parser.mechanism.reaction_block(MECH_STR)
    RXN_STRS = chemkin_io.parser.reaction.data_strings(RXN_BLOCK_STR)
    RCT_NAMES_LST = list(
        map(chemkin_io.parser.reaction.reactant_names, RXN_STRS))
    PRD_NAMES_LST = list(
        map(chemkin_io.parser.reaction.product_names, RXN_STRS))

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
        # print('formula test during append:', FORMULA_STR_LST)

    # print('formula test before sort:', FORMULA_STR_LST)
    RXN_INFO_LST = list(
        zip(FORMULA_STR_LST, RCT_NAMES_LST, PRD_NAMES_LST, RXN_NAME_LST))
    #PARAMS.SORT_RXNS = False
    print('SORT_RXNS test:', PARAMS.SORT_RXNS)
    if PARAMS.SORT_RXNS:
        RXN_INFO_LST.sort()
        FORMULA_STR_LST, RCT_NAMES_LST, PRD_NAMES_LST, RXN_NAME_LST = zip(*RXN_INFO_LST)

# Run EStokTPDriver for JSON inputs
    print('MECH_TYPE test:', MECH_TYPE)
    #print('RXN_INFO_LST test:', RXN_INFO_LST)
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

    # Convert essential pieces of json file to chemkin formatted data so
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
    # idxp = 0
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
                ich = rct['InChi']
                if PARAMS.CHECK_STEREO:
                    if not automol.inchi.is_complete(ich):
                        print('adding stereochemsiry for {}'.format(ich))
                        ich = automol.inchi.add_stereo(rct['InChi'])[0]
                rct_ichs.append(ich)
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
                ich = prd['InChi']
                if PARAMS.CHECK_STEREO:
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

            print('RAD_RAD_SORT test:', PARAMS.RAD_RAD_SORT)
            if PARAMS.RAD_RAD_SORT and not rad_rad_reac and not rad_rad_prod:
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

        formula_dct = ''
        for rct_ich in rct_ichs:
            formula_i_dct = automol.inchi.formula_dct(rct_ich)
            formula_dct = automol.formula._formula.join(formula_dct, formula_i_dct)
        FORMULA_STR = automol.formula._formula.string(formula_dct)
        FORMULA_STR_LST.append(FORMULA_STR)

    #print('FORMULA TEST:', FORMULA_STR_LST)
    print('RCT_MULS_LST:', RCT_MULS_LST)
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

    SORT_SMILES_PATH = os.path.join(MECH_PATH, 'smiles_sort.csv')
    with open(SORT_SMILES_PATH, 'w') as sorted_csv_file:
        sorted_csv_file.write(csv_str)

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
    CHG_DCT = {}
    MUL_DCT = {}
    SPC_DCT = {}
    TOT_LST = RCT_NAMES_LST + PRD_NAMES_LST
    for i, spc_names_lst in enumerate(TOT_LST):
        for j, spc_name in enumerate(spc_names_lst):
            chg = 0
            if spc_name not in SPC_NAMES:
                SPC_NAMES.append(spc_name)
                CHG_DCT[spc_name] = chg
                print('RCT_MULS TEST1:', spc_name, i, j)
                print('RCT_MULS TEST2:', RCT_MULS_LST[i])
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
                CHG_DCT[spc_name] = chg
                MUL_DCT[spc_name] = PRD_MULS_LST[i][j]
                SPC_DCT[spc_name] = {}
                SPC_DCT[spc_name]['chg'] = chg
                SPC_DCT[spc_name]['ich'] = PRD_ICHS_LST[i][j]
                SPC_DCT[spc_name]['mul'] = PRD_MULS_LST[i][j]
    RXN_INFO_LST = list(zip(FORMULA_STR_LST, RCT_NAMES_LST, PRD_NAMES_LST, RXN_NAME_LST))
    #print('RXN_INFO_LST TEST:', RXN_INFO_LST)

PES_LST = {}
CURRENT_FORMULA = ''
for fidx, formula in enumerate(FORMULA_STR_LST):
    #print('fidx test:', fidx, formula)
    if CURRENT_FORMULA == formula:
        PES_LST[formula]['RCT_NAMES_LST'].append(RCT_NAMES_LST[fidx])
        PES_LST[formula]['PRD_NAMES_LST'].append(PRD_NAMES_LST[fidx])
        PES_LST[formula]['RXN_NAME_LST'].append(RXN_NAME_LST[fidx])
        #print('old formula:')
        #print('NAMES_TST:', RCT_NAMES_LST[fidx], PRD_NAMES_LST[fidx], RXN_NAME_LST[fidx])
    else:
        CURRENT_FORMULA = formula
        PES_LST[formula] = {}
        PES_LST[formula]['RCT_NAMES_LST'] = [RCT_NAMES_LST[fidx]]
        PES_LST[formula]['PRD_NAMES_LST'] = [PRD_NAMES_LST[fidx]]
        PES_LST[formula]['RXN_NAME_LST'] = [RXN_NAME_LST[fidx]]
        #print('new formula:')
        #print('NAMES_TST:', RCT_NAMES_LST[fidx], PRD_NAMES_LST[fidx], RXN_NAME_LST[fidx])

for spc in SPC_DCT:
    ich = SPC_DCT[spc]['ich']
    mul = SPC_DCT[spc]['mul']
    if (ich, mul) in eleclvl.DCT:
        SPC_DCT[spc]['elec_levs'] = eleclvl.DCT[(ich, mul)]
    if (ich, mul) in symm.DCT:
        SPC_DCT[spc]['sym'] = symm.DCT[(ich, mul)]
    if ich in GEOM_DCT:
        SPC_DCT[spc]['geo_obj'] = GEOM_DCT[ich]
    SPC_DCT[spc]['hind_inc'] = PARAMS.HIND_INC * phycon.DEG2RAD

# Electronic structure parameters; code, method, basis, convergence control
ES_DCT = {
    'lvl_wbs': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'wb97xd',
        'basis': '6-31g*',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_wbm': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'wb97xd',
        'basis': '6-31+g*',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_wbt': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'wb97xd',
        'basis': 'cc-pvtz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_m06s': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'm062x',
        'basis': '6-31g*',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_m06m': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'm062x',
        'basis': '6-31+g*',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_m06t': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'm062x',
        'basis': 'cc-pvtz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_b2d': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b2plypd3',
        'basis': 'cc-pvdz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_b2t': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b2plypd3',
        'basis': 'cc-pvtz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_b2q': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b2plypd3',
        'basis': 'cc-pvqz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_b3s': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp',
        'basis': '6-31g*',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_b3t': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp',
        'basis': '6-31g*',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'cc_lvl_d': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)',
        'basis': 'cc-pvdz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'cc_lvl_t': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)',
        'basis': 'cc-pvtz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'cc_lvl_q': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)',
        'basis': 'cc-pvqz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'cc_lvl_df': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)-f12',
        'basis': 'cc-pvdz-f12',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'cc_lvl_tf': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)-f12',
        'basis': 'cc-pvtz-f12',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'cc_lvl_qf': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)-f12',
        'basis': 'cc-pvqz-f12',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'mlvl_cas_dz': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'caspt2',
        'basis': 'cc-pvdz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'mlvl_cas_tz': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'caspt2',
        'basis': 'cc-pvtz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        }
}

# The logic key in tsk_info_lst is for overwrite

if PARAMS.RUN_THERMO:
    SPC_QUEUE = list(SPC_NAMES)
    # thermodriver.driver.run(
    #     TSK_INFO_LST, ES_DCT, SPC_DCT, SPC_QUEUE, REF_MOLS, RUN_PREFIX,
    #     SAVE_PREFIX, options=OPTIONS)
    thermodriver.driver.run(
        PARAMS.TSK_INFO_LST, ES_DCT, SPC_DCT, SPC_QUEUE, PARAMS.REF_MOLS,
        PARAMS.RUN_PREFIX, PARAMS.SAVE_PREFIX,
        ene_coeff=PARAMS.ENE_COEFF, options=PARAMS.OPTIONS_THERMO)

if PARAMS.RUN_RATES:

    # print all the channels for all the PESs
    for pes_idx, PES in enumerate(PES_LST, start=1):
        print('PES test:', pes_idx, PES)
        PES_RXN_NAME_LST = PES_LST[PES]['RXN_NAME_LST']
        PES_RCT_NAMES_LST = PES_LST[PES]['RCT_NAMES_LST']
        PES_PRD_NAMES_LST = PES_LST[PES]['PRD_NAMES_LST']
        for chn_idx, _ in enumerate(PES_RXN_NAME_LST):
            print('channel {}: {} = {}'.format(
                chn_idx+1, ' + '.join(PES_RCT_NAMES_LST[chn_idx]),
                ' + '.join(PES_PRD_NAMES_LST[chn_idx])))

    if isinstance(PARAMS.PESNUMS, str):
        if PARAMS.PESNUMS == 'all':
            PARAMS.PESNUMS = numpy.arange(len(PES_LST)+1)
        elif '-' in PARAMS.PESNUMS:
            start, end = PARAMS.PESNUMS.split('-')
            PARAMS.PESNUMS = numpy.arange(int(start), int(end)+1)
        elif '[' in PARAMS.PESNUMS:
            nums = PARAMS.PESNUMS.replace('[', '').replace(']', '').split(',')
            PARAMS.PESNUMS = [int(num) for num in nums]
    # loop over PESs
    print('PARAMS.PESNUMS and PARAMS.CHANNELS:', PARAMS.PESNUMS, PARAMS.CHANNELS)
    for pes_idx, PES in enumerate(PES_LST, start=1):
        if pes_idx in PARAMS.PESNUMS:
            PES_RCT_NAMES_LST = PES_LST[PES]['RCT_NAMES_LST']
            PES_PRD_NAMES_LST = PES_LST[PES]['PRD_NAMES_LST']
            PES_RXN_NAME_LST = PES_LST[PES]['RXN_NAME_LST']
            if isinstance(PARAMS.CHANNELS, str):
                if PARAMS.CHANNELS == 'all':
                    print(len(PES_RXN_NAME_LST))
                    pes_chns = numpy.arange(len(PES_RXN_NAME_LST)+1)
                elif '-' in PARAMS.CHANNELS:
                    start, end = PARAMS.CHANNELS.split('-')
                    pes_chns = numpy.arange(int(start), int(end)+1)
                elif '[' in PARAMS.CHANNELS:
                    nums = PARAMS.CHANNELS.replace('[', '').replace(']', '').split(',')
                    pes_chns = [int(num) for num in nums]
            print('for pes:', pes_idx)
            RCT_NAMES_LST = []
            PRD_NAMES_LST = []
            RXN_NAME_LST = []

            # Split up sub-pes within a formula
            subpes_idx = 0
            conndct = {}
            connchnls = {}
            for chnl_idx, chnl_name in enumerate(PES_RXN_NAME_LST):
                #print('Channel name', chnl_name)
                connected_to = []
                chnl_species = [list(PES_RCT_NAMES_LST[chnl_idx]), list(PES_PRD_NAMES_LST[chnl_idx])]
                for conn_chnls_idx in conndct:
                    for spc_pair in chnl_species:
                        if len(spc_pair) == 1:
                            if spc_pair in conndct[conn_chnls_idx]:
                                if conn_chnls_idx not in connected_to:
                                    connected_to.append(conn_chnls_idx)
                            elif spc_pair[::-1] in conndct[conn_chnls_idx]:
                                if conn_chnls_idx not in connected_to:
                                    connected_to.append(conn_chnls_idx)
                if not connected_to:
                    conndct[subpes_idx] = chnl_species
                    connchnls[subpes_idx] = [chnl_idx]
                    subpes_idx += 1
                else:
                    conndct[connected_to[0]].extend(chnl_species)
                    connchnls[connected_to[0]].append(chnl_idx)
                    if len(connected_to) > 1:
                        for cidx, cval in enumerate(connected_to):
                            if cidx > 0:
                                conn_specs = conndct.pop(cval, None)
                                conn_chnls = connchnls.pop(cval, None)
                                conndct[connected_to[0]].extend(conn_specs)
                                connchnls[connected_to[0]].extend(conn_chnls)
                    for cidx in conndct:
                        conndct[cidx].sort()
                        conndct[cidx] = [conndct[cidx][i] for i in
                                         range(len(conndct[cidx])) if i == 0 or
                                         conndct[cidx][i] != conndct[cidx][i-1]]

            for cidx, cvals in enumerate(connchnls.values()):
                print('PES{}_{}: {} surface'.format(str(pes_idx), str(cidx+1), PES))
                for cval in cvals:
                    print(PES_RXN_NAME_LST[cval])


            # Loop ktp runs over the sub-pes
            for cidx, cvals in enumerate(connchnls.values()):
                print('ktp on PES{}_{}: {} for the following channels...'.format(str(pes_idx), str(cidx+1), PES))
                run_pes = False
                RCT_NAMES_LST = []
                PRD_NAMES_LST = []
                RXN_NAME_LST = []
                for chn_idx, _ in enumerate(PES_RXN_NAME_LST):
                    if chn_idx+1 in pes_chns and chn_idx in cvals:
                        run_pes = True
                        RCT_NAMES_LST.append(PES_RCT_NAMES_LST[chn_idx])
                        PRD_NAMES_LST.append(PES_PRD_NAMES_LST[chn_idx])
                        RXN_NAME_LST.append(PES_RXN_NAME_LST[chn_idx])
                        print('running channel {}: {} = {}'.format(
                            str(chn_idx+1),
                            ' + '.join(PES_RCT_NAMES_LST[chn_idx]),
                            ' + '.join(PES_PRD_NAMES_LST[chn_idx])))
                if run_pes:
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
                        SPC_DCT[tsname]['hind_inc'] = PARAMS.HIND_INC * phycon.DEG2RAD
                        ts_idx += 1
                        # if ts_mul_low != ts_mul_high and rad_rad:
                            #spc_dct = SPC_DCT[tsname].copy()
                            #tsname = 'ts_{:g}'.format(ts_idx)
                            #SPC_DCT[tsname] = spc_dct
                            #SPC_DCT[tsname]['mul'] = ts_mul_high
                            #ts_idx += 1

                    # Build the list for the energy transfer section to pass to ktpdriver
                    etrans_lst = [
                        PARAMS.EXP_FACTOR,
                        PARAMS.EXP_POWER,
                        PARAMS.EXP_CUTOFF,
                        PARAMS.EPS1,
                        PARAMS.EPS2,
                        PARAMS.SIG1,
                        PARAMS.SIG2,
                        PARAMS.MASS1
                    ]


                    print('RUNNING WITH MESS')
                    # run ktp for a given PES
                    ktpdriver.driver.run(
                        PARAMS.TSK_INFO_LST, ES_DCT, SPC_DCT, RCT_NAMES_LST, PRD_NAMES_LST,
                        PARAMS.RUN_PREFIX, PARAMS.SAVE_PREFIX,
                        ene_coeff=PARAMS.ENE_COEFF, options=PARAMS.OPTIONS_RATE, etrans=etrans_lst)
                        #'/lcrc/project/PACC/elliott/runhr', '/lcrc/project/PACC/elliott/savehr', options=OPTIONS)

# f. Partition function parameters determined internally
# TORS_MODEL can take values: 'RIGID', '1DHR', or 'TAU' and eventually 'MDHR'
# VIB_MODEL can take values: 'HARM', or 'VPT2' values.
