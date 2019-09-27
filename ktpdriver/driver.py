""" drivers for thermochemistry evaluations
"""
import os
from qcelemental import constants as qcc
import automol.inchi
import automol.geom
import chemkin_io
import scripts.es
import thermo.heatform
import esdriver.driver
import autofile.fs

TEMPS = [300., 500., 750., 1000., 1250., 1500., 1750., 2000.]
PRESS = [0.1, 1., 10., 100.]
EXP_FACTOR = 150.0
EXP_POWER = 0.85
EXP_CUTOFF = 15.
EPS1 = 100.0
EPS2 = 200.0
SIG1 = 6.
SIG2 = 6.
MASS1 = 15.0
PROJROT_SCRIPT_STR = ("#!/usr/bin/env bash\n"
                      "RPHt.exe >& /dev/null")

def run(
        tsk_info_lst, es_dct, spcdct, rct_names_lst, prd_names_lst, run_prefix,
        save_prefix, ene_coeff=[1.], vdw_params=[False, False, True],
        options=[True, True, True, False]):
    """ main driver for thermo run
    """

    #Determine options
    runes = options[0]  #run electronic structure theory (True/False)
    runspcfirst = options[1]
    runmess = options[2]  #run mess (True) / just make the mess input file (False)
    runrates = options[3]
    if not runmess:
        runrates = False

    spc_queue = []
    for rxn, _ in enumerate(rct_names_lst):
        rxn_spc = list(rct_names_lst[rxn])
        rxn_spc.extend(list(prd_names_lst[rxn]))
        for spc in rxn_spc:
            if spc not in spc_queue:
                spc_queue.append(spc)
    for spc in spc_queue:
        if not 'ich' in spcdct[spc]:
            spcdct[spc]['ich'] = automol.smiles.inchi(spcdct[spc]['smi'])

    rxn_lst = []
    for rxn, _ in enumerate(rct_names_lst):
        rxn_lst.append(
            {'species': [], 'reactants': list(rct_names_lst[rxn]), 'products':
             list(prd_names_lst[rxn])})

    #prepare filesystem
    if not os.path.exists(save_prefix):
        os.makedirs(save_prefix)
    if not os.path.exists(run_prefix):
        os.makedirs(run_prefix)

    #Run ESDriver
    if runes:
        if runspcfirst:
            rxn_lst[0]['species'] = spc_queue
        spc_dct = esdriver.driver.run(
            tsk_info_lst, es_dct, rxn_lst, spcdct.copy(), run_prefix, save_prefix, vdw_params)
    else:
        pass
    #BUT if we don't run ES I need to construct the following info right here for ts dict
    #ts_zma, rxn_fs,
    tsname_0 = 'ts_0'
    if 'original_zma' in spc_dct[tsname_0]:
        pes_formula = automol.geom.formula(automol.zmatrix.geometry(spc_dct[tsname_0]['original_zma']))
        print('Starting mess file preparation for {}:'.format(pes_formula))

    #Figure out the model and theory levels for the MESS files
    geo_lvl = ''
    harm_lvl = ''
    anharm_lvl = ''
    tors_lvl = ''
    sym_lvl = ''
    harm_lvl_ref = ''
    anharm_lvl_ref = ''
    tors_lvl_ref = ''
    sym_lvl_ref = ''

    ts_model = ['RIGID', 'HARM', '']
    for tsk in tsk_info_lst:
        if 'samp' in tsk[0] or 'geom' in tsk[0]:
            geo_lvl = tsk[1]
            geom = True
        if 'grad' in tsk[0] or 'hess' in tsk[0]:
            harm_lvl = tsk[1]
            harm_lvl_ref = tsk[2]
            if 'grad' in tsk[0]:
                grad = True
            if 'hess' in tsk[0]:
                hess = True
            if not geom:
                ene_lvl = tsk[1]
                geo_lvl = tsk[1]
        if 'hr' in tsk[0] or 'tau' in tsk[0]:
            tors_lvl = tsk[1]
            tors_lvl_ref = tsk[2]
            if 'md' in tsk[0]:
                ts_model[0] = 'MDHR'
            if 'tau' in tsk[0]:
                ts_model[0] = 'TAU'
            else:
                ts_model[0] = '1DHR'
        if 'anharm' in tsk[0] or 'vpt2' in tsk[0]:
            anharm_lvl = tsk[1]
            anharm_lvl_ref = tsk[2]
            ts_model[1] = 'ANHARM'
            if not hess:
                geo_lvl = tsk[1]
        if 'sym' in tsk[0]:
            sym_lvl = tsk[1]
            sym_lvl_ref = tsk[2]
            if 'samp' in tsk[0]:
                ts_model[2] = 'SAMPLING'
            if '1DHR' in tsk[0]:
                ts_model[2] = '1DHR'

    geo_thy_info = get_thy_info(es_dct, geo_lvl)
    harm_thy_info = get_thy_info(es_dct, harm_lvl)
    tors_thy_info = get_thy_info(es_dct, tors_lvl)
    anharm_thy_info = get_thy_info(es_dct, anharm_lvl)
    sym_thy_info = get_thy_info(es_dct, sym_lvl)
    harm_ref_thy_info = get_thy_info(es_dct, harm_lvl_ref)
    tors_ref_thy_info = get_thy_info(es_dct, tors_lvl_ref)
    anharm_ref_thy_info = get_thy_info(es_dct, anharm_lvl_ref)
    sym_ref_thy_info = get_thy_info(es_dct, sym_lvl_ref)
    pf_levels = [harm_thy_info, tors_thy_info, anharm_thy_info, sym_thy_info]
    ref_levels = [
        harm_ref_thy_info, tors_ref_thy_info, anharm_ref_thy_info, sym_ref_thy_info]

    #Collect energies for zero points
    spc_save_fs = autofile.fs.species(save_prefix)
    ts_queue = []
    for spc in spcdct:   #have to make sure you get them for the TS too
        if 'ts_' in spc:
            ts_queue.append(spc)
            if 'radical radical' in spcdct[spc]['class']:
                print('skipping rate for radical radical reaction: {}'.format(spc))
                continue
    for spc in spc_queue +  ts_queue:
        spc_info = (spcdct[spc]['ich'], spcdct[spc]['chg'], spcdct[spc]['mul'])
        if 'ts_' in spc:
            spc_save_path = spcdct[spc]['rxn_fs'][3]
            saddle = True
            save_path = spc_save_path
        else:
            spc_save_fs.leaf.create(spc_info)
            spc_save_path = spc_save_fs.leaf.path(spc_info)
            saddle = False
            save_path = save_prefix
        zpe, _ = scripts.thermo.get_zpe(
            spc, spcdct[spc], spc_save_path, pf_levels, ts_model)
        spcdct[spc]['zpe'] = zpe
        ene_strl = []
        ene_lvl = ''
        ene_lvl_ref = ''
        ene_idx = 0
        spcdct[spc]['ene'] = 0.
        ene_str = '! energy level:'
        for tsk in tsk_info_lst:
            if 'ene' in tsk[0]:
                if ene_idx > len(ene_coeff)-1:
                    print('Warning - an insufficient energy coefficient list was provided')
                    break
                ene_lvl = tsk[1]
                ene_lvl_ref = tsk[2]
                ene_ref_thy_info = scripts.es.get_thy_info(es_dct[ene_lvl_ref])
                ene_thy_info = scripts.es.get_thy_info(es_dct[ene_lvl])
                ene_strl.append(' {:.2f} x {}{}/{}//{}{}/{}\n'.format(
                    ene_coeff[ene_idx], ene_thy_info[3], ene_thy_info[1], ene_thy_info[2],
                    ene_ref_thy_info[3], ene_ref_thy_info[1], ene_ref_thy_info[2]))
                ene = scripts.thermo.get_electronic_energy(
                    spc_info, ene_ref_thy_info, ene_thy_info, save_path, saddle)
                spcdct[spc]['ene'] += ene*ene_coeff[ene_idx]
                ene_idx += 1
    ene_str += '!               '.join(ene_strl)

    #Collect formula and header string for the PES
    tsname_0 = 'ts_0'
    pes_formula = automol.geom.formula(automol.zmatrix.geometry(spc_dct[tsname_0]['original_zma']))
    rct_ichs = spc_dct[tsname_0]['rxn_ichs'][0]
    header_str, energy_trans_str = scripts.ktp.pf_headers(
        rct_ichs, TEMPS, PRESS, EXP_FACTOR, EXP_POWER, EXP_CUTOFF, EPS1, EPS2,
        SIG1, SIG2, MASS1)

    mess_strs = ['', '', '']
    idx_dct = {}
    first_ground_ene = 0.
    wells = scripts.ktp.make_all_well_data(
        rxn_lst, spcdct.copy(), save_prefix, ts_model, pf_levels, PROJROT_SCRIPT_STR)
    for idx, rxn in enumerate(rxn_lst):
        tsname = 'ts_{:g}'.format(idx)
        if 'radical radical' in spc_dct[tsname]['class']:
            print('skipping rate for radical radical reaction: {}'.format(tsname))
            continue
        tsform = automol.geom.formula(automol.zmatrix.geometry(spc_dct[tsname]['original_zma']))
        if tsform != pes_formula:
            print('Reaction list contains reactions on different potential energy surfaces: {} and {}'.format(
                tsform, pes_formula))
            print('Will proceed to construct only {}'.format(pes_formula))
            continue
        mess_strs, first_ground_ene = scripts.ktp.make_channel_pfs(
            tsname, rxn, wells, spcdct.copy(), idx_dct, mess_strs, first_ground_ene)
        print(idx_dct)
    well_str, bim_str, ts_str = mess_strs
    ts_str += '\nEnd\n'
    print(well_str)
    print(bim_str)
    print(ts_str)

    # run mess to produce rate output

    mess_path = scripts.ktp.run_rate(
        header_str, energy_trans_str, well_str, bim_str, ts_str,
        spcdct[tsname_0], geo_thy_info, spcdct[tsname_0]['rxn_fs'][3])

    # fit rate output to modified Arrhenius forms and print in ChemKin format
    pf_levels.append(ene_str)
    chemkin_header_str = scripts.thermo.run_ckin_header(pf_levels, ref_levels, ts_model)
    chemkin_str = chemkin_header_str
    starting_path = os.getcwd()
    labels = idx_dct.values()
    names = idx_dct.keys()
    err_thresh = 15.
    for lab_i, name_i in zip(labels, names):
        for lab_j, name_j in zip(labels, names):
            ene = 0.
            if lab_i != lab_j:
                for spc in name_i.split('+'):
                    ene += scripts.thermo.spc_energy(spcdct[spc]['ene'], spcdct[spc]['zpe'])
                for spc in name_j.split('+'):
                    ene -= scripts.thermo.spc_energy(spcdct[spc]['ene'], spcdct[spc]['zpe'])
                if ene > 0.:
                    reaction = name_i + '=' + name_j
                    sing_rate_params, sing_errs, doub_rate_params, doub_errs = scripts.ktp.mod_arr_fit(
                        lab_i, lab_j, mess_path)
                    max_err = max([vals[1] for vals in sing_errs.values()])
                    print('max_err test:', max_err, err_thresh)
                    print('sing err test:', sing_errs)
                    print('doub err test:', doub_errs)
                    print('sing_rate_params:', sing_rate_params)
                    if max_err < err_thresh:
                        chemkin_str += chemkin_io.mechwriter.reaction.plog(
                            reaction, sing_rate_params, sing_errs)
                    else:
                        chemkin_str += chemkin_io.mechwriter.reaction.plog(
                            reaction, doub_rate_params, doub_errs)

    print(chemkin_str)
    with open(starting_path+'/rates.ckin', 'w') as f:
        f.write(chemkin_str)


def get_thy_info(es_dct, key):
    """ setup theory info file from es dictionary
    """
    ret = []
    if key:
        ret = scripts.es.get_thy_info(es_dct[key])
    return ret

if __name__ == "__main__":

    MSG = """
           ================================================================
           ==                        AUTOMECHANIC                        ==
           ===         Andreas Copan, Sarah Elliott, Kevin Moore,       ===
           ===     Daniel Moberg, Carlo Cavallotti, Yuri Georgievski,   ===
           ==       Ahren Jasper, Murat Keceli, Stephen Klippenstein     ==
           ================================================================
           ==                         KTPDRIVER                          ==
           ===         Sarah Elliott, Kevin Moore, Andreas Copan,       ===
           ===      Daniel Moberg, Carlo Cavallotti, Yuri Georgievski,  ===
           ==            Ahren Jasper, Stephen Klippenstein              ==
           ================================================================\n"""
    print(MSG)
    #tsk_info_lst, es_dct, spcs = load_params()
    REF = 'cbh0'
    VDW_PARAMS = [True, True, False]  #[for reac (T/F), for reac (T/F), from reactants(T)/from TS (F)]
    TSK_INFO_LST = [
        ['conf_samp', 'mclev', 'mclev', False],
        ['find_ts', 'mclev', 'mclev', False],
        #['find_vdw', 'mclev', 'mclev', False],
        ['conf_samp', 'mclev', 'mclev', False],
        ['find_geom', 'optlev', 'mclev', False],
        ['conf_hess', 'optlev', 'optlev', False],
        # ['geom', 'b2tz', 'optlev', False],
        # ['conf_hess', 'b2tz', 'b2tz', False],
        # ['hr_scan', 'cheap', 'optlev', False],
        # ['hr_scan', 'cheap', 'optlev', False],
        ['conf_energy', 'optlev', 'optlev', False]
        # [ 'hr', 'hrlev', 'optlev', False]
        # [ 'sp', 'splev', 'optlev', False]]
        ]

    ES_DCT = {'mclev': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp', 'basis': '6-31g*',
        'ncycles': 60, 'mem': 32, 'nprocs': 8, 'econv': '1.e-8', 'gconv':
        '1.e-4', 'mc_nsamp': [True, 3, 1, 3 , 100, 5]},
              'optlev': {
                  'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp',
                  'basis': 'cc-pvdz', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
                  'econv': '1.e-8', 'gconv': '1.e-4', 'mc_nsamp': [True, 3, 1, 3 , 100, 5]},
              'hrlev':  {
                  'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp',
                  'basis': '6-31g*', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
                  'econv': '1.e-8', 'gconv': '1.e-4', 'mc_nsamp': [True, 3, 1, 3 , 100, 5]},
              'anlev':  {
                  'orb_res': 'RU', 'program': 'psi4', 'method': 'b3lyp',
                  'basis': 'cc-pvdz', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
                  'econv': '1.e-8', 'gconv': '1.e-4', 'mc_nsamp': [True, 3, 1, 3 , 100, 5]},
              '2':      {
                  'orb_res': 'RU', 'program': 'molpro', 'method': 'ccsd(t)',
                  'basis': 'cc-pvtz', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
                  'econv': '1.e-8', 'gconv': '1.e-4', 'mc_nsamp': [True, 3, 1, 3 , 100, 5]},
              'splev':  {
                  'orb_res': 'RU', 'program': 'molpro', 'method': 'b3lyp',
                  'basis': 'cc-pvtz', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
                  'econv': '1.e-8', 'gconv': '1.e-4', 'mc_nsamp': [True, 3, 1, 3 , 100, 5]},
              'cheap': {'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp', 'basis': 'sto-3g',  'mc_nsamp': [True, 3, 1, 3 , 100, 5]},
              'b2tz': {'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b2plypd3', 'basis': 'cc-pvtz',  'mc_nsamp': [True, 3, 1, 3 , 100, 5]}
              }


    #RCT_NAME_LST = [['nh3', 'oh']]
    #PRD_NAME_LST = [['nh2','water']]
    #RCT_NAME_LST = [['ch3oh', 'h'], ['ch3oh', 'h']]
    #PRD_NAME_LST = [['ch2oh', 'h2'], ['ch3o', 'h2']]
    RCT_NAME_LST = [['ch3oh', 'h']]
    PRD_NAME_LST = [['ch3o', 'h2']]

    #RCT_NAME_LST = [['methane', 'h'], ['methane', 'oh']]
    #PRD_NAME_LST = [['methyl','h2'], ['methyl','water']]

    SPCDCT = {
            'methane': {'smi': 'C', 'mul': 1, 'chg': 0},
            'h': {'smi': '[H]', 'mul': 2, 'chg': 0},
            'h2': {'smi': '[H][H]', 'mul': 1, 'chg': 0},
            'oh': {'smi': '[OH]', 'mul': 2, 'chg': 0},
            'methyl': {'smi': '[CH3]', 'mul': 2, 'chg': 0},
            'nh3': {'smi': '[NH3]', 'mul': 1, 'chg': 0},
            'nh2': {'smi': '[NH2]', 'mul': 2, 'chg': 0},
            'water': {'smi': 'O', 'mul': 1, 'chg': 0},
            'ch3oh': {'smi': 'CO', 'mul': 1, 'chg': 0},
            'ch2oh': {'smi': '[CH2]O', 'mul': 2, 'chg': 0},
            'ch3o': {'smi': 'C[O]', 'mul': 2, 'chg': 0}
             }
    #run(TSK_INFO_LST, ES_DCT, SPCDCT, RCT_NAME_LST, PRD_NAME_LST, '/lcrc/project/PACC/run', '/lcrc/project/PACC/save')
    run(
        TSK_INFO_LST, ES_DCT, SPCDCT, RCT_NAME_LST, PRD_NAME_LST,
        '/lcrc/project/PACC/elliott/run2', '/lcrc/project/PACC/elliott/save2',
        vdw_params=VDW_PARAMS)
    #run(tsk_info_lst, es_dct, spcdct, spcs, ref, 'runtest', 'savetest')
