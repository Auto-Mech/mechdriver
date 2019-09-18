""" drivers for thermochemistry evaluations
"""
import os
from qcelemental import constants as qcc
import automol.inchi
import automol.geom
import scripts.es
import thermo.heatform
import esdriver.driver
import autofile.fs


def run(tsk_info_lst, es_dct, spcdct, rct_names_lst, prd_names_lst, run_prefix, save_prefix, options=[True, True, True, False]):
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
        rxn_spc = rct_names_lst[rxn]
        rxn_spc.extend(prd_names_lst[rxn])
        for spc in rxn_spc:
            if spc not in spc_queue:
                spc_queue.append(spc)

    #Add reference molecules
    for spc in spc_queue:
        if not 'ich' in spcdct[spc]:
            spcdct[spc]['ich'] = automol.smiles.inchi(spcdct[spc]['smi'])

    #prepare filesystem
    if not os.path.exists(save_prefix):
        os.makedirs(save_prefix)
    if not os.path.exists(run_prefix):
        os.makedirs(run_prefix)

    #Run ESDriver
    if runes:
        rxn_lst = []
        for rxn, _ in enumerate(rct_names_lst):
            rxn_lst.append({'species': [], 'reactants': rct_names_lst[rxn], 'products': prd_names_lst[rxn]})
        if runspcfirst:
            rxn_lst[0]['species'] = spc_queue
        esdriver.driver.run(tsk_info_lst, es_dct, rxn_lst, spcdct, run_prefix, save_prefix)

    geo_lvl = ''
    harm_lvl = ''
    anharm_lvl = ''
    tors_lvl = ''

if __name__ == "__main__":

    MSG = """
           ================================================================
           ==                          AUTOMECHANIC                      ==
           ===         Andreas Copan, Sarah Elliott, Kevin Moore,       ===
           ===  Carlo Cavallotti, Yuri Georgievski, Murat Keceli,       === 
           ==                   Stephen Klippenstein                     ==
           ================================================================
           ==                          KTPDRIVER                      ==
           ===         Sarah Elliott, Kevin Moore, Andreas Copan,       ===
           ==   Carlo Cavallotti, Yuri Georgievski, Stephen Klippenstein ==
           ================================================================\n"""
    print(MSG)
    #tsk_info_lst, es_dct, spcs = load_params()
    REF = 'cbh0'

    TSK_INFO_LST = [
        ['conf_samp', 'mclev', 'mclev', False],
        ['find_ts', 'mclev', 'mclev', False],
        ['conf_samp', 'mclev', 'mclev', False],
        ['geom', 'optlev', 'optlev', False],
        ['conf_hess', 'optlev', 'optlev', False]
        # [ 'hr', 'hrlev', 'optlev', False]
        # [ 'sp', 'splev', 'optlev', False]]
        ]

    ES_DCT = {'mclev': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp', 'basis': '6-31g*',
        'ncycles': 60, 'mem': 32, 'nprocs': 8, 'econv': '1.e-8', 'gconv':
        '1.e-4'},
#   es_dct = {'mclev': {
#       'orb_res': 'RU', 'program': 'psi4', 'method': 'b3lyp', 'basis': '6-31g*',
#       'ncycles': 60, 'mem': 32, 'nprocs': 8, 'econv': '1.e-8', 'gconv': '1.e-4'},
#             'optlev': {
#                 'orb_res': 'RU', 'program': 'psi4', 'method': 'b3lyp',
#                 'basis': 'cc-pvdz', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
#                 'econv': '1.e-8', 'gconv': '1.e-4'},
#             'testlvl': {
#                 'orb_res': 'RU', 'program': 'psi4', 'method': 'b3lyp',
#                 'basis': 'cc-pvdz', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
#                 'econv': '1.e-8', 'gconv': '1.e-4'},
              'optlev': {
                  'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp',
                  'basis': 'cc-pvdz', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
                  'econv': '1.e-8', 'gconv': '1.e-4'},
              'hrlev':  {
                  'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp',
                  'basis': '6-31g*', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
                  'econv': '1.e-8', 'gconv': '1.e-4'},
              'anlev':  {
                  'orb_res': 'RU', 'program': 'psi4', 'method': 'b3lyp',
                  'basis': 'cc-pvdz', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
                  'econv': '1.e-8', 'gconv': '1.e-4'},
              '2':      {
                  'orb_res': 'RU', 'program': 'molpro', 'method': 'ccsd(t)',
                  'basis': 'cc-pvtz', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
                  'econv': '1.e-8', 'gconv': '1.e-4'},
              'splev':  {
                  'orb_res': 'RU', 'program': 'molpro', 'method': 'b3lyp',
                  'basis': 'cc-pvtz', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
                  'econv': '1.e-8', 'gconv': '1.e-4'}}


    # SPCS = ['prod1']
    RCT_NAME_LST = [['methane', 'oh']]
    PRD_NAME_LST = [['methyl','water']]

#    SPCS = ['reac1', 'prod1', 'methyl']
    SPCDCT = {
            'methane': {'smi': 'C', 'mul': 1, 'chg': 0},
            'oh': {'smi': '[OH]', 'mul': 2, 'chg': 0},
            'methyl': {'smi': '[C]', 'mul': 1, 'chg': 0},
            'water': {'smi': 'O', 'mul': 1, 'chg': 0}
             }
    run(TSK_INFO_LST, ES_DCT, SPCDCT, RCT_NAME_LST, PRD_NAME_LST, '/lcrc/project/PACC/run', '/lcrc/project/PACC/save')
    #run(tsk_info_lst, es_dct, spcdct, spcs, ref, 'runtest', 'savetest')
