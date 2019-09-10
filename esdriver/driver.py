from esdriver.load import load_logger
from esdriver.load import load_params
import automol.inchi
import automol.geom
import scripts.es
import thermo.heatform

import os
from qcelemental import constants as qcc
import logging
log   = logging.getLogger(__name__)



def run(tsk_info_lst, es_dct, rxn_lst, spcdct, run_prefix, save_prefix):

    #prepare species queue 
    spc_queue  = []
    ts_dct     = {}
    for i, rxn in enumerate(rxn_lst): 
        reacs = rxn['reactants']
        prods = rxn[ 'products']

        spc_queue.extend(rxn['species'])
        spc_queue.extend(reacs)
        spc_queue.extend(prods)
        if reacs and prods:
            ts_dct['ts_{:g}'.format(i)] = {'reacs': reacs, 'prods': prods}

    spc_queue = list(dict.fromkeys(spc_queue)) #removes duplicates
    
       # #Print message
       # msg  = '||| Reaction {} will compute {} for '.format(rxn, ' and '.join(rxns[rxn]['job']))
       # if rxns[rxn]['products']:
       #     msg += 'the reaction: {} --> {}'.format(
       #              ' + '.join([automol.inchi.smiles(spcs[x]['inchi']) for x in rxns[rxn]['reactants']]), 
       #              ' + '.join([automol.inchi.smiles(spcs[x]['inchi']) for x in rxns[rxn]['products']]))
       # else:
       #     msg += 'the species: '.format(rxns[rxn]['reactants'])
       # msg += '\n||| using the following rxnules: {} '.format(', '.join(tsk_info_lst[rxn]))
       # logging.info(msg)
        
    #Prepare filesystem
    if not os.path.exists(save_prefix):
        os.makedirs(save_prefix)
    if not os.path.exists(run_prefix):
        os.makedirs(run_prefix)
    
    #Loop over Tasks  
    for tsk_info in tsk_info_lst:
    
        #Task information
        tsk          = tsk_info[0]
        es_ini_key = tsk_info[1]
        es_run_key   = tsk_info[2]
        overwrite    = tsk_info[3]
    
        #Theory information
        ini_es_info = get_es_info(es_dct, es_ini_key)
        run_es_info = get_es_info(es_dct, es_run_key)
           
        #If task is to find the transition state, find all TSs for your reactionlist
        if tsk == 'tsfind':
            for ts in enumerate(ts_dct):
                log.info('  | Task {} \t\t\t'.format(tsk))
                geo, zma = scripts.es.find_ts(run_prefix, save_prefix, ts['reacs'],  ts['prods'], spcsdct, run_es_info, overwrite)
                if not isinstance(geo, str):
                    spcs[ts] = create_spec(geo)
                    spcs[ts]['zmatrix'] = zma
                    log.info('   | Success, transition state {} added to ES queue'.format(ts))
                    spc_queue.append(ts)
                    continue
    
        #Loop over all species
        for spc in spc_queue:
            msg = '\n  | Task {} \t\t\t Species {}: {}'.format(tsk, spc, automol.inchi.smiles(spcdct[spc]['inchi']))
            log.info(msg)
            
            #Get params
            spc_info = scripts.es.get_spc_info(spcdct[spc])
            
            #Run tasks
            if 'ts' in spc:
                #Check if the task has been completed at the requested running theory level
    
                #Every task starts with a geometry optimization at the running theory level
                scripts.es.ts_geo_ref(spc, spcdct, ts_dct[spc]['reacs'], ts_dct[spc]['prods'], ini_es_info, run_es_info, run_prefix, save_prefix, overwrite)
                #Run the requested task at the requested running theory level
                scripts.es.ts_run_task(tsk, spc, spcdct, ts_dct[spc]['reacs'], ts_dct[spc]['prods'], es_dct[run_es_key], ini_es_info, run_prefix, save_prefix, overwrite)
            else:
                #Check if the task has been completed at the requested running theory level
    
                #Every task starts with a geometry optimization at the running theory level
                geo, msg = scripts.es.geo_ref(      spcdct[spc], ini_es_info, run_es_info, run_prefix, save_prefix, overwrite)
                if not 'reference' in msg:   #newly computed geometry should be checked for imaginary frequencies
                    scripts.es.remove_imag(  spcdct[spc], es_dct[es_run_key], run_prefix, save_prefix, overwrite)
                #Run the requested task at the requested running theory level
                scripts.es.run_task(tsk, spcdct[spc], es_dct[es_run_key], ini_es_info, run_es_info, spc_info, run_prefix, save_prefix, overwrite)
    return
              

def create_spec(val, charge = 0, mc_nsamp = [True, 3, 1, 3, 100, 12], hind_inc = 360.):
    spec = {}
    if isinstance(val, str):
        ich = val
        geo = automol.ichi.geo(ich)
        zma = automol.geom.zmatrix(geo)
        spec['zmatrix'] = ich
    else:
        geo = val
        ich = automol.geom.inchi(geo)
    smi = automol.inchi.smiles(ich)
    mult = 1
    mult += smi.count('[')
    spec[   'inchi'] = ich
    spec[  'geoobj'] = geo
    spec[  'charge'] = charge
    spec[    'mult'] = mult
    spec['mc_nsamp'] = mc_nsamp
    spec['hind_inc'] = hind_inc * qcc.conversion_factor('degree', 'radian') 
    return spec

#def get_overwrite(array):
#    overwrite = False
#    if len(array) > 2: 
#         if array[2] == 'true':
#             overwrite = True
#    return overwrite 

def get_es_info(es_dct, key):
    if key == 'input':
        return  ['placeholder','placeholder','placeholder','placeholder']
    else:
        return scripts.es.get_thy_info(es_dct[key])

if __name__ == "__main__":

    #load_logger('output.dat')
    load_logger()
    msg = """
           ================================================================
           ==                            ESDRIVER                        ==
           ====        Andreas Copan, Kevin Moore, Sarah Elliott,      ====
           ==   Carlo Cavolotti, Yuri Georgievski, Stephen Klippenstein  ==
           ================================================================\n"""
    log.info(msg)
    
    tsk_info_lst, rxns, es_dct, spcs = load_params()

    run(tsk_info_lst, rxns, es_dct, spcs) 
