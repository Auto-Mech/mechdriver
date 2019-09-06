from rcdriver.load import load_logger
from rcdriver.load import load_params
import automol.inchi
import automol.geom
import scripts.es
import thermo.heatform

import os
from qcelemental import constants as qcc
import logging
log   = logging.getLogger(__name__)



def run(tsks, rxns, lvls, spcs):

    ######Loop over all of the rxnels specified (this type of thing should be part of top level driver)

    for rxn in rxns:
        #Species in this rxn
        species_queue = rxns[rxn]['reactants'] + rxns[rxn]['products'] 

        #Print message
        msg  = '||| Reaction {} will compute {} for '.format(rxn, ' and '.join(rxns[rxn]['job']))
        if rxns[rxn]['products']:
            msg += 'the reaction: {} --> {}'.format(
                     ' + '.join([automol.inchi.smiles(spcs[x]['inchi']) for x in rxns[rxn]['reactants']]), 
                     ' + '.join([automol.inchi.smiles(spcs[x]['inchi']) for x in rxns[rxn]['products']]))
        else:
            msg += 'the species: '.format(rxns[rxn]['reactants'])
        msg += '\n||| using the following rxnules: {} '.format(', '.join(tsks[rxn]))
        logging.info(msg)
        
        #Add reference molecules if thermo or rates calculations are requested
        if rxns[rxn]['job']:
            prepare_refs(rxns[rxn], spcs, species_queue)
            species_queue += rxns[rxn]['references']                  
            species_queue  = list(dict.fromkeys(species_queue))
  
        #Begin tasks  
        for tsk in tsks[rxn]:

            #prepare filesystem
            run_prefix, save_prefix = rxns[rxn]['paths']
            if not os.path.exists(save_prefix):
                os.makedirs(save_prefix)
            if not os.path.exists(run_prefix):
                os.makedirs(run_prefix)

            #Theory information
            overwrite        = get_overwrite(tsks[rxn][tsk])
            initial_thy_info = get_initial_thy_info_(lvls, tsks[rxn][tsk][1])
            running_thy_info = scripts.es.get_thy_info(lvls[tsks[rxn][tsk][0]])
               
            #Find TS at whatever point its supposed to 
            if tsk == 'tsfind':
                log.info('  | Task {} \t\t\t'.format(tsk))
                geo, zma = scripts.es.find_ts(run_prefix, save_prefix, rxns[rxn]['reactants'],  rxns[rxn]['products'], spcs, running_thy_info, overwrite)
                if not isinstance(geo, str):
                    spcs['ts_' + rxn] = create_spec(geo)
                    spcs['ts_' + rxn]['zmatrix'] = zma
                    log.info('   | Success, ts_{} added to queue'.format(tsk))
                    species_queue.append('ts_' + rxn)
                continue

            for spc in species_queue:
                msg = '\n  | Task {} \t\t\t Species {}: {}'.format(tsk, spc, automol.inchi.smiles(spcs[spc]['inchi']))
                log.info(msg)
                
                #Get params
                spc_info = scripts.es.get_spc_info(spcs[spc])
                
                #Run tasks
                if not 'ts' in spc:
                    scripts.es.geo_ref(spcs[spc], initial_thy_info, running_thy_info, run_prefix, save_prefix, overwrite)
                    scripts.es.remove_imag(spcs[spc], lvls[tsks[rxn][tsk][0]], run_prefix, save_prefix, overwrite)
                    scripts.es.run_task(tsk, spcs[spc], lvls[tsks[rxn][tsk][0]], initial_thy_info, run_prefix, save_prefix, overwrite)
                else:
                    scripts.es.ts_geo_ref(spc, spcs, rxns[rxn]['reactants'],  rxns[rxn]['products'], initial_thy_info, running_thy_info, run_prefix, save_prefix, overwrite)
                    scripts.es.ts_run_task(tsk, spc, spcs, rxns[rxn]['reactants'],  rxns[rxn]['products'], lvls[tsks[rxn][tsk][0]], initial_thy_info, run_prefix, save_prefix, overwrite)

    return
              

#####THESE FUNCTIONS BECOME A PART OF THERMO DRIVER
def is_scheme(entry):
    calls = {"basic": "get_basis",
             "cbh0":  "get_cbhzed",
             "cbh1":  "get_cbhone",
             "cbh2":  "get_cbhtwo"}
    return entry in calls.keys()

def get_ref(species, spcs, scheme):
    calls = {"basic": "get_basis",
             "cbh0":  "get_cbhzed",
             "cbh1":  "get_cbhone",
             "cbh2":  "get_cbhtwo"}
    call  = getattr(thermo.heatform, calls[scheme])
    ref   = []
    for spec in species:
        ref.extend(call(spcs[spec]['inchi']))
    return list(dict.fromkeys(ref))

def prepare_refs(rxn, spcs, species):
     if not rxn['references']:
         rxns['references'] =  'basic'
     if is_scheme(rxn['references'][0]):
         refs = get_ref(species, spcs, rxn['references'][0])
         rxn['references'] = []
         for ref in refs:
             needtoadd = True
             for spc in spcs:
                 if spcs[spc]['inchi'] == ref:
                     rxn['references'].append(spc)
                     needtoadd = False
                     break
             if needtoadd:
                 msg = ' || Adding reference species ref_{}'.format(ref)
                 logging.info(msg)
                 spcs['ref_' + ref] = create_spec(ref)
                 rxns[rxn]['references'].append('ref_' + ref)
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
    mult = 0
    mult += smi.count('[')
    spec[   'inchi'] = ich
    spec[  'geoobj'] = geo
    spec[  'charge'] = charge
    spec[    'mult'] = mult
    spec['mc_nsamp'] = mc_nsamp
    spec['hind_inc'] = hind_inc * qcc.conversion_factor('degree', 'radian') 
    return spec

def get_overwrite(array):
    overwrite = False
    if len(array) > 2: 
         if array[2] == 'true':
             overwrite = True
    return overwrite 

def get_initial_thy_info_(lvls, key):
    if key == 'input':
        return  ['placeholder','placeholder','placeholder','placeholder']
    else:
        return scripts.es.get_thy_info(lvls[key])

if __name__ == "__main__":

    load_logger('output.dat')
    #load_logger()
    msg = """
           ================================================================
           ==                          ES to k(T,P)                      ==
           ====        Andreas Copan, Kevin Moore, Sarah Elliott,      ====
           ==   Carlo Cavolotti, Yuri Georgievski, Stephen Klippenstein  ==
           ================================================================\n"""
    log.info(msg)
    
    tsks, rxns, lvls, spcs = load_params()

    run(tsks, rxns, lvls, spcs) 
