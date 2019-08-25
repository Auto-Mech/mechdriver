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

def prepare_refs(mods, spcs):
     if not mods[mod]['references']:
         mods[mod]['references'] =  'basic'
     if is_scheme(mods[mod]['references'][0]):
         refs = get_ref(species, spcs, mods[mod]['references'][0])
         mods[mod]['references'] = []
         for ref in refs:
             needtoadd = True
             for spc in spcs:
                 if spcs[spc]['inchi'] == ref:
                     mods[mod]['references'].append(spc)
                     needtoadd = False
                     break
             if needtoadd:
                 msg = ' || Adding reference species ref_{}'.format(ref)
                 logging.info(msg)
                 spcs['ref_' + ref] = create_spec(ref)
                 mods[mod]['references'].append('ref_' + ref)
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
        print(ich)
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
        return  ['fake','fake','fake','fake']
    else:
        return scripts.es.get_thy_info(lvls[key])

if __name__ == "__main__":

    load_logger('output.dat')
    msg = """
           ================================================================
           ==                ES to k(T,P) and RCDriver                   ==
           ====        Andreas Copan, Kevin Moore, Sarah Elliott,      ====
           ==   Carlo Cavolotti, Yuri Georgievski, Stephen Klippenstein  ==
           ================================================================\n"""
    log.info(msg)
    
    tsks, mods, lvls, spcs = load_params()
    
    for mod in mods:
        #Species in this model
        species = mods[mod]['reactants'] + mods[mod]['products'] 

        #Print message
        msg  = '||| Model {} will compute {} for '.format(mod, ' and '.join(mods[mod]['job']))
        if mods[mod]['products']:
            msg += 'the reaction: {} --> {}'.format(
                     ' + '.join([automol.inchi.smiles(spcs[x]['inchi']) for x in mods[mod]['reactants']]), 
                     ' + '.join([automol.inchi.smiles(spcs[x]['inchi']) for x in mods[mod]['products']]))
        else:
            msg += 'the species: '.format(mods[mod]['reactants'])
        msg += '\n||| using the following modules: {} '.format(', '.join(tsks[mod]))
        logging.info(msg)
        
        #Add reference molecules if thermo or rates calculations are requested
        if mods[mod]['job']:
            prepare_refs(mods, spcs)
            species += mods[mod]['references']                  
            species  = list(dict.fromkeys(species))
   
        #Begin tasks  
        for tsk in tsks[mod]:

            #prepare filesystem
            run_prefix, save_prefix = mods[mod]['paths']
            if not os.path.exists(save_prefix):
                os.makedirs(save_prefix)
            if not os.path.exists(run_prefix):
                os.makedirs(run_prefix)

            #Theory information
            overwrite = get_overwrite(tsks[mod][tsk])
            initial_thy_info = get_initial_thy_info_(lvls, tsks[mod][tsk][1])
            running_thy_info = scripts.es.get_thy_info(lvls[tsks[mod][tsk][0]])
                
            if tsk == 'tsfind':
                log.info('  | Task {} \t\t\t'.format(tsk))
                geo, zma = scripts.es.find_ts(run_prefix, save_prefix, mods[mod]['reactants'],  mods[mod]['products'], spcs, running_thy_info, overwrite)
                if not isinstance(geo, str):
                    spcs['ts_' + tsk] = create_spec(geo)
                    spcs['ts_' + tsk]['zmatrix'] = zma
                    log.info('   | Success, ts_{} added to queue'.format(tsk))
                    species.append('ts_' + tsk)
                continue

            for spc in species:
                msg = '\n  | Task {} \t\t\t Species {}: {}'.format(tsk, spc, automol.inchi.smiles(spcs[spc]['inchi']))
                log.info(msg)
                
                #Get params
                spc_info = scripts.es.get_spc_info(spcs[spc])
                
                #Make sure OPT exists for species
                if tsk not in ['sp']:
                    scripts.es.geo_at_lvl(spcs[spc], initial_thy_info, running_thy_info, run_prefix, save_prefix, overwrite, 'ts' in spc)
                    if not 'ts' in spc:
                        scripts.es.remove_imag(spcs[spc], lvls[tsks[mod][tsk][0]], run_prefix, save_prefix, overwrite)

                #Run tasks
                scripts.es.run_task(tsk, spcs[spc], lvls[tsks[mod][tsk][0]], initial_thy_info, run_prefix, save_prefix, overwrite, 'ts' in spc)
              

