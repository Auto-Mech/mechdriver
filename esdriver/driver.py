""" electronic structure drivers
"""
import os
from qcelemental import constants as qcc
import automol.inchi
import automol.geom
import moldr
import autofile.fs
import scripts.es
from esdriver.load import load_params

KICKOFF_SIZE = 0.1
KICKOFF_BACKWARD = False
PROJROT_SCRIPT_STR = ("#!/usr/bin/env bash\n"
                      "RPHt.exe >& /dev/null")

def run(tsk_info_lst, es_dct, rxn_lst, spcdct, run_prefix, save_prefix):
    """ driver for all electronic structure tasks
    """

    #prepare species queue
    spc_queue = []
    ts_dct = {}
    for i, rxn in enumerate(rxn_lst):
        reacs = rxn['reactants']
        prods = rxn['products']

        spc_queue.extend(rxn['species'])
        spc_queue.extend(reacs)
        spc_queue.extend(prods)
        if reacs and prods:
            ts_dct['ts_{:g}'.format(i)] = {'reacs': reacs, 'prods': prods}

    spc_queue = list(dict.fromkeys(spc_queue))
    # removes duplicates

    #Prepare filesystem
    if not os.path.exists(save_prefix):
        os.makedirs(save_prefix)
    if not os.path.exists(run_prefix):
        os.makedirs(run_prefix)

    #Loop over Tasks
    for tsk_info in tsk_info_lst:

        #Task information
        tsk = tsk_info[0]
        es_ini_key = tsk_info[2]
        es_run_key = tsk_info[1]
        overwrite = tsk_info[3]
        #Theory information
        ini_thy_info = get_es_info(es_dct, es_ini_key)
        thy_info = get_es_info(es_dct, es_run_key)

        #If task is to find the transition state, find all TSs for your reactionlist
        if tsk == 'find_ts':
            for ts in ts_dct:
                print('Task {} \t\t\t'.format(tsk))
                spcdct[ts] = create_spec(ts, ts_dct, spcdct)
                ts_info = (spcdct[ts]['ich'], spcdct[ts]['chg'], spcdct[ts]['mul'])
                spcdct[ts] = scripts.es.rxn_info(run_prefix, save_prefix, ts, spcdct, thy_info, ini_thy_info)
                rxn_run_fs, rxn_save_fs, rxn_run_path, rxn_save_path = scripts.es.get_rxn_fs(run_prefix, save_prefix, spcdct[ts])
                spcdct[ts]['rxn_fs'] = [rxn_run_fs, rxn_save_fs, rxn_run_path, rxn_save_path]
                rct_zmas, prd_zmas = scripts.es.get_zmas(spcdct[ts]['reacs'], spcdct[ts]['prods'], spcdct, ini_thy_info, save_prefix, run_prefix, KICKOFF_SIZE, KICKOFF_BACKWARD, PROJROT_SCRIPT_STR)
                rxn_class, ts_zma, dist_name, grid, tors_names = scripts.es.ts_params(rct_zmas, prd_zmas)
                spcdct[ts]['class']= rxn_class
                spcdct[ts]['tors_names'] = tors_names
                spcdct[ts]['original_zma']= ts_zma
                geo, zma = scripts.es.find_ts(spcdct[ts], ts_info, ts_zma, rxn_class, dist_name, grid, thy_info, rxn_run_path, rxn_save_path, overwrite)
                if not isinstance(geo, str):
                    print('Success, transition state {} added to species queue'.format(ts))
                    spc_queue.append(ts)
            continue

        #Loop over all species
        for spc in spc_queue:
            if 'ts_' in spc:
                print('\nTask {} \t\t\t Species {}'.format(tsk, spc))
                spc_run_fs, spc_save_fs, spc_run_path, spc_save_path = spcdct[spc]['rxn_fs']
                spc_info = scripts.es.get_spc_info(spcdct[spc])
            else:
                print('\nTask {} \t Species {}: {}'.format(
                    tsk, spc, automol.inchi.smiles(spcdct[spc]['ich'])))
                spc_info = scripts.es.get_spc_info(spcdct[spc])
                spc_run_fs = autofile.fs.species(run_prefix)
                spc_run_fs.leaf.create(spc_info)
                spc_run_path = spc_run_fs.leaf.path(spc_info)

                spc_save_fs = autofile.fs.species(save_prefix)
                spc_save_fs.leaf.create(spc_info)
                spc_save_path = spc_save_fs.leaf.path(spc_info)

            #Get params

            orb_restr = moldr.util.orbital_restriction(
                spc_info, thy_info)
            thy_level = thy_info[0:3]
            thy_level.append(orb_restr)
           
            thy_run_fs = autofile.fs.theory(spc_run_path)
            thy_run_fs.leaf.create(thy_level[1:4])
            thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])

            thy_save_fs = autofile.fs.theory(spc_save_path)
            thy_save_fs.leaf.create(thy_level[1:4])
            thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])

            if ini_thy_info[0] != 'input_geom':
                orb_restr = moldr.util.orbital_restriction(
                    spc_info, ini_thy_info)
                ini_thy_level = ini_thy_info[0:3]
                ini_thy_level.append(orb_restr)

                ini_thy_run_fs = autofile.fs.theory(spc_run_path)
                ini_thy_run_fs.leaf.create(ini_thy_level[1:4])
                ini_thy_run_path = ini_thy_run_fs.leaf.path(ini_thy_level[1:4])

                ini_thy_save_fs = autofile.fs.theory(spc_save_path)
                ini_thy_save_fs.leaf.create(ini_thy_level[1:4])
                ini_thy_save_path = ini_thy_save_fs.leaf.path(ini_thy_level[1:4])

                ini_cnf_run_fs = autofile.fs.conformer(ini_thy_run_path)
                ini_cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)

                ini_tau_run_fs = autofile.fs.tau(ini_thy_run_path)
                ini_tau_save_fs = autofile.fs.tau(ini_thy_save_path)
                min_cnf_locs = moldr.util.min_energy_conformer_locators(ini_cnf_save_fs)
                if min_cnf_locs:
                    min_cnf_run_path = ini_cnf_run_fs.leaf.path(min_cnf_locs)
                    min_cnf_save_path = ini_cnf_save_fs.leaf.path(min_cnf_locs)
                    ini_scn_run_fs = autofile.fs.conformer(min_cnf_run_path)
                    ini_scn_save_fs = autofile.fs.conformer(min_cnf_save_path)
                else:
                    ini_scn_run_fs = None
                    ini_scn_save_fs = None
                if  'ts' in spc:
                    ini_ts_run_fs = autofile.fs.ts(ini_thy_run_path)
                    ini_ts_run_fs.trunk.create()
                    ini_ts_save_fs = autofile.fs.ts(ini_thy_save_path)
                    ini_ts_save_fs.trunk.create()
                else:
                    ini_ts_run_fs = None
                    ini_ts_save_fs = None

            else:
                ini_thy_run_fs = None
                ini_thy_run_path = None
                ini_thy_save_fs = None
                ini_thy_save_path = None
                ini_cnf_run_fs = None
                ini_cnf_save_fs = None
                ini_tau_run_fs = None
                ini_tau_save_fs = None
                ini_thy_level = ini_thy_info
                ini_scn_run_fs = None
                ini_scn_save_fs = None
                ini_ts_run_fs = None
                ini_ts_run_path = None
                ini_ts_save_fs = None

            if 'ts' in spc:
                ts_run_fs = autofile.fs.ts(thy_run_path)
                ts_run_fs.trunk.create()
                ts_run_path = ts_run_fs.trunk.path()
                run_fs = autofile.fs.run(ts_run_path)
                ts_save_fs = autofile.fs.ts(thy_save_path)
                ts_save_fs.trunk.create()
            else:
                run_fs = autofile.fs.run(thy_run_path)
                run_fs.trunk.create()
                ts_run_fs = None
                ts_save_fs = None

            run_fs = autofile.fs.run(thy_run_path)
            run_fs.trunk.create()

            cnf_run_fs = autofile.fs.conformer(thy_run_path)
            cnf_save_fs = autofile.fs.conformer(thy_save_path)

            tau_run_fs = autofile.fs.tau(thy_run_path)
            tau_save_fs = autofile.fs.tau(thy_save_path)

            min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
            if min_cnf_locs:
                min_cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
                min_cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
                scn_run_fs = autofile.fs.conformer(min_cnf_run_path)
                scn_save_fs = autofile.fs.conformer(min_cnf_save_path)
            else:
                scn_run_fs = None
                scn_save_fs = None

            fs = [spc_run_fs, spc_save_fs, thy_run_fs, thy_save_fs, cnf_run_fs,
                  cnf_save_fs, tau_run_fs, tau_save_fs, scn_run_fs, scn_save_fs, ts_run_fs, ts_save_fs, run_fs]

            ini_fs = [ini_thy_run_fs, ini_thy_save_fs, ini_cnf_run_fs,
                    ini_cnf_save_fs, ini_tau_run_fs, ini_tau_save_fs,
                    ini_scn_run_fs, ini_scn_save_fs, ini_ts_run_fs, ini_ts_save_fs]

            #Run tasks
            if 'ts_' in spc:
                #Check if the task has been completed at the requested running theory level

                #Every task starts with a geometry optimization at the running theory level
                if 'samp' in tsk or 'scan' in tsk or 'geom' in tsk:
                    geo = moldr.ts.reference_geometry(spcdct[spc], thy_level, ini_thy_level, fs, ini_fs, overwrite)
                #Run the requested task at the requested running theory level
                    if geo:
                        scripts.es.ts_geometry_generation(tsk, spcdct[spc], es_dct[es_run_key],
                                                   thy_level, fs, spc_info, overwrite)
                else:
                    selection = 'min'
                    scripts.es.ts_geometry_analysis(tsk, thy_level, ini_fs,
                                                 selection, spc_info, overwrite)
            else:
                #Check if the task has been completed at the requested running theory level

                #Every task starts with a geometry optimization at the running theory level
                if 'samp' in tsk or 'scan' in tsk or 'geom' in tsk:
                # if not tsk == 'sp' or tsk == 'energy':
                    geo = moldr.geom.reference_geometry(
                        spcdct[spc], thy_level, ini_thy_level, fs, ini_fs,
                        kickoff_size=KICKOFF_SIZE,
                        kickoff_backward=KICKOFF_BACKWARD,
                        projrot_script_str=PROJROT_SCRIPT_STR,
                        overwrite=overwrite)
                    if geo:
                    #Run the requested task at the requested running theory level
                        scripts.es.geometry_generation(tsk, spcdct[spc], es_dct[es_run_key],
                                                   thy_level, fs, spc_info, overwrite)
                else:
                    selection = 'min'
                    if 'conf' in tsk:
                        min_cnf_locs = moldr.util.min_energy_conformer_locators(ini_cnf_save_fs)
                        if not min_cnf_locs:
                            print('Initial level of theory for conformers must be run before {} '.format(tsk))
                            continue
                        elif not ini_cnf_save_fs.leaf.file.geometry.exists(min_cnf_locs):
                            print('Initial level of theory for conformers must be run before {} '.format(tsk))
                            continue
                    elif 'tau' in tsk:
                        tau_locs = ini_tau_save_fs.leaf.existing()
                        if not tau_locs:
                            print('Initial level of theory for tau must be run before {} '.format(tsk))
                            continue
                        elif not ini_tau_save_fs.leaf.file.geometry.exists([tau_locs[0]]):
                            print('Initial level of theory for tau must be run before {} '.format(tsk))
                            continue
                    elif 'scan' in tsk:
                        scn_locs = ini_scn_save_fs.leaf.existing()
                        if not scn_locs:
                            print('Initial level of theory for scn must be run before {} '.format(tsk))
                            continue
                        elif not ini_scn_save_fs.leaf.file.geometry.exists([scn_locs[0]]):
                            print('Initial level of theory for scn must be run before {} '.format(tsk))
                            continue
                    scripts.es.geometry_analysis(tsk, thy_level, ini_fs,
                                                 selection, spc_info, overwrite)
                                                       thy_level, fs, spc_info, overwrite)
    return spcdct


def create_spec(ts, ts_dct, spcs, charge=0, mc_nsamp=[True, 3, 1, 3, 100, 12], hind_inc=30.):
    """
    Create a transition state entry for the spcdct
    """
    spec = {'ich': ''}
    spec['reacs'] = ts_dct[ts]['reacs']
    spec['prods'] = ts_dct[ts]['prods']
    ts_mul = automol.mult.ts.low([spcs[spc]['mul'] for spc in spec['reacs']], [spcs[spc]['mul'] for spc in spec['prods']])
    ts_chg = sum([spcs[spc]['chg'] for spc in spec['reacs']])
    spec['chg'] = ts_chg
    spec['mul'] = ts_mul
    spec['mc_nsamp'] = mc_nsamp
    spec['hind_inc'] = hind_inc * qcc.conversion_factor('degree', 'radian')
    return spec


def get_es_info(es_dct, key):
    """
    Turn es dictionary in theory info array
    """
    if key == 'input':
        ret = ['input_geom', None, None, None]
    else:
        ret = scripts.es.get_thy_info(es_dct[key])
    return ret

if __name__ == "__main__":

    MSG = """
           ================================================================
           ==                            ESDRIVER                        ==
           ====        Andreas Copan, Kevin Moore, Sarah Elliott,      ====
           ==   Carlo Cavolotti, Yuri Georgievski, Stephen Klippenstein  ==
           ================================================================\n"""
    print(MSG)

    TSK_INFO_LST, RXNS, ES_DCT, SPCS = load_params()

    run(tsk_info_lst, rxns, es_dct, spcs) 

