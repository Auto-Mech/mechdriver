""" electronic structure drivers
"""
import os
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

def run(tsk_info_lst, es_dct, rxn_lst, spc_dct, run_prefix, save_prefix,
        vdw_params=[False, False, True]):
    """ driver for all electronic structure tasks
    """

    print(tsk_info_lst)
    #prepare species queue
    spc_queue = []
    for _, rxn in enumerate(rxn_lst):
        reacs = rxn['reacs']
        prods = rxn['prods']

        spc_queue.extend(rxn['species'])
        spc_queue.extend(reacs)
        spc_queue.extend(prods)

    spc_queue = list(dict.fromkeys(spc_queue))
    # removes duplicates

    #Prepare filesystem
    if not os.path.exists(save_prefix):
        os.makedirs(save_prefix)
    if not os.path.exists(run_prefix):
        os.makedirs(run_prefix)

    #Loop over Tasks
    ts_found = []
    for tsk_info in tsk_info_lst:

        #Task information
        tsk = tsk_info[0]
        es_ini_key = tsk_info[2]
        es_run_key = tsk_info[1]
        overwrite = tsk_info[3]
        #Theory information
        ini_thy_info = scripts.es.get_thy_info(es_dct[es_ini_key])
        thy_info = scripts.es.get_thy_info(es_dct[es_run_key])
        ini_thy_info = get_es_info(es_dct, es_ini_key)
        thy_info = get_es_info(es_dct, es_run_key)

        #If task is to find the transition state, find all TSs for your reactionlist
        if tsk in ('find_ts', 'find_vdw'):
            for ts in spc_dct:
                if 'ts_' in ts:
                    print('Task {} \t for {} \t {}//{} \t {} = {}'.format(
                        tsk, ts, '/'.join(thy_info), '/'.join(ini_thy_info),
                        '+'.join(spc_dct[ts]['reacs']), '+'.join(spc_dct[ts]['prods'])))
                    ts_info = (spc_dct[ts]['ich'], spc_dct[ts]['chg'], spc_dct[ts]['mul'])
                    rxn_class = spc_dct[ts]['class']
                    if not rxn_class:
                        print('skipping reaction because type =', rxn_class)
                        continue
                    elif 'radical radical' in rxn_class and not 'high spin' in rxn_class:
                        print('skipping reaction because type =', rxn_class)
                        continue
                    ts_zma = spc_dct[ts]['original_zma']
                    dist_info = spc_dct[ts]['dist_info']
                    grid = spc_dct[ts]['grid']
                    _, _, rxn_run_path, rxn_save_path = spc_dct[ts]['rxn_fs']
                    if 'ts' in tsk:
                        geo, _, final_dist = scripts.es.find_ts(
                            spc_dct[ts], ts_info, ts_zma, rxn_class, dist_info,
                            grid, thy_info, rxn_run_path, rxn_save_path, overwrite)
                        spc_dct[ts]['dist_info'][1] = final_dist
                        if not isinstance(geo, str):
                            print('Success, transition state {} added to species queue'.format(ts))
                            spc_queue.append(ts)
                            ts_found.append(ts)
                    elif 'vdw' in tsk:
                        vdws = scripts.es.find_vdw(
                            ts, spc_dct, thy_info, ini_thy_info, ts_info, vdw_params,
                            es_dct[es_run_key]['mc_nsamp'], run_prefix,
                            save_prefix, KICKOFF_SIZE, KICKOFF_BACKWARD,
                            PROJROT_SCRIPT_STR, overwrite)
                        spc_queue.extend(vdws)
            continue


        #Loop over all species
        for spc in spc_queue:
            if 'ts_' in spc:
                print('\nTask {} \t {}//{} \t Species {}'.format(
                    tsk, '/'.join(thy_info), '/'.join(ini_thy_info), spc))
                spc_run_fs, spc_save_fs, spc_run_path, spc_save_path = spc_dct[spc]['rxn_fs']
                spc_info = scripts.es.get_spc_info(spc_dct[spc])

            else:
                print('\nTask {} \t {}//{} \t Species {}: {}'.format(
                    tsk, '/'.join(thy_info), '/'.join(ini_thy_info), spc,
                    automol.inchi.smiles(spc_dct[spc]['ich'])))
                spc_info = scripts.es.get_spc_info(spc_dct[spc])
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
            thy_save_fs = autofile.fs.theory(spc_save_path)

            print(spc)
            if 'ene' not in tsk and 'hess' not in tsk:
                if 'ts_' in spc:
                    thy_run_fs.leaf.create(thy_level[1:4])
                    thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])
                    thy_save_fs.leaf.create(thy_level[1:4])
                    thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])

                    thy_run_fs = autofile.fs.ts(thy_run_path)
                    thy_run_fs.trunk.create()
                    thy_run_path = thy_run_fs.trunk.path()

                    thy_save_fs = autofile.fs.ts(thy_save_path)
                    thy_save_fs.trunk.create()
                    thy_save_path = thy_save_fs.trunk.path()

                else:
                    thy_run_fs.leaf.create(thy_level[1:4])
                    thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])
                    thy_save_fs.leaf.create(thy_level[1:4])
                    thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])

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
            else:
                cnf_run_fs = None
                cnf_save_fs = None
                tau_run_fs = None
                tau_save_fs = None
                scn_run_fs = None
                scn_save_fs = None

            if ini_thy_info[0] != 'input_geom':
                orb_restr = moldr.util.orbital_restriction(
                    spc_info, ini_thy_info)
                ini_thy_level = ini_thy_info[0:3]
                ini_thy_level.append(orb_restr)

                ini_thy_run_fs = autofile.fs.theory(spc_run_path)
                ini_thy_save_fs = autofile.fs.theory(spc_save_path)
                if 'ts_' in spc:
                    ini_thy_run_fs.leaf.create(ini_thy_level[1:4])
                    ini_thy_run_path = ini_thy_run_fs.leaf.path(ini_thy_level[1:4])
                    ini_thy_save_fs.leaf.create(ini_thy_level[1:4])
                    ini_thy_save_path = ini_thy_save_fs.leaf.path(ini_thy_level[1:4])

                    ini_thy_run_fs = autofile.fs.ts(ini_thy_run_path)
                    ini_thy_run_fs.trunk.create()
                    ini_thy_run_path = ini_thy_run_fs.trunk.path()

                    ini_thy_save_fs = autofile.fs.ts(ini_thy_save_path)
                    ini_thy_save_fs.trunk.create()
                    ini_thy_save_path = ini_thy_save_fs.trunk.path()

                else:
                    ini_thy_run_fs.leaf.create(ini_thy_level[1:4])
                    ini_thy_run_path = ini_thy_run_fs.leaf.path(ini_thy_level[1:4])
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

            run_fs = autofile.fs.run(thy_run_path)
            run_fs.trunk.create()

            fs = [spc_run_fs, spc_save_fs, thy_run_fs, thy_save_fs,
                  cnf_run_fs, cnf_save_fs, tau_run_fs, tau_save_fs,
                  scn_run_fs, scn_save_fs, run_fs]

            ini_fs = [ini_thy_run_fs, ini_thy_save_fs, ini_cnf_run_fs,
                      ini_cnf_save_fs, ini_tau_run_fs, ini_tau_save_fs,
                      ini_scn_run_fs, ini_scn_save_fs]

            #Run tasks
            if 'ts_' in spc:
                if 'samp' in tsk or 'scan' in tsk or 'geom' in tsk:
                    geo = moldr.ts.reference_geometry(
                        spc_dct[spc], thy_level, ini_thy_level, fs, ini_fs,
                        spc_dct[spc]['dist_info'], overwrite)
                    if geo:
                        scripts.es.ts_geometry_generation(
                            tsk, spc_dct[spc], es_dct[es_run_key],
                            thy_level, fs, spc_info, overwrite)
                else:
                    selection = 'min'
                    scripts.es.ts_geometry_analysis(
                        tsk, thy_level, ini_fs, selection, spc_info, overwrite)
            else:
                if 'samp' in tsk or 'scan' in tsk or 'geom' in tsk:
                    geo = moldr.geom.reference_geometry(
                        spc_dct[spc], thy_level, ini_thy_level, fs, ini_fs,
                        kickoff_size=KICKOFF_SIZE,
                        kickoff_backward=KICKOFF_BACKWARD,
                        projrot_script_str=PROJROT_SCRIPT_STR,
                        overwrite=overwrite)
                    if geo:
                        if not 'vdw_' in spc:
                            scripts.es.geometry_generation(
                                tsk, spc_dct[spc], es_dct[es_run_key], thy_level,
                                fs, spc_info, overwrite)
                        else:
                            scripts.es.fake_geo_gen(
                                tsk, spc_dct[spc], es_dct[es_run_key], thy_level,
                                fs, spc_info, overwrite)
                else:
                    selection = 'min'
                    if 'conf' in tsk:
                        min_cnf_locs = moldr.util.min_energy_conformer_locators(ini_cnf_save_fs)
                        if not min_cnf_locs:
                            print(
                                'Initial level of theory for conformers must be ',
                                'run before {} '.format(tsk))
                            continue
                        elif not ini_cnf_save_fs.leaf.file.geometry.exists(min_cnf_locs):
                            print(
                                'Initial level of theory for conformers must be ',
                                'run before {} '.format(tsk))
                            continue
                    elif 'tau' in tsk:
                        tau_locs = ini_tau_save_fs.leaf.existing()
                        if not tau_locs:
                            print(
                                'Initial level of theory for tau must be ', 
                                'run before {} '.format(tsk))
                            continue
                        elif not ini_tau_save_fs.leaf.file.geometry.exists([tau_locs[0]]):
                            print(
                                'Initial level of theory for tau must be ',
                                'run before {} '.format(tsk))
                            continue
                    elif 'scan' in tsk:
                        scn_locs = ini_scn_save_fs.leaf.existing()
                        if not scn_locs:
                            print(
                                'Initial level of theory for scn must be run ',
                                'before {} '.format(tsk))
                            continue
                        elif not ini_scn_save_fs.leaf.file.geometry.exists([scn_locs[0]]):
                            print(
                                'Initial level of theory for scn must be run ',
                                'before {} '.format(tsk))
                            continue
                    scripts.es.geometry_analysis(tsk, thy_level, ini_fs,
                            selection, spc_info, overwrite)
    return ts_found


# def create_ts_spec(ts, ts_dct, spcs, charge=0, hind_inc=30.):
    # """
    # Create a transition state entry for the spc_dct
    # """
    # ts_spec = {'ich': ''}
    # ts_spec['reacs'] = ts_dct[ts]['reacs']
    # ts_spec['prods'] = ts_dct[ts]['prods']
    # ts_mul = ts_dct[ts]['mul']
    # ts_chg = sum([spcs[spc]['chg'] for spc in ts_spec['reacs']])
    # ts_spec['chg'] = ts_chg
    # ts_spec['mul'] = ts_mul
    # ts_spec['hind_inc'] = hind_inc * qcc.conversion_factor('degree', 'radian')
    # return ts_spec


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
           ==                        AUTOMECHANIC                        ==
           ===         Andreas Copan, Sarah Elliott, Kevin Moore,       ===
           ===     Daniel Moberg, Carlo Cavallotti, Yuri Georgievski,   ===
           ==       Ahren Jasper, Murat Keceli, Stephen Klippenstein     ==
           ================================================================
           ==                          ESDRIVER                          ==
           ====        Sarah Elliott, Andreas Copan, Kevin Moore,      ====
           ==            Carlo Cavolotti, Stephen Klippenstein           ==
           ================================================================\n"""
    print(MSG)

    # this line needs to be updated 
    TSK_INFO_LST, RXN_LST, ES_DCT, SPCDCT, RUN_PREFIX, SAVE_PREFIX = load_params()

    run(TSK_INFO_LST, ES_DCT, RXN_LST, SPCDCT, RUN_PREFIX, SAVE_PREFIX)

