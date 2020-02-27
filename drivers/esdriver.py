""" electronic structure drivers
"""

import routines
from lib.load import run as loadrun
from lib.load import species as loadspc
from lib.load import mechanism as loadmech
from lib.filesystem import check as fcheck
from lib.filesystem import path as fpath
from lib import printmsg


def run(rxn_lst,
        spc_dct,
        es_tsk_str,
        model_dct,
        thy_dct,
        run_options_dct,
        run_inp_dct):
    """ driver for all electronic structure tasks
    """

    # Print the header message for the driver
    printmsg.program_header('es')

    # Pull stuff from dcts for now
    run_prefix = run_inp_dct['run_prefix']
    save_prefix = run_inp_dct['save_prefix']
    # vdw_params = model_dct['options']['vdw_params']
    #freeze_all_tors = model_dct[model]['options']['freeze_all_tors']
    #ndim_tors = model_dct[model]['pf']['tors']
    # rad_rad_ts = model_dct[model]['pf']['ts_barrierless']
    freeze_all_tors = True
    ndim_tors = '1dhr'
    rad_rad_ts = 'pst'
    mc_nsamp = run_options_dct['mc_nsamp']
    kickoff = run_options_dct['kickoff']

    # Do some extra work to prepare the info to pass to the drivers
    es_tsk_lst = loadrun.build_run_es_tsks_lst(
        es_tsk_str, model_dct, thy_dct)

    # Species queue
    spc_queue = loadmech.build_spc_queue(rxn_lst)

    # Loop over Tasks
    for tsk_info in es_tsk_lst:

        # Task and theory information
        [tsk, thy_info, ini_thy_info, overwrite] = tsk_info

        # If task is to find the transition state, find all TSs for rxn lst
        if tsk in ('find_ts', 'find_vdw'):
            # Get info for the transition states
            ts_dct = loadspc.build_sadpt_dct(
                rxn_lst, model_dct, thy_dct, es_tsk_str,
                run_inp_dct, run_options_dct, spc_dct, {})
            spc_dct.update(ts_dct)
            for spc in spc_dct:
                if 'ts_' in spc:
                    # printmsg.sadpt_tsk_printmsg(
                    #     tsk, sadpt, spc_dct, thy_info, ini_thy_info)
                    if not spc_dct[spc]['class']:
                        print('skipping reaction because type =',
                              spc_dct[spc]['class'])
                        continue

                    # Find the transition state
                    if 'ts' in tsk:
                        geo, _, _ = routines.es.find.find_ts(
                            spc_dct, spc_dct[spc],
                            spc_dct[spc]['zma'],
                            # spc_dct[spc]['original_zma'],
                            ini_thy_info, thy_info,
                            run_prefix, save_prefix,
                            overwrite,
                            rad_rad_ts=rad_rad_ts)

                        # Add TS to species queue if TS is found
                        if not isinstance(geo, str):
                            print('Success, transition state',
                                  '{} added to species queue'.format(spc))
                            spc_queue.append((spc, ''))
                    elif 'vdw' in tsk:
                        pass
                        # vdws = routines.es.wells.find_vdw(
                        #     spc, spc_dct, thy_info, ini_thy_info,
                        #     vdw_params,
                        #     thy_dct[es_run_key]['mc_nsamp'], run_prefix,
                        #     save_prefix, 0.1, False,
                        #     overwrite)
                        # spc_queue.extend(vdws)
            continue

        # Loop over all species
        for spc in spc_queue:

            spc_name, _ = spc

            # Build the input and main run filesystem objects
            filesys, thy_level = fpath.set_fs(
                spc_dct, spc_name, thy_info,
                run_prefix, save_prefix,
                setfs_chk=True, ini_fs=False)
            ini_filesys, ini_thy_level = fpath.set_fs(
                spc_dct, spc_name, ini_thy_info,
                run_prefix, save_prefix,
                setfs_chk=bool(ini_thy_info[0] != 'input_geom'),
                ini_fs=True)

            # Run tasks
            saddle = bool('ts_' in spc_name)
            vdw = bool('vdw' in spc_name)
            avail = True
            if any(string in tsk for string in ('samp', 'scan', 'geom')):
                if not vdw:
                    routines.es.geometry_generation(
                        tsk, spc_name, spc_dct[spc_name], mc_nsamp,
                        ini_thy_level, thy_level, ini_filesys, filesys,
                        overwrite, saddle=saddle, kickoff=kickoff,
                        tors_model=(ndim_tors, freeze_all_tors))
                else:
                    routines.es.wells.fake_geo_gen(tsk, thy_level, filesys)
            else:
                if 'conf' in tsk and not saddle:
                    ini_cnf_save_fs = ini_filesys[3]
                    avail = fcheck.check_save(ini_cnf_save_fs, tsk, 'conf')
                    selection = 'min'
                elif 'tau' in tsk and not saddle:
                    ini_tau_save_fs = ini_filesys[5]
                    avail = fcheck.check_save(ini_tau_save_fs, tsk, 'tau')
                    selection = 'all'
                elif 'scan' in tsk and not saddle:
                    ini_scn_save_fs = ini_filesys[7]
                    avail = fcheck.check_save(ini_scn_save_fs, tsk, 'scan')
                    selection = 'all'
                if avail:
                    routines.es.geometry_analysis(
                        tsk, spc_name, thy_level, ini_filesys,
                        spc_dct[spc_name], overwrite,
                        saddle=saddle, selection=selection)
