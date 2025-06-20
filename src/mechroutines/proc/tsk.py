""" ProcessDriver Routines for Parsing the save filesystem
    and collating useful information in two a presentable
    format
"""

from mechlib.amech_io import printer as ioprinter
from mechroutines.models import typ
from mechroutines.proc import _util as util
from mechroutines.proc import _collect as collect
from autorun import execute_function_in_parallel

def run_tsk(tsk, obj_queue,
            proc_keyword_dct,
            spc_dct, thy_dct,
            spc_mod_dct, pes_mod_dct,
            run_prefix, save_prefix, mdriver_path):
    """ run a proc tess task
    for generating a list of conformer or tau sampling geometries
    """

    # Update spc_mod_dct from targetted levels in run.dat
    spc_mod_dct_i, pes_mod_dct_i = util.reconcile_spc_mod_and_proc_keyword_dcts(
        tsk, spc_mod_dct, proc_keyword_dct, pes_mod_dct, thy_dct)
    # Print the head of the task
    ioprinter.output_task_header(tsk)
    ioprinter.output_keyword_list(proc_keyword_dct, thy_dct)

    # Setup csv data dictionary for specific task
    csv_data = util.set_csv_data(tsk)
    filelabel, thylabel = util.get_file_label(
        tsk, pes_mod_dct_i, proc_keyword_dct, spc_mod_dct_i)
    chn_basis_ene_dct = {}
    col_array = [] 
    spc_array = []

    # Exclude unstable species
    # These species break certain checks (e.g. no ene exists for geo collect)
    obj_queue = util.remove_unstable(
        obj_queue, spc_dct, thy_dct, spc_mod_dct_i,
        proc_keyword_dct, save_prefix)
    obj_queue, ts_miss_data = util.remove_ts_missing(
        obj_queue, spc_dct)
    # obj_queue = util.remove_radrad_ts(
    #     obj_queue, spc_dct)


    # Set up lists for reporting missing data
    miss_data = ()
    # disp_dct = {}
    # Begin the loop over the species
    for spc_name in obj_queue:

        # Species Header in logfile
        ioprinter.obj('line_dash')
        ioprinter.info_message("Species: ", spc_name)

        # Heat of formation basis molecules and coefficients
        # does not require filesystem information
        if 'coeffs' in tsk:
            label = spc_name
            csv_data_i, spc_array = collect.coeffs(
                spc_name, spc_dct, pes_mod_dct_i, spc_array)
            csv_data[label] = csv_data_i
            col_array = spc_array
        # All other tasks require filesystem information
        else:
            # unpack spc and level info for conformers
            spc_dct_i = spc_dct[spc_name]
            ret = util.choose_conformers(
                spc_name, proc_keyword_dct, spc_mod_dct_i,
                save_prefix, run_prefix, spc_dct, thy_dct)
            cnf_info, mod_thy_info, symid_dct = ret
            cnf_fs, rng_cnf_locs_lst, rng_cnf_locs_path = cnf_info
            # Add geo to missing data task if locs absent
            if not rng_cnf_locs_lst:
                miss_data += ((spc_name, mod_thy_info, 'geometry'),)
            # Loop over conformers

            args = (
                tsk, spc_name, spc_dct, spc_dct_i, spc_mod_dct_i,
                proc_keyword_dct, thy_dct, mod_thy_info,
                pes_mod_dct_i, chn_basis_ene_dct, spc_array,
                rng_cnf_locs_lst, rng_cnf_locs_path, cnf_fs,
                symid_dct,
                run_prefix, save_prefix, col_array, miss_data)
            ret_lst = execute_function_in_parallel(
                _run_task_for_locs_lst, list(rng_cnf_locs_lst), args,
                nprocs=proc_keyword_dct['nprocs'])
            species_csv_data = {}
            species_miss_data = ()
            for (csv_data_j, col_array, miss_data_j) in ret_lst:
                if 'weight' in tsk:
                    if not 'hf_array' in species_csv_data:
                        species_csv_data['hf_array'] = []
                        species_csv_data['locs_lst'] = []
                        species_csv_data['pf_array'] = []
                    species_csv_data['hf_array'] += csv_data_j['hf_array']
                    species_csv_data['locs_lst'] += csv_data_j['locs_lst']
                    species_csv_data['pf_array'] += csv_data_j['pf_array']
                else:
                    species_csv_data.update(csv_data_j)
                species_miss_data += miss_data_j
            miss_data += species_miss_data
            if 'weight' in tsk:
                locs_lst, weight_lst = collect.pf_weights(
                    species_csv_data['locs_lst'], 
                    species_csv_data['hf_array'],
                    species_csv_data['pf_array'])
                for locs, weight in zip(locs_lst, weight_lst):
                    label = spc_name + ':' + ':'.join(locs)
                    csv_data[label] = [weight]
            elif 'freqs' in tsk:
                csv_data['freq'].update(species_csv_data['freq'])
                csv_data['tfreq'].update(species_csv_data['tfreq'])
                csv_data['allfreq'].update(species_csv_data['allfreq'])
                csv_data['scalefactor'].update(species_csv_data['scalefactor'])
                print("species: ",species_csv_data)
                print("total: ", csv_data)
            else:
                csv_data.update(species_csv_data)
    # Write a report that details what data is missing
    missing_data = miss_data + ts_miss_data

    # Write the csv data into the appropriate file
    util.write_csv_data(tsk, csv_data, filelabel, col_array, mdriver_path)

    # # Gather data that is provided for each species in files in a dir
    # data_dirs = (('displacements_'+thylabel, disp_dct),)
    # util.write_data_dirs(data_dirs, mdriver_path)

    return missing_data


# Task manager/skipper functions
def _skip(spc_name, spc_dct_i):
    """ check if frequencies should be skipped
    """
    skip = False
    if 'ts' not in spc_name:
        if typ.is_atom(spc_dct_i):
            skip = True
    return skip


def _run_task_for_locs_lst(
        tsk, spc_name, spc_dct, spc_dct_i, spc_mod_dct_i,
        proc_keyword_dct, thy_dct, mod_thy_info,
        pes_mod_dct_i, chn_basis_ene_dct, spc_array,
        locs_lst, locs_path_lst, cnf_fs, symid_dct,
        run_prefix, save_prefix,
        col_array, miss_data, proc_locs_lst, output_queue=None):
    """
    """
    csv_data = util.set_csv_data(tsk)
    for locs, locs_path in zip(locs_lst, locs_path_lst):
        if not any(locs[1] == sel_locs[1] for sel_locs in proc_locs_lst):
            continue
        miss_data_i = None
        label = spc_name + ':' + ':'.join(locs)
        print(label)

        if 'freq' in tsk and not _skip(spc_name, spc_dct_i):
            _dat, miss_data_i = collect.frequencies(
                spc_name, spc_dct_i, spc_mod_dct_i,
                proc_keyword_dct, thy_dct, spc_dct,
                cnf_fs, locs, locs_path, run_prefix, save_prefix)
            if _dat is not None:
                csv_data_i, csv_data_j, disp_str = _dat
                csv_data['freq'][label] = csv_data_i
                tors_freqs, all_freqs, sfactor = csv_data_j
                if tors_freqs is not None:
                    csv_data['tfreq'][label] = tors_freqs
                    csv_data['allfreq'][label] = all_freqs
                    csv_data['scalefactor'][label] = [sfactor]
                # if disp_str is not None:
                #     disp_dct.update({spc_name: disp_str})

        elif 'geo' in tsk:
            csv_data_i, miss_data_i = collect.geometry(
                spc_name, locs, locs_path, cnf_fs, mod_thy_info)
            print(csv_data_i)
            csv_data[label] = csv_data_i
        elif 'date' in tsk:
            csv_data_i, date_headers, miss_data_i = collect.time_stamp(
                spc_name, locs, locs_path, cnf_fs, mod_thy_info)
            csv_data[label] = csv_data_i
            col_array = date_headers
        elif 'molden' in tsk:
            csv_data_i, miss_data_i = collect.molden(
                spc_name, locs, locs_path, cnf_fs, mod_thy_info)
            print(csv_data_i)
            csv_data[label] = csv_data_i

        elif 'zma' in tsk:
            csv_data_i, miss_data_i = collect.zmatrix(
                spc_name, locs, locs_path, cnf_fs, mod_thy_info)
            csv_data[label] = csv_data_i

        elif 'torsion' in tsk and not _skip(spc_name, spc_dct_i):
            csv_data_i, miss_data_i = collect.torsions(
                spc_name, locs, locs_path, spc_dct_i, spc_mod_dct_i,
                mod_thy_info, run_prefix, save_prefix)
            print(csv_data_i)
            csv_data[label] = csv_data_i

        elif 'hess_json' in tsk:
            csv_data_i, miss_data_i = collect.hess_json(
                spc_name, spc_dct_i, spc_mod_dct_i,
                proc_keyword_dct, thy_dct,
                cnf_fs, locs, locs_path, run_prefix, save_prefix, mod_thy_info)
            csv_data[label] = csv_data_i

        elif 'ene' in tsk:
            csv_data_i, miss_data_i = collect.energy(
                spc_name, spc_dct_i, spc_mod_dct_i,
                proc_keyword_dct, thy_dct, locs, locs_path,
                cnf_fs, run_prefix, save_prefix)
            csv_data[label] = csv_data_i
            if ':'.join(locs) in symid_dct:
                rid_label = spc_name + ':' + locs[0] + ':'
                for symid in symid_dct[':'.join(locs)]:
                    csv_data[rid_label + ':'.join(symid)] = csv_data_i 
                

        elif 'enthalpy' in tsk or 'weight' in tsk:
            ret = collect.enthalpy(
                spc_name, spc_dct, spc_dct_i, spc_mod_dct_i,
                pes_mod_dct_i, chn_basis_ene_dct, spc_array,
                locs, locs_path, cnf_fs, run_prefix, save_prefix)
            csv_data_i, chn_basis_ene_dct, spc_array, miss_data_i = ret
            if 'weight' in tsk:
                csv_data['hf_array'].append(csv_data_i[2])
                csv_data['locs_lst'].append(locs)
            else:
                csv_data[label] = csv_data_i
                col_array = spc_array
                if ':'.join(locs) in symid_dct:
                    rid_label = spc_name + ':' + locs[0] + ':'
                    for symid in symid_dct[':'.join(locs)]:
                        csv_data[rid_label + ':'.join(symid)] = csv_data_i 

        elif 'entropy' in tsk:
            # this is enthalpy not entropy
            ret = collect.enthalpy(
                spc_name, spc_dct, spc_dct_i, spc_mod_dct_i,
                pes_mod_dct_i, chn_basis_ene_dct, spc_array,
                locs, locs_path, cnf_fs, run_prefix, save_prefix)
            csv_data_i, chn_basis_ene_dct, spc_array, miss_data_i = ret
            csv_data[label] = csv_data_i

        elif 'gibbs' in tsk:
            ret = collect.relative_gibbs(
                cnf_fs, mod_thy_info, locs, locs_path,
                locs_lst[0], locs_path_lst[0])
            csv_data[label] = [ret]
            if ':'.join(locs) in symid_dct:
                rid_label = spc_name + ':' + locs[0] + ':'
                for symid in symid_dct[':'.join(locs)]:
                    csv_data[rid_label + ':'.join(symid)] = [ret]
            
        elif 'heat' in tsk:
            ret = collect.enthalpy(
                spc_name, spc_dct, spc_dct_i, spc_mod_dct_i,
                pes_mod_dct_i, chn_basis_ene_dct, spc_array,
                locs, locs_path, cnf_fs, run_prefix, save_prefix)
            csv_data_i, chn_basis_ene_dct, spc_array, miss_data_i = ret
            csv_data[label] = csv_data_i
            if ':'.join(locs) in symid_dct:
                rid_label = spc_name + ':' + locs[0] + ':'
                for symid in symid_dct[':'.join(locs)]:
                    csv_data[rid_label + ':'.join(symid)] = csv_data_i

        elif 'messpf_inp' in tsk:
            ret = collect.messpf_input(
                spc_name, spc_dct_i, spc_mod_dct_i,
                pes_mod_dct_i, locs, locs_path,
                cnf_fs, run_prefix, save_prefix)
            csv_data_i, _, miss_data_i = ret
            print(csv_data_i)
            csv_data[label] = csv_data_i

        if 'pf' in tsk or 'weight' in tsk:
            ret = collect.partition_function(
                spc_name, spc_dct_i, spc_mod_dct_i,
                pes_mod_dct_i, locs, locs_path,
                cnf_fs, run_prefix, save_prefix)
            csv_data_i, miss_data_i = ret
            if 'weight' in tsk:
                csv_data['pf_array'].append(csv_data_i)
            else:
                csv_data[label] = csv_data_i[1]
                col_array = csv_data_i[0]
        elif 'si' in tsk:
            csv_data_i, miss_data_i = collect.sidata(
                spc_name, spc_dct_i, spc_mod_dct_i,
                proc_keyword_dct, thy_dct,
                cnf_fs, locs, locs_path, run_prefix, save_prefix, mod_thy_info)
            print(csv_data_i)
            csv_data[label] = csv_data_i

        if miss_data_i is not None:
            miss_data += (miss_data_i,)
    output_queue.put(((csv_data, col_array, miss_data,),))
    return output_queue
