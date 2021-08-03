""" ProcessDriver Routines for Parsing the save filesystem
    and collating useful information in two a presentable
    format
"""

from mechanalyzer.inf import thy as tinfo
from mechlib import filesys
from mechlib.amech_io import printer as ioprinter
from mechroutines.proc import _util as util
from mechroutines.proc import _collect as collect


def run_tsk(tsk, obj_queue,
            proc_keyword_dct,
            spc_dct, thy_dct,
            spc_mod_dct_i, model_dct,
            run_prefix, save_prefix):
    """ run a proc tess task
    for generating a list of conformer or tau sampling geometries
    """

    # Print the head of the task
    ioprinter.output_task_header(tsk)
    ioprinter.obj('line_dash')
    ioprinter.output_keyword_list(proc_keyword_dct, thy_dct)

    # Setup csv data dictionary for specific task
    csv_data = util.set_csv_data(tsk)
    filelabel = util.get_file_label(
        tsk, model_dct, proc_keyword_dct, spc_mod_dct_i)
    chn_basis_ene_dct = {}
    spc_array = []

    # Set up lists for reporting missing data
    miss_data = ()

<<<<<<< HEAD
    # Exclude unstable species
    # These species break certain checks (e.g. no ene exists for geo collect)
    obj_queue, ts_miss_data = _remove_ts_missing(
        obj_queue, spc_dct)
    obj_queue = _remove_unstable(
        obj_queue, spc_dct, thy_dct, proc_keyword_dct, save_prefix)

=======
>>>>>>> ts for procdriver
    # Begin the loop over the species
    for spc_name in obj_queue:

        # info printed to output file
        ioprinter.obj('line_dash')
        ioprinter.info_message("Species: ", spc_name)

        # Heat of formation basis molecules and coefficients
        # does not require filesystem information
        if 'coeffs' in tsk:
            label = spc_name
            csv_data_i, spc_array = collect.coeffs(
                spc_name, spc_dct, model_dct, spc_array)
            csv_data[label] = csv_data_i

        # All other tasks require filesystem information
        else:
            # unpack spc and level info for conformers
            spc_dct_i = spc_dct[spc_name]
            spc_mod_dct_i = util.choose_theory(
                proc_keyword_dct, spc_mod_dct_i)
            ret = util.choose_conformers(
                spc_name, proc_keyword_dct, spc_mod_dct_i,
                save_prefix, run_prefix, spc_dct_i, thy_dct)
            cnf_fs, rng_cnf_locs_lst, rng_cnf_locs_path, mod_thy_info = ret

            # Add geo to missing data task if locs absent
            if not rng_cnf_locs_lst:
                miss_data += ((spc_name, mod_thy_info, 'geometry'),)

            # Loop over conformers
            for locs, locs_path in zip(rng_cnf_locs_lst, rng_cnf_locs_path):
                label = spc_name + '_' + '_'.join(locs)
                if 'freq' in tsk:
                    csv_data_i, csv_data_j, miss_data_i = collect.frequencies(
                        spc_dct_i, spc_mod_dct_i, proc_keyword_dct, thy_dct,
                        cnf_fs, locs, locs_path, run_prefix, save_prefix)
                    csv_data['freq'][label] = csv_data_i
                    tors_freqs, all_freqs, sfactor = csv_data_j
                    if tors_freqs is not None:
                        csv_data['tfreq'][label] = tors_freqs
                        csv_data['allfreq'][label] = all_freqs
                        csv_data['scalefactor'][label] = [sfactor]
                    if miss_data_i is not None:
                        miss_data += (miss_data_i,)

                elif 'geo' in tsk:
                    csv_data_i, miss_data_i = collect.geometry(
                        spc_name, locs, locs_path, cnf_fs, mod_thy_info)
                    csv_data[label] = csv_data_i
                    if miss_data_i is not None:
                        miss_data += (miss_data_i,)

                elif 'molden' in tsk:
                    csv_data_i = collect.molden(
                        spc_name, locs, locs_path, cnf_fs, mod_thy_info)
                    csv_data[label] = csv_data_i

                elif 'zma' in tsk:
                    csv_data_i = collect.zmatrix(
                        spc_name, locs, locs_path, cnf_fs, mod_thy_info)
                    csv_data[label] = csv_data_i

                elif 'torsion' in tsk:
                    csv_data_i = collect.torsions(
                        spc_name, spc_dct_i, spc_mod_dct_i,
                        mod_thy_info,
                        run_prefix, save_prefix)
                    csv_data[label] = csv_data_i

                elif 'ene' in tsk:
                    csv_data_i = collect.energy(
                        spc_name, spc_dct_i, spc_mod_dct_i,
                        proc_keyword_dct, thy_dct, locs, locs_path,
                        cnf_fs, run_prefix, save_prefix)
                    csv_data[label] = csv_data_i

                elif 'enthalpy' in tsk:
                    ret = collect.enthalpy(
                        spc_name, spc_dct, spc_dct_i, spc_mod_dct_i,
                        model_dct, chn_basis_ene_dct, spc_array,
                        locs, locs_path, cnf_fs, run_prefix, save_prefix)
                    csv_data_i, chn_basis_ene_dct, spc_array = ret
                    csv_data[label] = csv_data_i

                elif 'entropy' in tsk:
                    ret = collect.enthalpy(
                        spc_name, spc_dct, spc_dct_i, spc_mod_dct_i,
                        model_dct, chn_basis_ene_dct, spc_array,
                        locs, locs_path, cnf_fs, run_prefix, save_prefix)
                    csv_data_i, chn_basis_ene_dct, spc_array = ret
                    csv_data[label] = csv_data_i

                elif 'heat' in tsk:
                    ret = collect.enthalpy(
                        spc_name, spc_dct, spc_dct_i, spc_mod_dct_i,
                        model_dct, chn_basis_ene_dct, spc_array,
                        locs, locs_path, cnf_fs, run_prefix, save_prefix)
                    csv_data_i, chn_basis_ene_dct, spc_array = ret
                    csv_data[label] = csv_data_i

    # Write a report that details what data is missing
<<<<<<< HEAD
    missing_data = miss_data + ts_miss_data
    util.write_missing_data_report(missing_data)
=======
    util.write_missing_data_report(miss_data)
>>>>>>> ts for procdriver

    # Write the csv data into the appropriate file
    util.write_csv_data(tsk, csv_data, filelabel, spc_array)


def _remove_unstable(spc_queue, spc_dct, thy_dct, proc_key_dct, save_prefix):
    """ For each species in the queue see if there are files
        in the save filesystem denoting they are unstable. If so,
        that species is removed from the queue for collection tasks.
    """

    thy_info = tinfo.from_dct(thy_dct.get(proc_key_dct.get('geolvl')))

    stable_queue = ()
    for spc_name in spc_queue:
        if 'ts_' in spc_name:
            stable_queue += (spc_name,)
        else:
            instab, path = filesys.read.instability_transformation(
                spc_dct, spc_name, thy_info, save_prefix)
            if instab is None:
                stable_queue += (spc_name,)
            else:
                ioprinter.info_message(
                    'Found instability file at path {}'.format(path),
                    newline=1)
                ioprinter.info_message(
                    'Removing {} from queue'.format(spc_name))

    return stable_queue


def _remove_ts_missing(obj_queue, spc_dct):
    """ Generate list of missing data, remove ts from queue
    """

    new_queue = ()
    ts_miss_data = ()
    for obj in obj_queue:
        if 'ts_' in obj:
            ts_dct = spc_dct[obj]
            miss = ts_dct.get('missdata')
            if miss is not None:
                ts_miss_data += ((obj, ts_dct['missdata'], 'geometry'),)
            else:
                new_queue += (obj,)
        else:
            new_queue += (obj,)

    return new_queue, ts_miss_data
