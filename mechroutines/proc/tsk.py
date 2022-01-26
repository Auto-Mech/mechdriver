""" ProcessDriver Routines for Parsing the save filesystem
    and collating useful information in two a presentable
    format
"""

from mechlib.amech_io import printer as ioprinter
from mechroutines.models import typ
from mechroutines.proc import _util as util
from mechroutines.proc import _collect as collect


def run_tsk(tsk, obj_queue,
            proc_keyword_dct,
            spc_dct, thy_dct,
            spc_mod_dct_i, pes_mod_dct_i,
            run_prefix, save_prefix, mdriver_path):
    """ run a proc tess task
    for generating a list of conformer or tau sampling geometries
    """

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
    # obj_queue, ts_miss_data = util.remove_ts_missing(
    #     obj_queue, spc_dct)
    # obj_queue = util.remove_radrad_ts(
    #     obj_queue, spc_dct)

    # Set up lists for reporting missing data
    miss_data = ()
    ts_miss_data = ()
    # Initialize dictionaries to carry strings for writing
    disp_dct = {}
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
                spc_name, spc_dct, pes_mod_dct_i, spc_array)
            csv_data[label] = csv_data_i
            col_array = spc_array
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

                miss_data_i = None
                label = spc_name + ':' + '_'.join(locs)
                print(label)

                if 'freq' in tsk and not _skip(spc_name, spc_dct_i):
                    _dat, miss_data_i = collect.frequencies(
                        spc_name, spc_dct_i, spc_mod_dct_i,
                        proc_keyword_dct, thy_dct,
                        cnf_fs, locs, locs_path, run_prefix, save_prefix)
                    if _dat is not None:
                        csv_data_i, csv_data_j, disp_str = _dat
                        csv_data['freq'][label] = csv_data_i
                        tors_freqs, all_freqs, sfactor = csv_data_j
                        if tors_freqs is not None:
                            csv_data['tfreq'][label] = tors_freqs
                            csv_data['allfreq'][label] = all_freqs
                            csv_data['scalefactor'][label] = [sfactor]
                        if disp_str is not None:
                            disp_dct.update({spc_name: disp_str})

                elif 'geo' in tsk:
                    csv_data_i, miss_data_i = collect.geometry(
                        spc_name, locs, locs_path, cnf_fs, mod_thy_info)
                    print(csv_data_i)
                    csv_data[label] = csv_data_i

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
                        spc_name, spc_dct_i, spc_mod_dct_i,
                        mod_thy_info, run_prefix, save_prefix)

                elif 'ene' in tsk:
                    csv_data_i, miss_data_i = collect.energy(
                        spc_name, spc_dct_i, spc_mod_dct_i,
                        proc_keyword_dct, thy_dct, locs, locs_path,
                        cnf_fs, run_prefix, save_prefix)
                    csv_data[label] = csv_data_i

                elif 'enthalpy' in tsk:
                    ret = collect.enthalpy(
                        spc_name, spc_dct, spc_dct_i, spc_mod_dct_i,
                        pes_mod_dct_i, chn_basis_ene_dct, spc_array,
                        locs, locs_path, cnf_fs, run_prefix, save_prefix)
                    csv_data_i, chn_basis_ene_dct, spc_array = ret
                    csv_data[label] = csv_data_i
                    col_array = spc_array

                elif 'entropy' in tsk:
                    ret = collect.enthalpy(
                        spc_name, spc_dct, spc_dct_i, spc_mod_dct_i,
                        pes_mod_dct_i, chn_basis_ene_dct, spc_array,
                        locs, locs_path, cnf_fs, run_prefix, save_prefix)
                    csv_data_i, chn_basis_ene_dct, spc_array = ret
                    csv_data[label] = csv_data_i

                elif 'heat' in tsk:
                    ret = collect.enthalpy(
                        spc_name, spc_dct, spc_dct_i, spc_mod_dct_i,
                        pes_mod_dct_i, chn_basis_ene_dct, spc_array,
                        locs, locs_path, cnf_fs, run_prefix, save_prefix)
                    csv_data_i, chn_basis_ene_dct, spc_array = ret
                    csv_data[label] = csv_data_i

                elif 'messpf_inp' in tsk:
                    ret = collect.messpf_input(
                        spc_name, spc_dct_i, spc_mod_dct_i,
                        pes_mod_dct_i, locs, locs_path,
                        cnf_fs, run_prefix, save_prefix)
                    csv_data_i, _, miss_data_i = ret
                    print(csv_data_i)
                    csv_data[label] = csv_data_i

                elif 'pf' in tsk:
                    ret = collect.partition_function(
                        spc_name, spc_dct_i, spc_mod_dct_i,
                        pes_mod_dct_i, locs, locs_path,
                        cnf_fs, run_prefix, save_prefix)
                    csv_data_i, miss_data_i = ret
                    csv_data[label] = csv_data_i[1]
                    col_array = csv_data_i[0]

                if miss_data_i is not None:
                    miss_data += (miss_data_i,)


    # Write a report that details what data is missing
    missing_data = miss_data + ts_miss_data

    # Write the csv data into the appropriate file
    util.write_csv_data(tsk, csv_data, filelabel, col_array, mdriver_path)

    # Gather data that is provided for each species in files in a dir
    data_dirs = (('displacements_'+thylabel, disp_dct),)
    util.write_data_dirs(data_dirs, mdriver_path)

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
