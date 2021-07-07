""" ProcessDriver Routines for Parsing the save filesystem
    and collating useful information in two a presentable
    format
"""

from mechlib.amech_io import printer as ioprinter
from mechlib import filesys
from mechroutines.proc import _util as util
# from mechroutines.proc import _collect as collect


def run_tsk(tsk, spc_dct, run_lst,
            thy_dct, proc_keyword_dct,
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

    # Begin the loop over the species
    for spc_name in run_lst:

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
                proc_keyword_dct, spc_mod_dct_i,
                save_prefix, run_prefix, spc_dct_i, thy_dct)
            cnf_fs, rng_cnf_locs_lst, rng_cnf_locs_path, mod_thy_info = ret

            # Loop over conformers
            for locs, locs_path in zip(rng_cnf_locs_lst, rng_cnf_locs_path):
                label = spc_name + '_' + '_'.join(locs)
                if 'freq' in tsk:
                    csv_data_i, csv_data_j = collect.freqs(
                        spc_dct_i, spc_mod_dct_i, proc_keyword_dct, thy_dct,
                        cnf_fs, locs, locs_path, run_prefix, save_prefix)
                    csv_data['freq'][label] = csv_data_i
                    tors_freqs, all_freqs, sfactor = csv_data_j
                    if tors_freqs is not None:
                        csv_data['tfreq'][label] = tors_freqs
                        csv_data['allfreq'][label] = all_freqs
                        csv_data['scalefactor'][label] = [sfactor]

                elif 'geo' in tsk:
                    csv_data_i = collect.geometry(
                        spc_name, locs, locs_path, cnf_fs, mod_thy_info)
                    csv_data[label] = csv_data_i

                elif 'molden' in tsk:
                    csv_data_i = collect.molden(
                        spc_name, locs, locs_path, cnf_fs, mod_thy_info)
                    csv_data[label] = csv_data_i

                elif 'zma' in tsk:
                    csv_data_i = collect.zmatrix(
                        spc_name, locs, locs_path, cnf_fs, mod_thy_info)
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

    # write the csv data into the appropriate file
    util.write_csv_data(tsk, csv_data, filelabel, spc_array)
