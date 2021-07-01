""" ProcessDriver Routines for Parsing the save filesystem
    and collating useful information in two a presentable
    format
"""

import os
import autofile
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo
from mechlib.amech_io import printer as ioprinter
from mechlib import filesys
from mechroutines.proc import _util as util
from mechroutines.proc import _collect as collect


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

    # print species
    for spc_name in run_lst:
        ioprinter.obj('line_dash')
        ioprinter.info_message("Species: ", spc_name)

        # Heat of formation basis molecules and coefficients
        # is not conformer specific but for the rest of the tasks
        # we have to generate a list of conformers to get properties for
        if 'coeffs' in tsk:
            label = spc_name
            csv_data_i, spc_array = collect.coeffs(
                spc_name, spc_dct, model_dct, spc_array)
            csv_data[label] = csv_data_i

        else:
            # unpack spc and level info for conformers
            spc_dct_i = spc_dct[spc_name]
            thy_info, spc_mod_dct_i = util.choose_theory(
                proc_keyword_dct, spc_mod_dct_i, thy_dct)
            rng_cnf_locs_lst, rng_cnf_locs_path = util.choose_conformers(
                proc_keyword_dct, spc_mod_dct_i,
                save_prefix, run_prefix, spc_dct_i, thy_dct)

            # Loop over conformers
            for locs, locs_path in zip(rng_cnf_locs_lst, rng_cnf_locs_path):

                label = spc_name + '_' + '_'.join(locs)
                _, cnf_fs = filesys.build_fs(
                    run_prefix, save_prefix, 'CONFORMER')
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
                        spc_name, locs, locs_path, cnf_fs)
                    csv_data[label] = csv_data_i

                elif 'zma' in tsk:
                    csv_data_i = collect.zmatrix(
                        spc_name, locs, locs_path, cnf_fs)
                    csv_data[label] = csv_data_i

                elif 'ene' in tsk:
                    energy = None
                    if spc_mod_dct_i:
                        pf_filesystems = filesys.models.pf_filesys(
                            spc_dct_i, spc_mod_dct_i,
                            run_prefix, save_prefix, saddle=False)
                        energy = ene.electronic_energy(
                            spc_dct_i, pf_filesystems, spc_mod_dct_i,
                            conf=(locs, locs_path, cnf_fs))
                    else:
                        spc_info = sinfo.from_dct(spc_dct_i)
                        thy_info = tinfo.from_dct(thy_dct.get(
                            proc_keyword_dct['proplvl']))
                        mod_thy_info = tinfo.modify_orb_label(
                            thy_info, spc_info)
                        sp_save_fs = autofile.fs.single_point(locs_path)
                        sp_save_fs[-1].create(mod_thy_info[1:4])
                        # Read the energy
                        sp_path = sp_save_fs[-1].path(mod_thy_info[1:4])
                        if os.path.exists(sp_path):
                            if sp_save_fs[-1].file.energy.exists(
                                    mod_thy_info[1:4]):
                                ioprinter.reading('Energy', sp_path)
                                energy = sp_save_fs[-1].file.energy.read(
                                    mod_thy_info[1:4])
                    csv_data[label] = [locs_path, energy]

                elif 'enthalpy' in tsk:
                    energy = None
                    pf_filesystems = filesys.models.pf_filesys(
                        spc_dct_i, spc_mod_dct_i,
                        run_prefix, save_prefix, saddle=False)
                    ene_abs = ene.read_energy(
                        spc_dct_i, pf_filesystems, spc_mod_dct_i,
                        run_prefix, conf=(locs, locs_path, cnf_fs),
                        read_ene=True, read_zpe=True, saddle=False)
                    hf0k, _, chn_basis_ene_dct, hbasis = basis.enthalpy_calculation(
                        spc_dct, spc_name, ene_abs,
                        chn_basis_ene_dct, model_dct, spc_mod_dct_i,
                        run_prefix, save_prefix, pforktp='pf', zrxn=None)
                    spc_basis, coeff_basis = hbasis[spc_name]
                    coeff_array = []
                    for spc_i in spc_basis:
                        if spc_i not in spc_array:
                            spc_array.append(spc_i)
                    for spc_i in spc_array:
                        if spc_i in spc_basis:
                            coeff_array.append(
                                coeff_basis[spc_basis.index(spc_i)])
                        else:
                            coeff_array.append(0)
                    csv_data[label] = [locs_path, ene_abs, hf0k, *coeff_array]

    util.write_csv_data(tsk, csv_data, filelabel, spc_array)
