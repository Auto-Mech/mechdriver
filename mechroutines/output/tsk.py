""" eletronic structure routines modules
"""

import sys
import importlib
import pandas

import autofile
import automol
from phydat import phycon
# from mechroutines.es.newts import run as runts
from mechlib.amech_io import printer as ioprinter
from mechlib.amech_io import parser
from mechlib import filesys
from mechlib.structure import instab
from mechroutines.pf.models import _vib as vib
from mechroutines.pf.models import ene
from mechroutines.pf.thermo import basis
from mechroutines.pf.models.inf import set_pf_info
from . import _util as util
from autofile import io_ as io

# Dictionary of Electronic Structure Calculators


def run_tsk(tsk, spc_dct, rxn_lst,
            thy_dct, print_keyword_dct,
            model_dct,
            run_prefix, save_prefix):
    """ run an electronic structure task
    for generating a list of conformer or tau sampling geometries
    """

    # Print the head of the task
    ioprinter.output_task_header(tsk)
    ioprinter.obj('line_dash')
    ioprinter.output_keyword_list(print_keyword_dct, thy_dct)

    # If species is unstable, set task to 'none'
    if 'enthalpy' in tsk:
        chn_basis_ene_dct = {}
        spc_array = []

    csv_data = {}
    spc_queue = parser.species.build_spc_queue(rxn_lst)
    for spc_name, (pf_model, spc_model) in spc_queue:

        if print_keyword_dct['geolvl']:
            thy_info = filesys.inf.get_es_info(
                print_keyword_dct['geolvl'], thy_dct)
        else:
            pf_levels, pf_models, _, _ = set_pf_info(
                model_dct, thy_dct, pf_model, pf_model)
            ioprinter.debug_message('test', pf_levels)
            ioprinter.debug_message('test2', thy_dct)
            thy_info = pf_levels['geo'][1]
            #    filesys.inf.get_es_info(
            #    pf_levels['geo'][1], thy_dct)

        ioprinter.obj('line_dash')
        ioprinter.info_message("Species: ", spc_name)
        stable = instab.check_unstable_species(
            tsk, spc_dct, spc_name, thy_info, save_prefix)
        if not stable:
            ioprinter.info_message(
                'Properties for species {}'.format(spc_name),
                'will not be included in output',
                'because species is unstable')
            continue

        spc_dct_i = spc_dct[spc_name]

        # # Set the filesystem objects
        # thy_run_fs, thy_run_path = filesys.build.spc_thy_fs_from_root(
        #     run_prefix, spc_info, mod_thy_info)
        # thy_save_fs, thy_save_path = filesys.build.spc_thy_fs_from_root(
        #     save_prefix, spc_info, mod_thy_info)
        # _, ini_thy_save_path = filesys.build.spc_thy_fs_from_root(
        #     save_prefix, spc_info, mod_ini_thy_info)
        # cnf_save_fs = autofile.fs.conformer(thy_save_path)
        # run_fs = filesys.build.run_fs_from_prefix(thy_run_path)

        # Loop over conformers
        if print_keyword_dct['geolvl']:
            cnf_fs, cnf_locs_lst, cnf_locs_paths = util.conformer_list(
                print_keyword_dct, save_prefix, run_prefix,
                spc_dct_i, thy_dct)
            pf_levels, pf_models = None, None
        else:
            pf_levels, pf_models, _, _ = set_pf_info(
                model_dct, thy_dct, spc_model, spc_model)
            ret = util.conformer_list_from_models(
                print_keyword_dct, save_prefix, run_prefix,
                spc_dct_i, thy_dct, pf_levels, pf_models)
            cnf_fs, cnf_locs_lst, cnf_locs_paths = ret
            # ref_scheme = pf_model['ref_scheme']
            # ref_enes = pf_models['ref_enes']
        for locs, locs_path in zip(cnf_locs_lst, cnf_locs_paths):
            label = spc_name + '_' + '_'.join(locs)
            if 'freq' in tsk:

                es_model = util.freq_es_levels(print_keyword_dct)
                pf_levels = parser.model.pf_level_info(
                    es_model, thy_dct)

                freqs, _, zpe = vib.read_locs_harmonic_freqs(
                    cnf_fs, locs_path, locs, saddle=False)
                #if pf_models:
                #    pf_filesystems = filesys.models.pf_filesys(
                #        spc_dct_i, pf_levels,
                #        run_prefix, save_prefix, saddle=False)
                #    rotors = tors.build_rotors(
                #        spc_dct_i, pf_filesystems, pf_models, pf_levels)
                #    #    rxn_class=rxn_class,
                #    #    frm_bnd_keys=frm_bnd_keys, brk_bnd_keys=brk_bnd_keys)
                #    if typ.nonrigid_tors(pf_models, rotors):
                #        freqs, imag, tors_zpe, pot_scalef = ret
                #        ret = vib.tors_projected_freqs_zpe(
                #            pf_filesystems, hr_str, prot_str,
                #            run_prefix, saddle=False) # saddle)
                #        if typ.scale_1d(pf_models):
                #            tors_strs = tors.make_hr_strings(
                #                rotors, run_path, pf_models['tors'],
                #                scale_factor=pot_scalef)
                #            [allr_str, hr_str, _, prot_str, mdhr_dat] = tors_strs
                #            _, _, tors_zpe, _ = vib.tors_projected_freqs_zpe(
                #                pf_filesystems, hr_str,
                #                prot_str, run_prefix, saddle=False) # saddle)
                #            # Calculate current zpe assuming no freq scaling: tors+projfreq
 #              #         zpe = tors_zpe + (sum(freqs) / 2.0) * phycon.WAVEN2EH

                tors_zpe = 0.0
                spc_data = []
                if freqs:
                    freqs, zpe = vib.scale_frequencies(
                        freqs, tors_zpe, pf_levels, scale_method='3c')
                    spc_data = [locs_path, zpe, *freqs]
                csv_data[label] = spc_data

            elif 'geo' in tsk:

                if cnf_fs[-1].file.geometry.exists(locs):
                    geo = cnf_fs[-1].file.geometry.read(locs)
                    energy = cnf_fs[-1].file.energy.read(locs)
                    comment = 'energy: {0:>15.10f}'.format(energy)
                    xyz_str = automol.geom.xyz_string(geo, comment=comment)
                else:
                    xyz_str = '\t -- Missing --'
                spc_data = '\n\nSPC: {}\tConf: {}\tPath: {}\n'.format(
                    spc_name, locs, locs_path) + xyz_str
                csv_data[label] = spc_data

            elif 'zma' in tsk:
                geo = cnf_fs[-1].file.geometry.read(locs)
                zma = automol.geom.zmatrix(geo)
                energy = cnf_fs[-1].file.energy.read(locs)
                comment = 'energy: {0:>15.10f}\n'.format(energy)
                zma_str = automol.zmatrix.string(zma)
                spc_data = '\n\nSPC: {}\tConf: {}\tPath: {}\n'.format(
                    spc_name, locs, locs_path) + comment + zma_str
                csv_data[label] = spc_data

            elif 'ene' in tsk:
                es_model = util.ene_es_levels(print_keyword_dct)
                energy = cnf_fs[-1].file.energy.read(locs)
                csv_data[label] = [locs_path, energy]

            elif 'enthalpy' in tsk:
                pf_filesystems = filesys.models.pf_filesys(
                    spc_dct_i, pf_levels,
                    run_prefix, save_prefix, saddle=False)
                ene_abs = ene.read_energy(
                    spc_dct_i, pf_filesystems, pf_models, pf_levels,
                    run_prefix, read_ene=True, read_zpe=True, saddle=False)
                ts_geom = None
                hf0k, _, chn_basis_ene_dct, hbasis = basis.enthalpy_calculation(
                    spc_dct, spc_name, ts_geom, ene_abs,
                    chn_basis_ene_dct, pf_levels, pf_models,
                    run_prefix, save_prefix, pforktp='pf', saddle=False)
                hf0k = hf0k * phycon.KCAL2EH
                spc_basis, coeff_basis = hbasis[spc_name]
                coeff_array = []
                for spc_i in spc_basis:
                    if spc_i not in spc_array:
                        spc_array.append(spc_i)
                for spc_i in spc_array:
                    if spc_i in spc_basis:
                        coeff_array.append(coeff_basis[spc_basis.index(spc_i)])
                    else:
                        coeff_array.append(0)
                csv_data[label] = [locs_path, ene_abs, hf0k, *coeff_array]

    if 'freq' in tsk:
        ncols = max([len(x) for x in csv_data.values()])
        df = pandas.DataFrame.from_dict(
            csv_data, orient='index',
            columns=['Path', 'ZPVE', *[''] * (ncols-2)])
        df.to_csv('freqs.csv', float_format='%.3f')
    if 'geo' in tsk:
        all_data = '\n'.join([spc_data for spc_data in csv_data.values()])
        io.write_file('geometry.txt', all_data)
    if 'zma' in tsk:
        all_data = '\n'.join([spc_data for spc_data in csv_data.values()])
        io.write_file('zmatrix.txt', all_data)
    if 'ene' in tsk:
        df = pandas.DataFrame.from_dict(
            csv_data, orient='index',
            columns=['Path', 'Energy [A.U.]'])
        df.to_csv('ene.csv', float_format='%.8f')
    if 'enthalpy' in tsk:
        df = pandas.DataFrame.from_dict(
            csv_data, orient='index',
            columns=[
                'Path', 'Energy [A.U.]', 'Hf (0 K) [kcal/mol]', *spc_array])
        df.to_csv('enthalpy.csv', float_format='%.6f')
