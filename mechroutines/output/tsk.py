""" eletronic structure routines modules
"""

import os
import sys
import importlib

import autofile
import automol
from phydat import phycon
# from mechroutines.es.newts import run as runts
from mechlib.amech_io import printer as ioprinter
from mechlib.amech_io import parser
from mechlib import filesys
# from mechlib.structure import instab
from mechroutines.pf.models import _vib as vib
from mechroutines.pf.models import _tors as tors
from mechroutines.pf.models import typ
from mechroutines.pf.models import ene
from mechroutines.pf.thermo import basis
from mechroutines.pf.models.inf import set_pf_info
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo
from . import _util as util

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

    # Setup csv data dictionary for specific task
    csv_data = util.set_csv_data(tsk)
    chn_basis_ene_dct = {}
    spc_array = []

    # Loop over species to collect task info
    spc_queue = parser.species.build_spc_queue(rxn_lst)
    for spc_name, (pf_model, spc_model) in spc_queue:

        # print species
        ioprinter.obj('line_dash')
        ioprinter.info_message("Species: ", spc_name)

        # If species is unstable, set task to 'none'
        # stable = instab.check_unstable_species(
        #     tsk, spc_dct, spc_name, thy_info, save_prefix)
        # if not stable:
        #     ioprinter.info_message(
        #         'Properties for species {}'.format(spc_name),
        #         'will not be included in output',
        #         'because species is unstable')
        #     continue

        # Heat of formation basis molecules and coefficients
        # is not conformer specific
        if 'coeffs' in tsk:
            pf_levels, pf_models, _, _ = set_pf_info(
                model_dct, thy_dct, pf_model, pf_model)
            thy_info = pf_levels['geo'][1]
            filelabel = 'coeffs'
            filelabel += '_{}'.format(pf_models['ref_scheme'])
            filelabel += '.csv'
            label = spc_name
            basis_dct, uniref_dct = basis.prepare_refs(
                pf_models['ref_scheme'], spc_dct, [[spc_name, None]])
            # Get the basis info for the spc of interest
            spc_basis, coeff_basis = basis_dct[spc_name]
            coeff_array = []
            for spc_i in spc_basis:
                if spc_i not in spc_array:
                    spc_array.append(spc_i)
            for spc_i in spc_array:
                if spc_i in spc_basis:
                    coeff_array.append(coeff_basis[spc_basis.index(spc_i)])
                else:
                    coeff_array.append(0)
            csv_data[label] = [*coeff_array]

        else:
            # unpack spc and level info
            spc_dct_i = spc_dct[spc_name]
            if print_keyword_dct['geolvl']:
                thy_info = tinfo.from_dct(thy_dct.get(
                    print_keyword_dct['geolvl']))
            else:
                pf_levels, pf_models, _, _ = set_pf_info(
                    model_dct, thy_dct, pf_model, pf_model)
                thy_info = pf_levels['geo'][1]

                # Loop over conformers
            if print_keyword_dct['geolvl']:
                _, rng_cnf_locs_lst, rng_cnf_locs_path = util.conformer_list(
                    print_keyword_dct, save_prefix, run_prefix,
                    spc_dct_i, thy_dct)
                pf_levels, pf_models = None, None
            else:
                pf_levels, pf_models, _, _ = set_pf_info(
                    model_dct, thy_dct, spc_model, spc_model)
                ret = util.conformer_list_from_models(
                    print_keyword_dct, save_prefix, run_prefix,
                    spc_dct_i, thy_dct, pf_levels, pf_models)
                _, rng_cnf_locs_lst, rng_cnf_locs_path = ret
            for locs, locs_path in zip(rng_cnf_locs_lst, rng_cnf_locs_path):

                label = spc_name + '_' + '_'.join(locs)
                _, cnf_fs = filesys.build_fs(
                    run_prefix, save_prefix, 'CONFORMER')
                if 'freq' in tsk:

                    filelabel = 'freq'
                    if pf_levels:
                        filelabel += '_m{}'.format(pf_levels['harm'][0])
                    else:
                        filelabel += '_{}'.format(print_keyword_dct['geolvl'])
                    filelabel += '.csv'

                    if pf_models:
                        pf_filesystems = filesys.models.pf_filesys(
                            spc_dct_i, pf_levels,
                            run_prefix, save_prefix, saddle=False)
                        ret = vib.full_vib_analysis(
                            spc_dct_i, pf_filesystems, pf_models, pf_levels,
                            run_prefix, saddle=False)
                        freqs, imag, tors_zpe, scale_fact, tors_freqs, all_freqs = ret
                        csv_data['tfreq'][label] = tors_freqs
                        csv_data['allfreq'][label] = all_freqs
                        csv_data['scalefactor'][label] = [scale_fact]
                    else:
                        es_model = util.freq_es_levels(print_keyword_dct)
                        pf_levels = parser.model.pf_level_info(
                            es_model, thy_dct)
                        try:
                            freqs, _, zpe = vib.read_locs_harmonic_freqs(
                                cnf_fs, locs_path, locs, saddle=False)
                        except:
                            freqs = []
                            zpe = 0

                    tors_zpe = 0.0
                    spc_data = []
                    zpe = tors_zpe + (sum(freqs) / 2.0) * phycon.WAVEN2EH
                    if freqs and print_keyword_dct['scale'] is not None:
                        freqs, zpe = vib.scale_frequencies(
                            freqs, tors_zpe, pf_levels, scale_method='3c')
                    spc_data = [locs_path, zpe, *freqs]
                    csv_data['freq'][label] = spc_data
                elif 'geo' in tsk:

                    filelabel = 'geo'
                    if pf_levels:
                        filelabel += '_{}'.format(pf_levels['harm'])
                    else:
                        filelabel += '_{}'.format(print_keyword_dct['geolvl'])
                    filelabel += '.txt'

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
                    
                    filelabel = 'zmat'
                    if pf_levels:
                        filelabel += '_{}'.format(pf_levels['harm'])
                    else:
                        filelabel += '_{}'.format(print_keyword_dct['geolvl'])
                    filelabel += '.txt'

                    geo = cnf_fs[-1].file.geometry.read(locs)
                    zma = automol.geom.zmatrix(geo)
                    energy = cnf_fs[-1].file.energy.read(locs)
                    comment = 'energy: {0:>15.10f}\n'.format(energy)
                    zma_str = automol.zmat.string(zma)
                    spc_data = '\n\nSPC: {}\tConf: {}\tPath: {}\n'.format(
                        spc_name, locs, locs_path) + comment + zma_str
                    csv_data[label] = spc_data

                elif 'ene' in tsk:
                    
                    filelabel = 'ene'
                    if pf_levels:
                        filelabel += '_{}'.format(pf_levels['harm'])
                        filelabel += '_{}'.format(pf_levels['ene'])
                    else:
                        filelabel += '_{}'.format(print_keyword_dct['geolvl'])
                        filelabel += '_{}'.format(print_keyword_dct['proplvl'])
                    filelabel += '.csv'

                    energy = None
                    if pf_levels:
                        pf_filesystems = filesys.models.pf_filesys(
                            spc_dct_i, pf_levels,
                            run_prefix, save_prefix, saddle=False)
                        energy = ene.electronic_energy(
                            spc_dct_i, pf_filesystems, pf_levels,
                            conf=(locs, locs_path, cnf_fs))
                    else:
                        spc_info = sinfo.from_dct(spc_dct_i)
                        thy_info = tinfo.from_dct(thy_dct.get(
                            print_keyword_dct['proplvl']))
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
                    filelabel = 'enthalpy'
                    if pf_levels:
                        filelabel += '_{}'.format(pf_levels['harm'])
                        filelabel += '_{}'.format(pf_levels['ene'])
                    else:
                        filelabel += '_{}'.format(print_keyword_dct['geolvl'])
                        filelabel += '_{}'.format(print_keyword_dct['proplvl'])
                    filelabel = '.csv'

                    energy = None
                    pf_filesystems = filesys.models.pf_filesys(
                        spc_dct_i, pf_levels,
                        run_prefix, save_prefix, saddle=False)
                    ene_abs = ene.read_energy(
                        spc_dct_i, pf_filesystems, pf_models, pf_levels,
                        run_prefix, conf=(locs, locs_path, cnf_fs),
                        read_ene=True, read_zpe=True, saddle=False)
                    ts_geom = None
                    hf0k, _, chn_basis_ene_dct, hbasis = basis.enthalpy_calculation(
                        spc_dct, spc_name, ts_geom, ene_abs,
                        chn_basis_ene_dct, pf_levels, pf_models,
                        run_prefix, save_prefix, pforktp='pf', saddle=False)
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
