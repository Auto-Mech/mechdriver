""" ProcessDriver Routines for Parsing the save filesystem
    and collating useful information in two a presentable
    format
"""

import os
import autofile
import automol
from phydat import phycon
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo
from mechlib.amech_io import printer as ioprinter
from mechlib.amech_io import parser
from mechlib import filesys
from mechroutines.pf.models import _vib as vib
from mechroutines.pf.models import ene
from mechroutines.pf.thermo import basis
from mechroutines.output import _util as util


def run_tsk(tsk, spc_dct, spc_name,
            thy_dct, proc_keyword_dct,
            pes_mod_dct_i, spc_mod_dct_i,
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
    chn_basis_ene_dct = {}
    spc_array = []

    # print species
    ioprinter.obj('line_dash')
    ioprinter.info_message("Species: ", spc_name)

    # Heat of formation basis molecules and coefficients
    # is not conformer specific
    if 'coeffs' in tsk:
        thy_info = spc_mod_dct_i['geo'][1]
        filelabel = 'coeffs'
        filelabel += '_{}'.format(pes_mod_dct_i['thermfit']['ref_scheme'])
        filelabel += '.csv'
        label = spc_name
        basis_dct, _ = basis.prepare_refs(
            pes_mod_dct_i['thermfit']['ref_scheme'],
            spc_dct, [[spc_name, None]],
            run_prefix, save_prefix)
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
        if proc_keyword_dct['geolvl']:
            thy_info = tinfo.from_dct(thy_dct.get(
                proc_keyword_dct['geolvl']))
        else:
            thy_info = spc_mod_dct_i['geo'][1]

            # Loop over conformers
        if proc_keyword_dct['geolvl']:
            _, rng_cnf_locs_lst, rng_cnf_locs_path = util.conformer_list(
                proc_keyword_dct, save_prefix, run_prefix,
                spc_dct_i, thy_dct)
            spc_mod_dct_i, pf_models = None, None
        else:
            ret = util.conformer_list_from_models(
                proc_keyword_dct, save_prefix, run_prefix,
                spc_dct_i, thy_dct, spc_mod_dct_i, pf_models)
            _, rng_cnf_locs_lst, rng_cnf_locs_path = ret
        for locs, locs_path in zip(rng_cnf_locs_lst, rng_cnf_locs_path):

            label = spc_name + '_' + '_'.join(locs)
            _, cnf_fs = filesys.build_fs(
                run_prefix, save_prefix, 'CONFORMER')
            if 'freq' in tsk:

                filelabel = 'freq'
                if spc_mod_dct_i:
                    filelabel += '_m{}'.format(spc_mod_dct_i['harm'][0])
                else:
                    filelabel += '_{}'.format(proc_keyword_dct['geolvl'])
                filelabel += '.csv'

                if pf_models:
                    pf_filesystems = filesys.models.pf_filesys(
                        spc_dct_i, spc_mod_dct_i,
                        run_prefix, save_prefix, saddle=False)
                    ret = vib.full_vib_analysis(
                        spc_dct_i, pf_filesystems, spc_mod_dct_i,
                        run_prefix, zrxn=None)
                    freqs, _, tors_zpe, sfactor, torsfreqs, all_freqs = ret
                    csv_data['tfreq'][label] = torsfreqs
                    csv_data['allfreq'][label] = all_freqs
                    csv_data['scalefactor'][label] = [sfactor]
                else:
                    es_model = util.freq_es_levels(proc_keyword_dct)
                    spc_mod_dct_i = parser.model.pf_level_info(
                        es_model, thy_dct)
                    try:
                        freqs, _, zpe = vib.read_locs_harmonic_freqs(
                            cnf_fs, locs, run_prefix, zrxn=None)
                    except:
                        freqs = []
                        zpe = 0

                tors_zpe = 0.0
                spc_data = []
                zpe = tors_zpe + (sum(freqs) / 2.0) * phycon.WAVEN2EH
                if freqs and proc_keyword_dct['scale'] is not None:
                    freqs, zpe = vib.scale_frequencies(
                        freqs, tors_zpe, spc_mod_dct_i, scale_method='3c')
                spc_data = [locs_path, zpe, *freqs]
                csv_data['freq'][label] = spc_data
            elif 'geo' in tsk:

                filelabel = 'geo'
                if spc_mod_dct_i:
                    filelabel += '_{}'.format(spc_mod_dct_i['harm'])
                else:
                    filelabel += '_{}'.format(proc_keyword_dct['geolvl'])
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
                if spc_mod_dct_i:
                    filelabel += '_{}'.format(spc_mod_dct_i['harm'])
                else:
                    filelabel += '_{}'.format(proc_keyword_dct['geolvl'])
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
                if spc_mod_dct_i:
                    filelabel += '_{}'.format(spc_mod_dct_i['harm'])
                    filelabel += '_{}'.format(spc_mod_dct_i['ene'])
                else:
                    filelabel += '_{}'.format(proc_keyword_dct['geolvl'])
                    filelabel += '_{}'.format(proc_keyword_dct['proplvl'])
                filelabel += '.csv'

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
                filelabel = 'enthalpy'
                if spc_mod_dct_i:
                    filelabel += '_{}'.format(spc_mod_dct_i['harm'])
                    filelabel += '_{}'.format(spc_mod_dct_i['ene'])
                else:
                    filelabel += '_{}'.format(proc_keyword_dct['geolvl'])
                    filelabel += '_{}'.format(proc_keyword_dct['proplvl'])
                filelabel = '.csv'

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
                    chn_basis_ene_dct, pes_mod_dct_i, spc_mod_dct_i,
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
