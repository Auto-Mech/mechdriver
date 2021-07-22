"""Collects the target info
"""

import os

import thermfit
import automol
import autofile
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo

from mechlib import filesys
from mechroutines.models import _vib as vib
from mechroutines.models import ene
from mechroutines.thermo import basis
from mechroutines.proc import _util as util
from mechlib.amech_io import printer as ioprinter


def zmatrix(spc_name, locs, locs_path, cnf_fs, mod_thy_info):
    """collect a zmatrix
    """
    if cnf_fs[-1].file.geometry.exists(locs):
        geo = cnf_fs[-1].file.geometry.read(locs)
        zma = automol.geom.zmatrix(geo)
        sp_fs = autofile.fs.single_point(locs_path)
        energy = sp_fs[-1].file.energy.read(mod_thy_info[1:4])
        comment = 'energy: {0:>15.10f}\n'.format(energy)
        zma_str = automol.zmat.string(zma)
    else:
        zma_str = '\t -- Missing --'
    spc_data = '\n\nSPC: {}\tConf: {}\tPath: {}\n'.format(
        spc_name, locs, locs_path) + comment + zma_str
    return spc_data


def molden(spc_name, locs, locs_path, cnf_fs, mod_thy_info):
    """collect a geometry
    """
    if cnf_fs[-1].file.geometry.exists(locs):
        geo = cnf_fs[-1].file.geometry.read(locs)
        sp_fs = autofile.fs.single_point(locs_path)
        energy = sp_fs[-1].file.energy.read(mod_thy_info[1:4])
        comment = 'energy: {0:>15.10f}'.format(energy)
        comment += 'SPC: {}\tConf: {}\tPath: {}'.format(
        spc_name, locs, locs_path)
        xyz_str = automol.geom.xyz_string(geo, comment=comment)
    else:
        xyz_str = '\t -- Missing --'
    spc_data = xyz_str
    return spc_data


def geometry(spc_name, locs, locs_path, cnf_fs, mod_thy_info):
    """collect a geometry
    """
    if cnf_fs[-1].file.geometry.exists(locs):
        geo = cnf_fs[-1].file.geometry.read(locs)
        sp_fs = autofile.fs.single_point(locs_path)
        energy = sp_fs[-1].file.energy.read(mod_thy_info[1:4])
        comment = 'energy: {0:>15.10f}'.format(energy)
        xyz_str = automol.geom.xyz_string(geo, comment=comment)
    else:
        xyz_str = '\t -- Missing --'
    spc_data = '\n\nSPC: {}\tConf: {}\tPath: {}\n'.format(
        spc_name, locs, locs_path) + xyz_str
    return spc_data


def freqs(
        spc_dct_i, spc_mod_dct_i, proc_keyword_dct, thy_dct,
        cnf_fs, locs, locs_path, run_prefix, save_prefix):
    """collect frequencies
    """
    if spc_mod_dct_i is not None:
        pf_filesystems = filesys.models.pf_filesys(
            spc_dct_i, spc_mod_dct_i,
            run_prefix, save_prefix, saddle=False)
        ret = vib.full_vib_analysis(
            spc_dct_i, pf_filesystems, spc_mod_dct_i,
            run_prefix, zrxn=None)
        freqs, _, zpe, sfactor, tors_strs, torsfreqs, all_freqs = ret
    else:
        es_levels = util.freq_es_levels(proc_keyword_dct)
        spc_mod_dct_i = util.generate_spc_model_dct(es_levels, thy_dct)
        freqs, _, zpe = vib.read_locs_harmonic_freqs(
            cnf_fs, locs, run_prefix, zrxn=None)
        if freqs and proc_keyword_dct['scale'] is not None:
            freqs, zpe = vib.scale_frequencies(
                freqs, 0.0, spc_mod_dct_i, scale_method='3c')
        torsfreqs = None
        all_freqs = None
        sfactor = None
    spc_data = [locs_path, zpe, *freqs]
    return spc_data, (torsfreqs, all_freqs, sfactor)


def coeffs(spc_name, spc_dct, model_dct, spc_array):
    """ get the heat of formation reference molecules for one species.
    """
    basis_dct, _ = thermfit.prepare_refs(
        model_dct['therm_fit']['ref_scheme'],
        spc_dct, (spc_name,))
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
    return [*coeff_array], spc_array


def energy(
    spc_name, spc_dct_i, spc_mod_dct_i, proc_keyword_dct,
    thy_dct, locs, locs_path,
    cnf_fs, run_prefix, save_prefix):
    """ collect energy
    """
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
    return [locs_path, energy]


def enthalpy(
        spc_name, spc_dct, spc_dct_i, spc_mod_dct_i, model_dct,
        chn_basis_ene_dct, spc_array, locs, locs_path,
        cnf_fs, run_prefix, save_prefix):
    """ collect enthalpies
    """
    print('spc_mod_dct_i', spc_mod_dct_i)
    pf_filesystems = filesys.models.pf_filesys(
        spc_dct_i, spc_mod_dct_i,
        run_prefix, save_prefix, saddle=False)
    print(pf_filesystems)
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
    return [locs_path, ene_abs, hf0k, *coeff_array], chn_basis_ene_dct, spc_array
