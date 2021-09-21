"""Collects the target info
"""

import os

import thermfit
import automol
import autofile
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo

from mechlib import filesys
from mechlib.amech_io.printer import reading
from mechroutines.models import _vib as vib
from mechroutines.models import _tors as tors
from mechroutines.models import ene
from mechroutines.thermo import basis
from mechroutines.proc import _util as util


def zmatrix(spc_name, locs, locs_path, cnf_fs, mod_thy_info):
    """collect a zmatrix
    """

    if cnf_fs[-1].file.geometry.exists(locs):
        geo = cnf_fs[-1].file.geometry.read(locs)
        zma = automol.geom.zmatrix(geo)
        sp_fs = autofile.fs.single_point(locs_path)
        _ene = sp_fs[-1].file.energy.read(mod_thy_info[1:4])
        comment = 'energy: {0:>15.10f}\n'.format(_ene)
        zma_str = automol.zmat.string(zma)
        miss_data = None
    else:
        zma_str = '\t -- Missing --'
        miss_data = (spc_name, mod_thy_info, 'zmatrix')

    spc_data = '\n\nSPC: {}\tConf: {}\tPath: {}\n'.format(
        spc_name, locs, locs_path) + comment + zma_str

    return spc_data, miss_data


def molden(spc_name, locs, locs_path, cnf_fs, mod_thy_info):
    """collect a geometry
    """
    if cnf_fs[-1].file.geometry.exists(locs):
        geo = cnf_fs[-1].file.geometry.read(locs)
        sp_fs = autofile.fs.single_point(locs_path)
        _ene = sp_fs[-1].file.energy.read(mod_thy_info[1:4])
        comment = 'energy: {0:>15.10f}'.format(_ene)
        comment += 'SPC: {}\tConf: {}\tPath: {}'.format(
            spc_name, locs, locs_path)
        xyz_str = automol.geom.xyz_string(geo, comment=comment)
        miss_data = None
    else:
        xyz_str = '\t -- Missing --'
        miss_data = (spc_name, mod_thy_info, 'geometry')

    spc_data = xyz_str
    return spc_data, miss_data


def geometry(spc_name, locs, locs_path, cnf_fs, mod_thy_info):
    """collect a geometry
    """
    if cnf_fs[-1].file.geometry.exists(locs):
        geo = cnf_fs[-1].file.geometry.read(locs)
        sp_fs = autofile.fs.single_point(locs_path)
        if sp_fs[-1].file.energy.exists(mod_thy_info[1:4]):
            _ene = sp_fs[-1].file.energy.read(mod_thy_info[1:4])
            comment = 'energy: {0:>15.10f}'.format(_ene)
            xyz_str = automol.geom.xyz_string(geo, comment=comment)
            miss_data = None
        else:
            comment = 'no energy found'
            xyz_str = automol.geom.xyz_string(geo, comment=comment)
            miss_data = (spc_name, mod_thy_info, 'energy')
    else:
        xyz_str = '\t -- Missing --'
        miss_data = (spc_name, mod_thy_info, 'geometry')

    spc_data = '\n\nSPC: {}\tConf: {}\tPath: {}\n'.format(
        spc_name, locs, locs_path) + xyz_str

    return spc_data, miss_data


def frequencies(
        spc_name, spc_dct_i, spc_mod_dct_i,
        proc_keyword_dct, thy_dct,
        cnf_fs, locs, locs_path, run_prefix, save_prefix):
    """collect frequencies
    """

    # Initialize the data objects to None
    freqs = None
    imag = None
    zpe = None
    sfactor = None
    torsfreqs = None
    all_freqs = None
    disps = None

    # Initialize a miss_data object that will be overwritten if data found
    if spc_mod_dct_i is not None:
        mod_thy_info = spc_mod_dct_i['vib']['geolvl'][1][1]
    else:
        mod_thy_info = tinfo.from_dct(thy_dct.get(
            proc_keyword_dct['proplvl']))

    miss_data = (spc_name, mod_thy_info, 'frequencies')

    # Get flags to to ID spc as a transiion state
    zrxn = spc_dct_i.get('zrxn', None)
    saddle = bool(zrxn)

    # Get vibrational frequencies
    if spc_mod_dct_i is not None:
        pf_filesystems = filesys.models.pf_filesys(
            spc_dct_i, spc_mod_dct_i,
            run_prefix, save_prefix,
            name=spc_name, saddle=saddle)

        ret = vib.full_vib_analysis(
            spc_dct_i, pf_filesystems, spc_mod_dct_i,
            run_prefix, zrxn=zrxn)
        if ret is not None:
            freqs, imag, zpe, sfactor, _, torsfreqs, all_freqs, disps = ret
            if saddle:
                print('Imaginary Frequencies[cm-1]: {}'.format(imag))
                freqs = (-1*imag,) + freqs
            miss_data = None
    else:
        es_levels = util.freq_es_levels(proc_keyword_dct)
        spc_mod_dct_i = util.generate_spc_model_dct(es_levels, thy_dct)
        ret = vib.read_locs_harmonic_freqs(
            cnf_fs, locs, run_prefix, zrxn=zrxn)
        if ret is not None:
            freqs, imag, zpe, disps = ret
            if freqs and proc_keyword_dct['scale'] is not None:
                freqs, zpe = vib.scale_frequencies(
                    freqs, 0.0, spc_mod_dct_i, scale_method='3c')
            if saddle:
                print('Imaginary Frequencies[cm-1]: {}'.format(imag))
                freqs = (-1*imag,) + freqs
            miss_data = None

    # Package up the frequencies data
    if freqs is not None:
        spc_data = [locs_path, zpe, *freqs]
        fxn_ret = spc_data, (torsfreqs, all_freqs, sfactor), disps
    else:
        fxn_ret = None

    return fxn_ret, miss_data


def torsions(spc_name, spc_dct_i, spc_mod_dct_i,
             run_prefix, save_prefix):
    """ get the torsion potentials
        currently just checks if there any non-empty potentials
    """

    mod_thy_info = spc_mod_dct_i['tors']['geolvl'][1][1]

    saddle = 'ts_' in spc_name
    pf_filesystems = filesys.models.pf_filesys(
        spc_dct_i, spc_mod_dct_i,
        run_prefix, save_prefix,
        name=spc_name, saddle=saddle)

    # Do initial check to see if a torsions file exists
    if pf_filesystems['tors'] is not None:
        [cnf_fs, _, min_cnf_locs, _, _] = pf_filesystems['tors']
        zma_fs = autofile.fs.zmatrix(cnf_fs[-1].path(min_cnf_locs))
        print('cnf_path', cnf_fs[-1].path(min_cnf_locs))
        if zma_fs[-1].file.torsions.exists([0]):
            rotors = tors.build_rotors(
                spc_dct_i, pf_filesystems, spc_mod_dct_i)
            pots = automol.rotor.potentials(rotors, flat=True)
            # Just check if there any potentials with
            if all(pot for pot in pots):
                miss_data = None
            else:
                miss_data = (spc_name, mod_thy_info, 'torsions')
        else:
            miss_data = None
            print('No rotors file in SAVE filesys. Assuming no rotors')
    else:
        miss_data = None

    return miss_data


def coeffs(spc_name, spc_dct, model_dct, spc_array):
    """ get the heat of formation reference molecules for one species.
    """

    spc_dct_i = spc_dct[spc_name]
    zrxn = spc_dct_i.get('zrxn')

    basis_dct = thermfit.prepare_basis(
        model_dct['therm_fit']['ref_scheme'],
        spc_dct, (spc_name,),
        zrxn=zrxn)
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


def energy(spc_name, spc_dct_i,
           spc_mod_dct_i, proc_keyword_dct,
           thy_dct, locs, locs_path,
           cnf_fs, run_prefix, save_prefix):
    """ collect energy
    """

    saddle = 'ts_' in spc_name

    _ene = None
    if spc_mod_dct_i:
        pf_filesystems = filesys.models.pf_filesys(
            spc_dct_i, spc_mod_dct_i,
            run_prefix, save_prefix,
            name=spc_name, saddle=saddle)
        _ene = ene.electronic_energy(
            spc_dct_i, pf_filesystems, spc_mod_dct_i,
            conf=(locs, locs_path, cnf_fs))
    else:
        spc_info = sinfo.from_dct(spc_dct_i)
        thy_info = tinfo.from_dct(thy_dct.get(
            proc_keyword_dct['proplvl']))
        mod_thy_info = tinfo.modify_orb_label(
            thy_info, spc_info)
        sp_save_fs = autofile.fs.single_point(locs_path)
        # Read the energy
        sp_path = sp_save_fs[-1].path(mod_thy_info[1:4])
        if os.path.exists(sp_path):
            if sp_save_fs[-1].file.energy.exists(
                    mod_thy_info[1:4]):
                reading('Energy', sp_path)
                _ene = sp_save_fs[-1].file.energy.read(
                    mod_thy_info[1:4])

    if _ene is not None:
        miss_data = None
    else:
        miss_data = (spc_name + '_'.join(locs), mod_thy_info, 'energy')

    return [locs_path, _ene], miss_data


def enthalpy(
        spc_name, spc_dct, spc_dct_i, spc_mod_dct_i, model_dct,
        chn_basis_ene_dct, spc_array, locs, locs_path,
        cnf_fs, run_prefix, save_prefix):
    """ collect enthalpies
    """

    zrxn = spc_dct_i.get('zrxn')
    saddle = bool(zrxn)

    # print('spc_mod_dct_i', spc_mod_dct_i)
    pf_filesystems = filesys.models.pf_filesys(
        spc_dct_i, spc_mod_dct_i,
        run_prefix, save_prefix,
        name=spc_name, saddle=saddle)
    # print(pf_filesystems)
    ene_abs = ene.read_energy(
        spc_dct_i, pf_filesystems, spc_mod_dct_i,
        run_prefix, conf=(locs, locs_path, cnf_fs),
        read_ene=True, read_zpe=True, saddle=saddle)
    hf0k, _, chn_basis_ene_dct, hbasis = basis.enthalpy_calculation(
        spc_dct, spc_name, ene_abs,
        chn_basis_ene_dct, model_dct, spc_mod_dct_i,
        run_prefix, save_prefix, pforktp='pf', zrxn=zrxn)
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
    return (
        [locs_path, ene_abs, hf0k, *coeff_array],
        chn_basis_ene_dct,
        spc_array
    )
