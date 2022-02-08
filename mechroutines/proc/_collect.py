"""Collects the target info
"""

import os
import tempfile

import thermfit
import automol
import autofile
import mess_io
import autorun
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo

from mechlib import filesys
import mechlib.amech_io.printer as ioprinter
from mechlib.amech_io import reader
from mechroutines.models import _rot as rot
from mechroutines.models import _vib as vib
from mechroutines.models import _tors as tors
from mechroutines.models import _symm as symm
from mechroutines.models import ene
from mechroutines.models import blocks
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
        comment = f'energy: {_ene:>15.10f}\n'
        zma_str = automol.zmat.string(zma)
        miss_data = None
    else:
        zma_str = '\t -- Missing --'
        miss_data = (spc_name, mod_thy_info, 'zmatrix')

    spc_data = f'\n\nSPC: {spc_name}\tConf: {locs}\tPath: {locs_path}\n'
    spc_data += comment + zma_str

    return spc_data, miss_data


def molden(spc_name, locs, locs_path, cnf_fs, mod_thy_info):
    """collect a geometry
    """
    if cnf_fs[-1].file.geometry.exists(locs):
        geo = cnf_fs[-1].file.geometry.read(locs)
        sp_fs = autofile.fs.single_point(locs_path)
        _ene = sp_fs[-1].file.energy.read(mod_thy_info[1:4])
        # comment = f'energy: {_ene:>15.10f}\n'
        # comment = f'\n\nSPC: {spc_name}\tConf: {locs}\tPath: {locs_path}\n'
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
            comment = f'energy: {_ene:>15.10f}\n'
            xyz_str = automol.geom.xyz_string(geo, comment=comment)
            miss_data = None
        else:
            comment = 'no energy found'
            xyz_str = automol.geom.xyz_string(geo, comment=comment)
            miss_data = (spc_name, mod_thy_info, 'energy')
    else:
        xyz_str = '\t -- Missing --'
        miss_data = (spc_name, mod_thy_info, 'geometry')

    spc_data = f'\n\nSPC: {spc_name}\tConf: {locs}\tPath: {locs_path}\n'
    spc_data += xyz_str

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
                print(f'Imaginary Frequencies[cm-1]: {imag}')
                freqs = (-1*imag,) + freqs
            miss_data = None

        # Do a TED check
        if zrxn is not None:
            vib.ted(spc_dct_i, pf_filesystems, spc_mod_dct_i,
                    run_prefix, zrxn=zrxn)
    else:
        es_levels = util.freq_es_levels(proc_keyword_dct)
        spc_mod_dct_i = util.generate_spc_model_dct(es_levels, thy_dct)
        ret = vib.read_locs_harmonic_freqs(
            cnf_fs, locs, run_prefix, zrxn=zrxn)
        if ret is not None:
            freqs, imag, zpe, disps = ret
            if freqs and proc_keyword_dct['scale'] is not None:
                freqs, zpe = vib.scale_frequencies(
                    freqs, 0.0, spc_mod_dct_i,
                    scale_method=proc_keyword_dct['scale'])
            if saddle:
                print(f'Imaginary Frequencies[cm-1]: {imag}')
                freqs = (-1*imag,) + freqs
            miss_data = None

        pf_filesystems = filesys.models.pf_filesys(
            spc_dct_i, spc_mod_dct_i,
            run_prefix, save_prefix,
            name=spc_name, saddle=saddle)

        # Do a TED check
        if zrxn is not None:
            vib.ted(spc_dct_i, pf_filesystems, spc_mod_dct_i,
                    run_prefix, zrxn=zrxn)

    # Package up the frequencies data
    if freqs is not None:
        spc_data = [locs_path, zpe, *freqs]
        fxn_ret = spc_data, (torsfreqs, all_freqs, sfactor), disps
    else:
        fxn_ret = None

    return fxn_ret, miss_data


def torsions(spc_name, locs, locs_path, spc_dct_i, spc_mod_dct_i,
             mod_thy_info, run_prefix, save_prefix):
    """ get the torsion potentials
        currently just checks if there any non-empty potentials
    """

    if spc_mod_dct_i is not None:
        mod_thy_info = spc_mod_dct_i['tors']['geolvl'][1][1]

    saddle = 'ts_' in spc_name
    pf_filesystems = filesys.models.pf_filesys(
        spc_dct_i, spc_mod_dct_i,
        run_prefix, save_prefix,
        name=spc_name, saddle=saddle, spc_locs=locs)

    # Do initial check to see if a torsions file exists
    if pf_filesystems['tors'] is not None:
        [cnf_fs, _, min_cnf_locs, _, _] = pf_filesystems['tors']
        zma_fs = autofile.fs.zmatrix(cnf_fs[-1].path(min_cnf_locs))
        zma_path = zma_fs[-1].path([0])
        print(f'Checking for torsions at {zma_path}')
        if zma_fs[-1].file.torsions.exists([0]):
            rotors = tors.build_rotors(
                spc_dct_i, pf_filesystems, spc_mod_dct_i)
            names = automol.rotor.names(rotors, flat=True)
            pots = automol.rotor.potentials(rotors, flat=True)
            # Just check if there any potentials with
            miss_data = None
            for name, pot in zip(names, pots):
                if pot:
                    # pot_str = ' '.join((f'{0:.1f}' for val in pot))
                    # print(f'Rotor {name}: {pot_str}')
                    print(f'Rotor {name}: {pot}')
                else:
                    print(f'Rotor {name}: MISSING POTENIAL')
                    miss_data = (spc_name, mod_thy_info, 'torsions')
        else:
            miss_data = None
            print('No rotors file in SAVE filesys. Assuming no rotors')
    else:
        miss_data = None

    return None, miss_data


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
                ioprinter.reading('Energy', sp_path)
                _ene = sp_save_fs[-1].file.energy.read(
                    mod_thy_info[1:4])
        else:
            ioprinter.warning_message(f'No Energy found at {sp_path}')

    if _ene is not None:
        miss_data = None
    else:
        # miss_data = (spc_name + '_'.join(locs), mod_thy_info, 'energy')
        miss_data = (spc_name, mod_thy_info, 'energy')

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
        name=spc_name, saddle=saddle, spc_locs=locs)
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


def partition_function(
        spc_name, spc_dct_i, spc_mod_dct_i,
        pes_mod_dct_i,
        locs, locs_path,
        cnf_fs, run_prefix, save_prefix):
    """ collect enthalpies
    """

    messpf_inp_str, dat_str_dct, miss_data = messpf_input(
        spc_name, spc_dct_i, spc_mod_dct_i,
        pes_mod_dct_i,
        locs, locs_path,
        cnf_fs, run_prefix, save_prefix)

    # Combine the strings together to create full MESS input file string
    # tempfile.tempdir = "./messpf_temp"
    # file_path = (
    # '/home/elliott/projects/AutoMech/RO2QOOH/all_conformers/all/temp')
    with tempfile.TemporaryDirectory() as file_path:
        autorun.write_input(
            file_path,
            messpf_inp_str,
            aux_dct=dat_str_dct,
            input_name='pf.inp')
        autorun.run_script(
            autorun.SCRIPT_DCT['messpf'],
            file_path)
        pf_arrays = reader.mess.messpf(
            file_path)
    return (
        pf_arrays, miss_data
    )


def messpf_input(
        spc_name, spc_dct_i, spc_mod_dct_i,
        pes_mod_dct_i,
        locs, locs_path,
        cnf_fs, run_prefix, save_prefix):
    """ collect enthalpies
    """

    _, _ = locs_path, cnf_fs  # needed

    zrxn = spc_dct_i.get('zrxn')
    saddle = bool(zrxn)
    miss_data = None

    # print('spc_mod_dct_i', spc_mod_dct_i)
    pf_filesystems = filesys.models.pf_filesys(
        spc_dct_i, spc_mod_dct_i,
        run_prefix, save_prefix,
        name=spc_name, saddle=saddle, spc_locs=locs)
    geom = rot.read_geom(pf_filesystems)
    rotors = tors.build_rotors(
        spc_dct_i, pf_filesystems, spc_mod_dct_i)
    freqs, imag, zpe, _, tors_strs, _, _, _ = vib.full_vib_analysis(
        spc_dct_i, pf_filesystems, spc_mod_dct_i,
        run_prefix, zrxn=zrxn)
    allr_str = tors_strs[0]

    zma = None
    sym_factor = symm.symmetry_factor(
        pf_filesystems, spc_mod_dct_i, spc_dct_i, rotors, grxn=zrxn, zma=zma)
    elec_levels = spc_dct_i['elec_levels']

    keys = ['writer', 'geom', 'sym_factor', 'freqs', 'imag', 'elec_levels',
            'mess_hr_str', 'mdhr_dat',
            'xmat', 'rovib_coups', 'rot_dists',
            'ene_chnlvl', 'ene_reflvl', 'zpe_chnlvl', 'ene_tsref',
            'edown_str', 'collid_freq_str']
    vals = ['species_block', geom, sym_factor, freqs, imag, elec_levels,
            allr_str, None,
            None, None, None,
            None, None, zpe, None,
            None, None]
    inf_dct = dict(zip(keys, vals))

    temps = pes_mod_dct_i['therm_temps']
    globkey_str = mess_io.writer.global_pf_input(
        temperatures=temps,
        rel_temp_inc=0.001,
        atom_dist_min=0.6
    )
    mess_writer = getattr(blocks, inf_dct['writer'])
    mess_block, dat_str_dct = mess_writer(inf_dct)
    if inf_dct['writer'] == 'tau_block':
        zero_energy = None
    else:
        zero_energy = inf_dct['zpe_chnlvl']

    spc_str = mess_io.writer.species(
        spc_label=spc_name,
        spc_data=mess_block,
        zero_ene=zero_energy
    )
    messpf_inp_str = mess_io.writer.messpf_inp_str(globkey_str, spc_str)

    return (messpf_inp_str, dat_str_dct, miss_data)
