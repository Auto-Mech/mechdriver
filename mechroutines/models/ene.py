""" calculates certain quantities of interest using MESS+filesytem
"""

import os
import automol
import autofile
from mechanalyzer.inf import rxn as rinfo
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo
from mechlib.amech_io import printer as ioprinter
from mechroutines.models import typ
from mechroutines.models import _vib as vib


# Functions to hand reading and formatting energies of single species
def read_energy(spc_dct_i, pf_filesystems,
                spc_model_dct_i, run_prefix,
                read_ene=True, read_zpe=True, conf=None, saddle=False):
    """ Get the energy for a species on a channel
    """

    # Read the electronic energy and ZPVE
    e_elec = None
    if read_ene:
        e_elec = electronic_energy(
            spc_dct_i, pf_filesystems, spc_model_dct_i, conf=conf)
        ioprinter.debug_message('e_elec in models ene ', e_elec)

    e_zpe = None
    if read_zpe:
        e_zpe = zero_point_energy(
            spc_dct_i, pf_filesystems, spc_model_dct_i,
            run_prefix, saddle=saddle)
        ioprinter.debug_message('zpe in models ene ', e_zpe)

    # Return the total energy requested
    ene = None
    if read_ene and read_zpe:
        if e_elec is not None and e_zpe is not None:
            ene = e_elec + e_zpe
    elif read_ene and not read_zpe:
        ene = e_elec
    elif read_ene and not read_zpe:
        ene = e_zpe

    return ene


def electronic_energy(spc_dct_i, pf_filesystems, spc_model_dct_i, conf=None):
    """ get high level energy at low level optimized geometry
    """

    ioprinter.info_message('- Calculating electronic energy')

    # spc_dct_i = spc_dct[spc_name]
    rxn_info = spc_dct_i.get('rxn_info', None)
    if rxn_info is not None:
        spc_info = rinfo.ts_info(rxn_info)
    else:
        spc_info = sinfo.from_dct(spc_dct_i)

    # Get the harmonic filesys information
    if conf:
        cnf_path = conf[1]
    else:
        [_, cnf_path, _, _, _] = pf_filesystems['harm']

    # Get the electronic energy levels
    ene_levels = tuple(val[1] for key, val in spc_model_dct_i['ene'].items()
                       if 'lvl' in key)
    print('ene levels', ene_levels)

    # Read the energies from the filesystem
    e_elec = None
    if os.path.exists(cnf_path):

        e_elec = 0.0
        # ioprinter.info_message('lvls', ene_levels)
        for (coeff, level) in ene_levels:
            # Build SP filesys
            mod_thy_info = tinfo.modify_orb_label(level, spc_info)
            sp_save_fs = autofile.fs.single_point(cnf_path)
            sp_save_fs[-1].create(mod_thy_info[1:4])
            # Read the energy
            sp_path = sp_save_fs[-1].path(mod_thy_info[1:4])
            if os.path.exists(sp_path):
                ioprinter.reading('Energy', sp_path)
                ene = sp_save_fs[-1].file.energy.read(mod_thy_info[1:4])
                e_elec += (coeff * ene)
            else:
                ioprinter.warning_message('No energy at path')
                e_elec = None
                break
    else:
        ioprinter.warning_message('No conformer to calculate the energy')

    return e_elec


def zero_point_energy(spc_dct_i,
                      pf_filesystems, spc_model_dct_i,
                      run_prefix, saddle=False):
    """ compute the ZPE including torsional and anharmonic corrections
    """

    ioprinter.info_message('- Calculating zero-point energy')

    # Calculate ZPVE
    is_atom = False
    if not saddle:
        if typ.is_atom(spc_dct_i):
            is_atom = True
    if is_atom:
        zpe = 0.0
    else:
        _, _, zpe, _ = vib.vib_analysis(
            spc_dct_i, pf_filesystems, spc_model_dct_i,
            run_prefix, zrxn=(None if not saddle else 'placeholder'))

    return zpe


def rpath_ref_idx(ts_dct, scn_vals, coord_name, scn_prefix,
                  ene_info1, ene_info2):
    """ Get the reference energy along a reaction path
    """

    # Set up the filesystem
    zma_fs = autofile.fs.zmatrix(scn_prefix)
    zma_path = zma_fs[-1].path([0])
    scn_fs = autofile.fs.scan(zma_path)

    ene_info1 = ene_info1[1][0][1]
    ene_info2 = ene_info2[0]
    ioprinter.debug_message('mod_eneinf1', ene_info1)
    ioprinter.debug_message('mod_eneinf2', ene_info2)
    mod_ene_info1 = tinfo.modify_orb_label(
        sinfo.from_dct(ts_dct), ene_info1)
    mod_ene_info2 = tinfo.modify_orb_label(
        sinfo.from_dct(ts_dct), ene_info2)

    ene1, ene2, ref_val = None, None, None
    for val in reversed(scn_vals):
        locs = [[coord_name], [val]]
        path = scn_fs[-1].path(locs)
        hs_fs = autofile.fs.high_spin(path)
        if hs_fs[-1].file.energy.exists(mod_ene_info1[1:4]):
            ene1 = hs_fs[-1].file.energy.read(mod_ene_info1[1:4])
        if hs_fs[-1].file.energy.exists(mod_ene_info2[1:4]):
            ene2 = hs_fs[-1].file.energy.read(mod_ene_info2[1:4])
        if ene1 is not None and ene2 is not None:
            ref_val = val
            break

    if ref_val is not None:
        scn_idx = scn_vals.index(ref_val)

    return scn_idx, ene1, ene2


# Writer
def zpe_str(spc_dct, zpe):
    """ return the zpe for a given species according a specified set of
    partition function levels
    """
    if automol.geom.is_atom(automol.inchi.geometry(spc_dct['inchi'])):
        zero_energy_str = 'End'
    else:
        zero_energy_str = ' ZeroEnergy[kcal/mol] ' + str(zpe)
        zero_energy_str += '\nEnd'

    return zero_energy_str
