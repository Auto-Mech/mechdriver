""" calculates certain quantities of interest using MESS+filesytem
"""

import os
import automol
import autofile
from routines.pf.models import typ
from routines.pf.models import _tors as tors
from routines.pf.models import _vib as vib
from routines.pf.models import _fs as fs
from routines.pf.models import _util as util
from lib.phydat import phycon
from lib.filesys import inf as finf
from lib.amech_io import parser


# Functions to hand reading and formatting energies of single species
def read_energy(spc_dct_i, pf_filesystems, pf_models, pf_levels,
                read_ene=True, read_zpe=True):
    """ Get the energy for a species on a channel
    """

    # Read the electronic energy and ZPVE
    e_elec = None
    if read_ene:
        e_elec = electronic_energy(
            spc_dct_i, pf_filesystems, pf_levels)

    e_zpe = None
    if read_zpe:
        e_zpe = zero_point_energy(
            spc_dct_i, pf_filesystems, pf_models, saddle=False)

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


def electronic_energy(spc_dct_i, pf_filesystems, pf_levels):
    """ get high level energy at low level optimized geometry
    """

    print('- Calculating electronic energy')

    # spc_dct_i = spc_dct[spc_name]
    spc_info = finf.get_spc_info(spc_dct_i)

    # Get the harmonic filesys information
    [_, cnf_path, _, _, _] = pf_filesystems['harm']

    # Get the electronic energy levels
    ene_levels = pf_levels['ene'][1]

    # Read the energies from the filesystem
    e_elec = None
    if os.path.exists(cnf_path):

        e_elec = 0.0
        for (coeff, level) in ene_levels:
            # Build SP filesys
            mod_thy_info = finf.modify_orb_restrict(spc_info, level)
            sp_save_fs = autofile.fs.single_point(cnf_path)
            sp_save_fs[-1].create(mod_thy_info[1:4])
            # Read the energy
            sp_path = sp_save_fs[-1].path(mod_thy_info[1:4])
            if os.path.exists(sp_path):
                print('Energy read from path {}'.format(sp_path))
                ene = sp_save_fs[-1].file.energy.read(mod_thy_info[1:4])
                e_elec += (coeff * ene)
            else:
                print('No energy at path')
                e_elec = None
                break
    else:
        print('No conformer to calculate the energy')

    return e_elec


def zero_point_energy(spc_dct_i,
                      pf_filesystems, pf_models, saddle=False):
    """ compute the ZPE including torsional and anharmonic corrections
    """

    print('- Calculating zero-point energy')

    # spc_dct_i = spc_dct[spc_name]
    frm_bnd_key, brk_bnd_key = util.get_bnd_keys(spc_dct_i, saddle)
    ts_bnd = util.set_ts_bnd(spc_dct_i, saddle)
    rxn_class = util.set_rxn_class(spc_dct_i, saddle)

    # Calculate ZPVE
    if typ.is_atom(spc_dct_i):
        zpe = 0.0
    else:
        rtr_names, rtr_grids, rtr_syms, const_dct, ref_ene = tors.rotor_info(
            spc_dct_i, pf_filesystems, pf_models,
            frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key)

        if typ.nonrigid_tors(pf_models, rtr_names):
            mess_hr_str, prot_hr_str, _, _ = tors.make_hr_strings(
                rtr_names, rtr_grids, rtr_syms, const_dct,
                ref_ene, pf_filesystems, pf_models,
                rxn_class, ts_bnd,
                saddle=saddle, tors_wgeo=True)

        # Obtain vibration partition function information
        if typ.nonrigid_tors(pf_models, rtr_names):
            _, _, zpe = vib.tors_projected_freqs_zpe(
                pf_filesystems, mess_hr_str, prot_hr_str)
        else:
            _, _, zpe = vib.read_harmonic_freqs(pf_filesystems)

    return zpe


# Functions to handle energies for a channel
def set_reference_ene(rxn_lst, spc_dct, thy_dct, model_dct,
                      run_prefix, save_prefix, ref_idx=0):
    """ Sets the reference species for the PES for which all energies
        are scaled relative to.
    """

    # Set the index for the reference species, right now defualt to 1st spc
    ref_rgts = rxn_lst[ref_idx]['reacs']
    ref_model = rxn_lst[ref_idx]['model'][1]
    print('\nDetermining the reference energy for PES...')
    print(' - Reference species assumed to be the',
          ' first set of reactants on PES: {}'.format('+'.join(ref_rgts)))
    print(' - Model for Reference Species: {}'.format(ref_model))

    # Get the model for the first reference species
    pf_levels = parser.model.pf_level_info(
        model_dct[ref_model]['es'], thy_dct)
    pf_models = parser.model.pf_model_info(
        model_dct[ref_model]['pf'])
    ref_ene_level = pf_levels['ene'][0]
    print(' - Energy Level for Reference Species: {}'.format(ref_ene_level))

    # Get the elec+zpe energy for the reference species
    print('')
    ref_ene = 0.0
    for rgt in ref_rgts:

        print(' - Calculating energy for {}...'.format(rgt))

        # Build filesystem
        pf_filesystems = fs.pf_filesys(
            spc_dct[rgt], pf_levels, run_prefix, save_prefix, saddle=False)

        # Calcualte the total energy
        ref_ene += read_energy(
            spc_dct[rgt], pf_filesystems, pf_models, pf_levels,
            read_ene=True, read_zpe=True)

    return ref_ene, ref_model


def calc_channel_enes(channel_infs, ref_ene,
                      chn_model, first_ground_model):
    """ Get the energies for several points on the reaction channel.
        The energy is determined by two different methods:
            (1) Read from the file system if chn_model == first_ground_model
            (2) Shift ene for the channel if chn_model != first_ground_model
    """

    if chn_model == first_ground_model:
        chn_enes = sum_enes(channel_infs, ref_ene, ene_lvl='ene_chnlvl')
    else:
        chn_enes1 = sum_enes(channel_infs, ref_ene, ene_lvl='ene_reflvl')
        chn_enes2 = sum_enes(channel_infs, ref_ene, ene_lvl='ene_reflvl')
        chn_enes = shift_enes(chn_enes1, chn_enes2)

    return chn_enes


def sum_enes(channel_infs, ref_ene, ene_lvl='ene_chnlvl'):
    """ sum the energies
    """
    # Calculate energies for species
    reac_ene = 0.0
    for rct in channel_infs['reacs']:
        reac_ene += rct[ene_lvl]
    prod_ene = 0.0
    for prd in channel_infs['prods']:
        prod_ene += prd[ene_lvl]

    # Set energy for inner transition state
    ts_ene = channel_infs['ts'][ene_lvl]

    # Initialize energy dct
    sum_ene = {
        'reacs': reac_ene,
        'prods': prod_ene,
        'ts': ts_ene
    }

    # Calculate energies for fake entrance- and exit-channel wells
    if 'fake_vdwr' in channel_infs:
        vdwr_ene = reac_ene - (1.0 * phycon.KCAL2EH)
        sum_ene.update(
            {'fake_vdwr': vdwr_ene, 'fake_vdwr_ts': reac_ene}
        )
    if 'fake_vdwp' in channel_infs:
        vdwp_ene = prod_ene - (1.0 * phycon.KCAL2EH)
        sum_ene.update(
            {'fake_vdwp': vdwp_ene, 'fake_vdwp_ts': reac_ene}
        )

    # Scale all of the energies in the dict
    for spc, ene in sum_ene.items():
        sum_ene[spc] = (ene - ref_ene) * phycon.EH2KCAL

    return sum_ene


def shift_enes(chn_enes1, chn_enes2):
    """ When two channels dont match, the energies need to be shifted
        to bring them into alignment.
    """

    # Find a species that has enes with both methods to be used to scale
    # I don't think we need to use any species, so I will use the first
    for spc in chn_enes1:
        if chn_enes1[spc] is not None and chn_enes2[spc] is not None:
            scale_ref_spcs = spc
            break
    scale_ref_ene1 = chn_enes1[scale_ref_spcs]
    scale_ref_ene2 = chn_enes2[scale_ref_spcs]

    # Now return a dct with the lvl1 enes or the scaled lvl2 enes
    fin_enes = {}
    for spc in chn_enes1:
        if chn_enes1[spc] is not None:
            fin_enes[spc] = chn_enes1[spc]
        else:
            fin_enes[spc] = scale_ref_ene1 + (chn_enes2[spc] - scale_ref_ene2)

    return fin_enes


# Writer
def zpe_str(spc_dct, zpe):
    """ return the zpe for a given species according a specified set of
    partition function levels
    """
    if automol.geom.is_atom(automol.inchi.geometry(spc_dct['ich'])):
        zero_energy_str = 'End'
    else:
        zero_energy_str = ' ZeroEnergy[kcal/mol] ' + str(zpe)
        zero_energy_str += '\nEnd'

    return zero_energy_str


# OLD SHIFT FUNCTION
# def calc_shift_ene2(spc_dct, spc_tgt, rxn,
#                    thy_dct, model_dct,
#                    chn_model, first_ground_model,
#                    save_prefix, saddle=False):
#    """ Function to shift the energy of a species to allow mixing of
#        channels calculated with two different methods.
#        We consider two levels of theory:
#          (1) method for PES reference species and
#          (2) method the user requested for this point on the PES.
#        Shifted energy will be calculated as:
#          E_t[1] = {E_s[1] + (E_t[2] - E_s[2])} - E_r[1]
#        where E_t = ene of target species missing info at lvl 1,
#              E_s = ene of other species on channel, and
#              E_r = ene of species that is reference for whole PES.
#        Note that this function only returns the part in curly braces and the
#        final E_r[1] term is subtracted in some other part of the code.
#    """
#
#    # Read the energy for the target species from the filesystem
#    tgt_ene_lvl2 = get_fs_ene_zpe(spc_dct, spc_tgt,
#                                  thy_dct, model_dct, chn_model,
#                                  save_prefix, saddle=saddle)
#
#    # Loop over the species in the channel and find one species
#    # where the energy and ZPVE has been calculated at levels 1 and 2.
#    chn_enes = {}
#    for spc in rxn:
#        # Try and read the energies from the filesystem
#        chn_ene1 = get_fs_ene_zpe(spc_dct, spc,
#                                  thy_dct, model_dct, first_ground_model,
#                                  save_prefix, saddle=bool('ts_' in spc))
#        chn_ene2 = get_fs_ene_zpe(spc_dct, spc,
#                                  thy_dct, model_dct, chn_model,
#                                  save_prefix, saddle=bool('ts_' in spc))
#        # Only add the energies to both dcts if ene1 and ene2 were found
#        if chn_ene1 and chn_ene2:
#            chn_enes[spc] = (chn_ene1, chn_ene2)
#
#    # Calculate the shifted energy for the species at level 1, if possible
#    if chn_enes:
#        # Get lvl1 and lvl2 enes for 1st spc in dcts (shouldn't matter which)
#        for spc, enes in chn_enes.values():
#            chn_ene_lvl1, chn_ene_lvl2 = enes
#            break
#        # Calculate energy
#        tgt_ene_lvl1 = chn_ene_lvl1 + (chn_ene_lvl2 - tgt_ene_lvl2)
#    else:
#        print('No species on the channel with energies at methods 1 and 2')
#        tgt_ene_lvl1 = None
#
#    return tgt_ene_lvl1