""" calculates certain quantities of interest using MESS+filesytem
"""

import os
import automol
import elstruct
import autofile

# New Libs
from lib.phydat import phycon
from lib.filesystem import minc as fsmin
from lib.filesystem import orb as fsorb
from lib.filesystem import inf as finf
from lib.load import model as loadmodel
import routines.pf.messf.models as pfmodels
from routines.pf.messf.blocks import set_model_filesys
from routines.pf.messf import _sym as sym
from routines.pf.messf import _util as messfutil
from routines.pf.messf import _tors as tors


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


def get_zero_point_energy(spc, spc_dct_i, pf_levels, spc_model,
                          save_prefix, saddle=False):
    """ compute the ZPE including torsional and anharmonic corrections
    """
    
    # Calculate ZPVE
    if util.is_atom(spc_dct_i):
        zpve = 0.0
    else:
        # Build dct combinining various information from the filesys and MESS
        inf_dct = models.read_filesys_for_spc(
            spc_dct_i, rxn, spc_model, pf_levels, save_prefix, saddle=saddle)

    return zpve


def get_high_level_energy(
        spc_info, thy_low_level, thy_high_level, save_prefix, saddle=False):
    """ get high level energy at low level optimized geometry
    """
    if saddle:
        spc_save_path = save_prefix
    else:
        spc_save_fs = autofile.fs.species(save_prefix)
        spc_save_fs[-1].create(spc_info)
        spc_save_path = spc_save_fs[-1].path(spc_info)

    thy_low_level = fsorb.mod_orb_restrict(
        spc_info, thy_low_level)

    ll_save_fs = autofile.fs.theory(spc_save_path)
    ll_save_path = ll_save_fs[-1].path(thy_low_level[1:4])

    if os.path.exists(ll_save_path):
        if saddle:
            ll_save_fs = autofile.fs.ts(ll_save_path)
            ll_save_fs[0].create()
            ll_save_path = ll_save_fs[0].path()

        cnf_save_fs = autofile.fs.conformer(ll_save_path)
        min_cnf_locs = fsmin.min_energy_conformer_locators(
            cnf_save_fs)
        if not min_cnf_locs:
            print('ERROR: No minimum conformer geometry for ',
                  'this species {}'.format(spc_info[0]))
            return 0.0
        cnf_save_path = cnf_save_fs[-1].path(min_cnf_locs)

        thy_high_level = fsorb.mod_orb_restrict(
            spc_info, thy_high_level)

        sp_save_fs = autofile.fs.single_point(cnf_save_path)
        sp_save_fs[-1].create(thy_high_level[1:4])

        if os.path.exists(sp_save_fs[-1].path(thy_high_level[1:4])):
            min_ene = sp_save_fs[-1].file.energy.read(thy_high_level[1:4])
        else:
            # print('No energy at path')
            min_ene = None
    else:
        # print('No energy at path')
        min_ene = None

    return min_ene


def set_reference_ene(rxn_lst, spc_dct, thy_dct, model_dct, save_prefix):
    """ Sets the reference species for the PES for which all energies
        are scaled relative to.
    """
    # Get the model for the first reference species
    ref_model = rxn_lst[0]['model'][1]

    # Get the elec+zpe energy for the reference species
    ref_ene = 0.0
    first_spc = rxn_lst[0]['reacs']
    for rct in first_spc:
        ref_ene += get_fs_ene_zpe(
            spc_dct, rct,
            thy_dct, model_dct, ref_model,
            save_prefix, saddle=False)

    return ref_ene, ref_model


def calc_channel_enes(spc_dct, rxn, tsname,
                      thy_dct, model_dct,
                      chn_model, first_ground_model,
                      save_prefix):
    """ Get the energies for several points on the reaction channel.
        The energy is determined by two different methods:
            (1) Read from the file system if chn_model == first_ground_model
            (2) Shift ene for the channel if chn_model != first_ground_model
    """

    reacs = rxn['reacs']
    prods = rxn['prods']
    trans_st = tsname
    species = [reacs, prods, trans_st]

    if chn_model == first_ground_model:
        chn_enes = read_channel_energies(spc_dct, species,
                                         thy_dct, model_dct,
                                         first_ground_model,
                                         save_prefix,)
    else:
        chn_enes1 = read_channel_energies(spc_dct, species,
                                          thy_dct, model_dct,
                                          first_ground_model,
                                          save_prefix)
        chn_enes2 = read_channel_energies(spc_dct, species,
                                          thy_dct, model_dct,
                                          chn_model,
                                          save_prefix)
        chn_enes = calc_shift_ene(chn_enes1, chn_enes2)

    return chn_enes


def calc_shift_ene(chn_enes1, chn_enes2):
    """ Another function to do the shifted energies
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


def read_channel_energies(spc_dct, species,
                          thy_dct, model_dct, model,
                          save_prefix):
    """ Read the energies for the channel
    """

    # Unpack the species list
    [reacs, prods, trans_st] = species

    # Read the reactants energies
    reac_ene = 0.0
    for reac in reacs:
        ene = get_fs_ene_zpe(spc_dct, reac,
                             thy_dct, model_dct, model,
                             save_prefix, saddle=False,
                             read_ene=True, read_zpe=False)
        if ene is not None:
            reac_ene += ene
        else:
            reac_ene = None
            break

    # Read the products energies
    prod_ene = 0.0
    for prod in prods:
        ene = get_fs_ene_zpe(spc_dct, prod,
                             thy_dct, model_dct, model,
                             save_prefix, saddle=False,
                             read_ene=True, read_zpe=False)
        if ene is not None:
            prod_ene += ene
        else:
            prod_ene = None
            break

    # Read the transition state energy
    ts_ene = get_fs_ene_zpe(spc_dct, trans_st,
                            thy_dct, model_dct, model,
                            save_prefix, saddle=False,
                            read_ene=True, read_zpe=False)

    # Set energies into the dct
    energy_dct = {}
    energy_dct['reacs'] = reac_ene
    energy_dct['prods'] = prod_ene
    energy_dct['ts'] = ts_ene

    # Add ZPEs if they are available
    if zpe_dct is not None:
        energy_dct['reacs'] += zpe_dct['reacs']
        energy_dct['prods'] += zpe_dct['prods']
        energy_dct['ts'] += zpe_dct['ts']

    return energy_dct


def get_fs_ene_zpe(spc_dct, spc,
                   thy_dct, model_dct, model,
                   save_prefix, saddle=False,
                   read_ene=True, read_zpe=True):
    """ Get the energy for a species on a channel
    """

    # Set the species save filesystem
    spc_save_fs = autofile.fs.species(save_prefix)

    # Set the spc info object
    spc_info = (spc_dct[spc]['ich'],
                spc_dct[spc]['chg'],
                spc_dct[spc]['mul'])

    # Set the model and info for the reaction
    pf_levels = loadmodel.set_es_model_info(model_dct[model]['es'], thy_dct)
    pf_model = loadmodel.set_pf_model_info(model_dct[model]['pf'])

    # Set paths
    if saddle:
        spc_save_path = spc_dct[spc]['rxn_fs'][3]
        save_path = spc_save_path
    else:
        spc_save_fs[-1].create(spc_info)
        spc_save_path = spc_save_fs[-1].path(spc_info)
        save_path = save_prefix

    # Read the electronic energy and ZPVE
    thy_low_level = pf_levels[0]
    thy_high_levels = pf_levels[1]
    e_elec = None
    e_zpe = None
    if read_ene:
        e_elec = 0.0
        for (coeff, level) in thy_high_levels:
            high_ene = get_high_level_energy(
                spc_info=spc_info,
                thy_low_level=thy_low_level,
                thy_high_level=level,
                save_prefix=save_path,
                saddle=saddle)
            e_elec += coeff * high_ene
    # if e_elec is not None:
    if read_zpe:
        e_zpe = get_zero_point_energy(
            spc, spc_dct[spc],
            pf_levels, pf_model,
            save_prefix=spc_save_path)
        e_zpe /= phycon.EH2KCAL

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
