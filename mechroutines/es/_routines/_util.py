""" utilites
"""

import numpy
import automol
from mechlib.amech_io.printer import info_message
from automol.extern import Ring_Reconstruction as RR


def calc_nsamp(tors_names, nsamp_par, zma, zrxn=None):
    """ Determine the number of samples to od
    """
    ntaudof = None
    if not any(tors_names):
        info_message(
            " - No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1
        tors_range_dct = None
    else:
        if zrxn is None:
            gra = automol.zmat.graph(zma)
        else:
            gra = automol.reac.ts_graph(zrxn)

        ntaudof = len(automol.graph.rotational_bond_keys(
            gra, with_ch_rotors=False))
        nsamp = nsamp_init(nsamp_par, ntaudof)

        tors_ranges = tuple((0, 2*numpy.pi) for tors in tors_names)
        tors_range_dct = dict(zip(tors_names, tors_ranges))

    return nsamp, tors_range_dct


def calc_nsampd(cnf_save_fs, cnf_run_fs, rid=None):
    """ Determine the number of samples completed
    """

    if rid is None:
        cnf_save_fs[0].create()
        if cnf_save_fs[0].file.info.exists():
            inf_obj_s = cnf_save_fs[0].file.info.read()
            nsampd = inf_obj_s.nsamp
        elif cnf_run_fs[0].file.info.exists():
            inf_obj_r = cnf_run_fs[0].file.info.read()
            nsampd = inf_obj_r.nsamp
        else:
            nsampd = 0
    else:
        cnf_save_fs[1].create([rid])
        if cnf_save_fs[1].file.info.exists([rid]):
            inf_obj_s = cnf_save_fs[1].file.info.read([rid])
            nsampd = inf_obj_s.nsamp
        elif cnf_run_fs[1].file.info.exists([rid]):
            inf_obj_r = cnf_run_fs[1].file.info.read([rid])
            nsampd = inf_obj_r.nsamp
        else:
            nsampd = 0

    return nsampd


def nsamp_init(nsamp_par, ntaudof):
    """ determine nsamp for given species"""
    if nsamp_par[0]:
        nsamp = min(nsamp_par[1] + nsamp_par[2] * nsamp_par[3]**ntaudof,
                    nsamp_par[4])
        # print('Setting nsamp using formula: min(A+B*C**n')
    else:
        nsamp = nsamp_par[5]
    return nsamp


def ring_samp_zmas(ring_atoms, nsamp_par, n_rings=1):
    """ choose starting number of sample zmas
    """
    ntors = len(ring_atoms) - 3 - 2*(n_rings-1)
    apar, bpar, cpar = nsamp_par[1:4]
    #Set maximum number of initial sampling points per run
    return min(10000,30*(apar + bpar * cpar**ntors))


########### GET INFORMATION ON SUBS ON FIRST AND LAST RING ATOM #############
def subs_analysis(all_ring_atoms,all_ring_atoms_list, ngbs, geo, unconnected_keys):
    """ Generate dicts of neighbours for the unconnected (in the Z-Matrix) atoms of
        a ring structure
    """
    last_ring_at_sub_dct,first_ring_at_sub_dct = {}, {}

    for key_dct, ring_atoms in all_ring_atoms_list.items():
        # Gets position of first substituent on first and last atom of ring (which usually messes up)
        first_ring_at,last_ring_at = unconnected_keys[key_dct][0],unconnected_keys[key_dct][1]
        first_ring_at_ngbs,last_ring_at_ngbs = set(ngbs[first_ring_at]), set(ngbs[last_ring_at])
        first_ring_at_subs = first_ring_at_ngbs.difference(set(all_ring_atoms))
        last_ring_at_subs = last_ring_at_ngbs.difference(set(all_ring_atoms))
        # Get xyz of ring and sub
        coord_ring = [xyz for i,(_,xyz) in enumerate(geo) if i in ring_atoms]
        first_coord_sub = [xyz for i,(_,xyz) in enumerate(geo) if i in first_ring_at_subs]
        last_coord_sub = [xyz for i,(_,xyz) in enumerate(geo) if i in last_ring_at_subs]
        # Call alpha beta calculator
        first_sub_params = [tuple(RR.GetRingSubstituentPosition(coord_ring,coord_sub,-1)
                                  ) for coord_sub in first_coord_sub]
        last_sub_params = [tuple(RR.GetRingSubstituentPosition(coord_ring,coord_sub,-1)
                                 ) for coord_sub in last_coord_sub]
        # Update sub dictionaries
        last_ring_at_sub_dct[key_dct] = {key:value for key,value in zip(
                                        last_ring_at_subs,last_sub_params)}
        first_ring_at_sub_dct[key_dct] = {key:value for key,value in zip(
                                        first_ring_at_subs,first_sub_params)}

    return first_ring_at_sub_dct,last_ring_at_sub_dct


########### FIX POSITIONS OF SUBS OF FIRST AND LAST RING ATOM #############
def fixings_subs_positions(samp_zma, all_ring_atoms_list, coos, unconnected_keys,
                           first_ring_at_sub_dct, last_ring_at_sub_dct, dist_thresh=1.2):
    """ maniuplates the dihedrals of a Z-Matrix to assure that the substituents of the
        unconnected atoms of a ring are not close to ring atoms or to each other
    """
    samp_geo = automol.zmat.geometry(samp_zma)
    # I need to perform the substituents check here AFTER I have built the samp ZMat!
    for key_dct,ring_atoms in all_ring_atoms_list.items():
        non_bonded_first,non_bonded_last = tuple(unconnected_keys[key_dct])
        new_key_dct = {}
        for name, cord in coos.items():
            atm_idxs = cord[0]
            if len(atm_idxs) == 2:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, angstrom=True)
            elif len(atm_idxs) == 3:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, degree=True)
            elif len(atm_idxs) == 4:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, degree=True)
                # First sub of first ring atom
                # Check that I am working on the third ring atom
                if atm_idxs[0] == ring_atoms[2]:
                    change_dh = False
                    for atm in first_ring_at_sub_dct[key_dct]:
                        dist_sub_to_last = automol.geom.distance(
                                samp_geo, atm, non_bonded_last, angstrom=True)
                        if dist_sub_to_last < dist_thresh: 
                            change_dh = True

                        for atm2 in last_ring_at_sub_dct[key_dct]:
                            dist_sub_to_sub = automol.geom.distance(
                                samp_geo, atm,atm2, angstrom=True)
                            if dist_sub_to_sub < dist_thresh: 
                                change_dh = True
                    if change_dh: 
                        new_key_dct[name] -= 60.

        samp_zma = automol.zmat.set_values_by_name(samp_zma, new_key_dct)
    samp_geo = automol.zmat.geometry(samp_zma)

    # I need to perform the substituents check here AFTER AFTER I have built the samp ZMat!
    for key_dct,ring_atoms in all_ring_atoms_list.items():
        new_key_dct = {}
        for name, cord in coos.items():
            atm_idxs = cord[0]
            if len(atm_idxs) == 2:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, angstrom=True)
            elif len(atm_idxs) == 3:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, degree=True)
            elif len(atm_idxs) == 4:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, degree=True)
                # First sub of first ring atom
                # If previous -60 rotation didn't help, rotate 60 on the other direction
                if atm_idxs[0] == ring_atoms[2]:
                    change_dh = False
                    for atm in first_ring_at_sub_dct[key_dct]:
                        dist_sub_to_last = automol.geom.distance(
                                samp_geo, atm,non_bonded_last, angstrom=True)
                        if dist_sub_to_last < dist_thresh: change_dh = True

                        for atm2 in last_ring_at_sub_dct[key_dct]:
                            dist_sub_to_sub = automol.geom.distance(
                                samp_geo, atm,atm2, angstrom=True)
                            if dist_sub_to_sub < dist_thresh: change_dh = True
                    if change_dh: 
                        new_key_dct[name] += 120.
                           
        samp_zma = automol.zmat.set_values_by_name(samp_zma, new_key_dct)
    samp_geo = automol.zmat.geometry(samp_zma)

    # I need to perform the substituents check here AFTER AFTER AFTER I have built the samp ZMat!
    for key_dct,ring_atoms in all_ring_atoms_list.items():
        new_key_dct = {}
        for name, cord in coos.items():
            atm_idxs = cord[0]
            if len(atm_idxs) == 2:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, angstrom=True)
            elif len(atm_idxs) == 3:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, degree=True)
            elif len(atm_idxs) == 4:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, degree=True)

                # Last ring atom first sub. 
                if atm_idxs[0] in last_ring_at_sub_dct[key_dct]:
                    change_dh = False
                    if atm_idxs[0] == min(last_ring_at_sub_dct[key_dct]):
                        for atm in last_ring_at_sub_dct[key_dct]:
                            dist_sub_to_first = automol.geom.distance(samp_geo,
                                                atm, non_bonded_first, angstrom=True)
                            if dist_sub_to_first < dist_thresh: change_dh = True
                            for atm2 in first_ring_at_sub_dct[key_dct]:
                                dist_sub_to_sub = automol.geom.distance(
                                    samp_geo, atm,atm2, angstrom=True)
                                if dist_sub_to_sub < dist_thresh: change_dh = True
                    if change_dh:  
                        new_key_dct[name] -= 60.     
        samp_zma = automol.zmat.set_values_by_name(samp_zma, new_key_dct)
    samp_geo = automol.zmat.geometry(samp_zma)

    for key_dct,ring_atoms in all_ring_atoms_list.items():
        new_key_dct = {}
        for name, cord in coos.items():
            atm_idxs = cord[0]
            if len(atm_idxs) == 2:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, angstrom=True)
            elif len(atm_idxs) == 3:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, degree=True)
            elif len(atm_idxs) == 4:
                new_key_dct[name] = automol.zmat.value(samp_zma, name, degree=True)

                # Last ring atom first sub. 
                if atm_idxs[0] in last_ring_at_sub_dct[key_dct]: 
                    change_dh = False
                    if atm_idxs[0] == min(last_ring_at_sub_dct[key_dct]):
                        for atm in last_ring_at_sub_dct[key_dct]:
                            dist_sub_to_first = automol.geom.distance(samp_geo,
                                                atm, non_bonded_first, angstrom=True)
                            if dist_sub_to_first < dist_thresh: change_dh = True
                            for atm2 in first_ring_at_sub_dct[key_dct]:
                                dist_sub_to_sub = automol.geom.distance(
                                    samp_geo, atm,atm2, angstrom=True)
                                if dist_sub_to_sub < dist_thresh: change_dh = True
                    if change_dh:  
                        new_key_dct[name] += 120.      
        samp_zma = automol.zmat.set_values_by_name(samp_zma, new_key_dct)

    return samp_zma
