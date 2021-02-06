""" drivers for coordinate scans
"""

import os
import itertools
import random
import copy
import numpy
import automol
import autofile
from autofile import fs
import mess_io
from phydat import phycon
from ioformat import run_script
from mechlib.structure import vib as vibprep
from mechlib.submission import DEFAULT_SCRIPT_DCT


def read_hr_pot(names, grid_vals, cnf_save_path,
                mod_tors_ene_info, ref_ene,
                constraint_dct,
                read_geom=False, read_grad=False,
                read_hess=False, read_zma=False):
    """ Get the potential for a hindered rotor
    """

    # Build initial lists for storing potential energies and Hessians
    print('gridval test', grid_vals)
    grid_points = automol.pot.points(grid_vals)
    pot, geoms, grads, hessians, zmas, paths = {}, {}, {}, {}, {}, {}

    # Set up filesystem information
    zma_fs = fs.zmatrix(cnf_save_path)
    zma_path = zma_fs[-1].path([0])
    if constraint_dct is None:
        scn_fs = autofile.fs.scan(zma_path)
    else:
        scn_fs = autofile.fs.cscan(zma_path)

    # Read the energies and Hessians from the filesystem
    for point, vals in zip(grid_points, grid_vals):

        locs = [names, vals]
        if constraint_dct is not None:
            locs = [constraint_dct] + locs

        ene = read_tors_ene(scn_fs, locs, mod_tors_ene_info)
        if ene is not None:
            pot[point] = (ene - ref_ene) * phycon.EH2KCAL
        else:
            pot[point] = -10.0

        # print('path test in read_hr_pot:', scn_fs[-1].path(locs))
        if read_geom:
            if scn_fs[-1].file.geometry.exists(locs):
                geoms[point] = scn_fs[-1].file.geometry.read(locs)
            else:
                geoms[point] = None

        if read_grad:
            if scn_fs[-1].file.gradient.exists(locs):
                grads[point] = scn_fs[-1].file.gradient.read(locs)
            else:
                grads[point] = None

        if read_hess:
            if scn_fs[-1].file.hessian.exists(locs):
                hessians[point] = scn_fs[-1].file.hessian.read(locs)
            else:
                hessians[point] = None

        if read_zma:
            if scn_fs[-1].file.zmatrix.exists(locs):
                zmas[point] = scn_fs[-1].file.zmatrix.read(locs)
            else:
                zmas[point] = None

        paths[point] = scn_fs[-1].path(locs)

    return pot, geoms, grads, hessians, zmas, paths


def calc_hr_frequencies(geoms, grads, hessians, run_path):
    """ Calculate the frequencies
    """

    # Initialize hr freqs list
    hr_freqs = {}
    for point in geoms.keys():
        _, proj_freqs, _, _ = vibprep.projrot_freqs(
            [geoms[point]],
            [hessians[point]],
            run_path,
            grads=[grads[point]])
        hr_freqs[point] = proj_freqs

    return hr_freqs


def read_tors_ene(filesys, locs, mod_tors_ene_info):
    """ read the energy for torsions
    """

    if filesys[-1].exists(locs):
        path = filesys[-1].path(locs)
        sp_fs = autofile.fs.single_point(path)
        if sp_fs[-1].file.energy.exists(mod_tors_ene_info[1:4]):
            ene = sp_fs[-1].file.energy.read(mod_tors_ene_info[1:4])
        else:
            ene = None
    else:
        ene = None

    return ene


def print_hr_pot(tors_pots):
    """ Check hr pot to see if a new mimnimum is needed
    """

    print('\nHR potentials...')
    for name in tors_pots:

        print('- Rotor {}'.format(name))
        pot_str = ''
        for pot in tors_pots[name].values():
            pot_str += ' {0:.2f}'.format(pot)

        print('- Pot:{}'.format(pot_str))


def check_hr_pot(tors_pots, tors_zmas, tors_paths, emax=-0.5, emin=-10.0):
    """ Check hr pot to see if a new mimnimum is needed
    """

    new_min_zma = None

    print('\nAssessing the HR potential...')
    for name in tors_pots:

        print('- Rotor {}'.format(name))
        pots = tors_pots[name].values()
        zmas = tors_zmas[name].values()
        paths = tors_paths[name].values()
        for pot, zma, path in zip(pots, zmas, paths):
            if emin < pot < emax:
                new_min_zma = zma
                emin = pot
                print(' - New minimmum energy ZMA found for torsion')
                print(' - Ene = {}'.format(pot))
                print(' - Found at path: {}'.format(path))
                print(automol.zmat.string(zma))
        
    return new_min_zma


# Building constraints
def set_constraint_names(zma, tors_names, tors_model):
    """ Determine the names of constraints along a torsion scan
    """

    const_names = tuple()
    if tors_names and tors_model in ('1dhrf', '1dhrfa'):
        if tors_model == '1dhrf':
            const_names = tuple(
                itertools.chain(*tors_names))
        elif tors_model == '1dhrfa':
            coords = list(automol.zmat.coordinates(zma))
            const_names = tuple(coord for coord in coords)

    return const_names


# CALCULATE THE ZPES OF EACH TORSION USING MESS
def mess_tors_zpes(tors_geo, hind_rot_str, tors_save_path,
                   script_str=DEFAULT_SCRIPT_DCT['messpf']):
    """ Calculate the frequencies and ZPVES of the hindered rotors
        create a messpf input and run messpf to get tors_freqs and tors_zpes
    """

    # Set up the filesys
    bld_locs = ['PF', 0]
    bld_save_fs = autofile.fs.build(tors_save_path)
    bld_save_fs[-1].create(bld_locs)
    pf_path = bld_save_fs[-1].path(bld_locs)

    pf_path = os.path.join(pf_path, str(random.randint(0,1234567)))
    if not os.path.exists(pf_path):
        os.makedirs(pf_path)

    print('Run path for MESSPF:')
    print(pf_path)

    # Write the MESSPF input file
    global_pf_str = mess_io.writer.global_pf(
        temperatures=[100.0, 200.0, 300.0, 400.0, 500],
        rel_temp_inc=0.001,
        atom_dist_min=0.6)
    dat_str = mess_io.writer.molecule(
        core=mess_io.writer.core_rigidrotor(tors_geo, 1.0),
        freqs=[1000.0],
        elec_levels=[[0.0, 1.0]],
        hind_rot=hind_rot_str,
    )
    spc_str = mess_io.writer.species(
        spc_label='Tmp',
        spc_data=dat_str,
        zero_ene=0.0
    )
    pf_inp_str = '\n'.join([global_pf_str, spc_str]) + '\n'

    with open(os.path.join(pf_path, 'pf.inp'), 'w') as pf_file:
        pf_file.write(pf_inp_str)

    # Run MESSPF
    run_script(script_str, pf_path)

    # Obtain the torsional zpes and freqs from the MESS output
    with open(os.path.join(pf_path, 'pf.log'), 'r') as mess_file:
        output_string = mess_file.read()

    tors_zpes = mess_io.reader.tors.zero_point_vibrational_energies(
        output_string)
    # tors_freqs = mess_io.reader.tors.freqs(output_string)
    tors_freqs = mess_io.reader.tors.grid_minimum_frequencies(output_string)

    return tors_zpes, tors_freqs
