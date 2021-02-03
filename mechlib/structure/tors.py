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


# FUNCTIONS TO SET UP TORSION NAME LISTS
def tors_name_prep(spc_dct_i, cnf_fs, min_cnf_locs, tors_model):
    """ Obtain torsion names through various means
    """

    if tors_model in ('1dhr', '1dhrf', '1dhrfa', 'mdhr', 'mdhrv', 'tau'):
        run_tors_names = ()
        if 'tors_names' in spc_dct_i:
            run_tors_names = spc_dct_i['tors_names']
            if '1dhr' in tors_model:
                run_tors_names = tuple(itertools.chain(*run_tors_names))
                run_tors_names = tuple((x,) for x in run_tors_names)
            else:
                run_tors_names = tuple(tuple(x) for x in run_tors_names)
            tloc = 'user'
        if not run_tors_names:
            run_tors_names = names_from_filesys(
                cnf_fs, min_cnf_locs, tors_model)
            tloc = 'fs'
        if not run_tors_names:
            run_tors_names = names_from_dct(spc_dct_i, tors_model)
            tloc = 'dct'
        if not run_tors_names and tors_model == 'tau':
            geo = cnf_fs[-1].file.geometry.read(min_cnf_locs)
            run_tors_names = names_from_geo(
                geo, tors_model, saddle=False)
            tloc = 'geo'
        if not run_tors_names:
            tloc = None
    else:
        tloc = None
        run_tors_names = ()

    if tloc is not None:
        if tloc == 'user':
            print(' - Reading tors names from user input...')
        elif tloc == 'fs':
            print(' - Reading tors names from the filesystem...')
        elif tloc == 'geo':
            print(' - Obtaining spc tors names from the geometry...')
        elif tloc == 'dct':
            print(' - Obtaining ts tors names from the geometry...')
    else:
        print('Warning: Your file system did not have any torsions defined')

    return run_tors_names


def names_from_geo(geo, tors_model, saddle=False):
    """ Build the tors name list from a geom
    """
    if not saddle:
        tors_names = [
            [name]
            for name in automol.geom.zmatrix_torsion_coordinate_names(geo)
        ]
        if tors_model in ('mdhr', 'mdhrv'):
            tors_names = [[tors
                           for rotor in tors_names
                           for tors in rotor]]
        tors_names = tuple(tuple(x) for x in tors_names)
    else:
        tors_names = tuple()

    return tors_names


def names_from_dct(spc_dct_i, tors_model):
    """ Build the tors name list from a dictionary
    """

    # Read names from dct
    amech_ts_tors_names = ()
    if 'amech_ts_tors_names' in spc_dct_i:
        amech_ts_tors_names = spc_dct_i['amech_ts_tors_names']
        if tors_model in ('1dhr', '1dhrf', '1dhrfa'):
            amech_ts_tors_names = [[name] for name in amech_ts_tors_names]
        else:
            amech_ts_tors_names = [amech_ts_tors_names]
        amech_ts_tors_names = tuple(tuple(x) for x in amech_ts_tors_names)

    return amech_ts_tors_names


def names_from_filesys(tors_cnf_fs, tors_min_cnf_locs, tors_model):
    """ read names from filesystem
    """

    tors_names = None
    if tors_min_cnf_locs is not None:
        if tors_cnf_fs[0].file.info.exists():
            inf_obj = tors_cnf_fs[0].file.info.read()
            tors_range_dct = dict(inf_obj.tors_ranges)
            tors_names = list(tors_range_dct.keys())

    if tors_names is not None:
        if tors_model in ('1dhr', '1dhrf', '1dhrfa', 'tau'):
            tors_names = [[name] for name in tors_names]
        else:
            tors_names = [[name for name in tors_names]]
        tors_names = tuple(tuple(x) for x in tors_names)

    return tors_names


# FUNCTIONS USED TO BUILD LSTS OF TORSIONS OF ANY DIMENSIONALITY
def hr_prep(zma, tors_name_grps, scan_increment=30.0, tors_model='1dhr',
            frm_bnd_keys=(), brk_bnd_keys=()):
    """ set-up the hr for different rotor combinations
        tors_names = [ ['D1'], ['D2', 'D3'], ['D4'] ]
    """

    # Get the tors names if thery have not already been supplied
    val_dct = automol.zmat.value_dictionary(zma)

    # Deal with the dimensionality of the rotors
    if tors_model in ('mdhr', 'mdhrv'):
        tors_name_grps = mdhr_prep(zma, tors_name_grps)

    # Build the grids corresponding to the torsions
    tors_grids, tors_sym_nums = [], []
    for tors_names in tors_name_grps:
        tors_linspaces = automol.zmat.torsional_scan_linspaces(
            zma, tors_names, scan_increment, frm_bnd_key=frm_bnd_keys,
            brk_bnd_key=brk_bnd_keys)
        tors_grids.append(
            [numpy.linspace(*linspace) + val_dct[name]
             for name, linspace in zip(tors_names, tors_linspaces)]
        )
        tors_sym_nums.append(
            automol.zmat.torsional_symmetry_numbers(
                zma, tors_names,
                frm_bnd_key=frm_bnd_keys, brk_bnd_key=brk_bnd_keys)
        )

    return tors_name_grps, tors_grids, tors_sym_nums


def mdhr_prep(zma, run_tors_names):
    """ Handle cases where the MDHR
    """

    # Figure out set of torsions are to be used: defined or AMech generated
    rotor_lst = run_tors_names

    # Check the dimensionality of each rotor to see if they are greater than 4
    # Call a function to reduce large rotors
    final_rotor_lst = []
    for rotor in rotor_lst:
        if len(rotor) > 4:
            for reduced_rotor in reduce_rotor_dimensionality(zma, rotor):
                final_rotor_lst.append(reduced_rotor)
        else:
            final_rotor_lst.append(rotor)

    return final_rotor_lst


def reduce_rotor_dimensionality(zma, rotor):
    """ For rotors with a dimensionality greater than 4, try and take them out
    """

    # Find the methyl rotors for that are a part of the MDHR
    reduced_rotor_lst = []
    methyl_rotors = []
    for tors in rotor:
        # If a methyl rotor add to methyl rotor list, or add to reduced lst
        if is_methyl_rotor(zma, rotor):   # Add arguments when ID methyls
            methyl_rotors.append(zma, tors)
        else:
            reduced_rotor_lst.append(tors)

    # Add each of methyl rotors, if any exist
    if methyl_rotors:
        for methyl_rotor in methyl_rotors:
            reduced_rotor_lst.append(methyl_rotor)

    # Check new dimensionality of list; if still high, flatten to lst of 1DHRs
    if len(reduced_rotor_lst) > 4:
        reduced_rotor_lst = [tors
                             for rotor in reduced_rotor_lst
                             for tors in rotor]

    return reduced_rotor_lst


def is_methyl_rotor(zma, rotor):
    """ Check if methyl rotor
    """
    raise NotImplementedError(zma, rotor)


# Handle the potential
def set_scan_dims(tors_grids):
    """ Determine the dimensions of the grid
    """
    assert len(tors_grids) in (1, 2, 3, 4), 'Rotor must be 1-4 dimensions'

    grid_points = ((i for i in range(len(grid)))
                   for grid in tors_grids)
    grid_vals = ((x for x in grid)
                 for grid in tors_grids)
    grid_points = tuple(itertools.product(*grid_points))
    grid_vals = tuple(itertools.product(*grid_vals))

    return grid_points, grid_vals


def read_hr_pot(torsionss, cnf_save_path,
                mod_tors_ene_info, ref_ene,
                constraint_dct,
                scan_increment=0.523599,
                read_geom=False, read_grad=False,
                read_hess=False, read_zma=False):
    """ Get the potential for a hindered rotor
    """
    
    names = automol.rotor.names(torsions)
    grids = automol.rotor.grids(torsions, increment=scan_increment)

    # Build initial lists for storing potential energies and Hessians
    grid_points, grid_vals = set_scan_dims(tors_grids)
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

        locs = [tors_names, vals]
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
            # if saddle:
            #     const_names = tuple(
            #         itertools.chain(*amech_sadpt_tors_names))
            # else:
            #     const_names = tuple(
            #         itertools.chain(*amech_spc_tors_names))
        elif tors_model == '1dhrfa':
            coords = list(automol.zmat.coordinates(zma))
            const_names = tuple(coord for coord in coords)

    return const_names


def build_constraint_dct(zma, const_names, scan_names=()):
    """ Build a dictionary of constraints
    """

    # print('const_names', const_names)
    # print('scan_names', scan_names)

    # Get the list names sorted for dictionary
    rnames = (name for name in const_names if 'R' in name)
    anames = (name for name in const_names if 'A' in name)
    dnames = (name for name in const_names if 'D' in name)
    rnames = tuple(sorted(rnames, key=lambda x: int(x.split('R')[1])))
    anames = tuple(sorted(anames, key=lambda x: int(x.split('A')[1])))
    dnames = tuple(sorted(dnames, key=lambda x: int(x.split('D')[1])))
    constraint_names = rnames + anames + dnames

    # Remove the scan coordinates so they are not placed in the dict

    constraint_names = tuple(name for name in constraint_names
                             if name not in scan_names)

    # Build dictionary
    if constraint_names:
        zma_vals = automol.zmat.value_dictionary(zma)
        zma_coords = automol.zmat.coordinates(zma)
        assert set(constraint_names) <= set(zma_coords.keys()), (
            'Attempting to constrain coordinates not in zma:\n{}\n{}'.format(
                constraint_names, zma_coords)
        )

        constraint_dct = dict(zip(
            constraint_names,
            (round(zma_vals[name], 2) for name in constraint_names)
        ))
    else:
        constraint_dct = None

    # print('dct')
    # print(constraint_dct)
    # import sys
    # sys.exit()

    return constraint_dct


# Functions to handle setting up groups and axes used to define torstions
def set_tors_def_info(zma, tors_name, tors_sym, pot,
                      frm_bnd_keys, brk_bnd_keys,
                      rxn_class, saddle=False):
    """ set stuff
    """

    # Set the ts bnd to the break or form keys based on rxn class
    if frm_bnd_keys:
        ts_bnd = frm_bnd_keys
    else:
        ts_bnd = brk_bnd_keys

    # print('set_tors_def_info test:', automol.zmat.string(zma), tors_name, ts_bnd, frm_bnd_keys, brk_bnd_keys)
    group, axis, atm_key = _set_groups_ini(
        zma, tors_name, ts_bnd, saddle)
    if saddle:
        # print('pot test:', pot)
        group, axis, pot, chkd_sym_num = _check_saddle_groups(
            zma, rxn_class, group, axis,
            pot, ts_bnd, tors_sym)
    else:
        chkd_sym_num = tors_sym
    group = list(numpy.add(group, 1))
    axis = list(numpy.add(axis, 1))
    if (atm_key+1) != axis[1]:
        axis.reverse()

    return group, axis, pot, chkd_sym_num


def _set_groups_ini(zma, tors_name, ts_bnd, saddle):
    """ Set the initial set of groups
    """
    gra = automol.zmat.graph(zma, stereo=False)
    # ts_gra = automol.graph.add_ts_bonds(gra, keys=[ts_bnd])
    coo_dct = automol.zmat.coordinates(zma, multi=False)
    axis = coo_dct[tors_name][1:3]
    atm_key = axis[1]
    if ts_bnd:
        for atm in axis:
            if atm in ts_bnd:
                atm_key = atm
                break
    # print('tors_def_info test:')
    # print('zma test:', automol.zmat.string(zma))
    # print('tors_name test:', tors_name)
    # print('ts_bnd test:', ts_bnd)
    # print('saddle test:', saddle)
    # print('axis test:', axis)
    # print('atm key test:', atm_key)
    # print('gra test:', gra)

    if saddle:
        if ts_bnd not in automol.graph.bond_keys(gra):
            gra = automol.graph.add_ts_bonds(gra, keys=[ts_bnd])
    group = list(
        automol.graph.branch_atom_keys(
            gra, atm_key, axis) - set(axis))
    if not group:
        for atm in axis:
            if atm != atm_key:
                atm_key = atm
        group = list(
            automol.graph.branch_atom_keys(
                gra, atm_key, axis) - set(axis))

    return group, axis, atm_key


def _check_saddle_groups(zma, rxn_class, group, axis, pot, ts_bnd, sym_num):
    """ Assess the hindered rotor groups and axes
    """
    n_atm = automol.zmat.count(zma)
    if 'addition' in rxn_class or 'abstraction' in rxn_class:
        group2 = []
        ts_bnd1 = min(ts_bnd)
        ts_bnd2 = max(ts_bnd)
        for idx in range(ts_bnd2, n_atm):
            group2.append(idx)
        if ts_bnd1 in group:
            for atm in group2:
                if atm not in group:
                    group.append(atm)

    # Check to see if symmetry of XH3 rotor was missed
    new_pot = False
    if sym_num == 1:
        group2 = []
        for idx in range(n_atm):
            if idx not in group and idx not in axis:
                group2.append(idx)
        all_hyd = True
        symbols = automol.zmat.symbols(zma)
        hyd_count = 0
        for idx in group2:
            if symbols[idx] != 'H' and symbols[idx] != 'X':
                all_hyd = False
                break
            if symbols[idx] == 'H':
                hyd_count += 1
        if all_hyd and hyd_count == 3:
            sym_num = 3
            lpot = int(len(pot)/3)
            new_pot = True
            print('sym_num test 2:', sym_num, lpot, new_pot)
    if new_pot:
        potr = {}
        for pval in range(lpot):
            potr[(pval,)] = pot[(pval,)]
    else:
        potr = copy.deepcopy(pot)

    return group, axis, potr, sym_num


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
