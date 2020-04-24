""" drivers for coordinate scans
"""

import numpy
import automol


# FUNCTIONS USED TO BUILD LSTS OF TORSIONS OF ANY DIMENSIONALITY
def hr_prep(zma, geo, run_tors_names=(), scan_increment=30.0, ndim_tors='1dhr',
            saddle=False, frm_bnd_key=(), brk_bnd_key=()):
    """ set-up the hr for different rotor combinations
        tors_names = [ ['D1'], ['D2', 'D3'], ['D4'] ]
    """

    # Get the tors names if thery have not already been supplied
    val_dct = automol.zmatrix.values(zma)
    if not run_tors_names:
        if not saddle:
            run_tors_names = [
                [name]
                for name in automol.geom.zmatrix_torsion_coordinate_names(geo)
            ]
            if ndim_tors in ('mdhr', 'mdhrv'):
                run_tors_names = [[tors
                                   for rotor in run_tors_names
                                   for tors in rotor]]

    # Deal with the dimensionality of the rotors
    if ndim_tors in ('mdhr', 'mdhrv'):
        run_tors_names = mdhr_prep(zma, run_tors_names)

    # Build the grids corresponding to the torsions
    run_tors_grids, run_tors_syms = [], []
    for tors_names in run_tors_names:
        tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
            zma, tors_names, scan_increment, frm_bnd_key=frm_bnd_key,
            brk_bnd_key=brk_bnd_key)
        run_tors_grids.append(
            [numpy.linspace(*linspace) + val_dct[name]
             for name, linspace in zip(tors_names, tors_linspaces)]
        )
        # tors_sym_nums.append(list(automol.zmatrix.torsional_symmetry_numbers(
        #    zma, tors_names, frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key))

    return run_tors_names, run_tors_grids  # run_tors_syms


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
        # If a methyl rotor add to methyl rotor list
        if is_methyl_rotor():   # Add arguments when ID methyls
            methyl_rotors.append(tors)
        # Add to reduced rotor list
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


def is_methyl_rotor():
    """ Check if methyl rotor
    """
    return False


# Building constraints
def build_constraint_dct(zma, tors_names):
    """ Build a dictionary of constraints
    """
    constraint_names = [name
                        for name_lst in tors_names
                        for name in name_lst]
    constraint_names.sort(key=lambda x: int(x.split('D')[1]))
    zma_vals = automol.zmatrix.values(zma)
    constraint_dct = dict(zip(
        constraint_names,
        (round(zma_vals[name], 2) for name in constraint_names)
    ))

    return constraint_dct


# Functions to handle setting up torsional defintion and potentials properly
def set_groups_ini(zma, tors_name, ts_bnd, saddle):
    """ Set the initial set of groups
    """
    gra = automol.zmatrix.graph(zma, remove_stereo=True)
    coo_dct = automol.zmatrix.coordinates(zma, multi=False)
    axis = coo_dct[tors_name][1:3]
    atm_key = axis[1]
    if ts_bnd:
        for atm in axis:
            if atm in ts_bnd:
                atm_key = atm
                break
    group = list(
        automol.graph.branch_atom_keys(
            gra, atm_key, axis, saddle=saddle, ts_bnd=ts_bnd) - set(axis))
    if not group:
        for atm in axis:
            if atm != atm_key:
                atm_key = atm
        group = list(
            automol.graph.branch_atom_keys(
                gra, atm_key, axis, saddle=saddle, ts_bnd=ts_bnd) - set(axis))

    return group, axis, atm_key


def check_saddle_groups(zma, rxn_class, group, axis, pot, ts_bnd, sym_num):
    """ Assess that hindered rotor groups and axes
    """
    n_atm = automol.zmatrix.count(zma)
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
    if sym_num == 1:
        group2 = []
        for idx in range(n_atm):
            if idx not in group and idx not in axis:
                group2.append(idx)
        all_hyd = True
        symbols = automol.zmatrix.symbols(zma)
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
            potp = []
            potp[0:lpot] = pot[0:lpot]
            pot = potp

    return group, axis, pot


def check_dummy_trans(zma):
    """ check trans
    """
    atom_symbols = automol.zmatrix.symbols(zma)
    dummy_idx = []
    for atm_idx, atm in enumerate(atom_symbols):
        if atm == 'X':
            dummy_idx.append(atm_idx)
    remdummy = numpy.zeros(len(zma[0]))
    for dummy in dummy_idx:
        for idx, _ in enumerate(remdummy):
            if dummy < idx:
                remdummy[idx] += 1

    return remdummy
