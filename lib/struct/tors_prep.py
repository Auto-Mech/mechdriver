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
    run_tors_grids = []
    for tors_names in run_tors_names:
        tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
            zma, tors_names, scan_increment, frm_bnd_key=frm_bnd_key,
            brk_bnd_key=brk_bnd_key)
        run_tors_grids.append(
            [numpy.linspace(*linspace) + val_dct[name]
             for name, linspace in zip(tors_names, tors_linspaces)]
        )
        # tors_sym_nums = tors.get_tors_sym_nums(
        #     spc_dct_i, tors_min_cnf_locs, tors_cnf_save_fs,
        #     frm_bnd_key, brk_bnd_key, saddle=False)

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
