"""
  Functions handling hindered rotor model calculations
"""

import numpy
from scipy.interpolate import interp1d
import automol
import mess_io
import projrot_io
import autofile
from lib.phydat import phycon
from lib.structure import tors as torsprep


# Function to deal with setting up all of the torsion info since it is a pain
def rotor_info(spc_dct_i, pf_filesystems, pf_models,
               saddle=False, frm_bnd_key=(), brk_bnd_key=()):
    """ get tors stuff
    """

    # Set up tors level filesystem and model
    tors_model = pf_models['tors']
    [cnf_fs, cnf_path, min_cnf_locs, _] = pf_filesystems['tors']

    # Read the increments from the filesystem
    if 'hind_inc' in spc_dct_i:
        scan_increment = spc_dct_i['hind_inc'] * phycon.DEG2RAD
    else:
        scan_increment = 30. * phycon.DEG2RAD

    # Set up ndim tors
    ndim_tors = '1dhr' if '1dhr' in tors_model else 'mdhr'

    # Set up the names of all the torsions in the rotors through hierarchy
    run_tors_names = ()
    if 'tors_names' in spc_dct_i:
        run_tors_names = torsprep.names_from_dct(spc_dct_i, ndim_tors)
        print('Reading tors names from user input')
        print(run_tors_names)
    if not run_tors_names:
        run_tors_names = torsprep.names_from_filesys(min_cnf_locs, cnf_path)
        print('Reading tors names from the filesystem')
        print(run_tors_names)

    # Obtain the info for each of the rotors
    rotor_names, rotor_grids, rotor_syms = (), (), ()
    if min_cnf_locs is not None:

        print('tors stuff')

        # Get geometry and reference energy for the torsional minimum
        zma = cnf_fs[-1].file.zmatrix.read(min_cnf_locs)
        geo = cnf_fs[-1].file.geometry.read(min_cnf_locs)
        ref_ene = cnf_fs[-1].file.energy.read(min_cnf_locs)

        # Get the hr prep stuff
        rotor_names, rotor_grids, rotor_syms = torsprep.hr_prep(
            zma, geo, run_tors_names=run_tors_names,
            scan_increment=scan_increment, tors_model=tors_model,
            saddle=saddle, frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key)

        # Build constraint dct
        if tors_model in ('1dhrf'):
            constraint_dct = torsprep.build_constraint_dct(
                zma, rotor_names)
        else:
            constraint_dct = None
    else:
        print('No tors info at path')

    return rotor_names, rotor_grids, rotor_syms, constraint_dct, ref_ene


# Functions to write MESS strings
def make_hr_strings(rotor_names, rotor_grids, rotor_syms, constraint_dct,
                    ref_ene, pf_filesystems, pf_models,
                    rxn_class, ts_bnd,
                    saddle=False, tors_wgeo=False):
    """ Procedure for building the MESS strings
    """

    # Set up tors level filesystem and model
    tors_model = pf_models['tors']
    [cnf_fs, cnf_path, min_cnf_locs, _] = pf_filesystems['tors']

    # Grab the torsional geometry if needed
    if tors_wgeo:
        hind_rot_geo = cnf_fs[-1].file.geometry.read(min_cnf_locs)
    else:
        hind_rot_geo = None

    # Grab the zmatrix
    zma = cnf_fs[-1].file.zmatrix.read(min_cnf_locs)

    # Write strings containing rotor info for MESS and ProjRot
    if tors_model in ('1dhr', '1dhrf', '1dhrv', '1dhrfv'):
        mess_str, projrot_str, chkd_sym_nums = _make_1dhr_tors_strs(
            zma, rxn_class, ts_bnd, ref_ene,
            rotor_names, rotor_grids, rotor_syms,
            constraint_dct, cnf_path,
            saddle=saddle,
            read_freqs=bool('v' in tors_model),
            hind_rot_geo=hind_rot_geo)
        mdhr_dats = []
    elif tors_model in ('mdhr', 'mdhrv'):
        mess_str, projrot_str, mdhr_dats, chkd_sym_nums = _make_mdhr_tors_strs(
            zma, rxn_class, ts_bnd, ref_ene,
            rotor_names, rotor_grids, rotor_syms,
            constraint_dct, cnf_path,
            saddle=saddle,
            read_freqs=bool('v' in tors_model),
            hind_rot_geo=hind_rot_geo)

    return mess_str, projrot_str, mdhr_dats, chkd_sym_nums


def _make_1dhr_tors_strs(zma, rxn_class, ts_bnd, ref_ene,
                         rotor_names, rotor_grids, rotor_syms,
                         constraint_dct, cnf_save_path,
                         saddle=False, read_freqs=False, hind_rot_geo=None):
    """ Gather the 1DHR torsional data and gather them into a MESS file
    """

    # Initialize empty strings and torsional sym num lsts
    mess_hr_str, projrot_hr_str = '', ''
    chkd_sym_nums = []

    # Loop over the torsions
    tors_info = zip(rotor_names, rotor_grids, rotor_syms)
    for tors_names, tors_grids, tors_sym in tors_info:

        # Grab zero elment because of formatting (each list has one tors)
        tors_name = tors_names[0]
        tors_grid = tors_grids[0]

        # Read the hindered rotor potential
        pot, _ = _read_hr_pot(
            [tors_name], tors_grid,
            cnf_save_path, ref_ene,
            constraint_dct, read_freqs=read_freqs)

        # Use a spline fitter to deal with failures in the potential
        pot = _hrpot_spline_fitter(pot)

        # Get the HR groups and axis for the rotor
        group, axis, chkd_sym_num = torsprep.set_tors_def_info(
            zma, tors_name, tors_sym, pot,
            ts_bnd, rxn_class, saddle=saddle)

        # Check for dummy atoms
        remdummy = torsprep.check_dummy_trans(zma)

        # Write the MESS and ProjRot strings for the rotor
        mess_hr_str += mess_io.writer.rotor_hindered(
            group=group,
            axis=axis,
            symmetry=chkd_sym_num,
            potential=pot,
            remdummy=remdummy,
            geom=hind_rot_geo,
            use_quantum_weight=True)
        projrot_hr_str += projrot_io.writer.rotors(
            axis=axis,
            group=group,
            remdummy=remdummy)

        # Append the sym number to lst
        chkd_sym_nums.append(chkd_sym_num)

    return mess_hr_str, projrot_hr_str, chkd_sym_nums


def _make_mdhr_tors_strs(zma, rxn_class, ts_bnd, ref_ene,
                         rotor_names, rotor_grids, rotor_syms,
                         constraint_dct, cnf_save_path,
                         saddle=False, read_freqs=False, hind_rot_geo=None):
    """ Gather the MDHR torsional data and gather them into a MESS file
    """

    # Initialize empty strings and torsional sym num lsts
    mess_hr_str, projrot_hr_str = '', ''
    mdhr_dat_str_lst = []
    chkd_sym_nums = []

    # Loop over the torsions
    tors_info = zip(rotor_names, rotor_grids, rotor_syms)
    for tors_names, tors_grids, tors_syms in tors_info:

        # Read the hindered rotor potential and add to master list
        hr_pot, hr_freqs = _read_hr_pot(
            tors_names, tors_grids,
            cnf_save_path, ref_ene,
            constraint_dct, read_freqs=read_freqs)

        # Write the MDHR potential file for each rotor set; append to lst
        mdhr_dat_str_lst.append(mess_io.mdhr_data(hr_pot, hr_freqs))

        # Check for dummy transformations
        remdummy = torsprep.check_dummy_trans(zma)

        # Loop over the rotors in the group and write the internal rotor strs
        for tors_name, tors_sym in zip(tors_names, tors_syms):

            # Set pot to empty list (may need fix)
            pot = ()  # Change hpw set_tors_def_info works?

            # Get the HR groups and axis for the rotor (new check?)
            group, axis, chkd_sym_num = torsprep.set_tors_def_info(
                zma, tors_name, tors_sym, pot,
                ts_bnd, rxn_class, saddle=saddle)

            # Write the MESS and ProjRot strings for the rotor
            mess_hr_str += mess_io.writer.mol_data.rotor_internal(
                group=group,
                axis=axis,
                symmetry=chkd_sym_num,
                grid_size=100,
                mass_exp_size=5,
                pot_exp_size=5,
                hmin=13,
                hmax=101,
                remdummy=remdummy,
                geom=hind_rot_geo,
                rotor_id=tors_name)
            projrot_hr_str += projrot_io.writer.rotors(
                axis=axis,
                group=group,
                remdummy=remdummy)

    return mess_hr_str, projrot_hr_str, mdhr_dat_str_lst, chkd_sym_nums


def make_flux_str(tors_min_cnf_locs, tors_cnf_save_fs,
                  tors_names, rotor_syms):
    """ Write out the input string for tau samling
    """
    # Loop over the torsions to get the flux strings
    flux_mode_str = ''
    if tors_min_cnf_locs is not None:

        # Get geometry for the torsional minimum
        zma = tors_cnf_save_fs[-1].file.zmatrix.read(
            tors_min_cnf_locs)
        name_matrix = automol.zmatrix.name_matrix(zma)
        key_matrix = automol.zmatrix.key_matrix(zma)

        # Write the MESS flux strings for each of the modes
        for tors_name, tors_sym in zip(tors_names, rotor_syms):

            # Get the idxs for the atoms used to define the dihedrals
            # Move code at some point to automol
            mode_idxs = [0, 0, 0, 0]
            for i, name_row in enumerate(name_matrix):
                if tors_name in name_row:
                    mode_idxs[0] = i
                    mode_idxs[1], mode_idxs[2], mode_idxs[3] = key_matrix[i]
                    break
            mode_idxs = [idx+1 for idx in mode_idxs]

            # Determine the flux span using the symmetry number
            mode_span = 360.0 / tors_sym

            # Write flux string for mode_ to overall flux string
            flux_mode_str += mess_io.writer.fluxional_mode(
                mode_idxs, span=mode_span)

    return flux_mode_str


# Functions to obtain values of the HR potentials from the filesystem
def _read_hr_pot(tors_names, tors_grids, cnf_save_path, ref_ene,
                 constraint_dct, read_freqs=False):
    """ Get the potential for a hindered rotor
    """

    # Build template pot lst and freqs list into a list-of-lists if ndim > 1
    if len(tors_names) == 1:
        dims = (len(tors_grids),)
    elif len(tors_names) == 2:
        dims = (len(tors_grids[0]), len(tors_grids[1]))
    elif len(tors_names) == 3:
        dims = (len(tors_grids[0]), len(tors_grids[1]), len(tors_grids[2]))
    elif len(tors_names) == 4:
        dims = (len(tors_grids[0]), len(tors_grids[1]),
                len(tors_grids[2]), len(tors_grids[3]))
    pot = numpy.zeros(dims).tolist()
    if read_freqs:
        freqs = numpy.zeros(dims).tolist()
    else:
        freqs = []

    # Read the energies from the filesystem
    if constraint_dct is None:
        scn_fs = autofile.fs.scan(cnf_save_path)
    else:
        scn_fs = autofile.fs.cscan(cnf_save_path)
    if len(tors_names) == 1:
        for i, grid_val_i in enumerate(tors_grids):
            # Set locs
            locs = [tors_names, [grid_val_i]]
            if constraint_dct is not None:
                locs.append(constraint_dct)
            # Read filesys
            if scn_fs[-1].exists(locs):
                ene = scn_fs[-1].file.energy.read(locs)
                pot[i] = (ene - ref_ene) * phycon.EH2KCAL
            else:
                pot[i] = 10.0
            if read_freqs:
                freqs[i] = scn_fs[-1].file.harmonic_frequencies.read(locs)
    elif len(tors_names) == 2:
        for i, grid_val_i in enumerate(tors_grids[0]):
            for j, grid_val_j in enumerate(tors_grids[1]):
                # Set locs
                locs = [tors_names, [grid_val_i, grid_val_j]]
                if constraint_dct is not None:
                    locs.append(constraint_dct)
                # Read filesys
                if scn_fs[-1].exists(locs):
                    ene = scn_fs[-1].file.energy.read(locs)
                    pot[i][j] = (ene - ref_ene) * phycon.EH2KCAL
                else:
                    pot[i][j] = 10.0
                if read_freqs:
                    freqs = scn_fs[-1].file.harmonic_frequencies.read(locs)
                    freqs[i][j] = freqs
    elif len(tors_names) == 3:
        for i, grid_val_i in enumerate(tors_grids[0]):
            for j, grid_val_j in enumerate(tors_grids[1]):
                for k, grid_val_k in enumerate(tors_grids[2]):
                    # Set locs
                    locs = [tors_names, [grid_val_i, grid_val_j, grid_val_k]]
                    if constraint_dct is not None:
                        locs.append(constraint_dct)
                    # Read filesys
                    if scn_fs[-1].exists(locs):
                        ene = scn_fs[-1].file.energy.read(locs)
                        pot[i][j][k] = (ene - ref_ene) * phycon.EH2KCAL
                    else:
                        pot[i][j][k] = 10.0
                    if read_freqs:
                        freqs = scn_fs[-1].file.harmonic_frequencies.read(locs)
                        freqs[i][j][k] = freqs
    elif len(tors_names) == 4:
        for i, grid_val_i in enumerate(tors_grids[0]):
            for j, grid_val_j in enumerate(tors_grids[1]):
                for k, grid_val_k in enumerate(tors_grids[2]):
                    for lma, grid_val_l in enumerate(tors_grids[3]):
                        # Set locs
                        locs = [tors_names,
                                [grid_val_i, grid_val_j,
                                 grid_val_k, grid_val_l]]
                        if constraint_dct is not None:
                            locs.append(constraint_dct)
                        # Read filesys
                        if scn_fs[-1].exists(locs):
                            ene = scn_fs[-1].file.energy.read(locs)
                            pot[i][j][k][lma] = (
                                (ene - ref_ene) * phycon.EH2KCAL)
                        else:
                            pot[i][j][k][lma] = 10.0
                        if read_freqs:
                            freqs = scn_fs[-1].file.harmonic_frequencies.read(
                                locs)
                            freqs[i][j][k][lma] = freqs

    return pot, freqs


def _hrpot_spline_fitter(pot, min_thresh=-0.05, max_thresh=15.0):
    """ Get a physical hindered rotor potential via a series of spline fits
    """

    # Initialize a variable for the size of the potential
    lpot = len(pot)+1

    # Build a potential list from only successful calculations
    idx_success = []
    pot_success = [0.0]
    for idx in range(lpot):
        if pot[idx] < 600.:
            idx_success.append(idx)
            pot_success.append(pot[idx])

    # Print warning messages
    if any(val > max_thresh for val in pot):
        max_pot = max(pot)
        print('Warning: Found pot val of {0:.2f}'.format(max_pot),
              ' which is larger than',
              'the typical maximum for a torsional potential')
        print('Potential before spline:', pot)
    if any(val < min_thresh for val in pot):
        min_pot = min(pot)
        print('Warning: Found pot val of {0:.2f}'.format(min_pot),
              ' which is below',
              '{0} kcal. Refit w/ positives'.format(min_thresh))
        print('Potential before spline:', pot)

    # Build a new potential list using a spline fit of the HR potential
    pot_spl = interp1d(
        numpy.array(idx_success), numpy.array(pot_success), kind='cubic')
    for idx in range(lpot):
        pot[idx] = float(pot_spl(idx))

    # Do second spline fit of only positive values if any negative values found
    if any(val < min_thresh for val in pot):
        print('Still found negative potential values after first spline')
        print('Potential after spline:', pot)
        x_pos = numpy.array([i for i in range(lpot)
                             if pot[i] >= min_thresh])
        y_pos = numpy.array([pot[i] for i in range(lpot)
                             if pot[i] >= min_thresh])
        pos_pot_spl = interp1d(x_pos, y_pos, kind='cubic')
        pot_pos_fit = []
        for idx in range(lpot):
            pot_pos_fit.append(pos_pot_spl(idx))

        print('Potential after spline:', pot_pos_fit)
        # Perform second check to see if negative potentials have been fixed
        if any(val < min_thresh for val in pot_pos_fit):
            print('Still found negative potential values after second spline')
            print('Replace with linear interpolation of positive values')
            neg_idxs = [i for i in range(lpot) if pot_pos_fit[i] < min_thresh]
            clean_pot = []
            for i in range(lpot):
                if i in neg_idxs:
                    # Find the indices for positive vals around negative value
                    idx_0 = i - 1
                    while idx_0 in neg_idxs:
                        idx_0 = idx_0 - 1
                    for j in range(i, lpot):
                        if pot_pos_fit[j] >= min_thresh:
                            idx_1 = j
                            break
                    # Get a new value for this point on the potential by
                    # doing a linear interp of positives
                    interp_val = (
                        pot_pos_fit[idx_0] * (1.0-((i-idx_0)/(idx_1-idx_0))) +
                        pot_pos_fit[idx_1] * ((i-idx_0)/(idx_1-idx_0))
                    )
                    clean_pot.append(interp_val)
                else:
                    clean_pot.append(pot_pos_fit[i])
            final_potential = clean_pot.copy()

        else:
            final_potential = pot_pos_fit.copy()

    else:
        final_potential = pot.copy()

    final_potential = final_potential[:-1]

    return final_potential
