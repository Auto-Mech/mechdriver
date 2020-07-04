"""
  Functions handling hindered rotor model calculations
"""

import numpy
from scipy.interpolate import interp1d
import automol
import mess_io
import projrot_io
import autofile
from autofile import fs
from lib.phydat import phycon
from lib.structure import tors as torsprep
from lib.structure import geom as geomprep


# Function to deal with setting up all of the torsion info since it is a pain
def rotor_info(spc_dct_i, pf_filesystems, pf_models,
               frm_bnd_key=(), brk_bnd_key=()):
    """ get tors stuff
    """

    # Set up tors level filesystem and model
    tors_model = pf_models['tors']
    [cnf_fs, cnf_path, min_cnf_locs, _, _] = pf_filesystems['tors']

    # Read the increments from the filesystem
    if 'hind_inc' in spc_dct_i:
        scan_increment = spc_dct_i['hind_inc']
    else:
        scan_increment = 30. * phycon.DEG2RAD

    # Set up ndim tors
    if '1dhr' in tors_model:
        tors_model = '1dhr'
    elif 'mdhr' in tors_model:
        tors_model = 'mdhr'

    # Set up the names of all the torsions in the rotors through hierarchy
    if tors_model in ('1dhr', 'mdhr', 'tau'):
        run_tors_names = ()
        if 'tors_names' in spc_dct_i:
            run_tors_names = torsprep.names_from_dct(spc_dct_i, tors_model)
            tloc = 'dct'
        if not run_tors_names:
            run_tors_names = torsprep.names_from_filesys(
                cnf_fs, min_cnf_locs, tors_model)
            tloc = 'fs'
        if not run_tors_names and tors_model == 'tau':
            geo = cnf_fs[-1].file.geometry.read(min_cnf_locs)
            run_tors_names = torsprep.names_from_geo(geo, tors_model, saddle=False)
            tloc = 'geo'
        if not run_tors_names:
            tloc = None
    else:
        tloc = None

    # Obtain the info for each of the rotors
    rotor_names, rotor_grids, rotor_syms = (), (), ()
    constraint_dct = None
    ref_ene = None
    if tloc is not None and min_cnf_locs is not None:

        if tloc == 'dct':
            print(' - Reading tors names from user input...')
        elif tloc == 'fs':
            print(' - Reading tors names from the filesystem...')
        elif tloc == 'geo':
            print(' - Obtaining tors names from the geometry...')
        for idx, tors_names in enumerate(run_tors_names):
            print(' - Rotor {}: {}'.format(str(idx+1), '-'.join(tors_names)))

        # Get geometry and reference energy for the torsional minimum
        zma_fs = fs.manager(cnf_fs[-1].path(min_cnf_locs), 'ZMATRIX')
        zma = zma_fs[-1].file.zmatrix.read([0])
        ref_ene = cnf_fs[-1].file.energy.read(min_cnf_locs)

        # Get the hr prep stuff
        rotor_names, rotor_grids, rotor_syms = torsprep.hr_prep(
            zma=zma,
            tors_name_grps=run_tors_names,
            scan_increment=scan_increment,
            tors_model=tors_model,
            frm_bnd_key=frm_bnd_key,
            brk_bnd_key=brk_bnd_key)

        # Build constraint dct
        if tors_model in ['1dhrf']:
            constraint_dct = torsprep.build_constraint_dct(
                zma, rotor_names)
    else:
        print('No tors info at path')

    return rotor_names, rotor_grids, rotor_syms, constraint_dct, ref_ene


# Functions to write MESS strings
def make_hr_strings(rotor_names, rotor_grids, rotor_syms, constraint_dct,
                    ref_ene, pf_filesystems, pf_models,
                    rxn_class, ts_bnd,
                    saddle=False, tors_wgeo=False,
                    build_mess=True, build_projrot=True):
    """ Procedure for building the MESS strings
    """

    # Set up tors level filesystem and model
    tors_model = pf_models['tors']
    [cnf_fs, cnf_path, min_cnf_locs, _, _] = pf_filesystems['tors']

    # Grab the torsional geometry if needed
    if tors_wgeo:
        hind_rot_geo = cnf_fs[-1].file.geometry.read(min_cnf_locs)
    else:
        hind_rot_geo = None

    # Grab the zmatrix
    zma_fs = fs.manager(cnf_fs[-1].path(min_cnf_locs), 'ZMATRIX')
    zma = zma_fs[-1].file.zmatrix.read([0])
    # zma = cnf_fs[-1].file.zmatrix.read(min_cnf_locs)

    # Write strings containing rotor info for MESS and ProjRot
    if tors_model in ('1dhr', '1dhrf', '1dhrv', '1dhrfv', 'tau'):
        mess_str, projrot_str, chkd_sym_nums = _make_1dhr_tors_strs(
            zma, rxn_class, ts_bnd, ref_ene,
            rotor_names, rotor_grids, rotor_syms,
            constraint_dct, cnf_path,
            saddle=saddle,
            read_freqs=bool('v' in tors_model),
            hind_rot_geo=hind_rot_geo,
            build_mess=build_mess, build_projrot=build_projrot)
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
                         saddle=False, read_freqs=False, hind_rot_geo=None,
                         build_mess=True, build_projrot=True):
    """ Gather the 1DHR torsional data and gather them into a MESS file
    """

    # Initialize empty strings and torsional sym num lsts
    mess_hr_str, projrot_hr_str = '', ''
    chkd_sym_nums = []

    # Loop over the torsions
    tors_info = zip(rotor_names, rotor_grids, rotor_syms)
    for tors_names, tors_grids, tors_sym in tors_info:

        # Grab zero elment because of formatting (each list has one tors)
        # tors_name = tors_names[0]
        # tors_grid = tors_grids[0]

        # Read potential if building MESS string
        if build_mess:
            # Read the hindered rotor potential
            pot, _ = _read_hr_pot(
                tors_names, tors_grids,
                cnf_save_path, ref_ene,
                constraint_dct, read_freqs=read_freqs)

            # Use a spline fitter to deal with failures in the potential
            pot = _hrpot_spline_fitter(pot)
        else:
            pot = ()

        # Get the HR groups and axis for the rotor
        group, axis, chkd_sym_num = torsprep.set_tors_def_info(
            zma, tors_names[0], tors_sym, pot,
            ts_bnd, rxn_class, saddle=saddle)

        # Check for dummy atoms
        remdummy = geomprep.build_remdummy_shift_lst(zma)

        # Write the MESS string for the rotor
        if build_mess:
            # Write string
            mess_hr_str += mess_io.writer.rotor_hindered(
                group=group,
                axis=axis,
                symmetry=chkd_sym_num,
                potential=pot,
                remdummy=remdummy,
                geom=hind_rot_geo,
                use_quantum_weight=True,
                rotor_id=tors_names[0])

        # Write the ProjRot string for the rotor
        if build_projrot:
            projrot_hr_str += '\n' + projrot_io.writer.rotors(
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
        mdhr_dat_str_lst.append(
            mess_io.writer.mdhr_data(hr_pot, freqs=hr_freqs))

        # Check for dummy transformations
        remdummy = geomprep.build_remdummy_shift_lst(zma)

        if not isinstance(tors_syms, list):
            tors_sym_list = [tors_syms]
        else:
            tors_sym_list = tors_syms

        # Loop over the rotors in the group and write the internal rotor strs
        for tors_name, tors_sym in zip(tors_names, tors_sym_list):

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
            projrot_hr_str += '\n' + projrot_io.writer.rotors(
                axis=axis,
                group=group,
                remdummy=remdummy)

    return mess_hr_str, projrot_hr_str, mdhr_dat_str_lst, chkd_sym_nums


def _make_1dhrv_tors_strs(zma, rxn_class, ts_bnd, ref_ene):
    """ Write the strings for the 1DHRV 
    """
    pass


def make_flux_str(tors_min_cnf_locs, tors_cnf_save_fs,
                  rotor_names, rotor_syms):
    """ Write out the input string for tau sampling
    """
    # Loop over the torsions to get the flux strings
    flux_mode_str = ''
    if tors_min_cnf_locs is not None:

        # Get geometry for the torsional minimum
        zma_fs = fs.manager(
            tors_cnf_save_fs[-1].path(tors_min_cnf_locs), 'ZMATRIX')
        zma = zma_fs[-1].file.zmatrix.read([0])
        name_matrix = automol.zmatrix.name_matrix(zma)
        key_matrix = automol.zmatrix.key_matrix(zma)
    
        # Write the MESS flux strings for each of the modes
        for tors_names, tors_sym in zip(rotor_names, rotor_syms):

            tors_name = tors_names[0]

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
        dims = (len(tors_grids[0]),)
    elif len(tors_names) == 2:
        dims = (len(tors_grids[0]), len(tors_grids[1]))
    elif len(tors_names) == 3:
        dims = (len(tors_grids[0]), len(tors_grids[1]), len(tors_grids[2]))
    elif len(tors_names) == 4:
        dims = (len(tors_grids[0]), len(tors_grids[1]),
                len(tors_grids[2]), len(tors_grids[3]))
    pot = numpy.zeros(dims).tolist()

    # Initialize freqs list
    freqs = numpy.zeros(dims).tolist() if read_freqs else []

    # Read the energies from the filesystem
    zma_fs = fs.manager(cnf_save_path, 'ZMATRIX')
    zma_path = zma_fs[-1].path([0])
    if constraint_dct is None:
        scn_fs = autofile.fs.scan(zma_path)
    else:
        scn_fs = autofile.fs.cscan(zma_path)
    if len(tors_names) == 1:
        for i, grid_val_i in enumerate(tors_grids[0]):
            # Set locs
            locs = [tors_names, [grid_val_i]]
            if constraint_dct is not None:
                locs.append(constraint_dct)
            # Read filesys
            # print('locs test:', locs)
            # print('scn fs test:', scn_fs[-1].path(locs))
            if scn_fs[-1].exists(locs):
                ene = scn_fs[-1].file.energy.read(locs)
                pot[i] = (ene - ref_ene) * phycon.EH2KCAL
            else:
                pot[i] = -10.0
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
                    pot[i][j] = -10.0
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
                        pot[i][j][k] = -10.0
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
                            pot[i][j][k][lma] = -10.0
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
    pot.append(0.0)

    # Build a potential list from only successful calculations
    idx_success = []
    pot_success = []
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
