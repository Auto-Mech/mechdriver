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
from lib.structure import vib as vibprep
from lib.structure import geom as geomprep


# FUNCTIONS TO BUILD ROTOR OBJECTS CONTAINING ALL NEEDED INFO
def build_rotors(spc_dct_i, pf_filesystems, pf_models,
                 rxn_class='', frm_bnd_keys=(), brk_bnd_keys=(),
                 tors_geo=True):
    """ Add more rotor info
    """

    saddle = bool(rxn_class and (frm_bnd_keys or brk_bnd_keys))

    # Set up tors level filesystem and model
    tors_model = pf_models['tors']
    [cnf_fs, cnf_save_path, min_cnf_locs, _, _] = pf_filesystems['tors']

    # Grab the zmatrix
    if min_cnf_locs is not None:
        zma_fs = fs.manager(cnf_fs[-1].path(min_cnf_locs), 'ZMATRIX')
        zma = zma_fs[-1].file.zmatrix.read([0])
        remdummy = geomprep.build_remdummy_shift_lst(zma)
        ref_ene = cnf_fs[-1].file.energy.read(min_cnf_locs)
        geo = cnf_fs[-1].file.geometry.read(min_cnf_locs) if tors_geo else None

    # Set the tors names
    rotor_inf = _rotor_info(
        zma, spc_dct_i, cnf_fs, min_cnf_locs, tors_model,
        frm_bnd_keys=frm_bnd_keys, brk_bnd_keys=brk_bnd_keys)

    # Build constraint dct
    if tors_model in ['1dhrf']:
        constraint_dct = torsprep.build_constraint_dct(
            zma, rotor_inf[0])
    else:
        constraint_dct = None

    # Read the potential energy surface for the rotors
    num_rotors = len(rotor_inf[0])
    rotors = []
    for tors_names, tors_grids, tors_syms in zip(*rotor_inf):

        # Initialize dct to hold info for each torsion of rotor
        rotor_dct = {}

        # Read the potential along the rotors
        if tors_model == 'mdhr':

            # Set to read additional info for vibrational adiabaticity
            if tors_model == 'mdhrv':
                read_geom, read_grad, read_hess = True, True, True
            else:
                read_geom, read_grad, read_hess = False, False, False

            # Read and MDHR potential for single MDHR rotor
            # Could be MDHR mod for sys w/ 1 Rotor
            if ((num_rotors > 1 and len(tors_names) > 1) or num_rotors == 1):
                rotor_dct['mdhr_pot_data'] = _read_hr_pot(
                    tors_names, tors_grids,
                    cnf_save_path, ref_ene,
                    constraint_dct,
                    read_geom=read_geom,
                    read_grad=read_grad,
                    read_hess=read_hess)

        for tname, tgrid, tsym in zip(tors_names, tors_grids, tors_syms):
            # Call read pot for 1DHR
            pot, _, _, _ = _read_hr_pot(
                [tname], [tgrid],
                cnf_save_path, ref_ene,
                constraint_dct)

            # Get the HR groups and axis for the rotor
            group, axis, sym_num = torsprep.set_tors_def_info(
                zma, tname, tsym, pot,
                frm_bnd_keys, rxn_class, saddle=saddle)
            remdummy = geomprep.build_remdummy_shift_lst(zma)

            # Get the indices for the torsion
            mode_idxs = automol.zmatrix.coord_idxs(zma, tname)
            mode_idxs = tuple((idx+1 for idx in mode_idxs))

            # Determine the flux span using the symmetry number
            mode_span = 360.0 / tsym

            # Build dictionary for the torsion
            keys = ['group', 'axis', 'sym_num', 'remdummy', 'pot',
                    'atm_idxs', 'span', 'hrgeo']
            vals = [group, axis, sym_num, remdummy, pot,
                    mode_idxs, mode_span, geo]
            rotor_dct[tname] = dict(zip(keys, vals))

        # Append to lst
        rotors.append(rotor_dct)

    return rotors


def _rotor_info(zma, spc_dct_i, cnf_fs, min_cnf_locs, tors_model,
                frm_bnd_keys=(), brk_bnd_keys=()):
    """ get tors stuff
    """

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

    # Obtain the info for each of the rotors
    run_tors_names = torsprep.tors_name_prep(
        spc_dct_i, cnf_fs, min_cnf_locs, tors_model)
    if run_tors_names:
        rotor_names, rotor_grids, rotor_syms = torsprep.hr_prep(
            zma=zma,
            tors_name_grps=run_tors_names,
            scan_increment=scan_increment,
            tors_model=tors_model,
            frm_bnd_keys=frm_bnd_keys,
            brk_bnd_keys=brk_bnd_keys)

    else:
        rotor_names, rotor_grids, rotor_syms = (), (), ()
        print('No tors info could be found from spcdct, filesys, or geom.')

    return rotor_names, rotor_grids, rotor_syms


# FUNCTIONS TO WRITE STRINGS FOR THE ROTORS FOR VARIOUS SITUATION
def make_hr_strings(rotors, run_path, tors_model,
                    mess_hr=True, mess_ir=True,
                    mess_flux=True, projrot=True):
    """ Procedure for building the MESS strings
        :return mess_allrot_str: combination of intl and hr strs
        :return mess_hr_str: hr strs
    """

    mess_allr_str, mess_hr_str, mess_flux_str, projrot_str = '', '', '', ''
    mdhr_dat = ''
    numrotors = len(rotors)
    for rotor in rotors:

        # Set some options for writing
        # if len(rotor) == 1:

        # Write the strings for each torsion of the rotor
        for tors_name, tors_dct in rotor.items():
            if 'D' in tors_name:
                tors_strs = _rotor_tors_strs(
                    tors_name, tors_dct['group'], tors_dct['axis'],
                    tors_dct['sym_num'], tors_dct['pot'],
                    tors_dct['remdummy'], tors_dct['hrgeo'],
                    tors_dct['atm_idxs'], tors_dct['span'],
                    mess_hr=mess_hr, mess_ir=mess_ir,
                    mess_flux=mess_flux, projrot=projrot)
                mess_hr_str += tors_strs[0]
                mess_flux_str += tors_strs[2]
                projrot_str += tors_strs[3]

                # For MDHR, add the appropriate string
                if 'mdhr' in tors_model:
                    if ((numrotors > 1 and len(rotor) > 1) or numrotors == 1):
                        mess_allr_str += tors_strs[1]
                    else:
                        mess_allr_str += tors_strs[0]
                else:
                    mess_allr_str += tors_strs[0]

        # Write the MDHR potential file string for the one MDHR
        # if len(rotor) > 1 and 'mdhr_pot_data' in rotor:
        if 'mdhr_pot_data' in rotor:
            pot, geoms, grads, hessians = rotor['mdhr_pot_data']
            hr_freqs = _calc_hr_frequenices(geoms, grads, hessians, run_path)
            mdhr_dat = mess_io.writer.mdhr_data(pot, freqs=hr_freqs)

    return mess_allr_str, mess_hr_str, mess_flux_str, projrot_str, mdhr_dat


def _rotor_tors_strs(tors_name, group, axis,
                     sym_num, pot, remdummy,
                     hr_geo, mode_idxs, mode_span,
                     mess_hr=True, mess_ir=True,
                     mess_flux=True, projrot=True):
    """ Gather the 1DHR torsional data and gather them into a MESS file
    """

    mess_hr_str = ''
    if mess_hr:
        mess_hr_str = mess_io.writer.rotor_hindered(
            group=group,
            axis=axis,
            symmetry=sym_num,
            potential=pot,
            remdummy=remdummy,
            geom=hr_geo,
            use_quantum_weight=True,
            rotor_id=tors_name)

    mess_ir_str = ''
    if mess_ir:
        mess_ir_str = mess_io.writer.mol_data.rotor_internal(
            group=group,
            axis=axis,
            symmetry=sym_num,
            grid_size=100,
            mass_exp_size=5,
            pot_exp_size=5,
            hmin=13,
            hmax=101,
            remdummy=remdummy,
            geom=None,
            rotor_id=tors_name)

    mess_flux_str = ''
    if mess_flux:
        mess_flux_str = mess_io.writer.fluxional_mode(
            mode_idxs, span=mode_span)

    projrot_str = ''
    if projrot:
        projrot_str = projrot_io.writer.rotors(
            axis=axis,
            group=group,
            remdummy=remdummy)

    return mess_hr_str, mess_ir_str, mess_flux_str, projrot_str


# Functions to obtain values of the HR potentials from the filesystem
def _read_hr_pot(tors_names, tors_grids, cnf_save_path, ref_ene,
                 constraint_dct,
                 read_geom=False, read_grad=False, read_hess=False):
    """ Get the potential for a hindered rotor
    """

    # print('cscn_path', scn_run_fs[1].path([coo_names]))

    # Build initial lists for storing potential energies and Hessians
    grid_points, grid_vals = torsprep.set_hr_dims(tors_grids)
    pot, geoms, grads, hessians = {}, {}, {}, {}

    # Set up filesystem information
    zma_fs = fs.manager(cnf_save_path, 'ZMATRIX')
    zma_path = zma_fs[-1].path([0])
    if constraint_dct is None:
        scn_fs = autofile.fs.scan(zma_path)
    else:
        scn_fs = autofile.fs.cscan(zma_path)

    # Read the energies and Hessians from the filesystem
    for point, vals in zip(grid_points, grid_vals):

        locs = [tors_names, vals]
        if constraint_dct is not None:
            locs.append(constraint_dct)

        # print('path', scn_fs[-1].path(locs))

        if scn_fs[-1].exists(locs):
            ene = scn_fs[-1].file.energy.read(locs)
            pot[point] = (ene - ref_ene) * phycon.EH2KCAL
        else:
            pot[point] = -10.0

        if read_geom:
            geoms[point] = scn_fs[-1].file.geometry.read(locs)

        if read_grad:
            grads[point] = scn_fs[-1].file.gradient.read(locs)

        if read_hess:
            hessians[point] = scn_fs[-1].file.hessian.read(locs)

    return pot, geoms, grads, hessians


def _calc_hr_frequenices(geoms, grads, hessians, run_path):
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
