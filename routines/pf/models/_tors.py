"""
  Functions handling hindered rotor model calculations
"""

import itertools
import numpy
from scipy.interpolate import interp1d
import automol
import mess_io
import projrot_io
from autofile import fs
from lib import filesys
from lib.phydat import phycon
from lib.structure import tors as torsprep
from lib.structure import geom as geomprep


# FUNCTIONS TO BUILD ROTOR OBJECTS CONTAINING ALL NEEDED INFO
def build_rotors(spc_dct_i, pf_filesystems, pf_models, pf_levels,
                 rxn_class='', frm_bnd_keys=(), brk_bnd_keys=()):
    """ Add more rotor info
    """

    saddle = bool(rxn_class and (frm_bnd_keys or brk_bnd_keys))

    # Set up tors level filesystem and model and level
    tors_model = pf_models['tors']
    tors_ene_info = pf_levels['tors'][1][0]
    mod_tors_ene_info = filesys.inf.modify_orb_restrict(
        filesys.inf.get_spc_info(spc_dct_i), tors_ene_info)
    [cnf_fs, cnf_save_path, min_cnf_locs, _, _] = pf_filesystems['tors']

    # Grab the zmatrix
    if min_cnf_locs is not None:
        zma_fs = fs.manager(cnf_fs[-1].path(min_cnf_locs), 'ZMATRIX')
        zma = zma_fs[-1].file.zmatrix.read([0])
        remdummy = geomprep.build_remdummy_shift_lst(zma)
        geo = cnf_fs[-1].file.geometry.read(min_cnf_locs)

        # Read the reference energy
        ref_ene = torsprep.read_tors_ene(
            cnf_fs, min_cnf_locs, mod_tors_ene_info)

    # Set the tors names
    rotor_inf = _rotor_info(
        zma, spc_dct_i, cnf_fs, min_cnf_locs, tors_model,
        frm_bnd_keys=frm_bnd_keys, brk_bnd_keys=brk_bnd_keys)

    # Read the potential energy surface for the rotors
    num_rotors = len(rotor_inf[0])
    rotors = []
    for tors_names, tors_grids, tors_syms in zip(*rotor_inf):

        # Initialize dct to hold info for each torsion of rotor
        rotor_dct = {}

        # Read the potential along the rotors
        if tors_model in ('mdhr', 'mdhrv'):

            # Set to read additional info for vibrational adiabaticity
            if tors_model == 'mdhrv':
                read_geom, read_grad, read_hess = True, True, True
            else:
                read_geom, read_grad, read_hess = False, False, False

            # Read and MDHR potential for single MDHR rotor
            # Could be MDHR mod for sys w/ 1 Rotor
            if ((num_rotors > 1 and len(tors_names) > 1) or num_rotors == 1):
                pot, geoms, grads, hessians, _ = torsprep.read_hr_pot(
                    tors_names, tors_grids,
                    cnf_save_path,
                    mod_tors_ene_info, ref_ene,
                    constraint_dct=None,   # No extra frozen treatments
                    read_geom=read_geom,
                    read_grad=read_grad,
                    read_hess=read_hess)
                rotor_dct['mdhr_pot_data'] = (pot, geoms, grads, hessians)

        for tname, tgrid, tsym in zip(tors_names, tors_grids, tors_syms):

            # Build constraint dct
            if tors_model == '1dhrf':
                tname_tup = tuple([tname])
                const_names = tuple(itertools.chain(*rotor_inf[0]))
                constraint_dct = torsprep.build_constraint_dct(
                    zma, const_names, tname_tup)
            elif tors_model == '1dhrfa':
                coords = list(automol.zmatrix.coordinates(zma))
                const_names = tuple(coord for coord in coords)
                tname_tup = tuple([tname])
                constraint_dct = torsprep.build_constraint_dct(
                        E
                    zma, const_names, tname_tup)
            else:
                constraint_dct = None

            # Call read pot for 1DHR
            pot, _, _, _, _ = torsprep.read_hr_pot(
                [tname], [tgrid],
                cnf_save_path,
                mod_tors_ene_info, ref_ene,
                constraint_dct)
            pot = _hrpot_spline_fitter(pot, min_thresh=-0.10, max_thresh=50.0)

            # Get the HR groups and axis for the rotor
            group, axis, pot, sym_num = torsprep.set_tors_def_info(
                zma, tname, tsym, pot,
                frm_bnd_keys, brk_bnd_keys,
                rxn_class, saddle=saddle)
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

        # print('rotors pot test:', pot)
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
                    scale_factor=None,
                    mess_hr=True, mess_ir=True,
                    mess_flux=True, projrot=True):
    """ Procedure for building the MESS strings
        :return mess_allrot_str: combination of intl and hr strs
        :return mess_hr_str: hr strs
    """

    mess_allr_str, mess_allr_nogeo_str = '', ''
    mess_hr_str, mess_flux_str, projrot_str = '', '', ''
    mdhr_dat = ''
    numrotors = len(rotors)
    for rotor in rotors:
        # Set some options for writing
        # if len(rotor) == 1:

        # Write the strings for each torsion of the rotor
        for tors_name, tors_dct in rotor.items():
            if 'D' in tors_name:

                if scale is None:
                    pot = tors_dct['pot']
                else:
                    pot = _scale_pot(tors_dct['pot'], scale_factor) 

                tors_strs = _rotor_tors_strs(
                    tors_name, tors_dct['group'], tors_dct['axis'],
                    tors_dct['sym_num'], pot,
                    tors_dct['remdummy'], tors_dct['hrgeo'],
                    tors_dct['atm_idxs'], tors_dct['span'],
                    scale=None,
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
            hr_freqs = torsprep.calc_hr_frequencies(
                geoms, grads, hessians, run_path)
            mdhr_dat = mess_io.writer.mdhr_data(
                pot, freqs=hr_freqs, nrot=numrotors)

    return mess_allr_str, mess_hr_str, mess_flux_str, projrot_str, mdhr_dat


def _rotor_tors_strs(tors_name, group, axis,
                     sym_num, pot, remdummy,
                     hr_geo, mode_idxs, mode_span,
                     mess_hr=True, mess_ir=True,
                     mess_flux=True, projrot=True):
    """ Gather the 1DHR torsional data and gather them into a MESS file
    """

    mess_hr_str, mess_hr_nogeo_str = '', ''
    if mess_hr:
        mess_hr_str = mess_io.writer.rotor_hindered(
            group=group,
            axis=axis,
            symmetry=sym_num,
            potential=pot,
            hmin=None,
            hmax=None,
            lvl_ene_max=None,
            therm_pow_max=None,
            remdummy=remdummy,
            geom=hr_geo,
            rotor_id=tors_name)

        mess_hr_nogeo_str = mess_io.writer.rotor_hindered(
            group=group,
            axis=axis,
            symmetry=sym_num,
            potential=pot,
            hmin=None,
            hmax=None,
            lvl_ene_max=None,
            therm_pow_max=None,
            remdummy=remdummy,
            geom=None,
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

    return (mess_hr_str, mess_ir_str,
            mess_flux_str, projrot_str)
    # return (mess_hr_str, mess_hr_nogeo_str, mess_ir_str,
    #         mess_flux_str, projrot_str)


def _need_tors_geo(pf_levels):
    """ Determine if a torsional geometry is geometry if needed
    """
    print('pflvl', pf_levels)
    return bool(pf_levels['tors'][1] == pf_levels['harm'])


def _scale_pot(pot, scale_factor):
    """ Scale the potential
    """

    new_pot = {}
    for idx, val in pot.items():
        new_pot[(idx,)] = pot[(idx,)] * scale_factor
        
    return new_pot 


def _hrpot_spline_fitter(pot_dct, min_thresh=-0.10, max_thresh=50.0):
    """ Get a physical hindered rotor potential via a series of spline fits
    """

    pot = list(pot_dct.values())

    # Initialize a variable for the size of the potential
    lpot = len(pot)+1
    pot.append(0.0)

    # Print warning messages
    print_pot = False
    if any(val > max_thresh for val in pot):
        print_pot = True
        max_pot = max(pot)
        print('Warning: Found pot val of {0:.2f}'.format(max_pot),
              ' which is larger than',
              'the typical maximum for a torsional potential')
    if any(val < min_thresh for val in pot):
        print_pot = True
        min_pot = min(pot)
        print('Warning: Found pot val of {0:.2f}'.format(min_pot),
              ' which is below',
              '{0} kcal. Refit w/ positives'.format(min_thresh))
    if print_pot:
        print('Potential before spline:', pot)

    # Build a potential list from only successful calculations
    # First replace high potential values with max_thresh
    # Then replace any negative potential values cubic spline fit values
    idx_success = []
    pot_success = []
    for idx in range(lpot):
        if pot[idx] < 600. and pot[idx] > min_thresh:
            idx_success.append(idx)
            if pot[idx] < max_thresh:
                pot_success.append(pot[idx])
            else:
                pot_success.append(max_thresh)

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

    fin_dct = {}
    for i, val in enumerate(final_potential):
        fin_dct[(i,)] = val

    return fin_dct
