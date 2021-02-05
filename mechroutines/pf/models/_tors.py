"""
  Functions handling hindered rotor model calculations
"""

import itertools
import numpy
from scipy.interpolate import interp1d
import automol
import mess_io
import projrot_io
from phydat import phycon
from autofile import fs
from mechlib import filesys
from mechlib.structure import tors as torsprep
from mechlib.amech_io import printer as ioprinter


# FUNCTIONS TO BUILD ROTOR OBJECTS CONTAINING ALL NEEDED INFO
def build_rotors(spc_dct_i, pf_filesystems, pf_models, pf_levels):
    """ Add more rotor info
    """

    # Set up tors level filesystem and model and level
    ioprinter.debug_message('pflvls', pf_levels)
    tors_model = pf_models['tors']
    tors_ene_info = pf_levels['tors'][1][0]
    [cnf_fs, cnf_save_path, min_cnf_locs, _, _] = pf_filesystems['tors']

    # Build the rotors
    if min_cnf_locs is not None:
        zma_fs = fs.zmatrix(cnf_fs[-1].path(min_cnf_locs))
        zma = zma_fs[-1].file.zmatrix.read([0])
        if zma_fs[-1].file.torsions.exists([0]):
            tors_dct = zma_fs[-1].file.torsions.read([0])
            tors_names = spc_dct_i.get('tors_names', None)
            rotors = automol.rotor.from_data(
                zma, tors_dct,
                tors_names=tors_names, multi=bool('1d' in tors_model))
        else:
            rotors = ()

    # Read the potential grids
    # rotors = _read_potentials

    return rotors


def _read_potentials(rotors, spc_dct_i, run_path, tors_model,
                     scale_factor=((), None)):
    """ read out the potentials
    """
    
    # Convert the rotor objects indexing to be in geoms
    increment = spc_dct_i.get('hind_inc', 30.0*phycon.DEG2RAD)
    rotor_zma = automol.rotor.zmatrix(rotors)
    tors_names = automol.rotor.names(rotors, flat=True)
    tors_grids = automol.rotor.grids(rotors, increment=increment, flat=True)

    # Count numbers
    numrotors = len(rotors)
    numtors = 0
    for rotor in rotors:
        numtors += len(rotor)

    # Calculate the scaling factors
    scale_indcs, factor = scale_factor
    nscale = numtors - len(scale_indcs)
    sfactor = factor**(2.0/nscale)
    ioprinter.debug_message(
        'scale_coeff test:', factor, nscale, sfactor)

    for ridx, rotor in enumerate(rotors):
        multi_idx = ridx
        for tidx, torsion in enumerate(rotor):

            # Read and spline-fit potential
            const_names = torsprep.set_constraint_names(
                rotor_zma, tors_names, tors_model)
            constraint_dct = torsprep.build_constraint_dct(
                rotor_zma, const_names, tors_names)
            pot, _, _, _, _, _ = torsprep.read_hr_pot(
                tors_names[tidx], tors_grids[tidx],
                cnf_save_path,
                mod_tors_ene_info, ref_ene,
                constraint_dct)
            pot = _hrpot_spline_fitter(
                pot, min_thresh=-0.0001, max_thresh=50.0)
            # Add potential to potential
            torsion.pot = pot

    if multi_idx is not None:
        tors_name = automol.rotor.names(rotors)[multi_idx]
        tors_grid = automol.rotor.grids(rotors, increment=increment)[multi_idx]

        pot, geoms, grads, hessians, _, _ = torsprep.read_hr_pot(
            tors_names, tors_grid,
            cnf_save_path,
            mod_tors_ene_info, ref_ene,
            constraint_dct=None,   # No extra frozen treatments
            read_geom=bool('v' in tors_model),
            read_grad=bool('v' in tors_model),
            read_hess=bool('v' in tors_model))
        rotor_dct['mdhr_pot_data'] = (pot, geoms, grads, hessians)

        if 'mdhr_pot_data' in rotor:
            pot, geoms, grads, hessians = rotor['mdhr_pot_data']
            hr_freqs = torsprep.calc_hr_frequencies(
                geoms, grads, hessians, run_path)
    else:
        pass 

    return rotors, (pot, freqs)


# FUNCTIONS TO WRITE STRINGS FOR THE ROTORS FOR VARIOUS SITUATION
def make_hr_strings(rotors, spc_dct_i, run_path, tors_model,
                    scale_factor=((), None)):
    """ Procedure for building the MESS strings
        :return mess_allrot_str: combination of intl and hr strs
        :return mess_hr_str: hr strs
    """

    # Initialize empty strings
    mess_allr_str  = ''
    mess_hr_str, mess_flux_str, projrot_str = '', '', ''
    mdhr_dat = ''

    # Convert the rotor objects indexing to be in geoms
    increment = spc_dct_i.get('hind_inc', 30.0*phycon.DEG2RAD)
    rotors, geo = automol.rotor.relabel_for_geometry(rotors)

    for ridx, rotor in enumerate(rotors):
        multirotor = bool(len(rotor) > 1)
        multi_idx = ridx
        for tidx, torsion in enumerate(rotor):

            # Write the rotor strings
            hr_str, ir_str, flux_str, prot_str = _tors_strs(torsion, geo)
            mess_hr_str += hr_str
            mess_flux_str += flux_str
            projrot_str += prot_str

            # For MDHR, add the appropriate string
            if 'mdhr' in tors_model:
                if ((numrotors > 1 and multirotor) or numrotors == 1):
                    mess_allr_str += ir_str
                else:
                    mess_allr_str += hr_str
            else:
                mess_allr_str += hr_str

    # write the mdhr dat string
    if mdhr_idx is not None:
        mdhr_dat = _mdhr_dat_str(rotors, spc_dct_i, mdhr_idx)
    else:
        mdhr_dat = ''

    return mess_allr_str, mess_hr_str, mess_flux_str, projrot_str, mdhr_dat
# Build prgram strings containing the rotor info

def _tors_strs(torsion, geo):
    """ Gather the 1DHR torsional data and gather them into a MESS file
    """

    mess_hr_str = mess_io.writer.rotor_hindered(
        group=torsion.groups[0],
        axis=torsion.axis[0],
        symmetry=torsion.symmetry,
        potential=torsion.pot,
        hmin=None,
        hmax=None,
        lvl_ene_max=None,
        therm_pow_max=None,
        geo=geo,
        rotor_id=torsion.name)

    mess_ir_str = mess_io.writer.mol_data.rotor_internal(
        group=torsion.groups[0],
        axis=torsion.axis[0],
        symmetry=torsion.symmetry,
        grid_size=100,
        mass_exp_size=5,
        pot_exp_size=5,
        hmin=13,
        hmax=101,
        geo=None,
        rotor_id=torsion.name)

    mess_flux_str = mess_io.writer.fluxional_mode(
        torsion.atom_idxs,
        span=torsion.span)

    projrot_str = projrot_io.writer.rotors(
        axis=torsion.axis[0],
        group=torsion.groups[0])

    return mess_hr_str, mess_ir_str, mess_flux_str, projrot_str


def _mdhr_dat_str(rotors, spc_dct_i):
    """ Build the data string for the MDHR
    """

    increment = spc_dct_i.get('hind_inc', 30.0*phycon.DEG2RAD)
    rotor_dims = automol.rotor.dimensions(rotors)

    # Read the potential along the rotors
    if tors_model in ('mdhr', 'mdhrv'):
        if tors_model == 'mdhrv':
            read_geom, read_grad, read_hess = True, True, True
        else:
            read_geom, read_grad, read_hess = False, False, False

        # Write the MDHR potential file string for the one MDHR
        # if len(rotor) > 1 and 'mdhr_pot_data' in rotor:
        if 'mdhr_pot_data' in rotor:
            pot, geoms, grads, hessians = rotor['mdhr_pot_data']
            hr_freqs = torsprep.calc_hr_frequencies(
                geoms, grads, hessians, run_path)
            mdhr_dat = mess_io.writer.mdhr_data(
                pot, freqs=hr_freqs, nrot=numrotors)

    return mdhr_dat


def _need_tors_geo(pf_levels):
    """ Determine if a torsional geometry is geometry if needed
    """
    ioprinter.debug_message('pflvl', pf_levels)
    return bool(pf_levels['tors'][1] == pf_levels['harm'])


def _hrpot_spline_fitter(pot_dct, min_thresh=-0.0001, max_thresh=50.0):
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
        ioprinter.warning_message(
            'Found pot val of {0:.2f}'.format(max_pot),
            ' which is larger than',
            'the typical maximum for a torsional potential')
    # reset any negative values for the first grid point to 0.
    if pot[0] < 0.:
        ioprinter.error_message('The first potential value should be 0.')
        pot[0] = 0.
    if any(val < min_thresh for val in pot):
        print_pot = True
        min_pot = min(pot)
        ioprinter.warning_message(
            'Found pot val of {0:.2f}'.format(min_pot),
            ' which is below',
            '{0} kcal. Refit w/ positives'.format(min_thresh))

    if print_pot:
        ioprinter.debug_message('Potential before spline:', pot)

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

    if len(pot_success) > 3:
    # Build a new potential list using a spline fit of the HR potential
        pot_spl = interp1d(
            numpy.array(idx_success), numpy.array(pot_success), kind='cubic')
        for idx in range(lpot):
            pot[idx] = float(pot_spl(idx))

    # Do second spline fit of only positive values if any negative values found
    if any(val < min_thresh for val in pot):
        ioprinter.warning_message(
            'Still found negative potential values after first spline')
        ioprinter.debug_message('Potential after spline:', pot)
        if len(pot_success) > 3:
            x_pos = numpy.array([i for i in range(lpot)
                                 if pot[i] >= min_thresh])
            y_pos = numpy.array([pot[i] for i in range(lpot)
                                 if pot[i] >= min_thresh])
            pos_pot_spl = interp1d(x_pos, y_pos, kind='cubic')
            pot_pos_fit = []
            for idx in range(lpot):
                pot_pos_fit.append(pos_pot_spl(idx))
        else:
            pot_pos_fit = []
            for idx in range(lpot):
                pot_pos_fit.append(pot[idx])

        ioprinter.debug_message('Potential after spline:', pot_pos_fit)
        # Perform second check to see if negative potentials have been fixed
        if any(val < min_thresh for val in pot_pos_fit):
            ioprinter.warning_message('Still found negative potential values after second spline')
            ioprinter.info_message('Replace with linear interpolation of positive values')
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
        val_fin = min(val, max_thresh)
        fin_dct[(i,)] = val_fin

    return fin_dct
