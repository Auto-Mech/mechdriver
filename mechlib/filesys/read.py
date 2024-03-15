""" drivers for coordinate scans
"""
import itertools
import numpy
from scipy.interpolate import CubicSpline
from scipy.interpolate import Akima1DInterpolator

import automol
import autofile
from phydat import phycon
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo
from mechanalyzer.inf import rxn as rinfo
from mechlib.filesys._build import build_fs
from mechlib.filesys.mincnf import min_energy_conformer_locators


def potential(names, grid_vals, cnf_save_path,
              mod_tors_ene_info, ref_ene,
              constraint_dct,
              read_geom=False, read_grad=False,
              read_hess=False, read_zma=False,
              read_energy_backstep=True,
              remove_bad_points=True):
    """ Get the potential for a hindered rotor
    """

    print('potential test:')
    print('names', names)
    print('grids', grid_vals)

    # Build initial lists for storing potential energies and Hessians
    # grid_points = automol.pot.points(grid_vals)
    grid_coords = tuple(itertools.product(*grid_vals))
    back_coords = tuple(tuple(val + 4*numpy.pi for val in grid)
                        for grid in grid_coords)
    pot, geoms, grads, hessians, zmas, paths = {}, {}, {}, {}, {}, {}

    # Set up filesystem information
    zma_fs = autofile.fs.zmatrix(cnf_save_path)
    zma_path = zma_fs[-1].path([0])
    if constraint_dct is None:
        scn_fs = autofile.fs.scan(zma_path)
    else:
        scn_fs = autofile.fs.cscan(zma_path)
    # Read the filesystem
    for idx, vals in enumerate(grid_coords):

        # Get angles in degrees for potential for now
        vals_conv = tuple(val*phycon.RAD2DEG for val in vals)

        # Get locs for reading filesysten
        locs = [names, vals]
        back_locs = [names, back_coords[idx]]
        if constraint_dct is not None:
            locs = [constraint_dct] + locs
            back_locs = [constraint_dct] + back_locs
        # print('path', scn_fs[-1].path(locs))

        # Read values of interest
        ene = energy(scn_fs, locs, mod_tors_ene_info)
        print("@AVC Energy read:", ene)
        if read_energy_backstep:
            back_ene = energy(scn_fs, back_locs, mod_tors_ene_info)
            step_ene = None
            if ene is not None:
                if back_ene is not None:
                    step_ene = min(ene, back_ene)
                else:
                    step_ene = ene
            elif back_ene is not None:
                step_ene = back_ene

            if step_ene is not None:
                enediff = (step_ene - ref_ene) * phycon.EH2KCAL
                if idx == 0:
                    if enediff > 0.05:
                        print('Warning the first potential value does not',
                              f'match the reference energy {enediff:.2f}')
                    ref_ene = step_ene
                    enediff = 0
                pot[vals_conv] = enediff
            else:
                pot[vals_conv] = None
        else:
            pot[vals_conv] = (ene - ref_ene) * phycon.EH2KCAL

        if read_geom:
            if scn_fs[-1].file.geometry.exists(locs):
                geoms[vals_conv] = scn_fs[-1].file.geometry.read(locs)
            else:
                geoms[vals_conv] = None

        if read_grad:
            if scn_fs[-1].file.gradient.exists(locs):
                grads[vals_conv] = scn_fs[-1].file.gradient.read(locs)
            else:
                grads[vals_conv] = None

        if read_hess:
            if scn_fs[-1].file.hessian.exists(locs):
                hessians[vals_conv] = scn_fs[-1].file.hessian.read(locs)
            else:
                hessians[vals_conv] = None

        if read_zma:
            if scn_fs[-1].file.zmatrix.exists(locs):
                zmas[vals_conv] = scn_fs[-1].file.zmatrix.read(locs)
            else:
                zmas[vals_conv] = None

        paths[vals] = scn_fs[-1].path(locs)

    # If potential has any terms that are not None, ID and remove bad points
    if remove_bad_points and len(names) == 1:
        if any(val is not None for val in pot.values()):
            pot = {k: v for k, v in pot.items() if v is not None}
            bad_angle = identify_bad_point(pot)
            if bad_angle is not None:
                pot = remove_bad_point(pot, bad_angle)
                pot = {k: v for k, v in pot.items() if v is not None}
        else:
            pot = {}

    return pot, geoms, grads, hessians, zmas, paths


def identify_bad_point(pot, thresh=0.05):
    """ Identifies a single bad point in a torsional potential based on a
        comparison of Akima and cubic spline fits
    """

    vals_conv = pot.keys()
    step_enes = pot.values()

    # Get the sorted angles
    shifted_angles = []
    for idx, angle in enumerate(vals_conv):
        if len(angle) == 1:  # if a tuple, just take the first value
            angle = angle[0]
        if idx == 0:
            start_angle = angle
        angle = angle - start_angle
        if angle > 180:
            angle = angle - 360
        shifted_angles.append(angle)
    shifted_angles = numpy.array(shifted_angles)

    # For methyl rotors, double the threshold
    if len(shifted_angles) == 4:
        thresh *= 2

    # Get the potentials and then sort them according to increasing angle
    step_enes = numpy.array(list(step_enes))
    sorted_idxs = numpy.argsort(shifted_angles)
    sorted_angles = shifted_angles[sorted_idxs]
    sorted_potentials = step_enes[sorted_idxs]

    # Fit cubic and Akima splines
    cub_spline = CubicSpline(sorted_angles, sorted_potentials)
    akima_spline = Akima1DInterpolator(sorted_angles, sorted_potentials)

    # Evaluate the splines on a fine grid to check for ringing
    fine_grid = numpy.arange(min(sorted_angles), max(sorted_angles), 1)
    diff = cub_spline(fine_grid) - akima_spline(fine_grid)
    max_fine_angle = fine_grid[numpy.argmax(diff)]
    max_norm_diff = max(diff) / max(step_enes)  # normalized by max potential

    # Remove the bad point
    bad_angle = None
    print('max norm diff of splines ', max_norm_diff, thresh)
    if max_norm_diff > thresh:
        max_idxs = numpy.argsort(abs(shifted_angles - max_fine_angle))[:2]
        suspect_enes = [step_enes[idx] for idx in max_idxs]
        bad_angle = shifted_angles[max_idxs[numpy.argmax(suspect_enes)]]
        if bad_angle < 0:
            bad_angle += 360  # convert back to original angle
        bad_angle += start_angle

    return bad_angle


def remove_bad_point(pot, bad_angle):
    """ Remove a single bad angle from a potential
    """

    # Find angle in pot that is within 0.1 degrees of bad_angle
    # to read the dictionary
    # Only works for 1D
    bad_tuple = None
    for angle in pot.keys():
        if numpy.isclose(angle, bad_angle, atol=0.1):
            bad_tuple = (angle,)

    assert bad_tuple is not None, (
        f'Angle {bad_angle*phycon.DEG2RAD} does not exist in pot dictionary')

    print(f'Removing bad angle at {bad_angle} degrees')

    return pot


# Single data point readers
def geometry(cnf_save_fs, mod_thy_info, conf='sphere', hbond_cutoffs=None):
    """ get the geometry
    """

    assert conf in ('minimum', 'sphere')

    # Read the file system
    if conf == 'minimum':
        geom = _min_energy_conformer(
            cnf_save_fs, mod_thy_info, hbond_cutoffs=hbond_cutoffs)
    elif conf == 'sphere':
        geom = _spherical_conformer(cnf_save_fs)

    return geom


def _min_energy_conformer(cnf_save_fs, mod_thy_info, hbond_cutoffs=None):
    """ Reads the minimum-energy conformer from the save FileSystem
    """

    ini_loc_info = min_energy_conformer_locators(
        cnf_save_fs, mod_thy_info, hbond_cutoffs=hbond_cutoffs)
    locs, path = ini_loc_info
    if path:
        min_conf = cnf_save_fs[-1].file.geometry.read(locs)
    else:
        min_conf = None

    return min_conf


def _spherical_conformer(cnf_save_fs):
    """ Reads the conformer from the Save FileSystem that is the most
        spherical.
    """

    cnf_locs_lst = cnf_save_fs[-1].existing()
    if cnf_locs_lst:
        cnf_geoms = [cnf_save_fs[-1].file.geometry.read(locs)
                     for locs in cnf_locs_lst]
        round_geom = automol.geom.minimum_volume_geometry(cnf_geoms)
    else:
        round_geom = None

    return round_geom


def energy(filesys, locs, mod_thy_info):
    """ Read the energy from an SP filesystem that is located in some
        root 'filesys object'
    """

    if filesys[-1].exists(locs):
        path = filesys[-1].path(locs)
        sp_fs = autofile.fs.single_point(path)
        if sp_fs[-1].file.energy.exists(mod_thy_info[1:4]):
            ene = sp_fs[-1].file.energy.read(mod_thy_info[1:4])
            # ioprinter.debug_message('sp_path', sp_fs[-1].path(thy_info))
        else:
            ene = None
            # ioprinter.info_message('No scan energy')
    else:
        ene = None
        # ioprinter.info_message('No scan energy')

    return ene


def reactions(rxn_info, ini_thy_info, zma_locs, save_prefix):
    """ Check if reaction exists in the filesystem and has been identified
    """

    zrxns, zmas = (), ()

    sort_rxn_info = rinfo.sort(rxn_info, scheme='autofile')
    ts_info = rinfo.ts_info(rxn_info)
    mod_ini_thy_info = tinfo.modify_orb_label(ini_thy_info, ts_info)

    rxn_fs = autofile.fs.reaction(save_prefix)
    if rxn_fs[-1].exists(sort_rxn_info):
        rxn_path = rxn_fs[-1].path(sort_rxn_info)
        _, ts_save_fs = build_fs(
            save_prefix, save_prefix, 'TRANSITION STATE',
            rxn_locs=sort_rxn_info,
            thy_locs=mod_ini_thy_info[1:])
        # print('ts save fs', ts_save_fs[0].path())
        ts_locs = ts_save_fs[-1].existing()
        if ts_locs:
            for locs in ts_save_fs[-1].existing():
                ts_path = ts_save_fs[-1].path(locs)
                print(f'    - Checking for TS info in CONFS/Z in {ts_path}')
                zrxn, zma = reaction(
                    rxn_info, ini_thy_info,
                    zma_locs, save_prefix, ts_locs=locs)
                if zrxn is not None:
                    zrxns += (zrxn,)
                    zmas += (zma,)
        else:
            print(f'No info at {rxn_path}')

    if not zrxns:
        zrxns, zmas = None, None

    return zrxns, zmas


def reaction(
        rxn_info, ini_thy_info, zma_locs, save_prefix, ts_locs=(0,),
        hbond_cutoffs=None):
    """ Check if reaction exists in the filesystem and has been identified
    """

    zrxn, zma = None, None

    sort_rxn_info = rinfo.sort(rxn_info, scheme='autofile')
    ts_info = rinfo.ts_info(rxn_info)
    mod_ini_thy_info = tinfo.modify_orb_label(ini_thy_info, ts_info)

    rxn_fs = autofile.fs.reaction(save_prefix)
    if rxn_fs[-1].exists(sort_rxn_info):
        _, cnf_save_fs = build_fs(
            save_prefix, save_prefix, 'CONFORMER',
            rxn_locs=sort_rxn_info,
            thy_locs=mod_ini_thy_info[1:],
            ts_locs=ts_locs)

        _, ini_min_cnf_path = min_energy_conformer_locators(
            cnf_save_fs, mod_ini_thy_info, hbond_cutoffs=hbond_cutoffs)
        if ini_min_cnf_path:
            zma_fs = autofile.fs.zmatrix(ini_min_cnf_path)
            if zma_fs[-1].file.reaction.exists(zma_locs):
                zrxn = zma_fs[-1].file.reaction.read(zma_locs)
                zma = zma_fs[-1].file.zmatrix.read(zma_locs)

        # For barrierless reactions with no conformer
        if zrxn is None:
            _, zma_fs = build_fs(
                save_prefix, save_prefix, 'ZMATRIX',
                rxn_locs=sort_rxn_info, ts_locs=ts_locs,
                thy_locs=mod_ini_thy_info[1:])

            if zma_fs[-1].file.reaction.exists(zma_locs):
                zrxn = zma_fs[-1].file.reaction.read(zma_locs)
                zma = zma_fs[-1].file.zmatrix.read(zma_locs)

    return zrxn, zma


def instability_transformation(spc_dct, spc_name, thy_info, save_prefix,
                               zma_locs=(0,), nprocs=1):
    """ see if a species and unstable and handle task management
    """

    spc_info = sinfo.from_dct(spc_dct[spc_name], canonical=True)
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)

    _, cnf_save_fs = build_fs(
        save_prefix, save_prefix, 'CONFORMER',
        spc_locs=spc_info,
        thy_locs=mod_thy_info[1:])

    # Check if any locs exist first?
    hbond_cutoffs = spc_dct[spc_name]['hbond_cutoffs']
    ini_loc_info = min_energy_conformer_locators(
        cnf_save_fs, mod_thy_info, hbond_cutoffs=hbond_cutoffs, nprocs=nprocs)
    _, min_cnf_path = ini_loc_info

    zma_save_fs = autofile.fs.zmatrix(min_cnf_path)

    # Check if the instability files exist
    if zma_save_fs[-1].file.instability.exists(zma_locs):
        instab_trans = zma_save_fs[-1].file.instability.read(zma_locs)
        zma = zma_save_fs[-1].file.zmatrix.read(zma_locs)
        _instab = (instab_trans, zma)
        path = zma_save_fs[-1].file.instability.path(zma_locs)
    else:
        _instab = None
        path = None

    return _instab, path


def energy_trans(etrans_save_fs, etrans_locs):
    """ Read out the the enery transfer parameters in the filesys
    """

    nsamp = etrans_save_fs[-1].file.read(etrans_locs)
    epsilon = etrans_save_fs[-1].file.read.epsilon(etrans_locs)
    sigma = etrans_save_fs[-1].file.read.sigma(etrans_locs)
    # alpha = etrans_fs[-1].file.alpha.read(etrans_locs)

    min_geo_traj = etrans_save_fs[-1].file.read.min_geos(etrans_locs)
    min_geos = automol.geom.from_xyz_trajectory_string(min_geo_traj)

    return nsamp, epsilon, sigma, min_geos
