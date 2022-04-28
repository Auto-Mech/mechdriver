"""
  Functions handling hindered rotor model calculations
"""

import automol
import autorun
import mess_io
import projrot_io
from phydat import phycon
from autofile import fs
from mechanalyzer.inf import thy as tinfo
from mechanalyzer.inf import spc as sinfo
from mechlib.amech_io import printer as ioprinter
from mechlib.amech_io import job_path
from mechlib import filesys


# FUNCTIONS TO BUILD ROTOR OBJECTS CONTAINING ALL NEEDED INFO
def build_rotors(spc_dct_i, pf_filesystems, spc_mod_dct_i,
                 read_potentials=True):
    """ Add more rotor info
    """

    run_prefix = pf_filesystems['run_prefix']
    spc_info = sinfo.from_dct(spc_dct_i)
    spc_fml = automol.inchi.formula_string(spc_info[0])
    if spc_fml is None:
        spc_fml = 'TS'
    run_path = job_path(run_prefix, 'PROJROT', 'FREQ', spc_fml, locs_id=None)

    # Set up tors level filesystem and model and level
    tors_model = spc_mod_dct_i['tors']['mod']
    tors_ene_info = spc_mod_dct_i['tors']['enelvl'][1][1]
    mod_tors_ene_info = tinfo.modify_orb_label(
        tors_ene_info, sinfo.from_dct(spc_dct_i))

    rotors, mdhr_dct = None, None
    print('tors_model', tors_model)
    if pf_filesystems['tors'] is not None:
        [cnf_fs, cnf_save_path, min_cnf_locs, _, _] = pf_filesystems['tors']

        # Build the rotors
        ref_ene = filesys.read.energy(cnf_fs, min_cnf_locs, mod_tors_ene_info)
        zma_fs = fs.zmatrix(cnf_fs[-1].path(min_cnf_locs))
        if (
            zma_fs[-1].file.torsions.exists([0]) and
            zma_fs[-1].file.zmatrix.exists([0]) and
            tors_model != 'rigid'
        ):
            rotors = automol.rotor.from_data(
                zma=zma_fs[-1].file.zmatrix.read([0]),
                tors_inf_dct=zma_fs[-1].file.torsions.read([0]),
                tors_names=(
                    spc_dct_i.get('tors_names', None)
                    if 'md' in tors_model else None),
                multi='md' in tors_model)
        # Read the potential grids
        if read_potentials and rotors is not None:
            rotors, mdhr_dct = _read_potentials(
                rotors, spc_dct_i, run_path, cnf_save_path,
                ref_ene, mod_tors_ene_info,
                tors_model)

    return rotors, mdhr_dct


def _read_potentials(rotors, spc_dct_i, run_path, cnf_save_path,
                     ref_ene, mod_tors_ene_info,
                     tors_model):
    """ read out the potentials
    """

    _ = run_path

    # Convert the rotor objects indexing to be in geoms
    increment = spc_dct_i.get('hind_inc', 30.0*phycon.DEG2RAD)
    rotor_zma = automol.rotor.zmatrix(rotors)

    # Determine base-line rotor non-specific info for constraints
    all_tors_names = automol.rotor.names(rotors)
    print('rotor names', all_tors_names)
    const_names = automol.zmat.set_constraint_names(
        rotor_zma, all_tors_names, tors_model)

    # Recalculate the rotor potential grids using desired increment
    rotor_grids = automol.rotor.grids(rotors, increment=increment)

    multi_idx = None
    for ridx, rotor in enumerate(rotors):
        if len(rotor) > 1:
            multi_idx = ridx

        for tidx, torsion in enumerate(rotor):
            # Read and spline-fit potential
            tors_grid = rotor_grids[ridx][tidx]
            constraint_dct = automol.zmat.constraint_dct(
                rotor_zma, const_names, (torsion.name,))
            pot, _, _, _, _, _ = filesys.read.potential(
                (torsion.name,), (tors_grid,),
                cnf_save_path,
                mod_tors_ene_info, ref_ene,
                constraint_dct,
                read_energy_backstep=True,
                remove_bad_points=True)
            if pot:
                # fit_pot = automol.pot.setup_1d_potential(
                fit_pot = automol.pot.fit_1d_potential(
                    pot, min_thresh=-0.0001, max_thresh=50.0)
                # Scale pot relative to first time
                # Code here dies if pot is empty
                final_pot = {}
                gridvals = tuple(pot.keys())
                ref_val = gridvals[0][0]
                # print(f'pot test: {torsion.name}')
                # print('gridvals', gridvals)
                # print('ref_val', ref_val)
                # print('fit_pot', fit_pot)
                for i, val in enumerate(gridvals):
                    final_pot[(val[0] - ref_val,)] = fit_pot[(i,)]
                # for idx, val in fit_pot.items():
                #     final_pot[(gridvals[idx[0]][0] - ref_val,)] = val
                torsion.pot = final_pot
            else:
                torsion.pot = pot

    if multi_idx is not None:
        is_mdhrv = 'v' in tors_model

        mdhr_name = automol.rotor.names(rotors)[multi_idx]
        mdhr_grid = automol.rotor.grids(rotors, increment=increment)[multi_idx]

        pot, geoms, grads, hessians, _, _ = filesys.read.potential(
            mdhr_name, mdhr_grid,
            cnf_save_path,
            mod_tors_ene_info, ref_ene,
            constraint_dct=None,   # No extra frozen treatments
            read_geom=is_mdhrv,
            read_grad=is_mdhrv,
            read_hess=is_mdhrv,
            read_energy_backstep=False,
            remove_bad_points=True)

        if is_mdhrv:
            script_str = autorun.SCRIPT_DCT['projrot']
            freqs = autorun.projrot.pot_frequencies(
                script_str, geoms, grads, hessians, run_path)
        else:
            freqs = None

        mdhr_dct = {'pot': pot, 'freqs': freqs}
    else:
        mdhr_dct = None

    return rotors, mdhr_dct


def scale_rotor_pots(rotors, scale_factor=((), None)):
    """ scale the pots
    """

    # Count numbers
    numtors = 0
    for rotor in rotors:
        numtors += len(rotor)

    # Calculate the scaling factors
    scale_indcs, factor = scale_factor
    nscale = numtors - len(scale_indcs)

    if nscale > 0:
        sfactor = factor**(2.0/nscale)
        ioprinter.debug_message(
            'scale_coeff test:', factor, nscale, sfactor)

        # test
        # sfactor = 1
        # test
        for tidx, rotor in enumerate(rotors):
            for torsion in rotor:
                if tidx not in scale_indcs and factor is not None:
                    torsion.pot = automol.pot.scale(torsion.pot, sfactor)
                    # following is being used in a test to see how effective
                    # a scaling of fixed scan torsional pots can be
                    # torsion.pot = automol.pot.relax_scale(torsion.pot)

    return rotors


# FUNCTIONS TO WRITE STRINGS FOR THE ROTORS FOR VARIOUS SITUATION
def make_hr_strings(rotors, mdhr_dct=None):
    """ Procedure for building the MESS strings
        :return mess_allrot_str: combination of intl and hr strs
        :return mess_hr_str: all 1dhr strs
    """

    # Initialize empty strings
    mess_allr_str = ''
    mess_hr_str, mess_flux_str, projrot_str = '', '', ''
    mdhr_dat = ''

    # Convert the rotor objects indexing to be in geoms
    geo, rotors = automol.rotor.relabel_for_geometry(rotors)

    # Get the number of rotors
    numrotors = len(rotors)
    for _, rotor in enumerate(rotors):
        multirotor = len(rotor) > 1

        for _, torsion in enumerate(rotor):

            # Write the rotor strings
            hr_str, ir_str, flux_str, prot_str = _tors_strs(torsion, geo)
            # mess_allr_str += hr_str
            mess_hr_str += hr_str
            mess_flux_str += flux_str
            projrot_str += prot_str

            # For MDHR, add the appropriate string
            if mdhr_dct is not None:
                if ((numrotors > 1 and multirotor) or numrotors == 1):
                    mess_allr_str += ir_str
                else:
                    mess_allr_str += hr_str
            else:
                mess_allr_str += hr_str

    # Write the mdhr dat string
    if mdhr_dct is not None:
        mdhr_dat = mess_io.writer.mdhr_data(
            mdhr_dct['pot'], freqs=mdhr_dct['freqs'], nrot=numrotors)

    return mess_allr_str, mess_hr_str, mess_flux_str, projrot_str, mdhr_dat


def _tors_strs(torsion, geo):
    """ Gather the 1DHR torsional data and gather them into a MESS file
    """

    mess_hr_str = mess_io.writer.rotor_hindered(
        group=torsion.groups[0],
        axis=torsion.axis,
        symmetry=torsion.symmetry,
        potential=torsion.pot,
        hmin=None,
        hmax=None,
        lvl_ene_max=None,
        therm_pow_max=None,
        geo=geo,
        rotor_id=torsion.name)

    mess_ir_str = mess_io.writer.rotor_internal(
        group=torsion.groups[0],
        axis=torsion.axis,
        symmetry=torsion.symmetry,
        grid_size=100,
        mass_exp_size=5,
        pot_exp_size=5,
        hmin=13,
        hmax=101,
        geo=None,
        rotor_id=torsion.name)

    mess_flux_str = mess_io.writer.fluxional_mode(
        torsion.indices,
        span=torsion.span)

    projrot_str = projrot_io.writer.rotors(
        axis=torsion.axis,
        group=torsion.groups[0])

    return mess_hr_str, mess_ir_str, mess_flux_str, projrot_str


def _need_tors_geo(pf_levels):
    """ Determine if a torsional geometry is geometry if needed
    """
    return bool(pf_levels['tors'][1] == pf_levels['harm'])
