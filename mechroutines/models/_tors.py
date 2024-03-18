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
from mechroutines.es.ts import ts_zma_locs
from mechroutines.models import typ


# FUNCTIONS TO BUILD ROTOR OBJECTS CONTAINING ALL NEEDED INFO
def build_rotors(spc_dct_i, pf_filesystems, spc_mod_dct_i,
                 read_potentials=True):
    """ Add more rotor info
    """

    run_prefix = pf_filesystems['run_prefix']
    spc_info = sinfo.from_dct(spc_dct_i, canonical=True)
    if 'rxn_info' not in spc_dct_i:
        spc_fml = automol.chi.formula_layer(spc_info[0])
    else:
        spc_fml = 'TS'
    run_path = job_path(run_prefix, 'PROJROT', 'FREQ', spc_fml, locs_id=None)

    # Set up tors level filesystem and model and level
    tors_model = spc_mod_dct_i['tors']['mod']
    tors_ene_info = spc_mod_dct_i['tors']['enelvl'][1][1]
    mod_tors_ene_info = tinfo.modify_orb_label(
        tors_ene_info, spc_info)

    rotors, mdhr_dct = None, None
    zma_locs = (0,)
    if pf_filesystems['tors'] is not None:
        [cnf_fs, cnf_save_path, min_cnf_locs, _, _] = pf_filesystems['tors']

        # Build the rotors
        ref_ene = filesys.read.energy(cnf_fs, min_cnf_locs, mod_tors_ene_info)
        zma_fs = fs.zmatrix(cnf_fs[-1].path(min_cnf_locs))
        zma_locs = (0,)
        if 'zrxn' in spc_dct_i:
            zma_locs = ts_zma_locs(None, None, zma_fs, spc_dct_i)
        if (
            zma_fs[-1].file.torsions.exists(zma_locs) and
            zma_fs[-1].file.zmatrix.exists(zma_locs) and
            tors_model != 'rigid'
        ):
            rotors = automol.data.rotor.rotors_from_data(
                zma=zma_fs[-1].file.zmatrix.read(zma_locs),
                tor_lst=zma_fs[-1].file.torsions.read(zma_locs),
                tor_names_lst=(
                    spc_dct_i.get('tors_names', None)
                    if 'md' in tors_model else None),
                multi='md' in tors_model)
        # Read the potential grids
        if read_potentials and rotors is not None:
            rotors, mdhr_dct = _read_potentials(
                rotors, spc_dct_i, run_path, cnf_save_path,
                ref_ene, mod_tors_ene_info,
                tors_model)

    # Squash the rotor potentials as necessary
    if rotors is not None:
        if typ.squash_tors_pot(spc_mod_dct_i):
            for rotor in rotors:
                for torsion in rotor:
                    torsion.pot = automol.pot.relax_scale(torsion.pot)
    return rotors, mdhr_dct, zma_locs


def _read_potentials(rotors, spc_dct_i, run_path, cnf_save_path,
                     ref_ene, mod_tors_ene_info,
                     tors_model):
    """ read out the potentials
    """

    _ = run_path

    # Convert the rotor objects indexing to be in geoms
    increment = spc_dct_i.get('hind_inc', 30.0*phycon.DEG2RAD)
    rotor_zma = automol.data.rotor.rotors_zmatrix(rotors)

    # Determine base-line rotor non-specific info for constraints
    all_tors_names = automol.data.rotor.rotors_torsion_names(rotors)
    print('rotor names', all_tors_names)
    const_names = automol.zmat.set_constraint_names(
        rotor_zma, all_tors_names, tors_model)

    # Recalculate the rotor potential grids using desired increment
    rotor_grids = automol.data.rotor.rotors_torsion_grids(rotors, increment=increment)

    multi_idx = None
    for ridx, rotor in enumerate(rotors):
        tor_lst = automol.data.rotor.torsions(rotor)
        tor_grids = automol.data.rotor.torsion_grids(rotor)

        if len(tor_lst) > 1:
            multi_idx = ridx

        for tidx, (tor, tor_grid) in enumerate(zip(tor_lst, tor_grids)):
            # Read and spline-fit potential
            tor_name = automol.data.tors.name(tor)
            constraint_dct = automol.zmat.constraint_dict(
                rotor_zma, const_names, (tor_name,))
            pot, scan_geos, _, _, _, _ = filesys.read.potential(
                (tor_name,), (tor_grid,),
                cnf_save_path,
                mod_tors_ene_info, ref_ene,
                constraint_dct,
                read_energy_backstep=True,
                remove_bad_points=True,
                read_geom=True)

            if pot:
                pot_obj = automol.data.potent.from_dict(
                    pot, aux_dct_dct={"geom": scan_geos}
                )
                pot_obj = automol.data.potent.clean(pot_obj)
                pot_obj = automol.data.potent.zero_coordinates_values(pot_obj)
                automol.data.rotor.set_potential(rotor, pot_obj, in_place=True)

    if multi_idx is not None:
        is_mdhrv = 'v' in tors_model

        mdhr_name = automol.data.rotor.rotors_torsion_names(rotors)[multi_idx]
        mdhr_grid = automol.data.rotor.rotors_torsion_grids(rotors, increment=increment)[multi_idx]

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


def scale_rotor_pots(rotors, scale_factor=((), None), scale_override=None):
    """ scale the pots
    """

    # Count numbers
    numtors = 0
    for rotor in rotors:
        numtors += len(rotor.torsions)

    # Calculate the scaling factors
    scale_indcs, factor = scale_factor
    nscale = numtors - len(scale_indcs)

    
    if nscale > 0:
        if scale_override is not None:
            sfactor = factor**(2.0/nscale)
            ioprinter.debug_message(
                'scale_coeff override test:', factor, nscale, sfactor)
        else:
            sfactor = factor**(2.0/nscale)
            sfacmax = 1.3
            sfacmin = 0.8
            if nscale > 4:
                sfacmax = 0.1*(nscale-4) + 1.3
                # sfacmin = -0.04*(nscale-4) + 1.22
            ioprinter.debug_message(
                'scale_coeff test:', factor, nscale, sfactor, sfacmax, sfacmin)
            if sfactor > sfacmax:
                print ('value of sfactor is greater than sfacmax')
                sfactor = sfacmax
                #elif sfactor < sfacmin:
            elif sfactor < sfacmin:
                sfactor = sfacmin
                print ('value of sfactor is less than sfacmin')

        # test
        # sfactor = 1
        # test
        for tidx, rotor in enumerate(rotors):
            if tidx not in scale_indcs and factor is not None:
                pot = automol.data.rotor.potential(rotor)
                pot = automol.data.potent.scale(pot, sfactor)
                automol.data.rotor.set_potential(rotor, pot, in_place=True)

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
    zma = automol.data.rotor.rotors_zmatrix(rotors)
    geo = automol.zmat.geometry(zma)

    # Get the number of rotors
    numrotors = len(rotors)
    for rotor in rotors:
        multirotor = automol.data.rotor.dimension(rotor)

        for torsion in automol.data.rotor.torsions(rotor, key_typ="geom"):

            # Write the rotor strings
            hr_str, ir_str, flux_str, prot_str = _tors_strs(torsion, rotor, geo)
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


def _tors_strs(torsion, rotor, geo):
    """ Gather the 1DHR torsional data and gather them into a MESS file
    """
    pot = automol.data.rotor.potential(rotor)
    pot_dct = automol.data.potent.dict_(pot)

    mess_hr_str = mess_io.writer.rotor_hindered(
        group=automol.data.tors.groups(torsion)[0],
        axis=automol.data.tors.axis(torsion),
        symmetry=automol.data.tors.symmetry(torsion),
        potential=pot_dct,
        hmin=None,
        hmax=None,
        lvl_ene_max=None,
        therm_pow_max=None,
        geo=geo,
        rotor_id=automol.data.tors.name(torsion))

    mess_ir_str = mess_io.writer.rotor_internal(
        group=automol.data.tors.groups(torsion)[0],
        axis=automol.data.tors.axis(torsion),
        symmetry=automol.data.tors.symmetry(torsion),
        grid_size=50,
        mass_exp_size=5,
        pot_exp_size=11,
        hmin=13,
        hmax=101,
        geo=None,
        rotor_id=automol.data.tors.name(torsion))

    mess_flux_str = mess_io.writer.fluxional_mode(
        automol.data.tors.coordinate(torsion),
        span=automol.data.tors.span(torsion))

    projrot_str = projrot_io.writer.rotors(
        axis=automol.data.tors.axis(torsion),
        group=automol.data.tors.groups(torsion)[0])

    return mess_hr_str, mess_ir_str, mess_flux_str, projrot_str


def _need_tors_geo(pf_levels):
    """ Determine if a torsional geometry is geometry if needed
    """
    return bool(pf_levels['tors'][1] == pf_levels['harm'])
