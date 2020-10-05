""" Writes the MESS strings for various species and treatments
    using data read from the save filesystem.
"""

import mess_io
from routines.pf.models import _fake as fake
from routines.pf.models import _util as util
from lib.phydat import phycon


# SINGLE SPECIES BLOCKS
def atom_block(inf_dct):
    """ prepare the species input for messpf
    """

    # Initialize an empty dat_dct for the return
    dat_dct = {}

    # Build the appropriate MESS string
    spc_str = mess_io.writer.atom(
        mass=inf_dct['mass'],
        elec_levels=inf_dct['elec_levels'])

    return spc_str, dat_dct


def species_block(inf_dct):
    """ prepare the species input for messpf
    """

    # Build the data files dct
    dat_dct = {}

    # Build the appropriate core string
    if inf_dct['mdhr_dat']:
        core_str = mess_io.writer.core_multirotor(
            geom=inf_dct['geom'],
            sym_factor=inf_dct['sym_factor'],
            pot_surf_file='mdhr_pot.dat',
            int_rot_str=inf_dct['mess_hr_str'],
            interp_emax=100,
            quant_lvl_emax=9
        )
        hind_rot_str = ''
        dat_dct['mdhr_pot.dat'] = inf_dct['mdhr_dat']
    else:
        core_str = mess_io.writer.core_rigidrotor(
            geom=inf_dct['geom'],
            sym_factor=inf_dct['sym_factor'],
            interp_emax=None
        )
        hind_rot_str = inf_dct['mess_hr_str']

    # Build the appropriate MESS string
    spc_str = mess_io.writer.molecule(
        core=core_str,
        elec_levels=inf_dct['elec_levels'],
        freqs=inf_dct['freqs'],
        hind_rot=hind_rot_str,
        xmat=inf_dct['xmat'],
        rovib_coups=inf_dct['rovib_coups'],
        rot_dists=inf_dct['rot_dists'],
        inf_intens=(),
        freq_scale_factor=None,
        use_harmfreqs_key=False
    )

    return spc_str, dat_dct


def fake_species_block(inf_dct_i, inf_dct_j):
    """ prepare a fake species block corresponding to the
        van der Waals well between two fragments
    """

    # Combine the electronic structure information for the two species together
    geom = fake.combine_geos_in_fake_well(inf_dct_i['geom'], inf_dct_j['geom'])

    sym_factor = inf_dct_i['sym_factor'] * inf_dct_j['sym_factor']

    elec_levels = util.combine_elec_levels(
        inf_dct_i['elec_levels'], inf_dct_j['elec_levels'])

    fake_freqs = fake.set_fake_freqs(inf_dct_i['geom'], inf_dct_j['geom'])
    freqs = fake_freqs + inf_dct_i['freqs'] + inf_dct_j['freqs']

    mess_hr_str = inf_dct_i['mess_hr_str'] + inf_dct_j['mess_hr_str']

    # Write the MESS string for the fake molecule
    core_str = mess_io.writer.core_rigidrotor(
        geom=geom,
        sym_factor=sym_factor,
        interp_emax=None
    )
    spc_str = mess_io.writer.molecule(
        core=core_str,
        freqs=freqs,
        elec_levels=elec_levels,
        hind_rot=mess_hr_str,
        xmat=(),
        rovib_coups=(),
        rot_dists=(),
        inf_intens=(),
        freq_scale_factor=None,
        use_harmfreqs_key=False
    )

    # Fix? I don't think you can do multirotor for the phase space theory
    dat_dct = {}

    return spc_str, dat_dct


def pst_block(ts_inf_dct, inf_dct_i, inf_dct_j):
    """ prepare a Phase Space Theory species block
    """

    # Combine the electronic structure information for the two species together
    sym_factor = inf_dct_i['sym_factor'] * inf_dct_j['sym_factor']

    elec_levels = util.combine_elec_levels(
        inf_dct_i['elec_levels'], inf_dct_j['elec_levels'])

    freqs = inf_dct_i['freqs'] + inf_dct_j['freqs']

    mess_hr_str = inf_dct_i['mess_hr_str'] + inf_dct_j['mess_hr_str']

    # Get the total stoichiometry of the two species
    stoich = util.get_stoich(inf_dct_i['geom'], inf_dct_j['geom'])

    # Write the MESS string for the Phase Space Theory TS
    core_str = mess_io.writer.core_phasespace(
        geom1=inf_dct_i['geom'],
        geom2=inf_dct_j['geom'],
        sym_factor=sym_factor,
        stoich=stoich,
        pot_prefactor=ts_inf_dct['cn_pst'],
        pot_exp=ts_inf_dct['n_pst'],
        tstlvl='e'
    )
    spc_str = mess_io.writer.molecule(
        core=core_str,
        freqs=freqs,
        elec_levels=elec_levels,
        hind_rot=mess_hr_str,
        xmat=(),
        rovib_coups=(),
        rot_dists=(),
        inf_intens=(),
        freq_scale_factor=None,
        use_harmfreqs_key=False
    )

    # Need to fix
    dat_dct = {}

    return spc_str, dat_dct


def tau_block(inf_dct):
    """ write  MESS string when using the Tau MonteCarlo
    """

    # Write the data string and set its name
    dat_str = mess_io.writer.monte_carlo.mc_data(
        geos=inf_dct['samp_geoms'],
        enes=inf_dct['samp_enes'],
        grads=inf_dct['samp_grads'],
        hessians=inf_dct['samp_hessians']
    )

    # Set the name of the tau dat file and add to dct
    tau_dat_file_name = 'tau.dat'
    dat_dct = {tau_dat_file_name: dat_str}

    # Write additional reference configuration file if needed
    if inf_dct['ref_geom'] and inf_dct['ref_grad'] and inf_dct['ref_hessian']:
        ref_config_file_name = 'reftau.dat'
        ref_dat_str = mess_io.writer.monte_carlo.mc_data(
            geos=inf_dct['ref_geom'],
            enes=['0.00'],
            grads=inf_dct['ref_grad'],
            hessians=inf_dct['ref_hessian']
        )
        ref_dat_str = "\n".join(ref_dat_str.splitlines()[1:])
        dat_dct.update({ref_config_file_name: ref_dat_str})
    else:
        ref_config_file_name = ''

    # Write the core string (seperate energies?)
    spc_str = mess_io.writer.monte_carlo.mc_species(
        geom=inf_dct['geom'],
        sym_factor=inf_dct['sym_factor'],
        elec_levels=inf_dct['elec_levels'],
        flux_mode_str=inf_dct['flux_mode_str'],
        data_file_name=tau_dat_file_name,
        reference_energy=inf_dct['reference_energy'],
        ref_config_file_name=ref_config_file_name,
        ground_energy=0.0,
        freqs=inf_dct['freqs'],
        use_cm_shift=True
    )

    return spc_str, dat_dct


# TS BLOCKS FOR VARIATIONAL TREATMENTS
def vrctst_block(inf_dct_ts, inf_dct_i, inf_dct_j):
    """ write a VRCTST block
    """

    # Build the data files dct
    dat_dct = {}

    # Combine electronic structure information for the two species together
    sym_factor = inf_dct_i['sym_factor'] * inf_dct_j['sym_factor'] * 0.850

    elec_levels = util.combine_elec_levels(
        inf_dct_i['elec_levels'], inf_dct_j['elec_levels'])

    freqs = inf_dct_i['freqs'] + inf_dct_j['freqs']

    mess_hr_str = inf_dct_i['mess_hr_str'] + inf_dct_j['mess_hr_str']

    # Get the total stoichiometry of the two species
    stoich = util.get_stoich(inf_dct_i['geom'], inf_dct_j['geom'])

    # Set the auxiliary flux file information
    flux_file_name = '{}_flux.dat'.format('ts')
    dat_dct[flux_file_name] = inf_dct_ts['flux_str']

    # Write the MESS string for the VRCTST TS
    core_str = mess_io.writer.core_rotd(
        sym_factor=sym_factor,
        flux_file_name=flux_file_name,
        stoich=stoich
    )
    spc_str = mess_io.writer.molecule(
        core=core_str,
        freqs=freqs,
        elec_levels=elec_levels,
        hind_rot=mess_hr_str,
        xmat=(),
        rovib_coups=(),
        rot_dists=(),
        inf_intens=(),
        freq_scale_factor=None,
        use_harmfreqs_key=False
    )

    return spc_str, dat_dct


def rpvtst_block(ts_inf_dct, inf_dct_i, inf_dct_j):
    """ prepare the mess input string for a variational TS that does not have
    a saddle point. Do it by calling the species block for each grid point
    in the scan file system
    """
    
    # Combine electronic structure information for the two species together
    sym_factor = inf_dct_i['sym_factor'] * inf_dct_j['sym_factor']
    mess_hr_str = inf_dct_i['mess_hr_str'] + inf_dct_j['mess_hr_str']
    elec_levels = [[0.0, 1.0]]
    # elec_levels = util.combine_elec_levels(
    #     inf_dct_i['elec_levels'], inf_dct_j['elec_levels'])

    rpath_strs = []
    for idx, dct in enumerate(ts_inf_dct['rpath']):

        # Iniialize the header of the rxn path pt string
        rpath_str = '!-----------------------------------------------\n'
        rpath_str += '! Rxn Path Pt {0}: '.format(str(idx+1))
        rpath_str += 'R = {0:.2f} Ang'.format(dct['rval'] * phycon.BOHR2ANG)
        rpath_str += '\n\n'

        # Write MESS string for the rxn path pt; add to rxn path pt string
        core_str = mess_io.writer.core_rigidrotor(
            geom=dct['geom'],
            sym_factor=sym_factor,
            interp_emax=None
        )
        rpath_str += mess_io.writer.molecule(
            core=core_str,
            freqs=dct['freqs'],
            elec_levels=elec_levels,
            hind_rot=mess_hr_str,
            xmat=(),
            rovib_coups=(),
            rot_dists=()
        )

        # Append rxn path pt string to full list of rpath strings
        rpath_strs.append(rpath_str)

    return rpath_strs
#
#
# def vtst_energy():
#     """  Get the VTST energy
#     """
#     if not saddle:
#         # Calcuate infinite separation ZPVE
#         # Assumes the ZPVE = ZPVE(1st grid pt) as an approximation
#         if idx == 0:
#             rct_zpe = zpe
#
#         # Calculate the reference energies
#         ene_rel = (ene - inf_sep_ene) * phycon.EH2KCAL
#         zpe_rel = zpe - rct_zpe
#         eref_abs = ene_rel + zpe_rel + spc_ene
#     if saddle:
#         # Calculate the relative energy
#         erel = ene * phycon.EH2KCAL + zpe - first_ground_ene
#
#     return erel
