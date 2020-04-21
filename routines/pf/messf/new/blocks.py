""" New blocks python
"""

import mess_io
from routines.pf.messf import _models as models
from routines.pf.messf import _fake as fake
from routines.pf.messf import _util as util


# SINGLE SPECIES BLOCKS
def species_block(spc_dct_i, rxn, spc_model, pf_levels,
                  save_prefix, saddle=False):
    """ prepare the species input for messpf
    """

    # Write the MESS string for the atom or molecule molecule
    if util.is_atom(spc_dct_i):
        # Get info and build the appropriate MESS string
        mass = util.atom_mass(spc_dct_i)
        elec_levels = spc_dct_i['elec_levels']
        spc_str = mess_io.writer.atom(mass, elec_levels)
    else:
        # Build dct combinining various information from the filesys and MESS
        inf_dct = models.read_filesys_for_spc(
            spc_dct_i, rxn, spc_model, pf_levels, save_prefix, saddle=saddle)

        # Build the appropriate core string
        tors_model, _, _ = spc_model
        if inf_dct['has_tors'] and tors_model == 'mdhr':
            core_str = inf_dct['mdhr_str']
        else:
            core_str = mess_io.writer.core_rigidrotor(
                geom=inf_dct['geom'],
                sym_factor=inf_dct['sym_factor'],
                interp_emax=None
            )

        # Build the appropriate MESS string
        spc_str = mess_io.writer.molecule(
            core=core_str,
            freqs=inf_dct['freqs'],
            elec_levels=inf_dct['elec_levels'],
            hind_rot=inf_dct['hind_rot_str'],
            xmat=inf_dct['xmat'],
            rovib_coups=inf_dct['rovib_coups'],
            rot_dists=inf_dct['rot_dists']
        )

    return spc_str


def fake_species_block(spc_dct_i, spc_dct_j, spc_model, pf_levels,
                       save_prefix, rxn=()):
    """ prepare a fake species block corresponding to the
        van der Waals well between two fragments
    """

    # Build a dct combinining various information from the filesys and MESS
    inf_dct_i = models.read_filesys_for_spc(
        spc_dct_i, rxn, spc_model, pf_levels, save_prefix, saddle=False)
    inf_dct_j = models.read_filesys_for_spc(
        spc_dct_j, rxn, spc_model, pf_levels, save_prefix, saddle=False)

    # Combine the electronic structure information for the two species together
    geom = fake.combine_geos_in_fake_well(inf_dct_i['geom'], inf_dct_j['geom'])

    sym_factor = inf_dct_i['sym_factor'] * inf_dct_j['sym_factor']

    elec_levels = util.combine_elec_levels(
        inf_dct_i['elec_levels'], inf_dct_j['elec_levels'])

    fake_freqs = fake.set_fake_freqs(inf_dct_i['geom'], inf_dct_j['geom'])

    freqs = fake_freqs + inf_dct_i['freqs'] + inf_dct_j['freqs']

    hind_rot_str = inf_dct_i['hind_rot_str'] + inf_dct_j['hind_rot_str']

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
        hind_rot=hind_rot_str,
        xmat=(),
        rovib_coups=(),
        rot_dists=()
    )

    return spc_str


def pst_block(spc_dct_i, spc_dct_j, spc_model, pf_levels,
              save_prefix, rxn=(), pst_params=(1.0, 6)):
    """ prepare a Phase Space Theory species block
    """

    # Build a dct combinining various information from the filesys and MESS
    inf_dct_i = models.read_filesys_for_spc(
        spc_dct_i, rxn, spc_model, pf_levels, save_prefix, saddle=False)
    inf_dct_j = models.read_filesys_for_spc(
        spc_dct_j, rxn, spc_model, pf_levels, save_prefix, saddle=False)

    # Combine the electronic structure information for the two species together
    sym_factor = inf_dct_i['sym_factor'] * inf_dct_j['sym_factor']

    elec_levels = util.combine_elec_levels(
        inf_dct_i['elec_levels'], inf_dct_j['elec_levels'])

    freqs = inf_dct_i['freqs'] + inf_dct_j['freqs']

    hind_rot_str = inf_dct_i['hind_rot_str'] + inf_dct_j['hind_rot_str']

    # Get the total stoichiometry of the two species
    stoich = util.get_stoich(inf_dct_i['geom'], inf_dct_j['geom'])

    # Write the MESS string for the Phase Space Theory TS
    core_str = mess_io.writer.core_phasespace(
        geom1=inf_dct_i['geom'],
        geom2=inf_dct_j['geom'],
        sym_factor=sym_factor,
        stoich=stoich,
        pot_prefactor=pst_params[0],
        pot_power_exp=pst_params[1]
    )
    spc_str = mess_io.writer.molecule(
        core=core_str,
        freqs=freqs,
        elec_levels=elec_levels,
        hind_rot=hind_rot_str,
        xmat=(),
        rovib_coups=(),
        rot_dists=()
    )

    return spc_str


def tau_block(spc_dct_i, spc_model, pf_levels, save_prefix,
              rxn=(), saddle=False):
    """ write  MESS string when using the Tau MonteCarlo
    """

    # Build a dct combinining various information from the filesys and MESS
    inf_dct = models.read_filesys_for_tau(
        spc_dct_i, spc_model, pf_levels,
        save_prefix, rxn=rxn, saddle=saddle)

    # Write the data string and set its name
    dat_str = mess_io.writer.monte_carlo.mc_data(
        geos=inf_dct['samp_geoms'],
        enes=inf_dct['samp_enes'],
        grads=inf_dct['samp_grads'],
        hessians=inf_dct['samp_hessians']
    )
    tau_dat_file_name = None

    # Write the core string (seperate energies?)
    spc_str = mess_io.writer.monte_carlo.mc_species(
        geom=inf_dct['ref_geo'],
        elec_levels=inf_dct['elec_levels'],
        flux_mode_str=inf_dct['flux_mode_str'],
        data_file_name=tau_dat_file_name,
        ground_energy=inf_dct['ground_ene'],
        reference_energy=inf_dct['reference_ene'],
        freqs=inf_dct['freqs'],
        no_qc_corr=True,
        use_cm_shift=True)

    return spc_str, dat_str


# TS BLOCKS FOR VARIATIONAL TREATMENTS
def vrctst_block(spc_dct_i, spc_dct_j, rxn, spc_model, pf_levels,
                 save_prefix):
    """ write a VRCTST block
    """

    # Build a dct combinining various information from the filesys and MESS
    inf_dct_i = models.read_filesys_for_spc(
        spc_dct_i, rxn, spc_model, pf_levels, save_prefix, saddle=False)
    inf_dct_j = models.read_filesys_for_spc(
        spc_dct_j, rxn, spc_model, pf_levels, save_prefix, saddle=False)

    # Get the flux file
    flux_file_name = '{}_flux.dat'.format('ts')  # fix
    flux_str = models.read_filesys_for_flux(rxn, spc_model)

    # Combine electronic structure information for the two species together
    sym_factor = inf_dct_i['sym_factor'] * inf_dct_j['sym_factor']

    elec_levels = util.combine_elec_levels(
        inf_dct_i['elec_levels'], inf_dct_j['elec_levels'])

    freqs = inf_dct_i['freqs'] + inf_dct_j['freqs']

    hind_rot_str = inf_dct_i['hind_rot_str'] + inf_dct_j['hind_rot_str']

    # Get the total stoichiometry of the two species
    stoich = util.get_stoich(inf_dct_i['geom'], inf_dct_j['geom'])

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
        hind_rot=hind_rot_str,
        xmat=(),
        rovib_coups=(),
        rot_dists=()
    )

    return spc_str, flux_str


def rpath_vtst_nosadpt_block(ts_dct, ts_label, reac_label, prod_label,
                             spc_ene, projrot_script_str, multi_info):
    """ prepare the mess input string for a variational TS that does not have
    a saddle point. Do it by calling the species block for each grid point
    in the scan file system
    """

    ts_info = ['', ts_dct['chg'], ts_dct['mul']]
    multi_level = fsorb.mod_orb_restrict(ts_info, multi_info)

    rxn_save_path = ts_dct['rxn_fs'][3]
    thy_save_fs = autofile.fs.theory(rxn_save_path)
    thy_save_fs[-1].create(multi_level[1:4])
    thy_save_path = thy_save_fs[-1].path(multi_level[1:4])
    scn_save_fs = autofile.fs.scan(thy_save_path)

    # Read the scan save filesystem to get the molecular info
    sym_factor = 1.
    irc_pt_strs = []
    proj_rotors_str = ''
    pot = []

    elec_levels = [[0., ts_dct['mul']]]
    grid = ts_dct['grid']
    grid = numpy.append(grid[0], grid[1])
    dist_name = ts_dct['dist_info'][0]

    # Read the infinite separation energy
    inf_locs = [[dist_name], [1000.]]
    inf_sep_ene = scn_save_fs[-1].file.energy.read(inf_locs)

    # Build dct of vtst info
    inf_dct = models.read_filesys_for_rpvtst()

    # Write the string
    variational_str = _write_mess_variational_str()

    return variational_str


def rpath_vtst_sadpt_block(ts_dct, ene_thy_level, geo_thy_level,
                           ts_label, reac_label, prod_label, first_ground_ene):
    """ prepare the mess input string for a variational TS where there is a
        saddle point on the MEP.
        In this case, there is limited torsional information.
    """
    [geo_thy_level, ene_thy_level, _, _, _, _] = pf_levels

    irc_idxs = ts_dct['irc_idxs']
    ts_info = ['', ts_dct['chg'], ts_dct['mul']]

    # Build the TS scan file system
    scn_save_fs, _, _, _ = irc.ts_scn_fs(
        ts_dct, ts_info, geo_thy_level)

    # Set the distance name for the reaction coordinate
    dist_name = 'RC'

    # Build dct of vtst info
    inf_dct_lst = models.read_filesys_for_rpvtst()

    # Get the saddle point

    # Write the string
    variational_str = _write_mess_variational_str(
        inf_dct_lst, irc_idxs, ts_idx=None)

    return variational_str


def _write_mess_variational_str(inf_dct_lst, rpath_idxs, ts_idx=None):
    """ write the variational string
    """

    rpath_pt_strs = []
    for idx, inf_dct in zip(rpath_idxs, inf_dct_lst):

        # Iniialize the header of the rxn path pt string
        rpath_pt_str = '!-----------------------------------------------\n'
        rpath_pt_str += '! RXN Path Point {0}'.format(str(int(idx+1)))
        if ts_idx is not None and ts_idx == idx:
            rpath_pt_str += '  (Saddle Point)'
        rpath_pt_str += '\n'

        # Write MESS string for the rxn path pt; add to rxn path pt string
        core_str = mess_io.writer.mol_data.core_rigidrotor(
            geom=inf_dct['geom'],
            sym_factor=inf_dct['sym_factor'],
            interp_emax=None
        )
        rpath_pt_str += mess_io.writer.species.molecule(
            core=core_str,
            freqs=inf_dct['freqs'],
            elec_levels=inf_dct['elec_levels'],
            hind_rot='',
            xmat=(),
            rovib_coups=(),
            rot_dists=()
        )

        # Add the ZPVE string for the rxn path point
        rpath_pt_str += (
            '  ZeroEnergy[kcal/mol]' +
            '      {0:<8.2f}'.format(inf_dct['zero_ene'])
        )

        # Append rxn path pt string to full list of rpath strings
        rpath_pt_strs.append(rpath_pt_str)

    return rpath_pt_strs


def vtst_energy():
    """  Get the VTST energy
    """
    if not saddle:
        # Calcuate infinite separation ZPVE
        # Assumes the ZPVE = ZPVE(1st grid pt) as an approximation
        if idx == 0:
            rct_zpe = zpe

        # Calculate the reference energies
        ene_rel = (ene - inf_sep_ene) * phycon.EH2KCAL
        zpe_rel = zpe - rct_zpe
        eref_abs = ene_rel + zpe_rel + spc_ene
    if saddle:
        # Calculate the relative energy
        erel = ene * phycon.EH2KCAL + zpe - first_ground_ene

    return erel
