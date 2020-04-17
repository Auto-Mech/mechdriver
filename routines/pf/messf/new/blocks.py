""" New blocks python
"""


# SINGLE SPECIES BLOCKS
def species_block(spc, spc_dct_i, spc_info, spc_model,
                  pf_levels, save_prefix):
    """ prepare the species input for messpf
    """

    # Build a dct combinining various information from the filesys and MESS
    inf_dct = read_species_filesys(spc, spc_dct_i, spc_info, spc_model,
                                   pf_levels, save_prefix)

    # Write the MESS string for the molecule
    if tors_model == 'tau':
        spc_str = inf_dct['mc_str']
    else:
        if tors_model == 'mdhr':
            core_str = inf_dct['mdhr_str']
        else:
            core_str = mess_io.writer.core_rigidrotor(
                geom=inf_dct['geom'],
                sym_factor=inf_dct['sym_factor'],
                interp_emax=None
            )
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


def fake_species_block(
        spc_dct_i, spc_dct_j, spc_info_i, spc_info_j,
        spc_model, pf_levels,
        save_prefix_i, save_prefix_j):
    """ prepare a fake species block corresponding to the
        van der Waals well between two fragments
    """

    # Build a dct combinining various information from the filesys and MESS
    inf_dct_i = read_species_filesys(spc, spc_dct_i, spc_info_i, spc_model,
                                     pf_levels, save_prefix)
    inf_dct_j = read_species_filesys(spc, spc_dct_j, spc_info_j, spc_model,
                                     pf_levels, save_prefix)

    # Combine the electronic structure information for the two species together
    geom = fake.combine_geos_in_fake_well(
        harm_min_cnf_locs_i, harm_min_cnf_locs_j,
        harm_cnf_save_fs_i, harm_cnf_save_fs_j)

    sym_factor = inf_dct_i['sym_factor'] * inf_dct_j['sym_factor']

    elec_levels = messfutil.combine_elec_levels(spc_dct_i, spc_dct_j)

    fake_freqs = fake.set_fake_freqs(
        harm_min_cnf_locs_i, harm_min_cnf_locs_j,
        harm_cnf_save_fs_i, harm_cnf_save_fs_j)
    freqs = fake_freqs + inf_dct_i['freqs'] * inf_dct_j['freqs']

    hind_rot_str = inf_dct_i['hind_rot_str'] * inf_dct_j['hind_rot_str']

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


def tau_block(spc_dct_i, spc_dct_j, spc_model, pf_levels,
    """ write  MESS string when using the Tau MonteCarlo
    """

    return tau_str


def pst_block(spc_dct_i, spc_dct_j, spc_model, pf_levels,
              spc_save_fs, pst_params=(1.0, 6)):
    """ prepare a Phase Space Theory species block
    """

    # Build a dct combinining various information from the filesys and MESS
    inf_dct_i = read_species_filesys(spc, spc_dct_i, spc_info_i, spc_model,
                                     pf_levels, save_prefix)
    inf_dct_j = read_species_filesys(spc, spc_dct_j, spc_info_j, spc_model,
                                     pf_levels, save_prefix)

    # Combine the electronic structure information for the two species together
    sym_factor = inf_dct_i['sym_factor'] * inf_dct_j['sym_factor']

    elec_levels = messfutil.combine_elec_levels(spc_dct_i, spc_dct_j)

    freqs = inf_dct_i['freqs'] * inf_dct_j['freqs']

    hind_rot_str = inf_dct_i['hind_rot_str'] * inf_dct_j['hind_rot_str']

    # Get the total stoichiometry of the two species
    stoich = messfutil.get_stoich(
        harm_min_cnf_locs_i, harm_min_cnf_locs_j,
        harm_cnf_save_fs_i, harm_cnf_save_fs_j)

    # Write the MESS string for the Phase Space Theory TS
    core_str = mess_io.writer.core_phasespace(
        geom1=geo_i,
        geom2=geo_j,
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


# TS BLOCKS FOR VARIATIONAL TREATMENTS
def vrctst_block():
    """ write a VRCTST block
    """

    # Build a dct combinining various information from the filesys and MESS
    # inf_dct_i = read_species_filesys(spc, spc_dct_i, spc_info_i, spc_model,
    #                                  pf_levels, save_prefix)
    # inf_dct_j = read_species_filesys(spc, spc_dct_j, spc_info_j, spc_model,
    #                                  pf_levels, save_prefix)

    # # Combine the electronic structure information for the two species together
    # sym_factor = inf_dct_i['sym_factor'] * inf_dct_j['sym_factor']

    # elec_levels = messfutil.combine_elec_levels(spc_dct_i, spc_dct_j)

    # freqs = inf_dct_i['freqs'] * inf_dct_j['freqs']

    # hind_rot_str = inf_dct_i['hind_rot_str'] * inf_dct_j['hind_rot_str']

    # Get the total stoichiometry of the two species
    stoich = messfutil.get_stoich(
        harm_min_cnf_locs_i, harm_min_cnf_locs_j,
        harm_cnf_save_fs_i, harm_cnf_save_fs_j)

    # Write the MESS string for the VRCTST TS
    core_str = mess_io.writer.core_rotd(
        sym_factor=sym_factor,
        flux_file=flux_file,
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

    # Write the string
    variational_str = _write_mess_variational_str()

    return variational_str


def rpath_vtst_wsadpt_block(ts_dct, ene_thy_level, geo_thy_level,
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

    # Write the string
    variational_str = _write_mess_variational_str()

    return variational_str


def _write_mess_variational_str():
    """ write the variational string
    """

    for x:
        # Iniialize the header of the rxn path pt string
        rpath_pt_str = '!-----------------------------------------------\n'
        rpath_pt_str += '! RXN Path Point {0}\n'.format(str(int(idx)))

        # Write MESS string for the rxn path pt; add to rxn path pt string
        core_str = mess_io.writer.mol_data.core_rigidrotor(
            geom=geom,
            sym_factor=sym_factor,
            interp_emax=None
        )
        irc_pt_str += mess_io.writer.species.molecule(
            core=core_str,
            freqs=freqs,
            elec_levels=elec_levels,
            hind_rot='',
            xmat=(),
            rovib_coups=(),
            rot_dists=()
        )

        # Add the ZPVE string for the rxn path point 
        rpat_pt_str += (
            '  ZeroEnergy[kcal/mol]      {0:<8.2f}'.format(erel)
        )

        # Append rxn path pt string to full list of rpath strings
        rpath_pt_strs.append(rpath_pt_str)

    # Write the MESS string for the variational section
    variational_str = mess_io.writer.rxnchan.ts_variational(
        ts_label=ts_label,
        reac_label=reac_label,
        prod_label=prod_label,
        rpath_pt_strs=rpath_pt_strs,
        tunnel=''
    )

    return variational_str


def vtst_energy():
    """  Get the VTST energy
    """
    if no saddle:
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
