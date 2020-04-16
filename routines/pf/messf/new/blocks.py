""" New blocks python
"""


# SINGLE SPECIES BLOCKS
def species_block(spc, spc_dct_i, spc_info, spc_model,
                  pf_levels, save_prefix):
    """ prepare the species input for messpf
    """

    # Build a dct combinining various information from the filesys and MESS
    inf_dct = read_species_filesys()

    # Write the species string for the molecule
    if tors_model == 'tau':
        spc_str = mc_str
    else:
        if tors_model == 'mdhr':
            core = mdhr_str
        else:
            core = mess_io.writer.core_rigidrotor(geo, symf)
        spc_str = mess_io.writer.molecule(
            core,
            freqs,
            elec_levels,
            hind_rot=hind_rot_str,
            xmat=xmat,
            rovib_coups=rovib_coups,
            rot_dists=rot_dists)

    return spc_str


def fake_species_block(
        spc_dct_i, spc_dct_j, spc_info_i, spc_info_j,
        spc_model, pf_levels,
        save_prefix_i, save_prefix_j):
    """ prepare a fake species block corresponding to the
        van der Waals well between two fragments
    """

    # Build a dct combinining various information from the filesys and MESS
    inf_dct_i = read_species_filesys()
    inf_dct_j = read_species_filesys()

    # Obtain combined values from both species
    elec_levels = messfutil.combine_elec_levels(spc_dct_i, spc_dct_j)
    freqs = fake.set_fake_freqs(
        harm_min_cnf_locs_i, harm_min_cnf_locs_j,
        harm_cnf_save_fs_i, harm_cnf_save_fs_j)
    freqs = freqs + freqs_i + freqs_j
    geo = fake.combine_geos_in_fake_well(
        harm_min_cnf_locs_i, harm_min_cnf_locs_j,
        harm_cnf_save_fs_i, harm_cnf_save_fs_j)
    sym_factor = sym_factor_i * sym_factor_j
    hind_rot_str = hr_str_i + hr_str_j

    # Write the MESS string
    core = mess_io.writer.core_rigidrotor(geo, sym_factor)
    spc_str = mess_io.writer.molecule(
        core, freqs, elec_levels,
        hind_rot=hind_rot_str)

    return spc_str


def pst_block(spc_dct_i, spc_dct_j, spc_model, pf_levels,
              spc_save_fs, pst_params=(1.0, 6)):
    """ prepare a Phase Space Theory species block
    """

    # Build a dct combinining various information from the filesys and MESS
    inf_dct_i = read_species_filesys()
    inf_dct_j = read_species_filesys()

    # Obtain combined values from both species
    elec_levels = messfutil.combine_elec_levels(spc_dct_i, spc_dct_j)
    sym_factor = sym_factor_i * sym_factor_j

    #

    # Get the stoichiometry
    stoich = messfutil.get_stoich(
        harm_min_cnf_locs_i, harm_min_cnf_locs_j,
        harm_cnf_save_fs_i, harm_cnf_save_fs_j)

    # Write the MESS input strings
    core = mess_io.writer.core_phasespace(
        geo_i, geo_j, sym_factor, stoich,
        pot_prefactor=pst_params[0], pot_power_exp=pst_params[1])
    spc_str = mess_io.writer.molecule(
        core, freqs, elec_levels,
        hind_rot=hind_rot_str)

    return spc_str


# TS BLOCKS FOR VTST TREATMENTS
def vtst_with_no_saddle_block(
        ts_dct, ts_label, reac_label, prod_label,
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

    for x:
        # Iniialize the header of the string
        irc_pt_str = '!-----------------------------------------------\n'
        irc_pt_str += '! IRC Point {0}\n'.format(str(int(idx)))

        # Write the molecule section for each irc point
        core = mess_io.writer.mol_data.core_rigidrotor(
            geom, sym_factor, interp_emax='')
        irc_pt_str += mess_io.writer.species.molecule(
            core, freqs, elec_levels,
            hind_rot='', xmat=None, rovib_coups='', rot_dists='')

        # Append the zero point energy string for the point
        irc_pt_str += ('  ZeroEnergy[kcal/mol]      '
                       '{0:<8.2f}'.format(erel))

        # Append string to list
        full_irc_str.append(irc_pt_str)

    # Write the MESS string for the variational sections
    variational_str = mess_io.writer.rxnchan.ts_variational(
        ts_label, reac_label, prod_label, full_irc_str)

    return vtst_str


def vtst_saddle_block(ts_dct, ene_thy_level, geo_thy_level,
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

    return


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


# FILE READERS
def read_file_sys_info():
    """ Pull all of the neccessary information from the filesystem for a species
    """

    # Unpack the models and levels
    [_, _, harm_level, vpt2_level, sym_level, tors_level] = pf_levels
    tors_model, vib_model, sym_model = spc_model

    # Lazy set for tors mod
    if tors_model == '1dhrf':
        tors_mod = ('1dhr', False)

    # Set theory filesystem used throughout
    thy_save_fs = autofile.fs.theory(save_prefix)

    # Set boolean to account for rad-rad reaction (not supported by vtst)
    rad_rad_ts = False
    if 'ts_' in spc:
        if spc_dct_i['rad_rad']:
            rad_rad_ts = True

    # Set the filesystem objects for various species models
    harmfs = set_model_filesys(
        thy_save_fs, spc_info, harm_level, saddle=('ts_' in spc))
    [harm_cnf_save_fs, _,
     harm_min_cnf_locs, harm_save_path] = harmfs
    if sym_level:
        symfs = set_model_filesys(
            thy_save_fs, spc_info, sym_level, saddle=('ts_' in spc))
        [sym_cnf_save_fs, _,
         sym_min_cnf_locs, _] = symfs
    if tors_level and not rad_rad_ts:
        torsfs = set_model_filesys(
            thy_save_fs, spc_info, tors_level[0], saddle=('ts_' in spc))
        [tors_cnf_save_fs, tors_cnf_save_path,
         tors_min_cnf_locs, tors_save_path] = torsfs
    if vpt2_level:
        vpt2fs = set_model_filesys(
            thy_save_fs, spc_info, vpt2_level, saddle=('ts_' in spc))
        [vpt2_cnf_save_fs, vpt2_cnf_save_path,
         vpt2_min_cnf_locs, vpt2_save_path] = vpt2fs

    # Check if any torsions to set the model

    # Build the species string and get the imaginary frequency
    # Initialize various auxiliary MESS data strings and names
    dat_str_dct = {}
    tau_dat_str = ''
    mdhr_dat_str = ''
    sct_dat_str = ''
    tau_dat_file_name = '{}_tau.dat'.format(spc[0])
    mdhr_dat_file_name = '{}_mdhr.dat'.format(spc[0])
    sct_dat_file_name = '{}_sct.dat'.format(spc[0])
    # Initialize various variables and strings
    symf = sym_factor
    imag = 0.0
    freqs = []
    hr_str = ""
    xmat = ()

    # Default strings
    hind_ot_str, proj_rot_str = '', ''

    # Pull all information from the filesys and write the MESS species block
    if messfutil.is_atom(harm_min_cnf_locs, harm_cnf_save_fs):
        mass = messfutil.atom_mass(harm_min_cnf_locs, harm_cnf_save_fs)
        spc_str = mess_io.writer.atom(
            mass, elec_levels)
    else:
        # Geo (Rotational)
        if nonrigid_rotations(rot_model):
            rovib_coup, rot_dists = rot.read_rotational_values(
                vpt2_save_fs, vpt2_min_cnf_locs)
        else:
            rovib_coup, rot_dists = (), ()

        # Torsions
        if nonrigid_tors(vib_model, tors_model, tors_names):
            mess_hr_str, proj_hr_str, mdhr_dat_str = tors.proc_mess_strings()
        else:
            hind_rot_str, proj_rotors_str, mdhr_dat_str = '', '', ''

        # Symmetry
        sym_factor = sym.symmetry_factor(
            sym_model, spc_dct_i, spc_info, dist_names,
            saddle, frm_bnd_key, brk_bnd_key, tors_names,
            tors_cnf_save_fs, tors_min_cnf_locs,
            sym_cnf_save_fs, sym_min_cnf_locs)
        if nonrigid_tors(vib_model, tors_model, tors_names):
            sym_nums = tors.get_tors_sym_nums(
                spc_dct_i, tors_min_cnf_locs, tors_cnf_save_fs,
                frm_bnd_key, brk_bnd_key, saddle=False)
            for num in sym_nums:
                sym_factor /= num

        # Vibrations
        harm_geo, freqs, imag_freq = harm_freqs()
        freqs, imag_freq, zpe = tors.tors_freqs_zpve()  # Could maybe split up
        if anharm_vib:
            xmat = vib.anharmonicity()
        else:
            xmat = ()

        # Set Tau if needed
        if need_tau:
            # Determine if hessians and grasdients are needed
            tau_freqs = () if vib_tau else freqs
            tau_gradient = vib_tau
            tau_hessian = vib_tau
            monte_carlo_str, tau_dat_str = tau.proc_tau(
                tors_min_cnf_locs, tors_cnf_save_fs,
                spc_dct_i,
                frm_bnd_key, brk_bnd_key,
                sym_factor, elec_levels,
                save_prefix,
                tau_dat_file_name,
                freqs=tau_freqs,
                saddle=False,
                gradient=tau_gradient,
                hessian=tau_hessian)
        else:
            monte_carlo_str, tau_dat_str = '', ''

        # Combine various data strings into a dct
        dat_str_dct = {
            'tau': (tau_dat_str, tau_dat_file_name),
            'mdhr': (mdhr_dat_str, mdhr_dat_file_name),
        }

    return spc_str, dat_str_dct, imag


# VTST
def vtst_read_file_sys_info():
    """ Pull all of the neccessary information from the filesystem for a species
    """
    for idx in irc_idxs:

        # Set the filesystem locators for each grid point
        locs = [[dist_name], [idx]]
        print(scn_save_fs[-1].path(locs))

        # Get geometry, energy, vibrational freqs, and zpe
        if scn_save_fs[-1].file.geometry.exists(locs):
            geom = scn_save_fs[-1].file.geometry.read(locs)
        else:
            print('no geom')
            continue
        if scn_save_fs[-1].file.energy.exists(locs):
            if ene_thy_level == geo_thy_level:
                ene = scn_save_fs[-1].file.energy.read(locs)
                print('ene', ene)
            else:
                scn_save_path = scn_save_fs[-1].path(locs)
                sp_save_fs = autofile.fs.single_point(scn_save_path)
                sp_level = fsorb.mod_orb_restrict(ts_info, ene_thy_level)
                if sp_save_fs[-1].file.energy.exists(sp_level[1:4]):
                    ene = sp_save_fs[-1].file.energy.read(sp_level[1:4])
                    print('ene-high', ene)
                else:
                    print('no energy')
                    continue
        else:
            print('no energy')
            continue
        if scn_save_fs[-1].file.hessian.exists(locs):
            proj_rotors_str = ''
            hess = scn_save_fs[-1].file.hessian.read(locs)
            scn_save_path = scn_save_fs[-1].path(locs)
            freqs, _, _ = vib.projrot_freqs_1(
                geom, hess,
                proj_rotors_str,
                scn_save_path, pot=False, saddle=True)
            zpe = sum(freqs)*phycon.WAVEN2KCAL/2.
        else:
            print('no hessian')
            continue


# Series of checks to determine what information is needed to be obtained
def nonrigid_rotations(rot_model):
    """ dtermine if a nonrigid rotation model is specified and further
        information is needed from the filesystem
    """
    return bool(rot_model in ('vpt2'))


def nonrigid_tors(vib_model, tors_model, tors_names):
    """ dtermine if a nonrigid torsional model is specified and further
        information is needed from the filesystem
    """
    has_tors = bool(tors_names)
    tors_hr_model = bool(tors_model in ('1dhr', 'mdhr', 'mdhrv'))
    tau_hr_model = bool(tors_model == 'tau' and vib_model != 'vib')
    return has_tors and (tors_hr_model or tau_hr_model)


def anharm_vib(vib_model):
    """ a
    """
    return bool(vib_model == 'vpt2')


def tau_pf(tors_model):
    """ determine if pf is done with tau
    """
    return bool(tors_model == 'tau')


def vib_tau(vib_model):
    """ determine if vibrations are treated via tau sampling
    """
    return bool(vib_model == 'tau')
