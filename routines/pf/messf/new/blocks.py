""" New blocks python
"""

def species_block():
    """ """

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
    hind_rot_str, proj_rot_str = '', ''

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
            hind_rot_str, proj_rotors_str, mdhr_dat_str = tors.proc_mess_strings()
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
            if vib_tau:
                tau_freqs = ()
                tau_gradient = True
                tau_hessian = True
            monte_carlo_str, tau_dat_str = tau.proc_tau(
                tors_min_cnf_locs, tors_cnf_save_fs,
                spc_dct_i,
                frm_bnd_key, brk_bnd_key,
                sym_factor, elec_levels,
                save_prefix,
                tau_dat_file_name,
                freqs=tau_freqs,
                saddle=False
                gradient=tau_gradient,
                hessian=tau_hessian)
        else:
            monte_carlo_str, tau_dat_str = '', ''

        # Write the species string for the molecule
        if tors_model == 'tau':
            spc_str = mc_str
        else:
            if tors_model == 'mdhr':
                core = mdhr_str
            else:
                core = mess_io.writer.core_rigidrotor(geo, symf)
            spc_str = mess_io.writer.molecule(
                core, freqs, elec_levels,
                hind_rot=hind_rot_str, xmat=xmat,
                rovib_coups=rovib_coups, rot_dists=rot_dists)

        # Combine various data strings into a dct
        dat_str_dct = {
            'tau': (tau_dat_str, tau_dat_file_name),
            'mdhr': (mdhr_dat_str, mdhr_dat_file_name),
        }

    return spc_str, dat_str_dct, imag


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
