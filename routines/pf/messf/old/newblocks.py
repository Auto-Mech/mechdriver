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
        # Geo

        # Torsions 
        if bld_tors():

        # Freqs

        # Symmetry

        # Set Tau if needed

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
                hind_rot=hr_str, xmat=xmat)

        # Combine various data strings into a dct
        dat_str_dct = {
            'tau': (tau_dat_str, tau_dat_file_name),
            'mdhr': (mdhr_dat_str, mdhr_dat_file_name),
        }

    return spc_str, dat_str_dct, imag


def bld_tors():
    """ decide to build torsions
    """
    has_tors = bool(tors_names)
    tors_hr_model = bool(tors in ('1dhr', 'mdhr')),
    tau_hr_model = bool(tors_model == 'tau' and vib_model != 'vib')
    return has_tors and (tors_hr_model or tau_hr_model)


