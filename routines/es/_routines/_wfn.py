""" libraries to build wavefunctions for Molpro multireference calculations
"""

import automol
import elstruct
from phydat import act_space


# BUILD THE MULTIREFERENCE WAVEFUNCTION TO PASS MOLPRO ELSTRUCT
def build_wfn(ref_zma, ts_info, ts_formula, high_mul,
              rct_ichs, rct_info,
              aspace, mod_var_scn_thy_info):
    """ Build wavefunction
    """

    if aspace is not None:

        num_act_orb, num_act_elc, num_states, guess_str = aspace
        guess_lines = guess_str.splitlines()
        casscf_options = cas_options(
            ts_info, ts_formula, num_act_elc, num_act_orb, num_states,
            add_two_closed=False)
        print('\nUsing wfn guess from file...')

    else:
        num_act_orb, num_act_elc, num_states = active_space(
            rct_ichs, rct_info)

        # Build the elstruct CASSCF options list used to build the wfn guess
        # (1) Build wfn with active space
        # (2) Build wfn with active space + 2 closed orbitals for stability
        cas_opt = []
        cas_opt.append(
            cas_options(
                ts_info, ts_formula, num_act_elc, num_act_orb, num_states,
                add_two_closed=False))
        cas_opt.append(
            cas_options(
                ts_info, ts_formula, num_act_elc, num_act_orb, num_states,
                add_two_closed=True))

        # Write string that has all the components for building the wfn guess
        guess_str = multiref_wavefunction_guess(
            high_mul, ref_zma, ts_info, mod_var_scn_thy_info, cas_opt)
        guess_lines = guess_str.splitlines()

        # Set casscf options
        casscf_options = cas_opt[0]

        print('\nGenerating wfn guess from internal options...')

    # Manipulate the opt kwargs to use the wavefunction guess
    cas_kwargs = {
        'casscf_options': casscf_options,
        'gen_lines': {1: guess_lines},
        'mol_options': ['nosym']  # Turn off symmetry
    }

    return cas_kwargs


# BUILD THE CASSCF OPTIONS
def active_space(rct_ichs, rct_info):
    """ Determine the active space for the multireference MEP scan
    """

    num_act_orb, num_act_elc, num_states = 0, 0, 1
    for idx, ich in enumerate(rct_ichs):
        if ich in act_space.DCT:
            num_act_orb += act_space.DCT[ich][0]
            num_act_elc += act_space.DCT[ich][1]
            num_states *= act_space.DCT[ich][2]
        else:
            rct_mult = rct_info[idx][2]
            num_act_orb += (rct_mult - 1)
            num_act_elc += (rct_mult - 1)
            num_states *= 1

    return num_act_orb, num_act_elc, num_states


def cas_options(spc_info, formula, num_act_elc, num_act_orb, num_states,
                add_two_closed=False):
    """ Prepare values prepare cas options for multireference wavefunctions
    """

    # Set the number of closed and occupied orbitals
    elec_cnt = automol.formula.electron_count(formula)
    closed_orb = (elec_cnt - num_act_elc) // 2
    occ_orb = closed_orb + num_act_orb
    if add_two_closed:
        closed_orb -= 2

    # Set the spin and charge values for the species
    spin = spc_info[2] - 1
    chg = spc_info[1]

    # Combine into a CASSCF options for elstruct
    cas_opt = [
        elstruct.option.specify(
            elstruct.Option.Scf.MAXITER_, 40),
        elstruct.option.specify(
            elstruct.Option.Casscf.OCC_, occ_orb),
        elstruct.option.specify(
            elstruct.Option.Casscf.CLOSED_, closed_orb),
        elstruct.option.specify(
            elstruct.Option.Casscf.WFN_, elec_cnt, 1, spin, chg, num_states)
        ]

    return cas_opt


# CONSTRUCT MULTIREFERENCE WAVEFUNCTION STRINGS FOR COMPLEX GUESSES
def multiref_wavefunction_guess(high_mul, zma,
                                spc_info, mod_thy_info,
                                casscf_options):
    """ Prepare wavefunction template for multireference electronic structure calcs
    """

    # Set variables for the programs
    [_, charge, mul] = spc_info
    [prog, _, basis, _] = mod_thy_info

    # Write a string to for high-spin UHF wfn calculation
    guess_str1 = elstruct.writer.energy(
        geom=zma,
        charge=charge,
        mult=high_mul,
        method='hf',
        basis=basis,
        prog=prog,
        orb_type='UU',
        mol_options=['nosym'],
        )
    guess_str1 += '\n\n'
    guess_str1 = '\n'.join(guess_str1.splitlines()[2:-6])

    # Write a string for low-spin CASSCF wfn calc
    guess_str2 = elstruct.writer.energy(
        geom=zma,
        charge=charge,
        mult=mul,
        method='casscf',
        basis=basis,
        prog=prog,
        orb_type='RR',
        casscf_options=casscf_options[0],
        mol_options=['nosym'],
        )
    guess_str2 += '\n\n'
    guess_str2 = '\n'.join(guess_str2.splitlines()[2:-6])

    # Combine two strings together
    guess_str = guess_str1 + '\n' + guess_str2 + '\n'

    # Write a second string for low-spin, lg active space CASSCF wfn calc
    if len(casscf_options) > 1:
        guess_str3 = elstruct.writer.energy(
            geom=zma,
            charge=charge,
            mult=mul,
            method='casscf',
            basis=basis,
            prog=prog,
            orb_type='RR',
            casscf_options=casscf_options[1],
            mol_options=['nosym'],
            )
        guess_str3 += '\n\n'
        guess_str3 = '\n'.join(guess_str3.splitlines()[2:])
        guess_str += guess_str3 + '\n'

    return guess_str


# STRING WRITER FOR VRC-TST TO REPLACE
def wfn_string(ts_info, mod_var_scn_thy_info, inf_sep_ene, cas_kwargs):
    """  Prepare a wfn string for VRC-TST
    """

    method_dct = {
        'caspt2': 'rs2',
        'caspt2c': 'rs2c',
    }
    method = method_dct[mod_var_scn_thy_info[1]]

    # Set the lines for methods
    method_lines = (
        "if (iterations.ge.0) then",
        "  {{{},shift=0.25}}".format(method),
        "  molpro_energy = energy + {}".format(inf_sep_ene),
        "else",
        "  molpro_energy = 10.0"
    )

    # Update the cas kwargs dct
    gen_lines_dct = cas_kwargs.get('gen_lines')
    gen_lines_dct.update({3: method_lines})
    cas_kwargs.update({'gen_lines': gen_lines_dct})

    inp_str = elstruct.writer.energy(
        geom='GEOMETRY_HERE',
        charge=ts_info[1],
        mult=ts_info[2],
        method='casscf',
        basis=mod_var_scn_thy_info[2],
        prog=mod_var_scn_thy_info[0],
        orb_type=mod_var_scn_thy_info[3],
        **cas_kwargs
        )

    return inp_str
