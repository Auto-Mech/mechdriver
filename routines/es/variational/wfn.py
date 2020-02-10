""" libraries to build wavefunctions for Molpro multireference calculations
"""

import automol
import elstruct
from lib.phydat import act_space


def cas_options(spc_info, formula, num_act_elc, num_act_orb,
                high_mul, add_two_closed=False):
    """ Prepare values prepare cas options for multireference wavefunctions
    """

    # Set electron counts and orbital indexes
    elec_count = automol.formula.electron_count(formula)
    closed_orb = (elec_count - num_act_elc) // 2
    occ_orb = closed_orb + num_act_orb
    if add_two_closed:
        closed_orb -= 2

    # Set the spin and charge values for the species
    two_spin = spc_info[2] - 1
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
            elstruct.Option.Casscf.WFN_, elec_count, 1, two_spin, chg)
        ]

    return cas_opt


def wfn_string(spc_info, formula, num_act_elc, num_act_orb,
               high_mul, add_two_closed=False):
    """ Prepare values prepare cas options for multireference wavefunctions
    """

    # Set electron counts and orbital indexes
    elec_count = automol.formula.electron_count(formula)
    closed_orb = (elec_count - num_act_elc) // 2
    occ_orb = closed_orb + num_act_orb
    if add_two_closed:
        closed_orb -= 2

    # Set the spin and charge values for the species
    two_spin = spc_info[2] - 1
    high_two_spin = high_mul - 1

    # Write a wavefunction card string
    wfn_str = (
        "{{uhf,maxit=300;wf,{0},1,{1}}}\n"
        "{{multi,maxit=40;closed,{2};occ,{3};wf,{0},1,{4};orbprint,3}}"
    ).format(elec_count, high_two_spin, closed_orb, occ_orb, two_spin)

    return wfn_str


def multiref_wavefunction_guess(high_mul, zma,
                                spc_info, thy_level,
                                casscf_options):
    """ Prepare wavefunction template for multireference electronic structure calcs
    """

    # Set variables for the programs
    charge = spc_info[1]
    mul = spc_info[2]
    basis = thy_level[2]
    prog = thy_level[0]
    prog = 'molpro2015'

    # Write a string to for high-spin UHF wfn calculation
    guess_str1 = elstruct.writer.energy(
        geom=zma,
        charge=charge,
        mult=high_mul,
        method='hf',
        basis=basis,
        prog=prog,
        orb_restricted=False,
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
        orb_restricted=True,
        casscf_options=casscf_options[0],
        mol_options=['nosym'],
        )
    guess_str2 += '\n\n'
    guess_str2 = '\n'.join(guess_str2.splitlines()[2:-6])

    # Combine two strings together
    guess_str = guess_str1 + guess_str2

    # Write a second string for low-spin, lg active space CASSCF wfn calc
    if len(casscf_options) > 1:
        guess_str3 = elstruct.writer.energy(
            geom=zma,
            charge=charge,
            mult=mul,
            method='casscf',
            basis=basis,
            prog=prog,
            orb_restricted=True,
            casscf_options=casscf_options[1],
            mol_options=['nosym'],
            )
        guess_str3 += '\n\n'
        guess_str3 = '\n'.join(guess_str3.splitlines()[2:])
        guess_str += guess_str3

    return guess_str


def active_space(ts_dct, spc_dct, ts_high_mul):
    """ Determine the active space for the multireference MEP scan
    """
    rcts = ts_dct['reacs']
    num_act_orb, num_act_elc = 0, 0
    for rct in rcts:
        rct_ich = spc_dct[rct]['ich']
        if rct_ich in act_space.DCT:
            num_act_orb += act_space.DCT[rct_ich][0]
            num_act_elc += act_space.DCT[rct_ich][1]
        else:
            num_act_orb += (ts_high_mul - 1)
            num_act_elc += (ts_high_mul - 1)

    return num_act_orb, num_act_elc
