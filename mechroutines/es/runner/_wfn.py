""" Set runtime parameters for multireference electronic structure calculations
"""

import copy
import automol
import elstruct
from phydat import act_space
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import rxn as rinfo
from mechlib.amech_io import printer as ioprinter


# BUILD THE MULTIREFERENCE WAVEFUNCTION TO PASS MOLPRO ELSTRUCT
def update_kwargs_for_multireference(kwargs, cas_kwargs):
    """ Update the elstruct kwargs dict to handle the multiref params
    """

    gen_lines = kwargs['gen_lines']
    cas_gen_lines = cas_kwargs['gen_lines']
    gen_lines = {
        1: cas_gen_lines[1] + gen_lines[1],
        2: gen_lines[2],
        3: gen_lines[3]
    }

    new_kwargs = copy.deepcopy(kwargs)
    new_kwargs.update(cas_kwargs)
    new_kwargs.update(gen_lines)

    return new_kwargs


def multireference_calculation_parameters(zma, spc_info, hs_spc_info,
                                          aspace, mod_thy_info, rxn_info=None):
    """ Prepares a keyword-argument dictionary that can be utilized by the
        elstruct library to perform multireference electronic structure
        calculations. The function is tasked with preparing two parts:

        (1) A string which is placed at the top of the input file that
            contains all neccessary commands for calculating a guess
            multireference wavefunction for a given program, and
        (2) Sets active-space casscf options used for the final part
            of the calculation that matches the user's desired electronic
            structure theory level

        Symmetry is also automatically turned off for all calculations

        Currently, this function only preps Molpro calculations
    """

    if aspace is not None:
        num_act_orb, num_act_elc, num_states, guess_str = aspace
        guess_lines = guess_str.splitlines()
        casscf_options = cas_options(
            zma, spc_info, num_act_elc, num_act_orb, num_states,
            add_two_closed=False)
        ioprinter.info_message('Using wfn guess from file...', newline=1)
    else:
        _inf = rxn_info if rxn_info is not None else spc_info
        typ = 'ts' if rxn_info is not None else 'spc'
        print('inf test', _inf, typ)
        num_act_orb, num_act_elc, num_states = active_space(
            _inf, typ=typ)

        # Build the elstruct CASSCF options list used to build the wfn guess
        # (1) Build wfn with active space
        # (2) Build wfn with active space + 2 closed orbitals for stability
        cas_opt = (
            cas_options(
                zma, spc_info, num_act_elc, num_act_orb, num_states,
                add_two_closed=False),
            cas_options(
                zma, spc_info, num_act_elc, num_act_orb, num_states,
                add_two_closed=True)
        )

        # Write string that has all the components for building the wfn guess
        guess_str = multiref_wavefunction_guess(
            zma, spc_info, hs_spc_info, mod_thy_info, cas_opt)
        guess_lines = guess_str.splitlines()

        # Set casscf options
        casscf_options = cas_opt[0]

        ioprinter.info_message(
                'Generating wfn guess from internal options...', newline=1)

    # Manipulate the opt kwargs to use the wavefunction guess
    cas_kwargs = {
        'casscf_options': casscf_options,
        'gen_lines': {1: guess_lines},
        'mol_options': ('nosym',)  # Turn off symmetry
    }

    return cas_kwargs


# BUILD THE CASSCF OPTIONS
def active_space(info_obj, typ='ts'):
    """ Determine the active space for the multireference MEP scan
    """

    def _active_space(ich, mul):
        """ Determine the active sapce for an InChI string
        """
        if ich in act_space.DCT:
            num_act_orb = act_space.DCT[ich][0]
            num_act_elc = act_space.DCT[ich][1]
            num_states = act_space.DCT[ich][2]
        else:
            num_act_orb = (mul - 1)
            num_act_elc = (mul - 1)
            num_states = 1

        print('active test')
        print(ich, mul)
        print(num_act_orb, num_act_elc, num_states)

        return num_act_orb, num_act_elc, num_states

    if typ == 'spc':
        ich = sinfo.value(info_obj, 'inchi')
        mul = sinfo.value(info_obj, 'mult')
        num_act_orb, num_act_elec, num_states = _active_space(ich, mul)
    elif typ == 'ts':
        rct_ichs = rinfo.value(info_obj, 'inchi')[0]
        rct_muls = rinfo.value(info_obj, 'mult')[0]

        print('rct test', rct_ichs, rct_muls)

        num_act_orb, num_act_elec, num_states = 0, 0, 1
        for ich, mul in zip(rct_ichs, rct_muls):
            norb, nelec, nstat = _active_space(ich, mul)
            num_act_orb += norb
            num_act_elec += nelec
            num_states *= nstat
    print('active test2')
    print(num_act_orb, num_act_elec, num_states)

    return num_act_orb, num_act_elec, num_states


def cas_options(zma, spc_info, num_act_elc, num_act_orb, num_states,
                add_two_closed=False):
    """ Prepare values prepare cas options for multireference wavefunctions
    """

    # Set the number of closed and occupied orbitals
    fml = automol.zmat.formula(zma)
    elec_cnt = automol.formula.electron_count(fml)

    print('count test', elec_cnt)

    closed_orb = (elec_cnt - num_act_elc) // 2
    occ_orb = closed_orb + num_act_orb
    if add_two_closed:
        closed_orb -= 2

    # Set the spin and charge values for the species
    spin = spc_info[2] - 1
    chg = spc_info[1]

    # Combine into a CASSCF options for elstruct
    cas_opt = (
        elstruct.option.specify(
            elstruct.Option.Scf.MAXITER_, 40),
        elstruct.option.specify(
            elstruct.Option.Casscf.OCC_, occ_orb),
        elstruct.option.specify(
            elstruct.Option.Casscf.CLOSED_, closed_orb),
        elstruct.option.specify(
            elstruct.Option.Casscf.WFN_, elec_cnt, 1, spin, chg, num_states)
    )

    return cas_opt


# CONSTRUCT MULTIREFERENCE WAVEFUNCTION STRINGS FOR COMPLEX GUESSES
def multiref_wavefunction_guess(zma, spc_info, hs_spc_info,
                                mod_thy_info, casscf_options):
    """ Prepare wavefunction template for multireference electronic structure calcs
    """

    # Set variables for the programs
    [_, charge, mul] = spc_info
    [_, _, high_mul] = hs_spc_info
    prog, _, basis, _ = mod_thy_info

    # Write a string to for high-spin UHF wfn calculation
    guess_str1 = elstruct.writer.energy(
        geo=zma,
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
        geo=zma,
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
            geo=zma,
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
        f"  {{{method},shift=0.25}}",
        f"  molpro_energy = energy + {inf_sep_ene}",
        "else",
        "  molpro_energy = 10.0"
    )

    # Update the cas kwargs dct
    gen_lines_dct = cas_kwargs.get('gen_lines')
    gen_lines_dct.update({3: method_lines})
    cas_kwargs.update({'gen_lines': gen_lines_dct})

    inp_str = elstruct.writer.energy(
        geo='GEOMETRY_HERE',
        charge=ts_info[1],
        mult=ts_info[2],
        method='casscf',
        basis=mod_var_scn_thy_info[2],
        prog=mod_var_scn_thy_info[0],
        orb_type=mod_var_scn_thy_info[3],
        **cas_kwargs
        )

    return inp_str
