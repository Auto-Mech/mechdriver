""" Electronic structure program input writing module.
"""

from elstruct import par
from elstruct.writer import program_modules as pm


# energy input writers
def programs():
    """ Constructs a list of available electronic structure programs.
        At minimum, each program must have an energy reader to be enumerated.
    """
    return pm.program_modules_with_function(pm.Job.ENERGY)


def energy(prog, geo, charge, mult, method, basis,
           # molecule options
           mol_options=(),
           # machine options
           memory=1, comment='', machine_options=(),
           # theory options
           orb_type='RU',
           scf_options=(), casscf_options=(), corr_options=(),
           # generic options
           gen_lines=None):
    """ Writes an input file string for an electronic energy calculation
        for a specified electronic structure program.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param geo: cartesian or z-matrix geometry
        :type geo: tuple
        :param charge: molecular charge
        :type charge: int
        :param mult: spin multiplicity
        :type mult: int
        :param method: electronic structure method
        :type method: str
        :param basis: basis set
        :type basis: str
        :param mol_options: options for the molecule block
        :type mol_options: tuple[str]
        ;param memory: memory in GB
        :type memory: int
        :param comment: a comment string to be placed at the top of the file
        :type comment: str
        :param machine_options: machine directives
            (num procs, num threads, etc.)
        :type machine_options: tuple[str]
        :param orb_type: 'R' indicates restricted orbitals, 'U' indicates
            unrestricted orbitals; can also be 'RR', 'RU', or 'UU'.
            Where first (second) character sets R/U for singlets (multiplets)
        :type orb_type: str
        :param scf_options: scf method directives
        :type scf_options: tuple[str]
        :param casscf_options: casscf method directives
        :type casscf_options: tuple[str]
        :param corr_options: correlation method directives
        :type corr_options: tuple[str]
        :param gen_lines: generic lines for the input file
        :type gen_lines: dict[idx:str]
    """
    prog, method, basis, orb_restricted = _process_theory_specifications(
        prog, method, basis, mult, orb_type)

    return pm.call_module_function(
        prog, pm.Job.ENERGY,
        # *args
        geo, charge, mult, method, basis, orb_restricted,
        # **kwargs
        # molecule options
        mol_options=mol_options,
        # machine options
        memory=memory, comment=comment, machine_options=machine_options,
        # theory options
        scf_options=scf_options, casscf_options=casscf_options,
        corr_options=corr_options,
        # generic options
        gen_lines=gen_lines,
        # job options
        job_options=(), frozen_coordinates=(),
        saddle=False)


# gradient input writers
def gradient_programs():
    """ Constructs a list of program modules implementing
        gradient input writers.
    """
    return pm.program_modules_with_function(pm.Job.GRADIENT)


def gradient(prog, geo, charge, mult, method, basis,
             # molecule options
             mol_options=(),
             # machine options
             memory=1, comment='', machine_options=(),
             # theory options
             orb_type='RU',
             scf_options=(), casscf_options=(), corr_options=(),
             # generic options
             gen_lines=None,
             # job options
             job_options=()):
    """ Writes an input file string for a gradient calculation
        for a specified electronic structure program.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param geo: cartesian or z-matrix geometry
        :type geo: tuple
        :param charge: molecular charge
        :type charge: int
        :param mult: spin multiplicity
        :type mult: int
        :param method: electronic structure method
        :type method: str
        :param basis: basis set
        :type basis: str
        :param mol_options: options for the molecule block
        :type mol_options: tuple[str]
        ;param memory: memory in GB
        :type memory: int
        :param comment: a comment string to be placed at the top of the file
        :type comment: str
        :param machine_options: machine directives
            (num procs, num threads, etc.)
        :type machine_options: tuple[str]
        :param orb_type: 'R' indicates restricted orbitals, 'U' indicates
            unrestricted orbitals; can also be 'RR', 'RU', or 'UU'.
            Where first (second) character sets R/U for singlets (multiplets)
        :type orb_type: str
        :param scf_options: scf method directives
        :type scf_options: tuple[str]
        :param casscf_options: casscf method directives
        :type casscf_options: tuple[str]
        :param corr_options: correlation method directives
        :type corr_options: tuple[str]
        :param gen_lines: generic lines for the input file
        :type gen_lines: dict[idx:str]
    """
    prog, method, basis, orb_restricted = _process_theory_specifications(
        prog, method, basis, mult, orb_type)

    return pm.call_module_function(
        prog, pm.Job.GRADIENT,
        # *args
        geo, charge, mult, method, basis, orb_restricted,
        # **kwargs
        # molecule options
        mol_options=mol_options,
        # machine options
        memory=memory, comment=comment, machine_options=machine_options,
        # theory options
        scf_options=scf_options, casscf_options=casscf_options,
        corr_options=corr_options,
        # generic options
        gen_lines=gen_lines,
        # job options
        job_options=job_options, frozen_coordinates=(),
        saddle=False)


# hessian input writers
def hessian_programs():
    """ Constructs a list of program modules implementing
        Hessian input writers.
    """
    return pm.program_modules_with_function(pm.Job.HESSIAN)


def hessian(prog, geo, charge, mult, method, basis,
            # molecule options
            mol_options=(),
            # machine options
            memory=1, comment='', machine_options=(),
            # theory options
            orb_type='RU',
            scf_options=(), casscf_options=(), corr_options=(),
            # generic options
            gen_lines=None,
            # job options
            job_options=()):
    """ Writes an input file string for a Hessian calculation
        for a specified electronic structure program.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param geo: cartesian or z-matrix geometry
        :type geo: tuple
        :param charge: molecular charge
        :type charge: int
        :param mult: spin multiplicity
        :type mult: int
        :param method: electronic structure method
        :type method: str
        :param basis: basis set
        :type basis: str
        :param mol_options: options for the molecule block
        :type mol_options: tuple[str]
        ;param memory: memory in GB
        :type memory: int
        :param comment: a comment string to be placed at the top of the file
        :type comment: str
        :param machine_options: machine directives
            (num procs, num threads, etc.)
        :type machine_options: tuple[str]
        :param orb_type: 'R' indicates restricted orbitals, 'U' indicates
            unrestricted orbitals; can also be 'RR', 'RU', or 'UU'.
            Where first (second) character sets R/U for singlets (multiplets)
        :type orb_type: str
        :param scf_options: scf method directives
        :type scf_options: tuple[str]
        :param casscf_options: casscf method directives
        :type casscf_options: tuple[str]
        :param corr_options: correlation method directives
        :type corr_options: tuple[str]
        :param gen_lines: generic lines for the input file
        :type gen_lines: dict[idx:str]
    """
    prog, method, basis, orb_restricted = _process_theory_specifications(
        prog, method, basis, mult, orb_type)

    return pm.call_module_function(
        prog, pm.Job.HESSIAN,
        # *args
        geo, charge, mult, method, basis, orb_restricted,
        # **kwargs
        # molecule options
        mol_options=mol_options,
        # machine options
        memory=memory, comment=comment, machine_options=machine_options,
        # theory options
        scf_options=scf_options, casscf_options=casscf_options,
        corr_options=corr_options,
        # generic options
        gen_lines=gen_lines,
        # job options
        job_options=job_options, frozen_coordinates=(),
        saddle=False)


# vpt2 input writers
def vpt2_programs():
    """ Constructs a list of program modules implementing
        2nd-order vibrational perturbation theory (VPT2) input writers.
    """
    return pm.program_modules_with_function(pm.Job.VPT2)


def vpt2(prog, geo, charge, mult, method, basis,
         # molecule options
         mol_options=(),
         # machine options
         memory=1, comment='', machine_options=(),
         # theory options
         orb_type='RU',
         scf_options=(), casscf_options=(), corr_options=(),
         # generic options
         gen_lines=None,
         # job options
         job_options=()):
    """ Writes an input file string for a
        2nd-order vibrational perturbation theory calculation
        for a specified electronic structure program.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param geo: cartesian or z-matrix geometry
        :type geo: tuple
        :param charge: molecular charge
        :type charge: int
        :param mult: spin multiplicity
        :type mult: int
        :param method: electronic structure method
        :type method: str
        :param basis: basis set
        :type basis: str
        :param mol_options: options for the molecule block
        :type mol_options: tuple[str]
        ;param memory: memory in GB
        :type memory: int
        :param comment: a comment string to be placed at the top of the file
        :type comment: str
        :param machine_options: machine directives
            (num procs, num threads, etc.)
        :type machine_options: tuple[str]
        :param orb_type: 'R' indicates restricted orbitals, 'U' indicates
            unrestricted orbitals; can also be 'RR', 'RU', or 'UU'.
            Where first (second) character sets R/U for singlets (multiplets)
        :type orb_type: str
        :param scf_options: scf method directives
        :type scf_options: tuple[str]
        :param casscf_options: casscf method directives
        :type casscf_options: tuple[str]
        :param corr_options: correlation method directives
        :type corr_options: tuple[str]
        :param gen_lines: generic lines for the input file
        :type gen_lines: dict[idx:str]
    """
    prog, method, basis, orb_restricted = _process_theory_specifications(
        prog, method, basis, mult, orb_type)

    return pm.call_module_function(
        prog, pm.Job.VPT2,
        # *args
        geo, charge, mult, method, basis, orb_restricted,
        # **kwargs
        # molecule options
        mol_options=mol_options,
        # machine options
        memory=memory, comment=comment, machine_options=machine_options,
        # theory options
        scf_options=scf_options, casscf_options=casscf_options,
        corr_options=corr_options,
        # generic options
        gen_lines=gen_lines,
        # job options
        job_options=job_options, frozen_coordinates=(),
        saddle=False)


# molec_properties input writers
def molecular_properties_programs():
    """ Constructs a list of program modules implementing
        molecular properties, including the
        dipole moment and polarizability, input writers.
    """
    return pm.program_modules_with_function(pm.Job.MOLPROP)


def molecular_properties(prog, geo, charge, mult, method, basis,
                         # molecule options
                         mol_options=(),
                         # machine options
                         memory=1, comment='', machine_options=(),
                         # theory options
                         orb_type='RU',
                         scf_options=(), casscf_options=(), corr_options=(),
                         # generic options
                         gen_lines=None,
                         # job options
                         job_options=()):
    """ Writes an input file string for molecular properties calculations,
        including the dipole moment and polarizability,
        for a specified electronic structure program.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param geo: cartesian or z-matrix geometry
        :type geo: tuple
        :param charge: molecular charge
        :type charge: int
        :param mult: spin multiplicity
        :type mult: int
        :param method: electronic structure method
        :type method: str
        :param basis: basis set
        :type basis: str
        :param mol_options: options for the molecule block
        :type mol_options: tuple[str]
        ;param memory: memory in GB
        :type memory: int
        :param comment: a comment string to be placed at the top of the file
        :type comment: str
        :param machine_options: machine directives
            (num procs, num threads, etc.)
        :type machine_options: tuple[str]
        :param orb_type: 'R' indicates restricted orbitals, 'U' indicates
            unrestricted orbitals; can also be 'RR', 'RU', or 'UU'.
            Where first (second) character sets R/U for singlets (multiplets)
        :type orb_type: str
        :param scf_options: scf method directives
        :type scf_options: tuple[str]
        :param casscf_options: casscf method directives
        :type casscf_options: tuple[str]
        :param corr_options: correlation method directives
        :type corr_options: tuple[str]
        :param gen_lines: generic lines for the input file
        :type gen_lines: dict[idx:str]
    """
    prog, method, basis, orb_restricted = _process_theory_specifications(
        prog, method, basis, mult, orb_type)

    return pm.call_module_function(
        prog, pm.Job.MOLPROP,
        # *args
        geo, charge, mult, method, basis, orb_restricted,
        # **kwargs
        # molecule options
        mol_options=mol_options,
        # machine options
        memory=memory, comment=comment, machine_options=machine_options,
        # theory options
        scf_options=scf_options, casscf_options=casscf_options,
        corr_options=corr_options,
        # generic options
        gen_lines=gen_lines,
        # job options
        job_options=job_options, frozen_coordinates=(),
        saddle=False)


# irc input writers
def irc_programs():
    """ Constructs a list of program modules implementing
        Intrinsic Reaction Coordinate input writers.
    """
    return pm.program_modules_with_function(pm.Job.IRC)


def irc(prog, geo, charge, mult, method, basis,
        # molecule options
        mol_options=(),
        # machine options
        memory=1, comment='', machine_options=(),
        # theory options
        orb_type='RU',
        scf_options=(), casscf_options=(), corr_options=(),
        # generic options
        gen_lines=None,
        # job options
        job_options=(), frozen_coordinates=()):
    """ Writes an input file string for an Intrinsic Reaction Coordinate
        calculation for a specified electronic structure program.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param geo: cartesian or z-matrix geometry
        :type geo: tuple
        :param charge: molecular charge
        :type charge: int
        :param mult: spin multiplicity
        :type mult: int
        :param method: electronic structure method
        :type method: str
        :param basis: basis set
        :type basis: str
        :param mol_options: options for the molecule block
        :type mol_options: tuple[str]
        ;param memory: memory in GB
        :type memory: int
        :param comment: a comment string to be placed at the top of the file
        :type comment: str
        :param machine_options: machine directives
            (num procs, num threads, etc.)
        :type machine_options: tuple[str]
        :param orb_type: 'R' indicates restricted orbitals, 'U' indicates
            unrestricted orbitals; can also be 'RR', 'RU', or 'UU'.
            Where first (second) character sets R/U for singlets (multiplets)
        :type orb_type: str
        :param scf_options: scf method directives
        :type scf_options: tuple[str]
        :param casscf_options: casscf method directives
        :type casscf_options: tuple[str]
        :param corr_options: correlation method directives
        :type corr_options: tuple[str]
        :param job_options: geometry optimization routine directives
        :type job_options: tuple[str]
        :param frozen_coordinates: only with z-matrix geometries; list of
            coordinate names to freeze
        :type fozen_coordinates: tuple[str]
        :param gen_lines: generic lines for the input file
        :type gen_lines: dict[idx:str]
    """
    prog, method, basis, orb_restricted = _process_theory_specifications(
        prog, method, basis, mult, orb_type)

    return pm.call_module_function(
        prog, pm.Job.IRC,
        # *args
        geo, charge, mult, method, basis, orb_restricted,
        # **kwargs
        # molecule options
        mol_options=mol_options,
        # machine options
        memory=memory, comment=comment, machine_options=machine_options,
        # theory options
        scf_options=scf_options, casscf_options=casscf_options,
        corr_options=corr_options,
        # generic options
        gen_lines=gen_lines,
        # job options
        job_options=job_options, frozen_coordinates=frozen_coordinates,
        saddle=False)


# optimization input writers
def optimization_programs():
    """ Constructs a list of program modules implementing
        geometry optimization input writers.
    """
    return pm.program_modules_with_function(pm.Job.OPTIMIZATION)


def optimization(prog, geo, charge, mult, method, basis,
                 # molecule options
                 mol_options=(),
                 # machine options
                 memory=1, comment='', machine_options=(),
                 # theory options
                 orb_type='RU',
                 scf_options=(), casscf_options=(), corr_options=(),
                 # generic options
                 gen_lines=None,
                 # job options
                 job_options=(), frozen_coordinates=(), saddle=False):
    """ Writes an input file string for a geometry optimization
        calculation for a specified electronic structure program.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param geo: cartesian or z-matrix geometry
        :type geo: tuple
        :param charge: molecular charge
        :type charge: int
        :param mult: spin multiplicity
        :type mult: int
        :param method: electronic structure method
        :type method: str
        :param basis: basis set
        :type basis: str
        :param mol_options: options for the molecule block
        :type mol_options: tuple[str]
        ;param memory: memory in GB
        :type memory: int
        :param comment: a comment string to be placed at the top of the file
        :type comment: str
        :param machine_options: machine directives
            (num procs, num threads, etc.)
        :type machine_options: tuple[str]
        :param orb_type: 'R' indicates restricted orbitals, 'U' indicates
            unrestricted orbitals; can also be 'RR', 'RU', or 'UU'.
            Where first (second) character sets R/U for singlets (multiplets)
        :type orb_type: str
        :param scf_options: scf method directives
        :type scf_options: tuple[str]
        :param casscf_options: casscf method directives
        :type casscf_options: tuple[str]
        :param corr_options: correlation method directives
        :type corr_options: tuple[str]
        :param job_options: geometry optimization routine directives
        :type job_options: tuple[str]
        :param frozen_coordinates: only with z-matrix geometries; list of
            coordinate names to freeze
        :type fozen_coordinates: tuple[str]
        :param saddle: optimize a saddle point?
        :type saddle: bool
        :param gen_lines: generic lines for the input file
        :type gen_lines: dict[idx:str]
    """
    prog, method, basis, orb_restricted = _process_theory_specifications(
        prog, method, basis, mult, orb_type)

    return pm.call_module_function(
        prog, pm.Job.OPTIMIZATION,
        # *args
        geo, charge, mult, method, basis, orb_restricted,
        # **kwargs
        # molecule options
        mol_options=mol_options,
        # machine options
        memory=memory, comment=comment, machine_options=machine_options,
        # theory options
        scf_options=scf_options, casscf_options=casscf_options,
        corr_options=corr_options,
        # generic options
        gen_lines=gen_lines,
        # job options
        job_options=job_options, frozen_coordinates=frozen_coordinates,
        saddle=saddle)


def _process_theory_specifications(prog, method, basis, mult, orb_type):
    """ Process the theory method including the orbital type conversion.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param method: electronic structure method
        :type method: str
        :param basis: basis set
        :type basis: str
        :param mult: spin multiplicity
        :type mult: int
        :param orb_type: 'R' indicates restricted orbitals, 'U' indicates
            unrestricted orbitals; can also be 'RR', 'RU', or 'UU'.
            Where first (second) character sets R/U for singlets (multiplets)
        :type orb_type: str
        :rtype: (str, str, str, str)
    """

    assert par.is_program(prog)

    # determine the orbital restriction
    singlet = (mult == 1)
    if len(orb_type) == 2:
        orb_type = orb_type[0] if singlet else orb_type[1]

    assert orb_type in ('R', 'U')
    orb_restricted = (orb_type == 'R')

    # for non-standard DFT/Basis, the user can input whatever they want
    if not par.Method.is_nonstandard_dft(method):
        core_method, _ = par.Method.evaluate_method_type(method)
        assert par.is_program_method(prog, core_method)
        assert par.is_program_method_orbital_type(
            prog, core_method, singlet, orb_type)

        prog = par.standard_case(prog)
        method = par.standard_case(method)

    if not par.Basis.is_nonstandard_basis(basis):
        assert par.is_program_basis(prog, basis)
        basis = par.standard_case(basis)

    return prog, method, basis, orb_restricted
