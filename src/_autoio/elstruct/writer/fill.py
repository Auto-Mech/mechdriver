"""
   Useful functions used for all the program writers
"""


import automol
import autowrite as aw
from elstruct.par import Reference, Program
from elstruct.option import is_valid as _option_is_valid
from elstruct.option import name as _option_name
from elstruct.pclass import values as _pclass_values
from elstruct.par import Option
from elstruct.par import program_basis_name
from elstruct.par import program_method_name
from elstruct.par import Method


class TemplateKey():
    """ mako template keys """
    # machine
    MEMORY = 'memory'
    MACHINE_OPTIONS = 'machine_options'
    NPROCS = 'nprocs'
    # theoretical method
    REFERENCE = 'reference'
    METHOD = 'method'
    BASIS = 'basis'
    SCF_OPTIONS = 'scf_options'
    SCF_GUESS_OPTIONS = 'scf_guess_options'
    XCGRID = 'xcgrid'
    # molecule / state
    MOL_OPTIONS = 'mol_options'
    CHARGE = 'charge'
    MULT = 'mult'
    SPIN = 'spin'
    GEOM = 'geom'
    ZMAT_VALS = 'zmat_vals'
    ZMAT_VAR_VALS = 'zmat_var_vals'
    ZMAT_CONST_VALS = 'zmat_const_vals'
    FROZEN_DIS_STRS = 'frozen_dis_strs'
    FROZEN_ANG_STRS = 'frozen_ang_strs'
    FROZEN_DIH_STRS = 'frozen_dih_strs'
    COORD_SYS = 'coord_sys'
    SADDLE = 'saddle'
    NUMERICAL = 'numerical'
    # job
    COMMENT = 'comment'
    JOB_KEY = 'job_key'
    JOB_OPTIONS = 'job_options'
    GEN_LINES = 'gen_lines'
    GEN_LINES_1 = 'gen_lines_1'
    GEN_LINES_2 = 'gen_lines_2'
    GEN_LINES_3 = 'gen_lines_3'
    # theoretical method
    BASIS = 'basis'
    SCF_METHOD = 'scf_method'
    SCF_OPTIONS = 'scf_options'
    ISMULTIREF = 'ismultiref'
    CASSCF_OPTIONS = 'casscf_options'
    CORR_METHOD = 'corr_method'
    CORR_OPTIONS = 'corr_options'


# Format strings
def geometry_strings(geo, frozen_coordinates, zma_sign='=', offset=False):
    """ Build the string for the input geometry

        :param geo: cartesian or z-matrix geometry
        :type geo: tuple
        :param frozen_coordinates: only with z-matrix geometries; list of
            coordinate names to freeze
        :type fozen_coordinates: tuple[str]
        :rtype: (str, str)
    """

    def _offset_planar_vals(val):
        """ Change 180.0 to 179.9 and 0.0 to 0.1 to deal with
            breaking symmetry.
        """
        offset_val_dct = {}
        for key, val in val_dct.items():
            _val = round(val, 1)
            if _val == 180.0:
                _val = 179.9
            elif _val == 0.0:
                _val = 0.1
            else:
                _val = val
            offset_val_dct[key] = _val

        return offset_val_dct

    if automol.geom.is_valid(geo):
        geo_str = automol.geom.string(geo)
        zmat_vval_str = ''
        zmat_cval_str = ''
    elif automol.zmat.is_valid(geo):
        zma = geo
        symbs = automol.zmat.symbols(zma)
        key_mat = automol.zmat.key_matrix(zma, shift=1)
        name_mat = automol.zmat.name_matrix(zma)
        val_dct = automol.zmat.value_dictionary(
            zma, angstrom=True, degree=True)
        geo_str = aw.zmat.matrix_block(symbs, key_mat, name_mat)

        vval_dct = {key: val for key, val in val_dct.items()
                    if key not in frozen_coordinates}
        cval_dct = {key: val for key, val in val_dct.items()
                    if key in frozen_coordinates}

        if offset:
            vval_dct = _offset_planar_vals(vval_dct)
            cval_dct = _offset_planar_vals(cval_dct)

        zmat_vval_str = aw.zmat.setval_block(
            vval_dct, setval_sign=zma_sign).strip()
        zmat_cval_str = aw.zmat.setval_block(
            cval_dct, setval_sign=zma_sign).strip()
    elif geo in ('GEOMETRY', 'GEOMETRY_HERE'):
        geo_str = geo
        zmat_cval_str = ''
        zmat_vval_str = ''
    else:
        raise ValueError(f"Invalid geometry value:\n{geo}")

    return geo_str, zmat_vval_str, zmat_cval_str


def build_gen_lines(gen_lines, line1=None, line2=None, line3=None):
    """ Set three lines for writing in various blocks of files.
        Function either grabs lines from the dictionary and if nothing
        present, then uses value provided by function
    """

    if gen_lines is not None:
        gen_lines_1 = '\n'.join(gen_lines[1]) if 1 in gen_lines else ''
        gen_lines_2 = '\n'.join(gen_lines[2]) if 2 in gen_lines else ''
        gen_lines_3 = '\n'.join(gen_lines[3]) if 3 in gen_lines else ''
    else:
        gen_lines_1 = ''
        gen_lines_2 = ''
        gen_lines_3 = ''

    if not gen_lines_1:
        gen_lines_1 = line1 if line1 is not None else ''
    if not gen_lines_2:
        gen_lines_2 = line2 if line2 is not None else ''
    if not gen_lines_3:
        gen_lines_3 = line3 if line3 is not None else ''

    return gen_lines_1, gen_lines_2, gen_lines_3


def update_gen_lines(gen_lines,
                     lines1=None, lines2=None, lines3=None):
    """ Update gen lines dictionary with new input lines
    """

    if gen_lines is not None:
        gen_lines_1 = gen_lines.get(1)
        gen_lines_2 = gen_lines.get(2)
        gen_lines_3 = gen_lines.get(3)

        if lines1 is not None:
            if gen_lines_1 is not None:
                gen_lines_1 += tuple(lines1)
            else:
                gen_lines_1 = tuple(lines1)
        if lines2 is not None:
            if gen_lines_2 is not None:
                gen_lines_2 += tuple(lines2)
            else:
                gen_lines_2 = tuple(lines2)
        if lines3 is not None:
            if gen_lines_3 is not None:
                gen_lines_3 += tuple(lines3)
            else:
                gen_lines_3 = tuple(lines3)

        gen_lines = {}
        for i, gen_line in enumerate([gen_lines_1, gen_lines_2, gen_lines_3]):
            if gen_line is not None:
                gen_lines[i+1] = gen_line

    else:
        gen_lines = {}
        if lines1:
            gen_lines.update({1: tuple(lines1)})
        if lines2:
            gen_lines.update({2: tuple(lines2)})
        if lines3:
            gen_lines.update({3: tuple(lines3)})

    return gen_lines


# Handle setting options for various programs
def evaluate_options(options, option_eval_dct):
    """ Build a list of program specific options.

        :param options: requested options.
        :type options: tuple(str)
        :param option_eval_dct: program specific values for an options
        :type option_eval_dct: dict[str: str]
        :type: dict[str: str]
    """

    options = list(options)
    option_names = tuple(sorted(option_eval_dct.keys()))

    for idx, option in enumerate(options):
        # Will evaluate option if possible, or just put in (very bad)
        try:
            name = _option_name(option)
            assert name in option_names
            options[idx] = option_eval_dct[name](option)
        except AssertionError:
            options[idx] = option

    return tuple(options)


def intercept_scf_guess_option(options, option_eval_dct):
    """ Set SCF guess options

        :param options: requested options.
        :type options: tuple(str)
        :param option_eval_dct: program specific values for an options
        :type option_eval_dct: dict[str: str]
        :rtype: (tuple(str), tuple(str))
    """

    guess_options = []
    scf_options = []
    for opt in options:
        if (_option_is_valid(opt) and opt in
                _pclass_values(Option.Scf.Guess)):
            guess_options.append(opt)
        else:
            scf_options.append(opt)
    scf_guess_options = evaluate_options(guess_options, option_eval_dct)
    scf_options = evaluate_options(scf_options, option_eval_dct)

    return scf_guess_options, scf_options


# Set the full description of the theoretical method with
def program_method_names(prog, method, basis, mult, orb_restricted):
    """ Sets all the names of all the components of a theoretical method
        to those specific to the program of interest so that a proper input
        file can be written.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param method: electronic structure method
        :type method: str
        :param basis: basis set
        :type basis: str
        :param mult: spin multiplicity
        :type mult: int
        :param orb_restricted: parameter designating if restriced refrence used
        :type orb_restricted: bool
        :rtype: (str, str, str)
    """

    # Set the singlet variable used by prog_method
    singlet = (mult == 1)

    # Determine the reference for the given method
    prog_reference = _reference(prog, method, mult, orb_restricted)

    # Determine the method
    if Method.is_casscf(method):
        prog_method = prog_reference
    elif method == Method.HF[0]:
        if prog in (Program.GAUSSIAN09, Program.GAUSSIAN03,
                    Program.GAUSSIAN16):
            prog_method = prog_reference if prog_reference else method
        else:
            prog_method = program_method_name(prog, method, singlet=singlet)
    else:
        prog_method = program_method_name(prog, method, singlet=singlet)

    # core_prog_method, mod = elstruct.Method.evaluate_method_type(prog_method)

    # Set the basis
    prog_basis = program_basis_name(prog, basis)

    return prog_method, prog_reference, prog_basis
    # return prog_method, prog_reference, prog_basis, mods


def _reference(prog, method, mult, orb_restricted):
    """ Determine the string for what the Hartree-Fock or Kohn-Sham
        reference should be based on the electronic structure method
        and electronic structure program is.

        :param prog: electronic structure program to use as a backend
        :type prog: str
        :param method: electronic structure method
        :type method: str
        :param mult: spin multiplicity
        :type mult: int
        :param orb_restricted: parameter designating if restriced refrence used
        :type orb_restricted: bool
        :rtype: str
    """
    # Need a multiref version
    if Method.is_dft(method) or Method.is_semi_empirical(method):
        reference = _dft_reference(prog, orb_restricted)
    elif method == Method.HF[0]:
        reference = _hf_reference(prog, mult, orb_restricted)
    else:
        reference = _corr_reference(prog, mult, orb_restricted)

    return reference


def _dft_reference(prog, orb_restricted):
    """ dft
    """
    if prog in (Program.GAUSSIAN09, Program.GAUSSIAN03, Program.GAUSSIAN16):
        reference = ''
    else:
        reference = (Reference.RKS if orb_restricted else
                     Reference.UKS)

    return reference


def _hf_reference(prog, mult, orb_restricted):
    """ hf
    """
    if prog in (Program.GAUSSIAN09, Program.GAUSSIAN03, Program.GAUSSIAN16):
        reference = ''
    else:
        if mult == 1:
            reference = Reference.RHF
        else:
            reference = (Reference.ROHF if orb_restricted else
                         Reference.UHF)

    if reference == Reference.ROHF:
        if prog == Program.MOLPRO2015:
            reference = Reference.RHF

    return reference


def _corr_reference(prog, mult, orb_restricted):
    """ correlated method reference
    """
    if mult == 1:
        reference = Reference.RHF
    else:
        reference = (Reference.ROHF if orb_restricted else
                     Reference.UHF)

    if reference == Reference.ROHF:
        if prog == Program.MOLPRO2015:
            reference = Reference.RHF

    return reference
