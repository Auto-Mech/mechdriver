""" cfour2 writer module """

import os
import automol
from ioformat import build_mako_str
import elstruct.par
from elstruct.writer import fill
from elstruct.writer._cfour2._par import OPTION_EVAL_DCT


PROG = elstruct.par.Program.CFOUR2

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_DIR = os.path.join(THIS_DIR, 'templates')


def write_input(job_key, geo, charge, mult, method, basis, orb_restricted,
                # molecule options
                mol_options=(),
                # machine options
                memory=1, comment='', machine_options=(),
                # theory options
                scf_options=(), casscf_options=(), corr_options=(),
                # generic options
                gen_lines=None,
                # job options
                job_options=(), frozen_coordinates=(), saddle=False):
    """ Write an input file string for an electronic structure calculation
        by processing all of the information and using it to fill in
        a Mako template of the input file.

        :param job_key: job contained in the input file
        :type job_key: str
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
        :param orb_restricted: parameter designating if restriced refrence used
        :type orb_restricted: bool
        :param mol_options: options for the molecule block
        :type mol_options: tuple[str]
        ;param memory: memory in GB
        :type memory: int
        :param comment: a comment string to be placed at the top of the file
        :type comment: str
        :param machine_options: machine directives
            (num procs, num threads, etc.)
        :type machine_options: tuple[str]
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
        :param saddle: parameter signifiying a saddle point calculation
        :type saddle: bool
        :param gen_lines: generic lines for the input file
        :type gen_lines: dict[idx:str]
    """

    # set correlated method; check if multiref
    prog_method, prog_reference, prog_basis = fill.program_method_names(
        PROG, method, basis, mult, orb_restricted)

    # Build the geometry object for the job
    geo_str, zmat_val_str, _ = fill.geometry_strings(
        geo, frozen_coordinates)

    geo_str = '\n'.join([' '.join(string.split())
                         for string in geo_str.splitlines()])
    zmat_val_str = '\n'.join([' '.join(string.split())
                              for string in zmat_val_str.splitlines()])
    zmat_val_str += '\n'

    # Set options for coupled cluster
    corr_options = _set_cc_prog(method, prog_reference)

    # Unused options
    _ = mol_options
    _ = machine_options

    scf_options = fill.evaluate_options(scf_options, OPTION_EVAL_DCT)
    casscf_options = fill.evaluate_options(casscf_options, OPTION_EVAL_DCT)
    job_options = fill.evaluate_options(job_options, OPTION_EVAL_DCT)

    numerical = None

    if automol.geom.is_valid(geo):
        coord_sys = 'CARTESIAN'
    elif automol.zmat.is_valid(geo):
        coord_sys = 'INTERNAL'

    # Set the gen lines blocks
    gen_lines_1, _, _ = fill.build_gen_lines(gen_lines)

    # Write the input file string
    fill_dct = {
        fill.TemplateKey.COMMENT: comment,
        fill.TemplateKey.MEMORY: memory,
        fill.TemplateKey.CHARGE: charge,
        fill.TemplateKey.MULT: mult,
        fill.TemplateKey.GEOM: geo_str,
        fill.TemplateKey.COORD_SYS: coord_sys,
        fill.TemplateKey.ZMAT_VAR_VALS: zmat_val_str,
        fill.TemplateKey.BASIS: prog_basis.upper(),
        fill.TemplateKey.METHOD: prog_method.upper(),
        fill.TemplateKey.REFERENCE: prog_reference.upper(),
        fill.TemplateKey.SCF_OPTIONS: '\n'.join(scf_options),
        fill.TemplateKey.CORR_OPTIONS: '\n'.join(corr_options),
        fill.TemplateKey.JOB_KEY: job_key,
        fill.TemplateKey.JOB_OPTIONS: '\n'.join(job_options),
        fill.TemplateKey.GEN_LINES: gen_lines_1,
        fill.TemplateKey.SADDLE: saddle,
        fill.TemplateKey.NUMERICAL: numerical
    }

    return build_mako_str(
        template_file_name='all.mako',
        template_src_path=TEMPLATE_DIR,
        template_keys=fill_dct,
        remove_whitespace=False)


# Helper functions for CFOUR
def _set_cc_prog(method, reference):
    """ Set the appropriate CC_PROG keyword based on method requested
    """

    if method in ('ccsd', 'ccsd(t)'):
        corr_options = (('ABCDTYPE=AOBASIS'),)
    if method in ('ccsd', 'ccsd(t)') and reference in ('rhf', 'uhf'):
        corr_options += (('CC_PROG=ECC'),)
    elif method in ('ccsd', 'ccsd(t)') and reference == 'rohf':
        corr_options += (('CC_PROG=VCC'),)
    elif method in ('ccsdt', 'ccsdt(q)') and reference == 'rhf':
        corr_options += (('CC_PROG=NCC'),)
    elif method in ('ccsdt', 'ccsdt(q)') and reference in ('uhf', 'rohf'):
        raise NotImplementedError("CFOUR ONLY ALLOWS CLOSED-SHELL")
    else:
        corr_options = ()

    return corr_options


# def core_writer():
#     """
#     """
#     keyword in input is FROZEN_CORE=ON which is the default
#     (i.e. always put this in
#     unless the all_electron option is given
