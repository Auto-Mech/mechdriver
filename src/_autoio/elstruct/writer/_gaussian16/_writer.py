""" gaussian09 writer module """

import os
from ioformat import build_mako_str
import elstruct.option
import elstruct.par
from elstruct.writer import fill
from elstruct.writer._gaussian16._par import OPTION_EVAL_DCT


PROG = elstruct.par.Program.GAUSSIAN16

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

    # Build the geometry object for the job
    geo_str, zmat_var_val_str, zmat_const_val_str = fill.geometry_strings(
        geo, frozen_coordinates)

    # Round the memory since Gaussian doesn't deal with floats
    memory = int(memory)

    # Set theory methods and options
    if elstruct.par.Method.is_correlated(method):
        assert not corr_options
    prog_method, prog_reference, prog_basis = fill.program_method_names(
        PROG, method, basis, mult, orb_restricted)

    if (prog_reference == elstruct.par.Reference.ROHF and
            job_key in (elstruct.par.Job.GRADIENT, elstruct.par.Job.HESSIAN)):
        job_options = list(job_options)
        job_options.insert(0, 'EnOnly')

    # Build various options
    scf_guess_options, scf_options = fill.intercept_scf_guess_option(
        scf_options, OPTION_EVAL_DCT)
    job_options = fill.evaluate_options(job_options, OPTION_EVAL_DCT)
    if saddle:
        job_options += ('CALCFC', 'TS', 'NOEIGEN', 'MAXCYCLES=60')

    if 'Tight' in job_options and job_key == 'optimization':
        job_key = 'tight_optimization'
    # Set the gen lines blocks
    gen_lines_1, _, _ = fill.build_gen_lines(gen_lines)

    # Build the dictionary to fill the Mako template
    fill_dct = {
        fill.TemplateKey.MEMORY: memory,
        fill.TemplateKey.MACHINE_OPTIONS: '\n'.join(machine_options),
        fill.TemplateKey.REFERENCE: prog_reference,
        fill.TemplateKey.METHOD: prog_method,
        fill.TemplateKey.BASIS: prog_basis,
        fill.TemplateKey.SCF_OPTIONS: ','.join(scf_options),
        fill.TemplateKey.SCF_GUESS_OPTIONS: ','.join(scf_guess_options),
        fill.TemplateKey.CASSCF_OPTIONS: ','.join(casscf_options),
        fill.TemplateKey.MOL_OPTIONS: ','.join(mol_options),
        fill.TemplateKey.COMMENT: comment,
        fill.TemplateKey.CHARGE: charge,
        fill.TemplateKey.MULT: mult,
        fill.TemplateKey.GEOM: geo_str,
        fill.TemplateKey.ZMAT_VAR_VALS: zmat_var_val_str,
        fill.TemplateKey.ZMAT_CONST_VALS: zmat_const_val_str,
        fill.TemplateKey.JOB_KEY: job_key,
        fill.TemplateKey.JOB_OPTIONS: ','.join(job_options),
        fill.TemplateKey.GEN_LINES: gen_lines_1,
    }

    return build_mako_str(
        template_file_name='all.mako',
        template_src_path=TEMPLATE_DIR,
        template_keys=fill_dct,
        remove_whitespace=False)
