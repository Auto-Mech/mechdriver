""" molpro2015 writer module """

import os
from ioformat import build_mako_str
import elstruct.par
import elstruct.option
from elstruct.writer import fill
from elstruct.writer._molpro2015._par import set_method_and_options
from elstruct.writer._molpro2015._par import OPTION_EVAL_DCT

PROG = elstruct.par.Program.MOLPRO2015

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

    # Set the spin
    spin = mult - 1

    # set correlated method; check if multiref
    core_method, method1, corr_options, gen_lines = set_method_and_options(
        method, corr_options, gen_lines)
    method2, prog_reference, prog_basis = fill.program_method_names(
        PROG, core_method, basis, mult, orb_restricted)

    # Set the program method to one of the following
    prog_method = method1 if method1 else method2
    
    # Reset the prog_method again to account for hf, cassf (this is dumb...)
    if prog_method not in ('hf', 'rhf', 'uhf', 'casscf'):
        corr_method = prog_method
    else:
        corr_method = ''
    # print('molpro method test since I cant write code...')
    # print(prog_reference, corr_method)

    # Set the geometry (combine opt+constrainted coordinates
    geo_str, zmat_vval_str, zmat_cval_str = fill.geometry_strings(
        geo, frozen_coordinates, offset=True)
    zmat_val_str = zmat_vval_str + '\n' + zmat_cval_str

    # Set the memory; convert from GB to MW
    memory_mw = int(memory * (1024.0 / 8.0))

    # Set the job directives and options
    # mol_options = fill.evaluate_options(mol_options, OPTION_EVAL_DCT)
    scf_options = fill.evaluate_options(scf_options, OPTION_EVAL_DCT)
    casscf_options = fill.evaluate_options(casscf_options, OPTION_EVAL_DCT)
    corr_options = fill.evaluate_options(corr_options, OPTION_EVAL_DCT)
    job_options = fill.evaluate_options(job_options, OPTION_EVAL_DCT)

    # Make no scf method if calling a multiref method
    ismultiref = elstruct.par.Method.is_multiref(method)
    if ismultiref:
        scf_options = []
        prog_reference = ''

    if saddle:
        job_options += ('root=2',)

    job_directives = [','.join(job_options)]
    if frozen_coordinates:
        job_directives.append('inactive,' + ','.join(frozen_coordinates))

    if job_key == 'hessian':
        job_directives.append('print,hessian,low=5')

    # Set the gen lines blocks
    gen_lines_1, gen_lines_2, gen_lines_3 = fill.build_gen_lines(
        gen_lines,
        line3='molpro_energy=energy\nshow[1,e25.15],molpro_energy')

    # Create a dictionary to fille the template
    fill_dct = {
        fill.TemplateKey.JOB_KEY: job_key,
        fill.TemplateKey.COMMENT: comment,
        fill.TemplateKey.MEMORY: memory_mw,
        fill.TemplateKey.MACHINE_OPTIONS: '\n'.join(machine_options),
        fill.TemplateKey.MOL_OPTIONS: '\n'.join(mol_options),
        fill.TemplateKey.GEOM: geo_str,
        fill.TemplateKey.ZMAT_VALS: zmat_val_str,
        fill.TemplateKey.CHARGE: charge,
        fill.TemplateKey.SPIN: spin,
        fill.TemplateKey.BASIS: prog_basis,
        fill.TemplateKey.SCF_METHOD: prog_reference,
        fill.TemplateKey.SCF_OPTIONS: ','.join(scf_options),
        fill.TemplateKey.ISMULTIREF: ismultiref,
        fill.TemplateKey.CASSCF_OPTIONS: '\n'.join(casscf_options),
        fill.TemplateKey.CORR_METHOD: corr_method,
        fill.TemplateKey.CORR_OPTIONS: ','.join(corr_options),
        fill.TemplateKey.JOB_OPTIONS: ';'.join(job_directives),
        fill.TemplateKey.GEN_LINES_1: gen_lines_1,
        fill.TemplateKey.GEN_LINES_2: gen_lines_2,
        fill.TemplateKey.GEN_LINES_3: gen_lines_3
    }

    return build_mako_str(
        template_file_name='all.mako',
        template_src_path=TEMPLATE_DIR,
        template_keys=fill_dct,
        remove_whitespace=False)
