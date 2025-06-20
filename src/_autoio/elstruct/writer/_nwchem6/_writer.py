""" NWChem 6.0 writer module """

import os
from ioformat import build_mako_str
import elstruct.par
import elstruct.option
from elstruct.writer import fill
from elstruct.writer._nwchem6._par import OPTION_EVAL_DCT


PROG = elstruct.par.Program.NWCHEM6

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

    _ = saddle

    prog_method, prog_reference, prog_basis = fill.program_method_names(
        PROG, method, basis, mult, orb_restricted)

    geo_str, zmat_val_str, _ = fill.geometry_strings(geo, frozen_coordinates)

    scf_options = fill.evaluate_options(scf_options, OPTION_EVAL_DCT)
    _ = casscf_options
    job_options = fill.evaluate_options(job_options, OPTION_EVAL_DCT)

    # Set the gen lines blocks
    gen_lines_1, _, _ = fill.build_gen_lines(gen_lines)

    fill_dct = {
        fill.TemplateKey.JOB_KEY: job_key,
        fill.TemplateKey.COMMENT: comment,
        fill.TemplateKey.MEMORY: memory,
        fill.TemplateKey.MACHINE_OPTIONS: '\n'.join(machine_options),
        fill.TemplateKey.MOL_OPTIONS: '\n'.join(mol_options),
        fill.TemplateKey.GEOM: geo_str,
        fill.TemplateKey.ZMAT_VALS: zmat_val_str,
        fill.TemplateKey.CHARGE: charge,
        fill.TemplateKey.MULT: mult,
        fill.TemplateKey.BASIS: prog_basis,
        fill.TemplateKey.SCF_METHOD: prog_reference,
        fill.TemplateKey.SCF_OPTIONS: ','.join(scf_options),
        fill.TemplateKey.CORR_METHOD: prog_method,
        fill.TemplateKey.CORR_OPTIONS: ','.join(corr_options),
        fill.TemplateKey.JOB_OPTIONS: ';'.join(job_options),
        fill.TemplateKey.GEN_LINES: gen_lines_1,
    }

    return build_mako_str(
        template_file_name='all.mako',
        template_src_path=TEMPLATE_DIR,
        template_keys=fill_dct,
        remove_whitespace=False)
