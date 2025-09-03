"""orca4 writer module"""

import os
import re

import automol
import pyparsing as pp
from ioformat import build_mako_str
from pyparsing import pyparsing_common as ppc

import elstruct.option
import elstruct.par
from elstruct.writer import fill

PROG = elstruct.par.Program.ORCA4

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_DIR = os.path.join(THIS_DIR, "templates")


def write_input(
    job_key,
    geo,
    charge,
    mult,
    method,
    basis,
    orb_restricted,
    # molecule options
    mol_options=(),
    # machine options
    memory=1,
    comment="",
    machine_options=(),
    # theory options
    scf_options=(),
    casscf_options=(),
    corr_options=(),
    # generic options
    gen_lines=None,
    # job options
    job_options=(),
    frozen_coordinates=(),
    saddle=False,
):
    """Write an input file string for an electronic structure calculation
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

    # Set the theoretical method
    prog_method, prog_reference, prog_basis = fill.program_method_names(
        PROG, method, basis, mult, orb_restricted
    )

    has_corr = True
    if automol.geom.electron_count(geo) <= 1:
        has_corr = False

    prog_bases = [prog_basis]
    # Add CABS basis for F12 methods
    if "f12" in prog_method.lower() and has_corr:
        cabs_basis = f"{prog_basis}-cabs"
        prog_bases.append(cabs_basis)

    # Add RI basis for RI methods
    if "ri" in prog_method.lower() and has_corr:
        # It is recommended to increase the basis set cardinality by 1
        bump_card_ = lambda m: m.group().translate(str.maketrans("dtq5", "tq56"))
        ri_basis = prog_basis.replace("-f12", "") + "/c"
        ri_basis = re.sub("[dtq5](?=z)", bump_card_, ri_basis)
        prog_bases.append(ri_basis)

    geo_str, zmat_var_val_str, zmat_const_val_str = fill.geometry_strings(
        geo, frozen_coordinates
    )

    _ = machine_options
    _ = casscf_options
    _ = job_options
    _ = zmat_var_val_str
    _ = zmat_const_val_str
    _ = gen_lines

    memory = memory * 1000.0

    if elstruct.par.Method.is_correlated(method):
        assert not corr_options

    numerical = False
    nprocs = 1
    coord_sys = "xyz"

    pal_option = next(
        (option for option in machine_options if "pal" in option.lower()), None
    )
    if pal_option is not None:
        nprocs_expr = (
            pp.CaselessLiteral("PAL")
            + pp.Opt(pp.CaselessLiteral("NPROCS"))
            + ppc.integer("nprocs")
        )
        expr = pp.SkipTo(nprocs_expr, include=True)
        result = expr.parse_string(pal_option)
        nprocs = int(result.get("nprocs", nprocs)) # pyright: ignore[reportArgumentType]

    # Orca won't parallelize over more than the number of electron pairs
    npairs = automol.geom.electron_count(geo) // 2
    nprocs = min(nprocs, npairs)

    if saddle:
        raise NotImplementedError

    fill_dct = {
        fill.TemplateKey.NPROCS: nprocs,
        fill.TemplateKey.NUMERICAL: numerical,
        fill.TemplateKey.COORD_SYS: coord_sys,
        fill.TemplateKey.MEMORY: memory,
        fill.TemplateKey.REFERENCE: prog_reference,
        fill.TemplateKey.METHOD: prog_method if has_corr else "",
        fill.TemplateKey.BASIS: " ".join(prog_bases),
        fill.TemplateKey.SCF_OPTIONS: ",".join(scf_options),
        fill.TemplateKey.MOL_OPTIONS: ",".join(mol_options),
        fill.TemplateKey.COMMENT: comment,
        fill.TemplateKey.CHARGE: charge,
        fill.TemplateKey.MULT: mult,
        fill.TemplateKey.GEOM: geo_str,
        fill.TemplateKey.JOB_KEY: job_key,
    }

    return build_mako_str(
        template_file_name="all.mako",
        template_src_path=TEMPLATE_DIR,
        template_keys=fill_dct,
        remove_whitespace=False,
    )
