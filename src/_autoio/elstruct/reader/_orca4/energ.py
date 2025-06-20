""" electronic energy readers
"""

import autoparse.pattern as app
import autoread as ar
import pyparsing as pp
from pyparsing import pyparsing_common as ppc

from elstruct.par import Method, Program, program_methods

PROG = Program.ORCA4


def _scf_energy(output_str):
    """Reads the SCF energy from the output file string.
    Returns the energy in Hartrees.

    :param output_str: string of the program's output file
    :type output_str: str
    :rtype: float
    """
    field = pp.Literal("Total Energy") + ":"
    # The "not followed by skipto(keyword)" makes sure we grab the last value, in case
    # this is an optimization
    parser = (
        pp.SkipTo(field, include=True)
        + ppc.number("ene")
        + ~pp.FollowedBy(pp.SkipTo(field))
    )
    return parser.parseString(output_str).get("ene")


def _one_electron_energy(output_str):
    """If there is only one electron, this returns the SCF energy"""
    parser = (
        ...
        + pp.Literal("Number of Electrons")
        + pp.Literal("NEL")
        + pp.Literal("....")
        + ppc.integer("nelec")
    )
    nelec = parser.parseString(output_str).get("nelec")
    if nelec == 1:
        return _scf_energy(output_str)
    return None


def _mp2_energy(output_str):
    """Reads the MP2 energy from the output file string.
    Returns the energy in Hartrees.

    :param output_str: string of the program's output file
    :type output_str: str
    :rtype: float
    """
    ene = _one_electron_energy(output_str)
    if ene is not None:
        return ene

    ene = ar.energy.read(
        output_str,
        app.one_of_these(
            [
                app.escape("MP2 TOTAL ENERGY:"),
                app.LINESPACES.join([app.escape("Initial E(tot)"), "..."]),
            ]
        ),
    )

    return ene


def _ccsd_energy(output_str):
    """Reads the CCSD energy from the output file string.
    Returns the energy in Hartrees.

    :param output_str: string of the program's output file
    :type output_str: str
    :rtype: float
    """
    ene = _one_electron_energy(output_str)
    if ene is not None:
        return ene

    ene = ar.energy.read(
        output_str, app.LINESPACES.join([app.escape("E(CCSD)"), "..."])
    )

    return ene


def _ccsd_t_energy(output_str):
    """Reads the CCSD(T) energy from the output file string.
    Returns the energy in Hartrees.

    :param output_str: string of the program's output file
    :type output_str: str
    :rtype: float
    """
    ene = _one_electron_energy(output_str)
    if ene is not None:
        return ene

    ene = ar.energy.read(
        output_str, app.LINESPACES.join([app.escape("E(CCSD(T))"), "..."])
    )

    return ene


def _ccsd_t_f12_ri_energy(output_str):
    """Reads the CCSD(T) energy from the output file string.
    Returns the energy in Hartrees.

    :param output_str: string of the program's output file
    :type output_str: str
    :rtype: float
    """
    ene = _one_electron_energy(output_str)
    if ene is not None:
        return ene

    label = pp.Literal("F12-ECCSD(T) with (T) scaled through CCSD")
    separator = pp.Literal("...")
    energy = ppc.number
    parser = ... + label + separator + energy("energy")
    ene = parser.parseString(output_str).get("energy")
    return ene


# A dictionary of functions for reading the energy from the output, by method
ENERGY_READER_DCT = {
    (Method.HF[0], frozenset({})): _scf_energy,
    (Method.Corr.MP2[0], frozenset({})): _mp2_energy,
    (Method.Corr.CCSD[0], frozenset({})): _ccsd_energy,
    (Method.Corr.CCSD_T[0], frozenset({})): _ccsd_t_energy,
    (Method.Corr.CCSD_T_F12_RI[0], frozenset({})): _ccsd_t_f12_ri_energy,
}

# Add DFT methods to the reader dictionary
METHODS = program_methods(PROG)
for METHOD in METHODS:
    if Method.is_standard_dft(METHOD):
        ENERGY_READER_DCT[(METHOD, frozenset({}))] = _scf_energy

# Add DFT methods to the reader dictionary
READ_METHODS = set(method[0] for method in ENERGY_READER_DCT)


def method_list():
    """Constructs a list of available electronic structure methods."""
    return tuple(sorted(ENERGY_READER_DCT.keys()))


def energy(method, output_str):
    """get total energy from output"""

    # Parse the method and lists
    core_method, pfxs = Method.evaluate_method_type(method)
    full_method = (core_method, frozenset(pfxs))
    assert full_method in method_list()

    # Get the appropriate reader and call it
    if Method.is_nonstandard_dft(method):
        energy_reader = _scf_energy
    else:
        energy_reader = ENERGY_READER_DCT[full_method]

    return energy_reader(output_str)
