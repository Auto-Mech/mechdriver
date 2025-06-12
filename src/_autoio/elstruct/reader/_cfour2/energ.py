""" electronic energy readers
"""

import autoread as ar
import autoparse.pattern as app
from elstruct.par import Program, Method, program_methods


PROG = Program.CFOUR2


def _hf_energy(output_str):
    """ Reads the Hartree-Fock energy from the output file string.
        Returns the energy in Hartrees.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    ene = ar.energy.read(
        output_str,
        app.one_of_these([
            app.escape('E(SCF)='),
            app.escape('E(ROHF)=')])
        )

    return ene


def _mp2_energy(output_str):
    """ Reads the MP2 energy from the output file string.
        Returns the energy in Hartrees.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    ene = ar.energy.read(
        output_str,
        app.one_of_these([
            app.escape('Total MP2 energy'),
            app.escape('MP2 energy')
        ]))

    return ene


def _ccsd_energy(output_str):
    """ Reads the CCSD energy from the output file string.
        Returns the energy in Hartrees.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    ene = ar.energy.read(
        output_str,
        app.escape('CCSD energy')
        )

    return ene


def _ccsd_t_energy(output_str):
    """ Reads the CCSD energy from the output file string.
        Returns the energy in Hartrees.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    ene = ar.energy.read(
        output_str,
        app.escape('CCSD(T) energy')
        )

    return ene


# a dictionary of functions for reading the energy from the output, by method
ENERGY_READER_DCT = {
    (Method.HF[0], frozenset({})): _hf_energy,
    (Method.Corr.MP2[0], frozenset({})): _mp2_energy,
    (Method.Corr.CCSD[0], frozenset({})): _ccsd_energy,
    (Method.Corr.CCSD_T[0], frozenset({})): _ccsd_t_energy,
}
METHODS = program_methods(PROG)

# Check if we have added any unsupported methods to the energy reader
READ_METHODS = set(method[0] for method in ENERGY_READER_DCT)
assert READ_METHODS <= set(METHODS)


def method_list():
    """ Constructs a list of available electronic structure methods.
    """
    return tuple(sorted(ENERGY_READER_DCT.keys()))


def energy(method, output_str):
    """ Reads the the total electronic energy from the output file string.
        Returns the energy in Hartrees.

        :param method: electronic structure method to read
        :type method: str
        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    # Parse the method and lists
    core_method, pfxs = Method.evaluate_method_type(method)
    full_method = (core_method, frozenset(pfxs))
    assert full_method in method_list()

    # Get the appropriate reader and call it
    energy_reader = ENERGY_READER_DCT[full_method]

    return energy_reader(output_str)
