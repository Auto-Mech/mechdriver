""" electronic energy readers
"""

import autoread as ar
import autoparse.pattern as app
from elstruct.par import Program, Method, program_methods


PROG = Program.GAUSSIAN09

DOUB_HYB_DFT = [
    'b2plypd3'
]


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
            app.escape('E(RHF) ='),
            app.escape('E(UHF) ='),
            app.escape('E(ROHF) =')]))

    return ene


def _mp2_energy(output_str):
    """ Reads the MP2 energy from the output file string.
        Returns the energy in Hartrees.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """
    return ar.energy.read(
        output_str,
        app.escape('EUMP2 = '),
        val_ptt=app.EXPONENTIAL_FLOAT_D)


def _dft_energy(output_str):
    """ Reads the energy from most density functional theory methods
        from the output file string. Returns the energy in Hartrees.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    e_pattern = app.escape('E(') + r'[^\s]+' + app.escape(')')
    ene = ar.energy.read(
        output_str,
        start_ptt=app.LINESPACES.join([
            'SCF Done:', e_pattern, '=']))

    return ene


def _doub_hyb_dft_energy(output_str):
    """ Read the energy from double-hybdrid density functional theory methods
        from the output file string. Return the energy in Hartrees.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    e_pattern = (
        app.escape('E') + app.maybe('2') + app.escape('(') +
        app.one_of_these([dft.upper() for dft in DOUB_HYB_DFT]) +
        app.escape(')')
    )
    dft_pattern = (
        e_pattern + app.SPACES + '=' + app.SPACES +
        app.EXPONENTIAL_FLOAT_D + app.SPACES +
        e_pattern + app.SPACES + '='
    )

    ene = ar.energy.read(
        output_str,
        start_ptt=dft_pattern
        )

    return ene


# A dictionary of functions for reading the energy from the output, by method
ENERGY_READER_DCT = {
    (Method.HF[0], frozenset({})): _hf_energy,
    (Method.Corr.MP2[0], frozenset({})): _mp2_energy,
}
METHODS = program_methods(PROG)

# Add DFT methods to the reader dictionary
for METHOD in METHODS:
    if Method.is_standard_dft(METHOD):
        if METHOD not in DOUB_HYB_DFT:
            ENERGY_READER_DCT[(METHOD, frozenset({}))] = _dft_energy
        else:
            ENERGY_READER_DCT[(METHOD, frozenset({}))] = _doub_hyb_dft_energy
    elif Method.is_semi_empirical(METHOD):
        ENERGY_READER_DCT[(METHOD, frozenset({}))] = _dft_energy

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
    if Method.is_nonstandard_dft(core_method):
        if core_method not in DOUB_HYB_DFT:
            energy_reader = _dft_energy
        else:
            energy_reader = _doub_hyb_dft_energy
    else:
        energy_reader = ENERGY_READER_DCT[full_method]

    return energy_reader(output_str)
