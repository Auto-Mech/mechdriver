""" electronic energy readers
"""

import autoread as ar
import autoparse.pattern as app
import autoparse.find as apf
from elstruct.par import Program, Method, program_methods


PROG = Program.MOLPRO2015

BASIS_PATTERN = app.one_or_more(
    app.one_of_these(
        [app.LETTER, app.NUMBER, app.escape('*'), app.escape('-'),
         app.escape('('), app.escape(')')])
    )


# Single reference wavefunction methods
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
            app.LINESPACES.join([
                app.escape('!RHF STATE'), app.FLOAT, app.escape('Energy')]),
            app.LINESPACES.join([
                app.escape('!UHF STATE'), app.FLOAT, app.escape('Energy')]),
        ]))

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
            app.escape('!MP2 total energy') + app.maybe(':'),
            app.escape('!RMP2 energy'), app.escape('!DF-RMP2 energy'), 
            app.escape('!MP2 energy'),
        ]))

    return ene


def _ccsd_energy(output_str):
    """ Reads the CCSD/UCCSD energy from the output file string.
        Returns the energy in Hartrees.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    ene = ar.energy.read(
        output_str,
        app.one_of_these([
            app.escape('!CCSD total energy') + app.maybe(':'),
            app.escape('!RHF-UCCSD energy'),
        ]))

    return ene


def _ccsd_t_energy(output_str):
    """ Reads the CCSD(T)/UCCSD(T) energy from the output file string.
        Returns the energy in Hartrees.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    ene = ar.energy.read(
        output_str,
        app.one_of_these([
            app.escape('!CCSD(T) total energy') + app.maybe(':'),
            app.escape('!RHF-UCCSD(T) energy'),
            app.LINESPACES.join([
                app.escape('!CCSD(T) STATE'),
                app.FLOAT,
                app.escape('Energy')]),
        ]))

    return ene


def _ccsdt_energy(output_str):
    """ Reads the CCSDT energy from the output file string.
        Returns the energy in Hartrees.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    ene = ar.energy.read(
        output_str,
        app.one_of_these([
            app.LINESPACES.join([
                app.escape('!CCSDT STATE'),
                app.FLOAT,
                app.escape('Energy')]),
        ]))

    return ene


def _ccsdt_q_energy(output_str):
    """ Reads the CCSDT(Q) energy from the output file string.
        Returns the energy in Hartrees.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    ene = ar.energy.read(
        output_str,
        app.one_of_these([
            app.LINESPACES.join([
                app.escape('!CCSDT(Q) STATE'),
                app.FLOAT,
                app.escape('Energy')]),
        ]))

    return ene


def _ccsd_t_f12_energy(output_str):
    """ Reads the CCSD(T)-F12/UCCSD(T)-F12 energy from the output file string.
        Currently, the function only reads the energy of the F12b variant.
        Returns the energy in Hartrees.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    ene = ar.energy.read(
        output_str,
        app.one_of_these([
            app.escape('!CCSD(T)-F12b total energy') + app.maybe(':'),
            app.escape('!RHF-UCCSD(T)-F12b energy'), #Open shell, Molpro old
            app.escape('!RHF-UCCSD(T)-F12 energy') + app.maybe(':'), #Molpro 2024
        ]))
    return ene


def _lmp2_f12_energy(output_str):
    """ Get the MP2 value

        Check if we need to grab _CORR value or first val
        Or if it is the later value given in the pattern...not clear..
    """
    ene = ar.energy.read(
        output_str, app.escape('PNO-LMP2-F12 total energy'))
    return ene


def _lccsd_f12_energy(output_str):
    """ Reads the LCCSD-F12 energy from the output file string.
        Currently, the function only reads the energy of the F12b variant.

        For local, open-shell R/L seem interchangable.

        I think it works with just PNO variant? (because of search line)

        I think it is the same for closed- and open-shell systems.

        Returns the energy in Hartrees.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    ptt = app.one_of_these([
        app.escape('!PNO-RCCSD-F12b total energy'),
        app.escape('!PNO-LCCSD-F12b total energy')])

    ene = ar.energy.read(output_str, ptt)
    return ene


def _lccsd_t_f12_energy(output_str):
    """ Reads the LCCSD(T)-F12 energy from the output file string.
        Currently, the function only reads the energy of the F12b variant.

        I think it works with the PNO variants of LCC as well as just LCC

        I think it is the same for closed- and open-shell systems.
        For local, open-shell R/L seem interchangable.

        Returns the energy in Hartrees.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """
    ptt = app.one_of_these([
        app.escape('!RCCSD(T)-F12b total energy'),
        app.escape('!LCCSD(T)-F12b total energy')])
    ene = ar.energy.read(output_str, ptt)
    return ene


# Multi reference wavefunction methods
def _casscf_energy(output_str):
    """ Reads the CASSCF energy from the output file string.
        Returns the energy in Hartrees.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    ene = ar.energy.read(
        output_str,
        app.LINESPACES.join([
            app.escape('!MCSCF STATE'), app.FLOAT, app.escape('Energy')]),
        )

    return ene


def _caspt2_energy(output_str):
    """ Reads the CASPT2 energy from the output file string.
        Returns the energy in Hartrees.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    ene = ar.energy.read(
        output_str,
        app.LINESPACES.join([
            app.escape('!RSPT2 STATE'), app.FLOAT, app.escape('Energy')]),
        )

    return ene


def _mrci_energy(output_str):
    """ Reads the MRCISD+Q energy from the output file string.
        Returns the energy in Hartrees.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: float
    """

    ene = ar.energy.read(
        output_str,
        app.LINESPACES.join([
            app.escape('!MRCI STATE'), app.FLOAT, app.escape('Energy')]),
        )

    return ene


# Relativistic methods
# def _dkh_energy(output_str):
#     """ Read the energy from a relativistic calculation using a
#         Douglas-Kroll Hamiltonian.
#     """
#     # SETTING E_DK           =     -1056.62103554  AU
# def _cowan_griffin_energy(output_str):
#     """ Read the energy from a relativistic calculation using a
#         the Cowan-Griffin operator, this includes the Darwin and
#         mass-velocity terms.
#     """
#
#     # Read the NREL and then the REL and add the two
#     # SETTING E_NREL         =     -1053.56240805  AU
#     # MASSV            =       -14.84964286 AU
#     # DARWIN           =        11.25455695 AU
#     # EREL             =        -3.59508592 AU


def _end_file_energy(output_str):
    """ Reads a user-defined electronic energy in the output string that has
        be given the variable name is MOLPRO_ENERGY.
        Returns the energy in Hartrees.
    """

    end_file_ptt = (
        'MOLPRO_ENERGY' + app.SPACES +
        app.escape('=') + app.SPACES +
        app.capturing(app.FLOAT) + app.SPACES +
        'AU'
    )
    ene = apf.last_capture(end_file_ptt, output_str)
    ene = float(ene) if ene is not None else ene

    return ene


# A dictionary of functions for reading the energy from the output, by method
ENERGY_READER_DCT = {
    (Method.HF[0], frozenset({})): _hf_energy,
    (Method.Corr.MP2[0], frozenset({})): _mp2_energy,
    (Method.Corr.CCSD[0], frozenset({})): _ccsd_energy,
    (Method.Corr.CCSD_T[0], frozenset({})): _ccsd_t_energy,
    (Method.Corr.CCSDT[0], frozenset({})): _ccsdt_energy,
    (Method.Corr.CCSDT_Q[0], frozenset({})): _ccsdt_q_energy,
    (Method.Corr.CCSD_T_F12[0], frozenset({})): _ccsd_t_f12_energy,
    (Method.Corr.MP2[0], frozenset({Method.ModPrefix.DF[0]})): _mp2_energy,
    (Method.Corr.CCSD[0], frozenset({Method.ModPrefix.ALL_ELEC[0],
                                     Method.ModPrefix.REL_DKH[0]})):
    _ccsd_energy,
    (Method.Corr.CCSD_T[0], frozenset({Method.ModPrefix.ALL_ELEC[0],
                                       Method.ModPrefix.REL_DKH[0]})):
    _ccsd_t_energy,
    (Method.Corr.CCSD_T_F12[0], frozenset({Method.ModPrefix.ALL_ELEC[0],
                                           Method.ModPrefix.REL_DKH[0]})):
    _ccsd_t_f12_energy,
    (Method.Corr.CCSD[0], frozenset({Method.ModPrefix.REL_DKH[0]})):
        _ccsd_energy,
    (Method.Corr.CCSD_T[0], frozenset({Method.ModPrefix.REL_DKH[0]})):
        _ccsd_t_energy,
    (Method.Corr.CCSD_T_F12[0], frozenset({Method.ModPrefix.REL_DKH[0]})):
        _ccsd_t_f12_energy,
    (Method.Corr.CCSD[0], frozenset({Method.ModPrefix.ALL_ELEC[0]})):
        _ccsd_energy,
    (Method.Corr.CCSD_T[0], frozenset({Method.ModPrefix.ALL_ELEC[0]})):
        _ccsd_t_energy,
    (Method.Corr.CCSD_T_F12[0], frozenset({Method.ModPrefix.ALL_ELEC[0]})):
        _ccsd_t_f12_energy,
    (Method.Corr.MP2_F12[0], frozenset({Method.ModPrefix.L_PNO[0]})):
        _lmp2_f12_energy,
    (Method.Corr.CCSD_F12[0], frozenset({Method.ModPrefix.L_PNO[0]})):
        _lccsd_f12_energy,
    (Method.Corr.CCSD_T_F12[0], frozenset({Method.ModPrefix.L_PNO[0]})):
        _lccsd_t_f12_energy,
    (Method.MultiRef.CASSCF[0], frozenset({})): _casscf_energy,
    (Method.MultiRef.CASPT2[0], frozenset({})): _caspt2_energy,
    (Method.MultiRef.CASPT2I[0], frozenset({})): _caspt2_energy,
    (Method.MultiRef.CASPT2C[0], frozenset({})): _caspt2_energy,
    (Method.MultiRef.MRCISDQ[0], frozenset({})): _mrci_energy,
}
METHODS = program_methods(PROG)
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
    core, pfxs = Method.evaluate_method_type(method)
    _method = (core, frozenset(pfxs))
    assert _method in method_list()

    # First try and grab the energy printed at the end of the file
    ene = _end_file_energy(output_str)

    # If no energy is found, get the appropriate reader and call it
    if ene is None:
        energy_reader = ENERGY_READER_DCT[_method]
        ene = energy_reader(output_str)

    return ene
