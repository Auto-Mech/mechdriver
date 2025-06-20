""" electronic energy readers
"""
import autoread as ar
import autoparse.pattern as app
from elstruct.par import Program, Method, program_methods


PROG = Program.MRCC2018


def _hf_energy(output_string):
    ene = ar.energy.read(
        output_string,
        app.one_of_these([
            app.escape('***FINAL HARTREE-FOCK ENERGY:'),
            app.escape('***SEMICANONICAL ROHF ENERGY:'),
        ]))
    return ene


def _mp2_energy(output_string):
    ene = ar.energy.read(
        output_string,
        app.one_of_these([
            app.escape('Total MP2 energy [au]:'),
            app.escape('MP2 energy [au]:')
        ]))
    return ene


def _ccsd_energy(output_string):
    ene = ar.energy.read(
        output_string,
        app.escape('Total CCSD energy [au]:')
        )
    return ene


def _ccsd_t_energy(output_string):
    ene = ar.energy.read(
        output_string,
        app.escape('CCSD(T) total energy [au]:')
        )
    return ene


def _ccsdt_energy(output_string):
    ene = ar.energy.read(
        output_string,
        app.escape('Total CCSDT energy [au]:')
        )
    return ene


def _ccsdt_q_energy(output_string):
    ene = ar.energy.read(
        output_string,
        app.one_of_these([
            app.escape('Total CCSDT(Q) energy [au]:'),
            app.escape('Total CCSDT(Q)/B energy [au]:'),
        ]))
    return ene


# a dictionary of functions for reading the energy from the output, by method
ENERGY_READER_DCT = {
    (Method.HF[0], frozenset({})): _hf_energy,
    (Method.Corr.MP2[0], frozenset({})): _mp2_energy,
    (Method.Corr.CCSD[0], frozenset({})): _ccsd_energy,
    (Method.Corr.CCSD_T[0], frozenset({})): _ccsd_t_energy,
    (Method.Corr.CCSDT[0], frozenset({})): _ccsdt_energy,
    (Method.Corr.CCSDT_Q[0], frozenset({})): _ccsdt_q_energy,
}
METHODS = program_methods(PROG)

# Check if we have added any unsupported methods to the energy reader
READ_METHODS = set(method[0] for method in ENERGY_READER_DCT)
assert READ_METHODS <= set(METHODS)


def method_list():
    """ list of available electronic structure methods
    """
    return tuple(sorted(ENERGY_READER_DCT.keys()))


def energy(method, output_string):
    """ get total energy from output
    """

    # Parse the method and lists
    core_method, pfxs = Method.evaluate_method_type(method)
    full_method = (core_method, frozenset(pfxs))
    assert full_method in method_list()

    # get the appropriate reader and call it
    energy_reader = ENERGY_READER_DCT[full_method]

    return energy_reader(output_string)
