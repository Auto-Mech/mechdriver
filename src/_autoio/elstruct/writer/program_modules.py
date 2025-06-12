""" helpers for importing and managing program modules
"""

import importlib
from elstruct import par
from elstruct import pclass


# Functions to import and call the appropriate writer function
def call_module_function(prog, function, *args, **kwargs):
    """ call the module implementation of a given function

        :param prog: the program
        :type prog: str
        :param function_template: a function with the desired signature
        :type function_template: function
    """

    def _rename_prog(prog):
        """ Rename a program if number does not match module name """
        if prog in ('molpro2021', 'molpro2021_mppx'):
            prog = 'molpro2015'
        elif prog in ('gaussian03'):
            prog = 'gaussian09'
        return prog

    new_name = _rename_prog(prog)
    assert new_name in pclass.values(par.Program)
    assert new_name in program_modules_with_function(function)

    name = f'_{_rename_prog(prog)}'
    module = importlib.import_module(f'elstruct.writer.{name:s}')
    writer = getattr(module, 'write_input')

    return writer(function, *args, **kwargs)


def program_modules_with_function(function):
    """
        :param function: a function with the desired signature
        :type function: function
    """

    progs = []
    for prog in pclass.values(par.Program):
        if function in WRITER_MODULE_DCT[prog]:
            progs.append(prog)

    return progs


# Information on what writers have been implemented
class Job():
    """ Names of electronic structure jobs to ne written
    """
    ENERGY = 'energy'
    GRADIENT = 'gradient'
    HESSIAN = 'hessian'
    VPT2 = 'vpt2'
    IRC = 'irc'
    MOLPROP = 'molecular_properties'
    OPTIMIZATION = 'optimization'


# Dictionaries that dictate what writer/reader functionality
WRITER_MODULE_DCT = {
    par.Program.CFOUR2: (
        Job.ENERGY, Job.GRADIENT, Job.HESSIAN, Job.OPTIMIZATION),
    par.Program.GAUSSIAN09: (
        Job.ENERGY, Job.GRADIENT, Job.HESSIAN, Job.OPTIMIZATION,
        Job.MOLPROP, Job.IRC, Job.VPT2),
    par.Program.GAUSSIAN03: (
        Job.ENERGY, Job.GRADIENT, Job.HESSIAN, Job.OPTIMIZATION,
        Job.MOLPROP, Job.IRC, Job.VPT2),
    par.Program.GAUSSIAN16: (
        Job.ENERGY, Job.GRADIENT, Job.HESSIAN, Job.OPTIMIZATION,
        Job.MOLPROP, Job.IRC, Job.VPT2),
    par.Program.MOLPRO2015: (
        Job.ENERGY, Job.GRADIENT, Job.HESSIAN, Job.OPTIMIZATION,
        Job.MOLPROP, Job.IRC, Job.VPT2),
    par.Program.MOLPRO2021: (
        Job.ENERGY, Job.GRADIENT, Job.HESSIAN, Job.OPTIMIZATION,
        Job.MOLPROP, Job.IRC, Job.VPT2),
    par.Program.MRCC2018: (
        Job.ENERGY, Job.HESSIAN, Job.OPTIMIZATION),
    par.Program.NWCHEM6: (),
    # par.Program.NWCHEM6: (
    #     Job.ENERGY, Job.OPTIMIZATION),
    par.Program.ORCA4: (
        Job.ENERGY, Job.GRADIENT, Job.HESSIAN, Job.OPTIMIZATION),
    par.Program.PSI4: (
        Job.ENERGY, Job.GRADIENT, Job.HESSIAN, Job.OPTIMIZATION,
        Job.MOLPROP, Job.IRC),
    par.Program.QCHEM5: (
        Job.ENERGY, Job.GRADIENT, Job.HESSIAN, Job.OPTIMIZATION
    )
}
