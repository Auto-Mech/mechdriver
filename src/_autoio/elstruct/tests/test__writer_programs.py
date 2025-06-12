""" test elstruct.writer
"""

from elstruct import writer


def test__programs():
    """ test writer.programs
    """
    assert set(writer.programs()) >= {
        'cfour2', 'gaussian09', 'gaussian03', 'gaussian16', 'molpro2015',
        'mrcc2018', 'psi4'}


def test__gradient_programs():
    """ test writer.gradient_programs
    """
    assert set(writer.gradient_programs()) >= {
        'cfour2', 'gaussian09', 'gaussian03', 'gaussian16', 'molpro2015',
        'orca4', 'psi4'}


def test__hessian_programs():
    """ test writer.hessian_programs
    """
    assert set(writer.hessian_programs()) >= {
        'cfour2', 'gaussian09', 'gaussian03', 'gaussian16', 'molpro2015',
        'mrcc2018', 'orca4', 'psi4'}


def test__vpt2_programs():
    """ test writer.vpt2_programs
    """
    assert set(writer.vpt2_programs()) >= {
        'gaussian09', 'gaussian03', 'gaussian16'}


def test__molecular_properties_programs():
    """ test writer.molecular_properties_programs
    """
    assert set(writer.molecular_properties_programs()) >= {
        'gaussian09', 'gaussian03', 'gaussian16', 'molpro2015'}


def test__irc_programs():
    """ test writer.irc_programs
    """
    assert set(writer.irc_programs()) >= {
        'gaussian09', 'gaussian03', 'gaussian16', 'psi4'}


def test__optimization_programs():
    """ test writer.optimization_programs
    """
    assert set(writer.optimization_programs()) >= {
        'cfour2', 'gaussian09', 'gaussian03', 'gaussian16', 'molpro2015',
        'mrcc2018', 'orca4', 'psi4'}


if __name__ == '__main__':
    test__programs()
