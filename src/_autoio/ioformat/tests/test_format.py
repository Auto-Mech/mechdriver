"""
    test ioformat._format functions
"""

import os
import ioformat


MAKO_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), 'data')
INI_STR = (
    '! From Theory\n'
    'A + B = C + D    1.00 0.00 0.00  \n'
    '  PLOG / 0.01  3.0E+8 2.0 100.0\n'
    '\n'
    'A + D = E        2.00 0.00 0.00  \n'
    '\n'
    '! From Experiment\n'
    'E = F            3.00 0.00 0.00  \n'
    '  PLOG / 0.01  6.0 2.5 105\n'
    '\n'
    '\n'
)
INI_STR2 = (
    'ChemicalEigenvalueMax                 0.2\n'
    'ChemicalEigenvalueMin                 1.e-6\n'
    'CalculationMethod                     direct\n'
    'EigenvalueOutput			  eigenvalue.out\n'
    'Model\n'
    '  EnergyRelaxation\n'
    '    Exponential\n'
    '      Factor[1/cm]                   234'
)


def test__string_build():
    """ test ioformat.build_mako_str
        test ioformat.indent
        test ioformat.addchar
    """

    mako_keys = {'param1': 'molecule', 'param2': 'atom', 'param3': 3}
    mako_str = ioformat.build_mako_str(
        'test.mako', MAKO_PATH, mako_keys)
    assert mako_str == (
        'param1 is molecule\n'
        'param2 is atom\n'
        'param3 is 3\n'
        'param3 is 3\n'
        'param3 is 3\n'
    )

    ini_string = 'molecule'
    assert ioformat.indent(ini_string, 4) == '    molecule'
    assert ioformat.addchar(ini_string, '- ', side='pre') == '- molecule'
    assert ioformat.addchar(ini_string, ' +++', side='post') == 'molecule +++'


def test__string_alter():
    """ test ioformat.headlined_sections
    """

    head_secs = ioformat.headlined_sections(INI_STR, '=')
    assert head_secs == [
        'A + B = C + D    1.00 0.00 0.00  \n  PLOG / 0.01  3.0E+8 2.0 100.0\n',
        'A + D = E        2.00 0.00 0.00  \n\n! From Experiment',
        'E = F            3.00 0.00 0.00  \n  PLOG / 0.01  6.0 2.5 105\n\n'
    ]


def test__string_remove():
    """ test ioformat.remove_whitespace
        test ioformat.remove_trail_whitespace
        test ioformat.remove_comment_lines
    """

    out_str = ioformat.remove_whitespace_from_string(INI_STR)
    assert out_str == (
        '! From Theory\n'
        'A + B = C + D    1.00 0.00 0.00\n'
        'PLOG / 0.01  3.0E+8 2.0 100.0\n'
        'A + D = E        2.00 0.00 0.00\n'
        '! From Experiment\n'
        'E = F            3.00 0.00 0.00\n'
        'PLOG / 0.01  6.0 2.5 105\n'
    )

    out_str = ioformat.remove_trail_whitespace(INI_STR)
    assert out_str == (
        '! From Theory\n'
        'A + B = C + D    1.00 0.00 0.00\n'
        '  PLOG / 0.01  3.0E+8 2.0 100.0\n'
        'A + D = E        2.00 0.00 0.00\n'
        '! From Experiment\n'
        'E = F            3.00 0.00 0.00\n'
        '  PLOG / 0.01  6.0 2.5 105\n'
    )

    out_str = ioformat.remove_comment_lines(INI_STR, '!')
    assert out_str == (
        '\n'
        'A + B = C + D    1.00 0.00 0.00  \n'
        '  PLOG / 0.01  3.0E+8 2.0 100.0\n'
        '\n'
        'A + D = E        2.00 0.00 0.00  \n'
        '\n'
        '\n'
        'E = F            3.00 0.00 0.00  \n'
        '  PLOG / 0.01  6.0 2.5 105\n'
        '\n'
        '\n'
    )


def test__string_add():
    """ test ioformat.add_line
    """

    out_str = ioformat.add_line(
        INI_STR2, 'WellExtension', 'Model', 'before')
    assert out_str == (
        'ChemicalEigenvalueMax                 0.2\n'
        'ChemicalEigenvalueMin                 1.e-6\n'
        'CalculationMethod                     direct\n'
        'EigenvalueOutput			  eigenvalue.out\n'
        'WellExtension\n'
        'Model\n'
        '  EnergyRelaxation\n'
        '    Exponential\n'
        '      Factor[1/cm]                   234'
    )

    out_str = ioformat.add_line(
        INI_STR2, '  LumpingScheme  W1+W2', 'Model', 'after')
    assert out_str == (
        'ChemicalEigenvalueMax                 0.2\n'
        'ChemicalEigenvalueMin                 1.e-6\n'
        'CalculationMethod                     direct\n'
        'EigenvalueOutput			  eigenvalue.out\n'
        'Model\n'
        '  LumpingScheme  W1+W2\n'
        '  EnergyRelaxation\n'
        '    Exponential\n'
        '      Factor[1/cm]                   234'
    )


def test__string_change():
    """ test ioformat.change_line
    """

    out_str = ioformat.change_line(
        INI_STR2,
        '      Factor[1/cm]                   550',
        '      Factor[1/cm]                   234')
    assert out_str == (
        'ChemicalEigenvalueMax                 0.2\n'
        'ChemicalEigenvalueMin                 1.e-6\n'
        'CalculationMethod                     direct\n'
        'EigenvalueOutput			  eigenvalue.out\n'
        'Model\n'
        '  EnergyRelaxation\n'
        '    Exponential\n'
        '      Factor[1/cm]                   550'
    )
