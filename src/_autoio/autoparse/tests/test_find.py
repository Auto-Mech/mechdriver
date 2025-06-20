""" test autoparse
"""

import numpy as np
import autoparse


PATTERNS = (
    autoparse.pattern.LOWERCASE_LETTER,
    autoparse.pattern.UPPERCASE_LETTER,
    autoparse.pattern.LETTER,
    autoparse.pattern.DIGIT,
)
STRING = ' ___A_a_ * & b ___c_ C 1d 2 __e_ D__'
STRING2 = (
    'A*_ 1 a z>\n'
    'B$- 2 b y~\n'
    'C%+ 3 c x,'
)

XYZ_STRING = """6
charge: 0, mult: 1
F    1.584823  -0.748487  -0.427122
C    0.619220   0.190166  -0.271639
C   -0.635731  -0.183914  -0.180364
Cl  -1.602333   0.736678  -0.026051
H    0.916321   1.229946  -0.227127
H   -0.882300  -1.224388  -0.229636
"""

ATOM_SYMBOL_PATTERN = (
    autoparse.pattern.LETTER +
    autoparse.pattern.maybe(autoparse.pattern.LETTER)
)
NUMBER_PATTERN = autoparse.pattern.FLOAT
XYZ_LINE_PATTERN = autoparse.pattern.LINESPACES.join([
    autoparse.pattern.capturing(ATOM_SYMBOL_PATTERN),
    autoparse.pattern.capturing(NUMBER_PATTERN),
    autoparse.pattern.capturing(NUMBER_PATTERN),
    autoparse.pattern.capturing(NUMBER_PATTERN),
])
BAD_XYZ_LINE_PATTERN = autoparse.pattern.LINESPACES.join([
    autoparse.pattern.capturing('BAD'),
    autoparse.pattern.capturing(NUMBER_PATTERN),
    autoparse.pattern.capturing(NUMBER_PATTERN),
    autoparse.pattern.capturing(NUMBER_PATTERN),
])


def test__first_capture():
    """ test autoparse.find.first_capture
    """
    cap = autoparse.find.first_capture('(cl)', XYZ_STRING)
    assert cap is None
    cap = autoparse.find.first_capture('(cl)', XYZ_STRING, case=False)
    assert cap == 'Cl'


def test__remove_empty_lines():
    """ test autoparse.find.remove_empty_lines
    """
    string = ' dkjsdf alsdh\n\n\n\n    asdf09j123\n\ndsaf09uasdf'
    assert autoparse.find.remove_empty_lines(string) == (
        ' dkjsdf alsdh\n    asdf09j123\ndsaf09uasdf'
    )


def test__is_number():
    """ test autoparse.find.is_number
    """
    assert autoparse.find.is_number(' 5 ') is True
    assert autoparse.find.is_number(' 1e-5 ') is True
    assert autoparse.find.is_number(' 1e- 5 ') is False
    assert autoparse.find.is_number(' .1e-200     \n \t \n ') is True


def test__split():
    """ test find.split
        test find.split_words
        test find.split_lines
    """

    pattern = '___c_'
    assert autoparse.find.split(pattern, STRING, case=True) == (
        ' ___A_a_ * & b ', ' C 1d 2 __e_ D__')
    assert autoparse.find.split_words(STRING) == (
        '___A_a_', '*', '&', 'b', '___c_', 'C', '1d', '2', '__e_', 'D__')
    assert autoparse.find.split_lines(STRING2) == (
        'A*_ 1 a z>', 'B$- 2 b y~', 'C%+ 3 c x,')


def test__simple_finders():
    """ test find.starts_with
    """

    pattern = autoparse.pattern.escape('A*_')
    assert autoparse.find.starts_with(pattern, STRING2, case=True)

    pattern = autoparse.pattern.escape('x,')
    assert autoparse.find.ends_with(pattern, STRING2, case=True)

    pattern = autoparse.pattern.escape('C%+')
    ptt_matcher = autoparse.find.matcher(pattern, case=True)
    assert ptt_matcher(STRING2)

    pattern = (
        'charge:' +
        autoparse.pattern.capturing(autoparse.pattern.NUMBER))
    assert autoparse.find.all_captures(pattern, XYZ_STRING, case=True) is None

    pattern = None
    assert autoparse.find.all_captures(pattern, XYZ_STRING, case=True) is None

    pattern = (
        'charge: ' +
        autoparse.pattern.capturing(autoparse.pattern.NUMBER))
    assert autoparse.find.all_captures(pattern, XYZ_STRING) == ('0',)

    pattern = autoparse.pattern.capturing(autoparse.pattern.LETTER)
    assert autoparse.find.all_captures_with_spans(pattern, 'abcd') == (
        ('a', (0, 1)), ('b', (1, 2)), ('c', (2, 3)), ('d', (3, 4)))


def test__advanced_finders():
    """ test find.first_matching_pattern
        test find.first_matching_pattern_all_captures
        test find.first_matching_pattern_first_capture
        test find.first_matching_pattern_last_capture
    """

    assert autoparse.find.first_matching_pattern(PATTERNS, 'A') == (
        autoparse.pattern.UPPERCASE_LETTER)
    assert autoparse.find.first_matching_pattern(PATTERNS, 'a') == (
        autoparse.pattern.LOWERCASE_LETTER)
    assert autoparse.find.first_matching_pattern(PATTERNS, '5') == (
        autoparse.pattern.DIGIT)

    patterns = list(map(autoparse.pattern.capturing, PATTERNS))
    assert (
        autoparse.find.first_matching_pattern_all_captures(patterns, STRING)
        == ('a', 'b', 'c', 'd', 'e')
    )

    patterns = list(map(autoparse.pattern.capturing, PATTERNS))
    assert (
        autoparse.find.first_matching_pattern_first_capture(patterns, STRING)
        == 'a'
    )

    patterns = list(map(autoparse.pattern.capturing, PATTERNS))
    assert (
        autoparse.find.first_matching_pattern_last_capture(patterns, STRING)
        == 'e'
    )


def test__single():
    """ test conv.single
    """
    bad_mult_pattern = autoparse.pattern.LINESPACES.join([
        autoparse.pattern.escape('multiplicity:'),
        autoparse.pattern.capturing(autoparse.pattern.UNSIGNED_INTEGER),
    ])
    mult_pattern = autoparse.pattern.LINESPACES.join([
        autoparse.pattern.escape('mult:'),
        autoparse.pattern.capturing(autoparse.pattern.UNSIGNED_INTEGER),
    ])
    cap = autoparse.find.first_capture(bad_mult_pattern, XYZ_STRING)
    val = autoparse.cast(cap)
    assert val is None

    mcap = autoparse.find.first_capture(mult_pattern, XYZ_STRING)
    mval = autoparse.cast(mcap)
    assert mval == 1


def test__singles():
    """ test conv.singles
    """
    pattern = autoparse.pattern.capturing(autoparse.pattern.FLOAT)
    caps = autoparse.find.all_captures(pattern, XYZ_STRING)
    vals = autoparse.cast(caps)
    assert vals == (
        1.584823, -0.748487, -0.427122, 0.61922, 0.190166, -0.271639,
        -0.635731, -0.183914, -0.180364, -1.602333, 0.736678, -0.026051,
        0.916321, 1.229946, -0.227127, -0.8823, -1.224388, -0.229636
    )


def test__multi():
    """ test conv.multi
    """
    mcap = autoparse.find.first_capture(XYZ_LINE_PATTERN, XYZ_STRING)
    mval = autoparse.cast(mcap)
    assert mval == ('F', 1.584823, -0.748487, -0.427122)

    mcap = autoparse.find.first_capture(BAD_XYZ_LINE_PATTERN, XYZ_STRING)
    mval = autoparse.cast(mcap)
    assert mval is None


def test__multis():
    """ test conv.multis
    """
    mcaps = autoparse.find.all_captures(XYZ_LINE_PATTERN, XYZ_STRING)
    mvals = autoparse.cast(mcaps)
    assert mvals == (('F', 1.584823, -0.748487, -0.427122),
                     ('C', 0.61922, 0.190166, -0.271639),
                     ('C', -0.635731, -0.183914, -0.180364),
                     ('Cl', -1.602333, 0.736678, -0.026051),
                     ('H', 0.916321, 1.229946, -0.227127),
                     ('H', -0.8823, -1.224388, -0.229636))


STRING_TESTWHERE = [
    'Species CH3C2CH3',
    '      RRHO',
    '        Geometry[angstrom]      10',
    '        C         0.000000        0.000000        0.000000',
    '        C         0.000000        0.000000        1.486233',
    '        C         1.282262        0.000000        2.237693',
    '        H        -0.917155       -0.221794        2.011789',
    '        H        -0.994756        0.165640       -0.409314',
    '        H         0.665488        0.771240       -0.396783',
    '        H         0.363930       -0.953867       -0.405851',
    '        H         1.132439        0.165640        3.302884',
    '        H         1.961070        0.771240        1.864155',
    '        H         1.816422       -0.953867        2.128913',
    '        Core    RigidRotor',
    '          SymmetryFactor        2',
    '        End',
    '        Rotor     Hindered',
    '          Group                  5 6 7          ',
    '          Axis                   2 1            ',
    '          Symmetry               3              ',
    '          Potential[kcal/mol]    6              ',
    '                0.00  0.06  0.18  0.47  0.15  0.04',
    '          End',
    '        Rotor     Hindered',
    '          Group                  8 9 10         ',
    '          Axis                   2 3            ',
    '          Symmetry               3              ',
    '          Potential[kcal/mol]    6              ',
    '                0.00  0.06  0.18  0.47  0.15  0.04',
    '          End',
    '        Frequencies[1/cm]       22',
    '         344.5',
    '         390.8   895.9   945.7',
    '         954.8  1037.6  1155.3',
    '        1191.5  1378.3  1422.2',
    '        1422.9  1482.7  1490.8',
    '        1493.7  1504.5  2970.6',
    '        2973.5  3043.5  3043.6',
    '        3116.5  3117.6  3200.6',
    '        ZeroEnergy[kcal/mol]    0.0',
    '        ElectronicLevels[1/cm]  1',
    '            0   2',
    '      End'
]


def test__where_is():
    """ test find.where_is
    """
    line = ('                0.00  0.06  0.18  0.47  0.15  0.04')
    assert (autoparse.find.where_is(line, STRING_TESTWHERE)
            == np.array([21, 28])).all()
    line_single = '        ElectronicLevels[1/cm]  1'
    assert autoparse.find.where_is(line_single, STRING_TESTWHERE)[0] == 40


def test__where_in():
    """ test find.where_in
    """
    line = 'End'
    assert (autoparse.find.where_in(line, STRING_TESTWHERE)
            == np.array([15, 22, 29, 42])).all()
    nonlines = ['ciao']
    assert len(autoparse.find.where_in(nonlines, STRING_TESTWHERE)) == 0


def test__where_in_any():
    """ test find.where_in_any
    """
    lines = ['End', 'RRHO']
    assert (autoparse.find.where_in_any(lines, STRING_TESTWHERE)
            == np.array([1, 15, 22, 29, 42])).all()
