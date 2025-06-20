""" test autoparse
"""

import autoparse


STRING1 = '''
name = "Bob"
Age = 54
has_W2 = True
print(name, Age, has_W2)
Bob 54 True
1099_filed = False
SyntaxError: invalid token'''
STRING2 = 'a1y'
STRING3 = '  AbbAbbA123245'


def test__variable_name():
    """ test autoparse.pattern.VARIABLE_NAME
    """

    pattern = autoparse.pattern.capturing(autoparse.pattern.VARIABLE_NAME)
    captures = autoparse.find.all_captures(pattern, STRING1)
    assert captures == (
        'name', 'Bob', 'Age', 'has_W2', 'True', 'print', 'name', 'Age',
        'has_W2', 'Bob', 'True', '_filed', 'False', 'SyntaxError', 'invalid',
        'token'
    )


def test__position_checks():
    """ test autoparse.pattern.preceded_by
        test autoparse.pattern.not_preceded_by
        test autoparse.pattern.followed_by
        test autoparse.pattern.not_followed_by
    """

    ptt1 = autoparse.pattern.preceded_by('1') + 'a'
    assert not autoparse.find.has_match(ptt1, STRING2)

    ptt2 = autoparse.pattern.preceded_by('1') + 'y'
    assert autoparse.find.has_match(ptt2, STRING2)

    ptt3 = autoparse.pattern.not_preceded_by('1') + 'a'
    assert autoparse.find.has_match(ptt3, STRING2)

    ptt4 = autoparse.pattern.not_preceded_by('1') + 'y'
    assert not autoparse.find.has_match(ptt4, STRING2)

    ptt5 = 'a' + autoparse.pattern.followed_by('1')
    assert autoparse.find.has_match(ptt5, STRING2)

    ptt6 = 'y' + autoparse.pattern.followed_by('1')
    assert not autoparse.find.has_match(ptt6, STRING2)

    ptt7 = 'a' + autoparse.pattern.not_followed_by('1')
    assert not autoparse.find.has_match(ptt7, STRING2)

    ptt8 = 'y' + autoparse.pattern.not_followed_by('1')
    assert autoparse.find.has_match(ptt8, STRING2)


def test__padded():
    """ test autoparse.pattern.lpadded
        test autoparse.pattern.rpadded
        test autoparse.pattern.padded
    """

    ptt = 'AAA'
    lmatch_ptt = '            AAAA'
    mmatch_ptt = '            AAAA            '
    rmatch_ptt = 'AAAA            '
    assert autoparse.find.has_match(
        autoparse.pattern.lpadded(ptt),
        lmatch_ptt
    )
    assert autoparse.find.has_match(
        autoparse.pattern.padded(ptt),
        mmatch_ptt
    )
    assert autoparse.find.has_match(
        autoparse.pattern.rpadded(ptt),
        rmatch_ptt
    )


def test__block_pattern():
    """ test autoparse.pattern.block_pattern
    """

    test_str = (
        'begin'
        '   3r230refkwepfkwpekfwekf2\n'
        'qmvnboefe02-0943-102d29  \n'
        'fkwef320f0wek   vdv2\n'
        'end'
    )
    begin_ptt = 'begin'
    end_ptt = 'end'

    block_ptt = autoparse.pattern.block_pattern(begin_ptt, end_ptt)
    block = autoparse.find.first_capture(block_ptt, test_str)
    assert block == (
        '   3r230refkwepfkwpekfwekf2\n'
        'qmvnboefe02-0943-102d29  \n'
        'fkwef320f0wek   vdv2\n'
    )


def test__named_capturing():
    """ test autoparse.pattern.named_capturing
    """

    ptt = autoparse.pattern.named_capturing('Abb', 'mol')
    assert autoparse.find.first_named_capture(ptt, STRING3) == {'mol': 'Abb'}


def test__series():
    """ test autoparse.pattern.series
    """

    ptt = autoparse.pattern.capturing(
        autoparse.pattern.series('A', 'bb'))
    assert autoparse.find.first_capture(ptt, STRING3) == 'AbbAbbA'
