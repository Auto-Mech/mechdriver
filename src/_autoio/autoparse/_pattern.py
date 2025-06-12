""" re pattern generators
"""
from re import escape as re_escape


def escape(pattern):
    """ escape special characters in pattern

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return re_escape(pattern)


def maybe(pattern):
    """ a pattern that may or may not be present

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return rf'(?:{pattern:s})?'


def preceded_by(pattern):
    """ matches if the current position is preceded by the pattern

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return rf'(?<={pattern:s})'


def not_preceded_by(pattern):
    """ matches if the current position is not preceded by the pattern

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return rf'(?<!{pattern:s})'


def followed_by(pattern):
    """ matches if the current position is followed by the pattern

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return rf'(?={pattern:s})'


def not_followed_by(pattern):
    """ matches if the current position is not followed by the pattern

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return rf'(?!{pattern:s})'


def zero_or_more(pattern, greedy=True):
    """ zero or more repeats of a pattern

    :param pattern: an `re` pattern
    :type pattern: str
    :param greedy: match as much as possible?
    :type greedy: bool

    :rtype: str
    """
    return (rf'(?:{pattern:s})*' if greedy else rf'(?:{pattern:s})*?')


def one_or_more(pattern, greedy=True):
    """ one or more repeats of a pattern

    :param pattern: an `re` pattern
    :type pattern: str
    :param greedy: match as much as possible?
    :type greedy: bool

    :rtype: str
    """
    return (rf'(?:{pattern:s})+' if greedy else rf'(?:{pattern:s})+?')


def one_of_these(patterns):
    """ any one of a series of patterns

    :param patterns: a series of `re` patterns
    :type patterns: list of strings

    :rtype: str
    """
    patterns_str = '|'.join(patterns)
    return rf'(?:{patterns_str:s})'


def capturing(pattern):
    """ generate a capturing pattern

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return rf'({pattern:s})'


def named_capturing(pattern, name):
    """ generate a named capturing pattern

    :param pattern: an `re` pattern
    :type pattern: str
    :param name: a name for the capture
    :type name: str

    :rtype: str
    """
    return rf'(?P<{name:s}>{pattern:s})'


def series(pattern, sep_pattern):
    """ repeated patters with an intervening separator pattern

    :param pattern: an `re` pattern
    :type pattern: str
    :param sep_pattern: an `re` pattern
    :type sep_pattern: str

    :rtype: str
    """
    return pattern + zero_or_more(sep_pattern + pattern)
