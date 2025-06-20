""" autoparse.find
**************

Extract information from a file using re patterns.
"""
import re
from functools import partial
import numpy as np
from autoparse._lib import STRING_START as _STRING_START
from autoparse._lib import STRING_END as _STRING_END
from autoparse._lib import LINE_START as _LINE_START
from autoparse._lib import NEWLINE as _NEWLINE
from autoparse._lib import SPACES as _SPACES
from autoparse._lib import LINESPACES as _LINESPACES
from autoparse._lib import NUMBER as _NUMBER
from autoparse._pattern import maybe as _maybe


def has_match(pattern, string, case=True):
    """ does this string have a pattern match?

    :param pattern: pattern to search for
    :type pattern: str
    :param string: string to search
    :type string: str
    :param case: if capitalization matters
    :type case: bool
    :return: is there a match
    :rtype: bool
    """
    match = _re_search(pattern, string, case=case)
    return match is not None


def full_match(pattern, string, case=True):
    """ does this pattern match this *entire* string?

    :param pattern: pattern to search for
    :type pattern: str
    :param string: string to search
    :type string: str
    :param case: if capitalization matters
    :type case: bool
    :return: does it fully match
    :rtype: bool
    """
    pattern_ = _STRING_START + pattern + _STRING_END
    return has_match(pattern_, string, case=case)


def starts_with(pattern, string, case=True):
    """ does the string start with this pattern

    :param pattern: pattern to search for
    :type pattern: str
    :param string: string to search
    :type string: str
    :param case: if capitalization matters
    :type case: bool
    :return: does it start with the pattern
    :rtype: bool
    """
    start_pattern = _STRING_START + pattern
    return has_match(start_pattern, string, case=case)


def ends_with(pattern, string, case=True):
    """ does the string end with this pattern

    :param pattern: pattern to search for
    :type pattern: str
    :param string: string to search
    :type string: str
    :param case: if capitalization matters
    :type case: bool
    :return: does it end with the pattern
    :rtype: bool
    """
    end_pattern = pattern + _STRING_END
    return has_match(end_pattern, string, case=case)


def matcher(pattern, case=True):
    """ return a boolean matching function
    """
    return partial(has_match, pattern, case=case)


def all_captures(pattern, string, case=True):
    """ capture(s) for all matches of a capturing pattern

    :param pattern: pattern to search for
    :type pattern: str
    :param string: string to search
    :type string: str
    :param case: if capitalization matters
    :type case: bool
    :return: all the instances of this pattern
    :rtype: tuple
    """
    caps = _re_findall(pattern, string, case=case)
    if caps is not None:
        cap_lst = tuple(caps)
    else:
        cap_lst = None
    return cap_lst


def all_captures_with_spans(pattern, string, case=True):
    """ capture(s) for all matches of a capturing pattern, with spans showing
        the start and end of the match

    :param pattern: pattern to search for
    :type pattern: str
    :param string: string to search
    :type string: str
    :param case: if capitalization matters
    :type case: bool
    :return: all the instances of this pattern, with spans
    :rtype: tuple
    """
    lst = []
    for match in _re_finditer(pattern, string, case=case):
        cap = match.groups()

        if len(cap) == 0:
            cap = None
        elif len(cap) == 1:
            cap = cap[0]

        span = match.span()
        lst.append((cap, span))

    lst = tuple(lst)
    return lst


def first_capture(pattern, string, case=True):
    """ capture(s) from first match for a capturing pattern

    :param pattern: pattern to search for
    :type pattern: str
    :param string: string to search
    :type string: str
    :param case: if capitalization matters
    :type case: bool
    :return: first instance of this pattern
    :rtype: str
    """
    match = _re_search(pattern, string, case=case)
    return (match.group(1) if match and len(match.groups()) == 1 else
            match.groups() if match else None)


def last_capture(pattern, string, case=True):
    """ capture(s) from first match for a capturing pattern

    :param pattern: pattern to search for
    :type pattern: str
    :param string: string to search
    :type string: str
    :param case: if capitalization matters
    :type case: bool
    :return: last instance of this pattern
    :rtype: str
    """
    caps_lst = all_captures(pattern, string, case=case)
    return caps_lst[-1] if caps_lst else None


def first_named_capture(pattern, string, case=True):
    """ capture dictionary from first match for a pattern with named captures
    """
    match = _re_search(pattern, string, case=case)
    return match.groupdict() if match and match.groupdict() else None


def split(pattern, string, case=True):
    """ split string at matches
    """
    return tuple(_re_split(pattern, string, case=case))


def split_words(string):
    """ split string at whitespaces
    """
    return split(_SPACES, strip_spaces(string))


def split_lines(string):
    """ split string at newlines
    """
    return split(_NEWLINE, string)


def remove(pattern, string, case=True):
    """ remove pattern matches
    """
    return replace(pattern, '', string, case=case)


def remove_empty_lines(string):
    """ remove empty lines from a string
    """
    pattern = _LINE_START + _maybe(_LINESPACES) + _NEWLINE
    return remove(pattern, string)


def strip_spaces(string):
    """ strip spaces from the string ends
    """
    lspaces = _STRING_START + _SPACES
    rspaces = _SPACES + _STRING_END
    rstrip_string = remove(rspaces, string)
    return remove(lspaces, rstrip_string)


def replace(pattern, repl, string, case=True):
    """ replace pattern matches
    """
    return _re_sub(pattern, repl, string, case=case)


# data type checkers
def is_number(string):
    """ does this string encode a (real) number?
    """
    return full_match(_NUMBER, strip_spaces(string))


# advanced finders
def first_matching_pattern(patterns, string, case=True):
    """ from a series of patterns, return the first one matching the string
    """
    _has_match = partial(has_match, string=string,
                         case=case)
    pattern = next(filter(_has_match, patterns), None)
    return pattern


def first_matching_pattern_all_captures(patterns, string, case=True):
    """ all captures from the first matching pattern
    """
    pattern = first_matching_pattern(patterns, string,
                                     case=case)
    return all_captures(pattern, string, case=case)


def first_matching_pattern_first_capture(patterns, string,
                                         case=True):
    """ first capture from the first matching pattern
    """
    pattern = first_matching_pattern(patterns, string,
                                     case=case)
    return first_capture(pattern, string, case=case)


def first_matching_pattern_last_capture(patterns, string, case=True):
    """ last capture from the first matching pattern
    """
    pattern = first_matching_pattern(patterns, string,
                                     case=case)
    return last_capture(pattern, string, case=case)


def _re_search(pattern, string, case=True):
    flags = _re_flags(case=case)
    return re.search(pattern, string, flags=flags)


def _re_findall(pattern, string, case=True):
    if pattern and string is not None:
        flags = _re_flags(case=case)
        ptt = re.findall(pattern, string, flags=flags)
        if ptt:
            ret = ptt
        else:
            ret = None
    else:
        ret = None
    return ret


def _re_finditer(pattern, string, case=True):
    if pattern and string is not None:
        flags = _re_flags(case=case)
        match_iter = re.finditer(pattern, string, flags=flags)
    else:
        match_iter = iter([])
    return match_iter


def _re_split(pattern, string, case=True):
    flags = _re_flags(case=case)
    return re.split(pattern, string, maxsplit=0, flags=flags)


def _re_sub(pattern, repl, string, case=True):
    flags = _re_flags(case=case)
    return re.sub(pattern, repl, string, count=0, flags=flags)


def _re_flags(case=True):
    flags = re.MULTILINE
    if not case:
        flags |= re.IGNORECASE
    return flags


def where_in(word, lines):
    """ Finds where word is in lines and returns array
        For multiple words: all words must be find in lines
        :param word: word/s to look for
        :type word: str/list for multiple words
        :param lines: lines to scan
        :type lines: list(str)
        :return where_array: array with the indices
        :rtype: numpy array
    """
    if isinstance(word, str):
        word = [word]

    _arr = np.array(
        [all(word_i in line for word_i in word) for line in lines], dtype=int)

    return np.where(_arr == 1)[0]


def where_in_any(word, lines):
    """ Finds where word is in lines and returns array
        For multiple words: any of the listed words may be found in line
        :param word: word/s to look for
        :type word: str/list for multiple words
        :param lines: list to scan
        :type lines: list(str)
        :return where_array: array with the indices
        :rtype: numpy array
    """
    if isinstance(word, str):
        word = [word]

    _arr = np.array(
        [any(word_i in line for word_i in word) for line in lines], dtype=int)

    return np.where(_arr == 1)[0]


def where_is(word, lines):
    """ Finds where list corresponds to the required word
        :param word: word to look for
        :type word: str
        :param lines: list to scan
        :type lines: list(str)
        :return where_array: array with the indices
        :rtype: numpy array
    """

    where_array = np.where(
        np.array([line == word for line in lines], dtype=int) == 1)[0]

    return where_array
