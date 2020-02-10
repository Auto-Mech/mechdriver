"""
Parser functions for the theory.dat file
"""

from autoparse.find import first_capture
from autoparse.pattern import capturing
from autoparse.pattern import zero_or_more
from autoparse.pattern import one_or_more
from autoparse.pattern import NONSPACE
from autoparse.pattern import SPACE
from autoparse.pattern import WILDCARD


ELSTRUCT_REQUIRED_KEYWORDS = [
    'program',
    'method',
    'basis',
    'orb_restrict'
]
ELSTRUCT_SUPPORTED_KEYWORDS = [
    'program',
    'method',
    'basis',
    'orb_restrict',
    'memory',
    'nprocs',
    'econv',
    'gconv'
]


def read_program(input_string, level):
    """ obtain the theory level
    """
    pattern = ('program' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_level_section(input_string, level)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_method(input_string, level):
    """ obtain the method level
    """
    pattern = ('method' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_level_section(input_string, level)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_basis(input_string, level):
    """ obtain the basis level
    """
    pattern = ('basis' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_level_section(input_string, level)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_orb_restrict(input_string, level):
    """ obtain the orb_restricted level
    """
    pattern = ('orb_restrict' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_level_section(input_string, level)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_memory(input_string, level):
    """ obtain the memory level
    """
    pattern = ('memory' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_level_section(input_string, level)

    keyword = first_capture(pattern, block)

    if keyword is None:
        keyword = 0.5
    else:
        keyword = float(keyword)

    return keyword


def read_nprocs(input_string, level):
    """ obtain the memory level
    """
    pattern = ('memory' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_level_section(input_string, level)

    keyword = first_capture(pattern, block)

    if keyword is None:
        keyword = 1
    else:
        keyword = int(keyword)

    assert keyword is not None

    return keyword


def check_defined_theory_level_keywords(input_string, level):
    """ obtains the keywords defined in the input by the user
    """

    section_string = _get_level_section(input_string, level)
    defined_keywords = _get_defined_keywords(section_string)

    # Check if keywords are supported
    if not all(keyword in ELSTRUCT_SUPPORTED_KEYWORDS
               for keyword in defined_keywords):
        raise NotImplementedError

    # Check if elements of keywords
    if not all(keyword in defined_keywords
               for keyword in ELSTRUCT_REQUIRED_KEYWORDS):
        raise NotImplementedError


def _get_level_section(input_string, level):
    """ species input
    """

    pattern = ('level' + one_or_more(SPACE) + level +
               capturing(one_or_more(WILDCARD, greedy=False)) +
               'end')
    section = first_capture(pattern, input_string)

    return section


def _get_defined_keywords(section_string):
    """ gets a list of all the keywords defined in a section
    """

    defined_keys = []
    for line in section_string.splitlines():
        if '=' in line:
            tmp = line.strip().split('=')[0]
            defined_keys.append(tmp.strip())

    return defined_keys
