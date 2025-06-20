""" re pattern generators that use _lib constants
"""


from autoparse._pattern import capturing as _capturing
from autoparse._pattern import one_or_more as _one_or_more
from autoparse._pattern import zero_or_more as _zero_or_more
from autoparse._lib import WILDCARD as _WILDCARD
from autoparse._lib import LINESPACE as _LINESPACE


def block_pattern(begin_pattern, end_pattern):
    """ a pattern that will grab all of the block of text
        that exists between the specified begin and end
        patterns
    """
    return (
        begin_pattern +
        _capturing(_one_or_more(_WILDCARD, greedy=False)) +
        end_pattern
    )


def lpadded(pattern, fill_pattern=_LINESPACE):
    """ a pattern allowing optional fill patterns to the left
    """
    return _zero_or_more(fill_pattern) + pattern


def rpadded(pattern, fill_pattern=_LINESPACE):
    """ a pattern allowing optional fill patterns to the right
    """
    return pattern + _zero_or_more(fill_pattern)


def padded(pattern, fill_pattern=_LINESPACE):
    """ a pattern allowing optional fill patterns to the right or left
    """
    return _zero_or_more(fill_pattern) + pattern + _zero_or_more(fill_pattern)
