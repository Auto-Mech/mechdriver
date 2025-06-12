""" convert string captures to various datatypes
"""

from collections.abc import Sequence as _Sequence


def cast(seq):
    """ cast each string in a nested sequence to int or float, if possible
    (recursive)
    """
    if _is_string(seq):
        ret = _cast_string(seq)
    elif _is_sequence(seq):
        ret = tuple(cast(obj) for obj in seq)
    else:
        ret = seq
    return ret


def _is_string(obj):
    return isinstance(obj, (str, bytes, bytearray))


def _is_sequence(obj):
    return isinstance(obj, _Sequence)


def _cast_string(string):
    """ cast an individual string to int or float, if possible
    """

    ret = string
    try:
        ret = int(string.replace('D+', 'E+').replace('D-', 'E-'))
    except ValueError:
        try:
            ret = float(string.replace('D+', 'E+').replace('D-', 'E-'))
        except ValueError:
            pass

    return ret
