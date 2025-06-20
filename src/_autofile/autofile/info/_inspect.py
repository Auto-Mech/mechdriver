""" handles argument specification
"""

from inspect import getfullargspec as function_argspec


def function_keys(function):
    """ returns the keys to a function's arguments
    """
    argspec = function_argspec(function)
    return frozenset(argspec.args)


__all__ = ['function_keys']
