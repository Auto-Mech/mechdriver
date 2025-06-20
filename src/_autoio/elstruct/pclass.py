""" functions for working with parameter classes
"""

import inspect
import itertools


def values(cls):
    """ List the values of a parameter class.

        :param cls: class object
        :type cls: obj
    """
    assert inspect.isclass(cls)
    vals = tuple(val for val in _public_attributes(cls)
                 if not inspect.isclass(val))
    return vals


def all_values(cls):
    """ Recursively list the values of a parameter class tree.

        :param cls: class object
        :type cls: obj
    """
    assert inspect.isclass(cls)
    vals = tuple(itertools.chain(*(
        [val] if not inspect.isclass(val) else all_values(val)
        for val in _public_attributes(cls))))
    return vals


def _public_attributes(cls):
    """ Recursively list the public values of a parameter class tree.

        :param cls: class object
        :type cls: obj
    """
    return tuple(val for name, val in
                 inspect.getmembers(cls, lambda x: not inspect.isroutine(x))
                 if not name.startswith('_') and not inspect.isfunction(val))
