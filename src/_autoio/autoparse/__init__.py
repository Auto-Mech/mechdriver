"""
autoparse
*********
"""
from ._conv import cast
#: pattern generators
from . import pattern
#: text parsers
from . import find

__all__ = ['pattern', 'find', 'cast']
