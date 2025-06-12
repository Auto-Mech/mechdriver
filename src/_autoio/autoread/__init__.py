""" Flexibly read data from a string
"""

from autoread import energy
from autoread import geom
from autoread import matrix
from autoread._zmat import zmat
from autoread._zmat import vmat
from autoread._zmat import setval

__all__ = [
    'energy',
    'geom',
    'matrix',
    'zmat',
    'vmat',
    'setval',
]
