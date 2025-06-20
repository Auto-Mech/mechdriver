"""
Calculates k(T, P) and using various functional forms
"""

from ratefit.calc._rates import single_arrhenius
from ratefit.calc._rates import double_arrhenius
from ratefit.calc._rates import arrhenius
from ratefit.calc._rates import lowp_limit
from ratefit.calc._rates import lindemann
from ratefit.calc._rates import troe
from ratefit.calc._rates import plog
from ratefit.calc._rates import cheb
from ratefit.calc._rates import p_to_m


__all__ = [
    'single_arrhenius',
    'double_arrhenius',
    'arrhenius',
    'lowp_limit',
    'lindemann',
    'troe',
    'plog',
    'cheb',
    'p_to_m',
]
