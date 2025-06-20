""" module for working with YAML-style information
"""
from autofile.info._info import object_
from autofile.info._info import string
from autofile.info._info import dict_
from autofile.info._info import from_string
from autofile.info._info import matches_function_signature
from autofile.info._info import Info

__all__ = [
    'object_',
    'string',
    'dict_',
    'from_string',
    'matches_function_signature',
    'Info',
]
