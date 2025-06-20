""" autofile: filesystem schema and interface
"""

from autofile import io_
from autofile import json_
from autofile import model
from autofile import info
from autofile import data_types
from autofile import schema
from autofile import fs
from autofile._conv import directory_to_dictionary
from autofile._safemode import turn_off_safemode
from autofile._safemode import turn_on_safemode
from autofile._safemode import safemode_is_on


__all__ = [
    'io_',
    'json_',
    'model',
    'info',
    'data_types',
    'schema',
    'fs',
    'directory_to_dictionary',
    'turn_off_safemode',
    'turn_on_safemode',
    'safemode_is_on'
]
