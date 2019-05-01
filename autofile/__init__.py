""" autofile: filesystem schema and interface
"""
from autofile import file
from autofile import info
from autofile import system
from autofile._fs import FS_ as fs

__all__ = [
    'file',
    'info',
    'system',
    'fs',
]
