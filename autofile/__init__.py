""" autofile: filesystem schema and interface
"""
from autofile import file
from autofile import info
from autofile import system
from autofile._fs import species_filesystem
from autofile._fs import reaction_filesystem

SFS = species_filesystem()
RFS = reaction_filesystem()

__all__ = [
    'file',
    'info',
    'system',
    'SFS',
    'RFS',
]
