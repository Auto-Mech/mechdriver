""" read/write quantum chemistry data in standard formats
"""
from autofile.file import name
from autofile.file import write
from autofile.file import read
from autofile.file._util import read_file
from autofile.file._util import write_file

__all__ = [
    'name',
    'write',
    'read',
    'write_file',
    'read_file',
]
