"""
  Handle moving about paths
"""

import os


def starting_path():
    """ get original working directory
    """
    return os.getcwd()


def go_to(path):
    """ change directory to path and return the original working directory
    """
    os.chdir(path)


def prepare_path(path, loc):
    """ change directory to starting path, return chemkin path
    """
    return os.path.join(path, loc)


def read_file(path, file_name=''):
    """ Open a file and read it as a string
    """
    fname = os.path.join(path, file_name)
    with open(fname, errors='ignore') as file_obj:
        file_str = file_obj.read()
    return file_str
