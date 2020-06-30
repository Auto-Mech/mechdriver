"""
  Handle moving about paths
"""

import os


def starting_path():
    """ get original working directory
    """
    starting_path = os.getcwd()
    return starting_path


def go_to(path):
    """ change directory to path and return the original working directory
    """
    os.chdir(path)


def prepare_path(path, loc):
    """ change directory to starting path, return chemkin path
    """
    new_path = os.path.join(path, loc)
    return new_path
