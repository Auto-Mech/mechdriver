"""
  Handle ckin paths
"""

import os


def ckin_path():
    """ Set path and make directories
    """
    starting_path = os.getcwd()
    path = os.path.join(starting_path, 'ckin')
    if not os.path.exists(path):
        os.makedirs(path)

    return path
