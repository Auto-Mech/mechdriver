"""
Library of functions to interact with the filesystem
"""

import os


def prefix_filesystem(run_prefix, save_prefix):
    """ Create the AutoMech run and save directory file systems
    """
    if not os.path.exists(run_prefix):
        os.mkdir(run_prefix)
    if not os.path.exists(save_prefix):
        os.mkdir(save_prefix)
