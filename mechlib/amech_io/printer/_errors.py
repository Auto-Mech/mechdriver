"""
ES routines prints
"""

# import logging
import automol


from mechlib.amech_io.printer import error_message


def missing_input(prop):
    """ a
    """
    if prop is not None:
        error_message('The {} needs to be specified in input'.format(prop))
