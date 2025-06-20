"""
ES routines prints
"""

from mechlib.amech_io.printer._print import error_message


def missing_input(prop):
    """ Print a message saying that some required information is
        missing from an input file.
    """
    if prop is not None:
        error_message(f'The {prop} needs to be specified in input')
