"""
 Program printing messages
"""

from mechlib.amech_io.printer._print import message


def program_run_message(prog, path):
    """ message for a program
    """

    message(f'Run path for {prog}: {path}')
