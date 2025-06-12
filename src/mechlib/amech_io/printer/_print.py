"""
  Various status messages
"""

from mechlib.amech_io.printer._format import format_message


def message(message_label, *args, newline=None, indent=None):
    """ Print a general message to output.
    """
    if len(args) > 0:
        print(format_message(message_label, newline, indent), *args)
    else:
        print(format_message(message_label, newline, indent))


def debug_message(message_label, *args,
                  newline=None, indent=None, print_debug=True):
    """ Print a debug message to output.
    """
    if print_debug:
        _msg = format_message(message_label, newline, indent)
        if len(args) > 0:
            print('Debug: ', _msg, *args)
        else:
            print('Debug: ', _msg)


def info_message(message_label, *args, newline=None, indent=None):
    """ Print an info message to output.
    """
    if len(args) > 0:
        print(format_message(message_label, newline, indent), *args)
    else:
        print(format_message(message_label, newline, indent))


def error_message(message_label, *args, newline=None, indent=None):
    """ Print an error message to output.
    """
    if len(args) > 0:
        print('ERROR: ', format_message(message_label, newline, indent), *args)
    else:
        print('ERROR: ', format_message(message_label, newline, indent))


def warning_message(message_label, *args, newline=None, indent=None):
    """ Print a warning message to output.
    """
    _msg = format_message(message_label, newline, indent)
    if len(args) > 0:
        print('WARNING: ', _msg, *args)
    else:
        print('WARNING: ', _msg)
