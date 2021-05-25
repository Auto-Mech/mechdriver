"""
Format a string
"""


def _indent(message, number):
    message = '\t'.expandtabs(int(number * 8)) + message
    return message


def _newline(message, number):
    message = '\n' * number + message
    return message


def format_message(message, newline, indent):
    """ Format a message string
    """
    if newline:
        message = _newline(message, newline)
    if indent:
        message = _indent(message, indent)
    return message
