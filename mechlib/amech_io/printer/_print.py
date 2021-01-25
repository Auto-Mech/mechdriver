"""
  Various status messages
"""

from mechlib.amech_io.printer._format import format_message


def message(message_label, *args, newline=None, indent=None):
    if len(args) > 0:
        print(format_message(message_label, newline, indent), *args)
    else:
        print(format_message(message_label, newline, indent))


def debug_message(message_label, *args, newline=None, indent=None):
    if len(args) > 0:
        print('Debug: ', format_message(message_label, newline, indent), *args)
    else:
        print('Debug: ', format_message(message_label, newline, indent))


def info_message(message_label, *args, newline=None, indent=None):
    if len(args) > 0:
        print(format_message(message_label, newline, indent), *args)
    else:
        print(format_message(message_label, newline, indent))


def error_message(message_label, *args, newline=None, indent=None):
    if len(args) > 0:
        print('**ERROR: ', format_message(message_label, newline, indent), *args)
    else:
        print('**ERROR: ', format_message(message_label, newline, indent))


def warning_message(message_label, *args, newline=None, indent=None):
    if len(args) > 0:
        print('*WARNING: ', format_message(message_label, newline, indent), *args)
    else:
        print('*WARNING: ', format_message(message_label, newline, indent))


