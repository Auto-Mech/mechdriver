"""
  Various status messages
"""

from mechlib.amech_io.printer._print import message


def running(label=None, path=None, newline=False, indent=False):
    """ Print a message that some task is running.
    """
    print_str = 'Running '
    status_message(print_str, label, path, newline, indent)


def writing(label=None, path=None, newline=False, indent=False):
    """ Print a message that some data is being written to a file.
    """
    print_str = 'Writing '
    status_message(print_str, label, path, newline, indent)


def saving(label=None, path=None, newline=False, indent=False):
    """ Print a message that some data is being saved to a file.
    """
    print_str = 'Writing '
    status_message(print_str, label, path, newline, indent)


def reading(label=None, path=None, newline=False, indent=False):
    """ Print a message that a file is being read.
    """
    print_str = 'Reading '
    status_message(print_str, label, path, newline, indent)


def checking(label=None, path=None, newline=False, indent=False):
    """ Print a message that a file/object is being checked..
    """
    print_str = 'Checking '
    status_message(print_str, label, path, newline, indent)


def generating(label=None, path=None, newline=False, indent=False):
    """ Print a message that a file/object is being generated.
    """
    print_str = 'Generating '
    status_message(print_str, label, path, newline, indent)


def results():
    """ Print static results method.
    """
    message('Results:')


def status_message(print_str, label=None, path=None,
                   newline=False, indent=False):
    """ Print a general status message.
    """
    if label:
        print_str += label
    if path:
        print_str += ' at {}'.format(path)
    message(print_str, newline=newline, indent=indent)
