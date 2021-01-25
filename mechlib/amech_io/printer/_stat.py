"""
  Various status messages
"""

from mechlib.amech_io.printer import message


def running(label=None, path=None, newline=False, indent=False):
    print_str = 'Running '    
    status_message(print_str, label, path, newline, indent)


def writing(label=None, path=None, newline=False, indent=False):
    print_str = 'Writing '    
    status_message(print_str, label, path, newline, indent)


def saving(label=None, path=None, newline=False, indent=False):
    print_str = 'Writing '    
    status_message(print_str, label, path, newline, indent)


def reading(label=None, path=None, newline=False, indent=False):
    print_str = 'Reading '    
    status_message(print_str, label, path, newline, indent)


def checking(label=None, path=None, newline=False, indent=False):
    print_str = 'Checking '    
    status_message(print_str, label, path, newline, indent)


def generating(label=None, path=None, newline=False, indent=False):
    print_str = 'Generating '    
    status_message(print_str, label, path, newline, indent)


def results():
    message('Results:')


def status_message(print_str, label=None, path=None, newline=False, indent=False):
    if label:    
        print_str += label
    if path:
        print_str += ' at {}'.format(path)
    message(print_str, newline=newline, indent=indent)
