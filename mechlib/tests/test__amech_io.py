"""
Tests amech_io functions
"""

from mechlib import amech_io

DRIVERS = ['amech', 'inp', 'ktp', 'thermo', 'trans', 'es']


def test__print_headers():
    """ prints header and exit messages for the various drivers
    """
    for driver in DRIVERS:
        print(amech_io.printer.program_header(driver))
        print(amech_io.printer.program_exit(driver))


if __name__ == '__main__':
    test__print_headers()
