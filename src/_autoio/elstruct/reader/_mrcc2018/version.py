"""
Reads the program name and version number from the output file
"""

import autoparse.pattern as app
import autoparse.find as apf


def program_name(output_string):
    """ reads the program name (here: MainName + MainVersion)
    """
    if _check_name_string(output_string) is not None:
        version_string = _get_version_string(output_string)
        year = version_string.split()[2].strip()
        prog_name = 'mrcc' + year

    return prog_name


def program_version(output_string):
    """ reads the program version number
    """
    version_string = _get_version_string(output_string)
    version_string = version_string.split()
    month = version_string[0].strip().lower()
    day = version_string[1].strip().replace(',', '')
    prog_version = month + day

    return prog_version


def _check_name_string(output_string):
    """ checks to see if the orca program string is in the output
    """

    pattern = 'MRCC program system'

    prog_string = apf.has_match(pattern, output_string)

    return prog_string


def _get_version_string(output_string):
    """ obtains the string containing the version number
    """

    pattern = ('Release date:' + app.SPACE +
               app.capturing(app.one_or_more(app.NONNEWLINE)))

    version_string = apf.first_capture(pattern, output_string)

    return version_string
