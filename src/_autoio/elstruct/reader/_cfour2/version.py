""" Reads the program name and version number from the output file
"""

import autoparse.pattern as app
import autoparse.find as apf


def program_name(output_str):
    """ Reads the program name from the output file string.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: str
    """

    if _check_name_string(output_str) is not None:
        version_string = _get_version_string(output_str)
        num = version_string.split('.')[0].strip()
        prog_name = 'cfour' + num
    else:
        prog_name = None

    return prog_name


def program_version(output_str):
    """ Reads the program version number from the output file string.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: str
    """

    version_string = _get_version_string(output_str)
    num = version_string.split('.')[0].strip()
    prog_version = version_string.split('.')[1].strip()
    prog_version = num + '.' + prog_version

    return prog_version


def _check_name_string(output_str):
    """ checks to see if the cfour program string is in the output
    """

    pattern = (app.escape('* CFOUR Coupled-Cluster techniques ') +
               app.escape('for Computational Chemistry *'))

    prog_string = apf.has_match(pattern, output_str)

    return prog_string


def _get_version_string(output_str):
    """ obtains the string containing the version number
    """

    pattern = ('Version ' +
               app.capturing(app.one_or_more(app.NONNEWLINE)))

    version_string = apf.first_capture(pattern, output_str)

    return version_string
