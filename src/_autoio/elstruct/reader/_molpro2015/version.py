"""
Reads the program name and version number from the output file
"""

import autoparse.pattern as app
import autoparse.find as apf


def program_name(output_str):
    """ Reads the program name from the output file string.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: str
    """

    prog_string = _get_prog_string(output_str)
    prog_name = prog_string[0].strip()
    prog_name = 'molpro' + prog_name

    return prog_name


def program_version(output_str):
    """ Reads the program version number from the output file string.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: str
    """

    prog_string = _get_prog_string(output_str)
    if prog_string:
        prog_version = prog_string[1].strip()
    else:
        prog_version = None
    return prog_version


def _get_prog_string(output_str):
    """ Parses out the program information from the output file string.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: str
    """

    pattern = ('Version' +
               app.SPACES +
               app.capturing(app.INTEGER) +
               app.escape('.') +
               app.capturing(app.INTEGER) +
               app.SPACES +
               'linked')

    prog_string = apf.first_capture(pattern, output_str)

    return prog_string
