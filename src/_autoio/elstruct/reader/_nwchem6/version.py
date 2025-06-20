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

    prog_string = _get_prog_string(output_str)
    num = prog_string.split()[-1]
    num = num.split('.')[0].strip()
    prog_name = prog_string.split('(')[1][:6]
    prog_name = prog_name.lower().strip()
    prog_name = prog_name + num

    return prog_name


def program_version(output_str):
    """ Reads the program version number from the output file string.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: str
    """

    prog_string = _get_prog_string(output_str)
    prog_version = prog_string.split()[-1].strip()

    return prog_version


def _get_prog_string(output_str):
    """ Parses out the program information from the output file string.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: str
    """

    pattern = app.capturing(
        'Northwest Computational Chemistry Package' + app.SPACE +
        app.escape('(NWChem)') + app.SPACE +
        app.DIGIT + '.' + app.DIGIT)

    prog_string = apf.first_capture(pattern, output_str)

    return prog_string
