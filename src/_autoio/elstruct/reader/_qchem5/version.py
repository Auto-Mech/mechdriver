""" Reads the program name and version number from the output file

    have program name and version together and return a tuple
"""


def program_name(output_str):
    """ Reads the program name from the output file string.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: str
    """

    prog_string = _get_prog_string(output_str)
    prog_name = prog_string.split()[0]
    prog_name = prog_name.replace('-', '').lower()

    return prog_name


def program_version(output_str):
    """ Reads the program version number from the output file string.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: str
    """

    prog_string = _get_prog_string(output_str)
    prog_version = prog_string.split()[1]
    prog_version = prog_version.replace(',', '')

    return prog_version


def _get_prog_string(output_str):
    """ Parses out the program information from the output file string.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: str
    """

    prog_string = None
    for line in output_str.splitlines():
        if 'Q-Chem' in line and 'Inc.' in line:
            prog_string = line
            break

    return prog_string
