"""
Reads the program name and version number from the output file
"""
import pyparsing as pp
import autoparse.pattern as app
import autoparse.find as apf


def program_name(output_string):
    """ reads the program name (here: MainName + MainVersion)
    """
    if _check_name_string(output_string) is not None:
        version_string = _get_version_string(output_string)
        num = version_string.split('.')[0].strip()
        prog_name = 'orca' + num

    return prog_name


def program_version(output_string):
    """ reads the program version number
    """
    version_string = _get_version_string(output_string)
    num = version_string.split('.')[0].strip()
    prog_version = version_string.split('.')[1:]
    prog_version = ('.'.join(prog_version)).strip()
    prog_version = num + '.' + prog_version

    return prog_version


def _check_name_string(output_string):
    """ checks to see if the orca program string is in the output
    """

    pattern = app.escape('* O   R   C   A *')

    prog_string = apf.has_match(pattern, output_string)

    return prog_string


def _get_version_string(output_string):
    """ obtains the string containing the version number
    """
    version = pp.DelimitedList(pp.Word(pp.nums), '.', combine=True)
    parser = ... + pp.Literal("Program Version") + version("version")

    return parser.parseString(output_string).get("version")
