""" Various formatting functions used by each I/O module
"""

import os
from mako.template import Template
import more_itertools as mit
import autoparse.pattern as app
import autoparse.find as apf


# Build formatted strings
def build_mako_str(template_file_name, template_src_path, template_keys,
                   remove_whitespace=True):
    """ Uses an input dictionary to fill in Mako template file containing the
        keys of the dictionary, then writes a string corresponding to the
        filled-in Mako template.

        :param template_file_name: Name of the Mako template file
        :type template_file_name: str
        :param template_src_path: Path where Mako template file resides
        :type template_str_path: str
        :param template_keys: keys and values used to fill Mako template
        :type template_keys: dict[template key: template value]
        :rtype: str
    """

    template_file_path = os.path.join(template_src_path, template_file_name)
    mako_str = Template(filename=template_file_path).render(**template_keys)

    if remove_whitespace:
        mako_str = remove_trail_whitespace(mako_str)

    return mako_str


def indent(string, nspaces):
    """ Indents each of the lines of a multiline string.

        :param string: Input string to indent
        :type string: str
        :param nspaces: number of spaces to indent the lines of the string
        :type nspaces: int
        :return indented_string: string indented by nspaces
        :rtype string
    """

    pad = nspaces * ' '
    indented_string = ''.join(pad+line for line in string.splitlines(True))

    return indented_string


def change_line(string, newline, searchline):
    """ Search for a line in a string and replace it with a new one.
    """

    string_lines = string.splitlines()

    # Find line with search string
    for i, line in enumerate(string_lines):
        if searchline.strip() == line.strip():
            linenum = i
            break

    # Change the line by changing list at index of search string
    string_lines[linenum] = newline

    return '\n'.join(string_lines)


def add_line(string, addline, searchline, position):
    """ Add a line to some string at some positin specified by a line
        currently present in the string.
    """

    string_lines = string.splitlines()

    # Find line with search string
    for i, line in enumerate(string_lines):
        if searchline.strip() == line.strip():
            linenum = i
            break

    # Add the line
    if position == 'before':
        string_lines.insert(linenum, addline)
    else:
        string_lines.insert(linenum+1, addline)

    return '\n'.join(string_lines)


def addchar(string, char, side='pre'):
    """ Pre- or Post-pends a character to a string.

        :param string: Input string to add a character to
        :type string: str
        :param char: character to to prepend to string
        :type char: str
        :return new_string: string with new character
        :rtype string
    """

    assert side in ('pre', 'post'), (
        'side must be pre or post'
    )
    char = str(char).strip()

    if side == 'pre':
        new_string = char + ' ' + string
    else:
        new_string = string + ' ' + char

    return new_string


# Parsing methods for generic strings
def headlined_sections(string, headline_pattern):
    """ Returns sections with headlines matching a pattern.

        :param string: string to return sections
        :type string: str
        :param headline_pattern: pattern to demarcate strings into sections
        :type headline_pattern: str
        :rtype: str
    """

    lines = string.splitlines()
    join_lines = '\n'.join
    pattern_matcher = apf.matcher(headline_pattern)
    lines = mit.lstrip(lines, pred=lambda line: not pattern_matcher(line))
    sections = list(map(join_lines, mit.split_before(lines, pattern_matcher)))

    return sections


# Clean up strings by removing unneccessary whitespace and comments
def remove_whitespace_from_string(string):
    """ Removes leading spaces, trailing spaces, and empty lines from a string.

        :param string: string to clean up
        :type string: str
        :rtype: str
    """

    empty_line = app.LINE_START + app.maybe(app.LINESPACES) + app.NEWLINE
    trailing_spaces = app.LINESPACES + app.LINE_END
    leading_spaces = app.LINE_START + app.LINESPACES
    pattern = app.one_of_these([empty_line, trailing_spaces, leading_spaces])

    return apf.remove(pattern, string)


def remove_trail_whitespace(string):
    """ Removes trailing spaces and empty lines from a input string.

        :param string: string to remove trailing whitespace
        :type string: str
        :rtype: str
    """

    empty_line = app.LINE_START + app.maybe(app.LINESPACES) + app.NEWLINE
    trailing_spaces = app.LINESPACES + app.LINE_END
    pattern = app.one_of_these([empty_line, trailing_spaces])

    return apf.remove(pattern, string)


def remove_comment_lines(string, delim_pattern):
    """ Remove comment lines marked by a delimiter pattern.

        :param string: string to remove comments
        :type string: str
        :param delim_pattern: pattern of delimiter to identify comment lines
        :type delim_pattern: str
        :rtype: str
    """

    pattern = delim_pattern + app.zero_or_more(app.NONNEWLINE)

    return apf.remove(pattern, string)


def remove_empty_lines(string):
    """ Removes empty lines from a input string.

        :param string: string to remove trailing whitespace
        :type string: str
        :rtype: str
    """

    empty_line = app.LINE_START + app.maybe(app.LINESPACES) + app.NEWLINE

    return apf.remove(empty_line, string)
