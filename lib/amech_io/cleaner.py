"""
 Functions to edit MESS files
"""

import autoparse.pattern as app
import autoparse.find as apf


def remove_whitespace(string):
    """ remove leading spaces, trailing spaces, and empty lines from a string
    """
    empty_line = app.LINE_START + app.maybe(app.LINESPACES) + app.NEWLINE
    trailing_spaces = app.LINESPACES + app.LINE_END
    leading_spaces = app.LINE_START + app.LINESPACES
    pattern = app.one_of_these([empty_line, trailing_spaces, leading_spaces])
    return apf.remove(pattern, string)


def remove_line_comments(string, delim_pattern):
    """ remove line comments marked by a delimiter pattern
    """
    pattern = delim_pattern + app.zero_or_more(app.NONNEWLINE)
    return apf.remove(pattern, string)
