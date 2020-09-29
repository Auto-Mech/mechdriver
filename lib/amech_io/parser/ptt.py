""" Library of patterns to simplify the parsing of input files
"""

import sys
import os
import ast
import ioformat
import autoparse.find as apf
import autoparse.pattern as app


# General Patterns that will be helpful
KEYWORD_KEYVALUE_PATTERN = (
    app.capturing(app.one_or_more(app.NONSPACE)) +
    app.zero_or_more(app.SPACE) +
    '=' +
    app.zero_or_more(app.SPACE) +
    app.capturing(app.LINE_FILL)
)


def read_inp_str(filepath, filename, remove_comments=None):
    """ read the run parameters from a file
    """
    input_file = os.path.join(filepath, filename)
    try:
        with open(input_file, 'r') as inp_file:
            inp_str = inp_file.read()
        if remove_comments:
            inp_str = ioformat.remove_comment_lines(inp_str, remove_comments)
    except FileNotFoundError:
        print('*ERROR: Input file does not exist: ', input_file)
        sys.exit()

    return inp_str


# Reading various section strings from the files #
def paren_section(string):
    """ Read the string that has the global model information
    """
    return (string + app.SPACE + app.escape('=') + app.SPACE +
            app.escape('(') +
            app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
            app.escape(')'))


def end_section(string):
    """ Read the string that has the global model information
    """
    return (string + app.one_or_more(app.SPACE) +
            app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
            'end')


def end_section_wname(string, name):
    """ Read the string that has the global model information
    """
    return (string + app.one_or_more(app.SPACE) + name +
            app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
            'end')


def end_section_wname2(string):
    """ Read the string that has the global model information
    """
    return (string + app.SPACES +
            app.capturing(app.one_or_more(app.NONSPACE)) +
            app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
            'end')


def keyword_pattern(string):
    """ Generates the key pattern string
    """
    value = (string +
             app.zero_or_more(app.SPACE) +
             '=' +
             app.zero_or_more(app.SPACE) +
             app.capturing(app.one_or_more(app.NONSPACE)))
    return set_value_type(value)


def parse_idx_inp(idx_str):

    """ parse idx string
    """
    idx_str = idx_str.strip()
    if idx_str.isdigit():
        idxs = [int(idx_str)]
    if '-' in idx_str:
        [idx_begin, idx_end] = idx_str.split('-')
        idxs = list(range(int(idx_begin), int(idx_end)+1))
    return idxs


# Build a keyword dictionary
def build_keyword_dct(section_str):
    """ Take a section with keywords defined and build
        a dictionary for the keywords
    """

    keyword_dct = {}
    section_str = ioformat.remove_whitespace(section_str)
    for line in section_str.splitlines():
        key_val = apf.first_capture(KEYWORD_KEYVALUE_PATTERN, line)
        formtd_key, formtd_val = format_param_vals(key_val)
        keyword_dct[formtd_key] = formtd_val
    return keyword_dct


def build_keyword_lst(section_str):
    """ build lst
    """

    keyword_lst = []
    section_str = ioformat.remove_whitespace(section_str)
    for line in section_str.splitlines():
        keyword_lst.append(line)
    return keyword_lst


def build_vals_lst(section_str):
    """ build lst
    """
    val_lst = []
    section_str = ioformat.remove_whitespace(section_str)
    for line in section_str.splitlines():
        val_lst.extend((float(val) for val in line.split()))

    return val_lst


# Functions from read_dat file
def format_param_vals(pvals):
    """ format param vals string
    """
    [keyword, value] = pvals

    frmtd_keyword = keyword.strip().lower()

    value = value.strip()
    # Format values if it is a list (of string(s), boolean(s), int(s))
    if all(sym in value for sym in ('[[', ']]')):
        value = value.replace('D', '').replace('d', '')
        value = ast.literal_eval(value)
        frmtd_value = []
        for sub_lst in value:
            assert all(isinstance(val, int) for val in sub_lst)
            frmtd_value.append(['D{}'.format(val) for val in sub_lst])
    elif all(sym in value for sym in ('[', ']')):
        value = value.replace('[', '').replace(']', '')
        value = value.split(',')
        frmtd_value = []
        # Set string in list to boolean or integer if needed
        for elm in value:
            elm = elm.strip()
            if ':' in elm:
                elm_lst = elm.split(':')
                frmtd_value.append([float(elm_lst[0]), elm_lst[1]])
            else:
                frmtd_value.append(set_value_type(elm))
    else:
        # Format values if it has singular value
        frmtd_value = set_value_type(value)

    return frmtd_keyword, frmtd_value


def set_value_type(value):
    """ set type of value
        right now we handle True/False boolean, int, float, and string
    """
    if value == 'True':
        frmtd_value = True
    elif value == 'False':
        frmtd_value = False
    elif value == 'None':
        frmtd_value = None
    elif value.isdigit():
        frmtd_value = int(value)
    elif 'e' in value:
        try:
            frmtd_value = float(value)
        except ValueError:
            frmtd_value = value
    elif '.' in value:
        if value.replace('.', '').replace('-', '').isdigit():
            frmtd_value = float(value)
    else:
        frmtd_value = value

    return frmtd_value
