""" Library of patterns to simplify the parsing of input files
"""

import ast
import autoparse.find as apf
import autoparse.pattern as app
from ioformat._format import remove_whitespace_from_string


# PATTERNS
KEYWORD_KEYVALUE_PATTERN = (
    app.capturing(app.one_or_more(app.NONSPACE)) +
    app.zero_or_more(app.SPACE) +
    '=' +
    app.zero_or_more(app.SPACE) +
    app.capturing(app.LINE_FILL)
)
FULL_BLOCK_PTT = app.capturing(app.one_or_more(app.WILDCARD, greedy=False))
# GEN_PTT = app.capturing(app.one_or_more(app.NONSPACE))
GEN_PTT = app.capturing(app.one_or_more(app.LINE_FILL))


def keyword_value_ptt(key=None):
    """ Build a regex pattern to search for a keyword and its value, defined as
            key = value

        :param key: string for a specific key to find
        :type key: str
        :rtype: str
    """
    keywrd = key if key is not None else app.one_or_more(app.NONSPACE)
    return (app.capturing(keywrd) +
            app.ZSPACES + app.not_preceded_by('ChI') +  app.escape('=') + app.ZSPACES +
            GEN_PTT)


def paren_block_ptt(key=None):
    """ Build a regex pattern to search for a keyword and its value, defined as
            key = (value)

        :param key: string for a specific key to find
        :type key: str
        :rtype: str
    """
    keywrd = key if key is not None else app.one_or_more(app.NONSPACE)
    return (app.capturing(keywrd) + app.SPACES + app.escape('=') + app.SPACES +
            app.escape('(') + FULL_BLOCK_PTT + app.escape(')'))


def named_end_block_ptt(header, footer=None):
    """ Build a regex pattern to search for blocks of end
            key = (value)
            {header} {name}
               **TEXT**
            end {footer}
    """

    # Set the top and bottom of the end block pattern
    top_ptt = (
        header + app.SPACES + app.capturing(app.one_or_more(app.NONSPACE)))

    bot_ptt = 'end'
    if footer is not None:
        bot_ptt += app.SPACES + footer

    return top_ptt + FULL_BLOCK_PTT + bot_ptt


def end_block_ptt(header, name=None, footer=None):
    """ Read the string that has the global model information
        can also capture with the name if needed...
        {header} {name}
          DATA
        end {footer}
    """

    # Set the top and bottom of the end block pattern
    top_ptt = header
    if name is not None:
        top_ptt += app.SPACES + app.capturing(name)

    bot_ptt = 'end'
    if footer is not None:
        bot_ptt += app.SPACES + footer

    return top_ptt + FULL_BLOCK_PTT + bot_ptt


def symb_block_ptt(symb, header):
    """ Read the string that has the global model information
        can also capture with the name if needed...
        {symb}{header}
          DATA
        {symb}end
    """

    # Set the top and bottom of the end block pattern
    top_ptt = app.escape(symb) + app.capturing(header)
    bot_ptt = app.escape(symb) + 'end'

    return top_ptt + FULL_BLOCK_PTT + bot_ptt


# Simple Block Parsers
def paren_blocks(string, key=None):
    """ A pattern for a certain block
    """
    return apf.all_captures(paren_block_ptt(key=key), string)


def keyword_value_blocks(string, key=None):
    """ A pattern for a certain block
    """
    return apf.all_captures(keyword_value_ptt(key=key), string)


def named_end_blocks(string, header, footer=None):
    """ A pattern for a certian block

        rtype: dict[str: str]
    """
    caps = apf.all_captures(
        named_end_block_ptt(header, footer=footer), string)
    if caps is not None:
        caps = dict(zip((cap[0] for cap in caps), (cap[1] for cap in caps)))
    return caps


# single blcok
def symb_block(string, symb, header):
    """ AA pattern for a certan block
    """
    return apf.first_capture(symb_block_ptt(symb, header), string)


def end_block(string, header, name=None, footer=None):
    """ A pattern for a certain block

        rtype: str
    """
    _ptt = end_block_ptt(header, name=name, footer=footer)
    cap = apf.first_capture(_ptt, string)

    return cap if cap is not None else None


# Build keyword dictiaonries
def keyword_dcts_from_blocks(block_dct):
    """ Build a dict of dicts from string blocks that are stored in
        a dictionary.

        :param block_dct: dictionary of strings with a name for a key
        :type block_dct: dict[str: str]
        :rtype: dict[str: dict[str:str]]
    """

    new_block_dct = {}
    if block_dct is not None:
        for key, block in block_dct.items():

            key_dct = keyword_dct_from_paren_blocks(block)
            if key_dct is None or key_dct == {}:
                key_dct = keyword_dct_from_block(block)

            new_block_dct[key] = key_dct

    return new_block_dct


def keyword_dct_from_paren_blocks(block):
    """ Obtains the keyword-value pairs that are defined in blocks from parentheses
        Could be just a bunch of values or a set of keyword-value pairs
            # First check for keyword list or just vals
    """

    ret = {}

    pblocks = paren_blocks(block)
    if pblocks is not None:
        for pblock in pblocks:
            key_dct = keyword_dct_from_block(pblock[1])
            if key_dct is not None:
                ret[pblock[0]] = key_dct
            else:
                vals = values_from_block(pblock[1], non_val_ptt=app.LETTER)
                if vals:
                    ret[pblock[0]] = vals
    else:
        ret = None

    return ret


def keyword_dct_from_block(block, formatvals=True):
    """ Take a section with keywords defined and build
        a dictionary for the keywords
        assumes a block that is a list of key-val pairs
    """

    key_dct = None

    if block is not None:
        block = remove_whitespace_from_string(block)
        key_val_blocks = keyword_value_blocks(block)
        if key_val_blocks is not None:
            key_dct = {}
            for key, val in key_val_blocks:
                if formatvals:
                    formtd_key, formtd_val = format_keyword_values(key, val)
                    key_dct[formtd_key] = formtd_val
                else:
                    key_dct[key] = val
    return key_dct


# Build various objects containing keyword and value information
def values_from_block(block, val_ptt=app.NUMBER, non_val_ptt=None):
    """ Takes a multiline string that consists solely of floats and
        converts this block into a list of numbers
        could call set_value_type for generality I guess
        prob just do a capture of nums (floats, int, etc)
    """
    caps = apf.all_captures(val_ptt, block)
    if non_val_ptt is not None:
        non_caps = apf.all_captures(non_val_ptt, block)
        if non_caps:
            caps = ()
    if caps:
        vals = tuple(set_value_type(cap) for cap in caps)
    else:
        vals = None

    return vals


def idx_lst_from_line(line):
    """ Build a list of indices from a block of tests
    """

    idxs = []
    for string in line.strip().split(','):
        _strip_string = string.strip()
        if _strip_string.isdigit():
            idxs.append(int(_strip_string))
        elif '-' in line:
            [idx_begin, idx_end] = _strip_string.split('-')
            idxs.extend(list(range(int(idx_begin), int(idx_end)+1)))

    return tuple(idxs)


# Formats the values associated with various keywords
def format_keyword_values(keyword, value):
    """ Takes a keyword-value pair in string formats and then returns
        the pair with their types matching the internal Python version.
        Convert string to string, boolean, int, float, etc

        :param key_val_pair:  keyword and its
        :type key_val_pair: (str, str)
        :rtype: (type(str), type(str))
    """

    # Format the keyword
    frmtd_keyword = set_value_type(keyword.strip().lower())
    # Format values if it is a list (of string(s), boolean(s), int(s))
    # Additional functionality is used to handle when values are lists
    value = value.strip()
    if all(sym in value for sym in ('[[', ']]')):
        # Need to fix to handle dihedrals
        # ast converson breaks for strings if they are not in quotes
        if frmtd_keyword == 'tors_names':
            value = value.replace('D', '').replace('d', '')
            value = ast.literal_eval(value)
            frmtd_value = tuple(
                tuple(f'D{val}' for val in sub_lst) for sub_lst in value
            )
        else:
            frmtd_value = ast.literal_eval(value)
            frmtd_value = tuple(tuple(x) for x in frmtd_value)
    elif all(sym in value for sym in ('[', ']')):
        value = value.replace('[', '').replace(']', '')
        value = value.split('"')
        flat_value = []
        for i, val in enumerate(value):
            if i%2 == 0 and val:
                flat_value.extend(
                    va for va in val.split(',') if va.replace(' ','') != '')
            elif val:
                flat_value.append(val)
        value = flat_value
        frmtd_value = ()
        # Set string in list to boolean or integer if needed
        for elm in value:
            elm = elm.strip()
            if ':' in elm:
                elm_lst = elm.split(':')
                if 'ene' in frmtd_keyword:
                    frmtd_value += ((float(elm_lst[0]), elm_lst[1]),)
                else:
                    frmtd_value += ((float(elm_lst[0]), float(elm_lst[1])),)
            else:
                frmtd_value += (set_value_type(elm),)
    else:
        # Format values if it has singular value
        frmtd_value = set_value_type(value)

    return frmtd_keyword, frmtd_value


def set_value_type(value):
    """ set type of value
        right now we handle True/False boolean, int, float, and string
    """

    if value.lower() == 'true':
        frmtd_value = True
    elif value.lower() == 'false':
        frmtd_value = False
    elif value.lower() == 'none':
        frmtd_value = None
    elif value.isdigit():
        frmtd_value = int(value)
    elif 'e' in value or '.' in value:
        try:
            frmtd_value = float(value)
        except ValueError:
            frmtd_value = value
    else:
        frmtd_value = value

    return frmtd_value
