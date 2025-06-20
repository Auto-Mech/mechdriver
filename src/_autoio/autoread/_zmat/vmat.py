""" z-matrix matrix block parsers
"""

import numpy
from autoparse import cast as _cast
import autoparse.find as apf
import autoparse.pattern as app
from autoread import par

KEY_PATTERN = app.UNSIGNED_INTEGER
NAME_PATTERN = app.VARIABLE_NAME
ENTRY_SEP_PATTERN = app.LINESPACE


def read(string,
         start_ptt=None,
         symb_ptt=par.Pattern.ATOM_SYMBOL,
         key_ptt=KEY_PATTERN,
         name_ptt=NAME_PATTERN,
         entry_start_ptt=None,
         entry_sep_ptt=ENTRY_SEP_PATTERN,
         entry_end_ptt=None,
         line_start_ptt=None,
         line_end_ptt=None,
         last=True,
         case=False):
    """ Read atom symbols and coordinate names and keys of a Z-matrix
        from a string by capturing symbols (from column 1),
        keys (from columns 2, 4, 6), and variable names/values
        (from columns 3, 5, 7).

        :param start_ptt: pattern before the start of the z-matrix
        :type start_ptt: str
        :param symb_ptt: matches atom symbol in first column of the z-matrix
        :type symb_ptt: str
        :param key_ptt: matches key/index in columns 2, 4, 6 of the z-matrix
        :type key_ptt: str
        :param name_ptt: matches z-matrix variable names in columns 3, 5, 7;
            can also be used to match numbers at these positions
        :type name_ptt: str
        :param entry_start_ptt: matches before key_ptt
        :type entry_start_ptt: str
        :param entry_sep_ptt: matches between key_ptt and name_ptt
        :type entry_sep_ptt: str
        :param entry_end_ptt: matches after name_ptt
        :type entry_end_ptt: str
        :param line_start_ptt: matches at the start of a z-matrix line
        :type line_start_ptt: str
        :param line_end_ptt: matches at the end of a z-matrix line
        :type line_end_ptt: str
        :param last: capture the last match, instead of the first?
        :type last: bool
        :param case: make the match case-sensitive?
        :type case: bool
        :return: symbols, key matrix, variable name matrix
        :rtype: tuple
    """

    line_ptts_ = [
        line_pattern(
            num,
            symb_ptt=app.capturing(symb_ptt),
            key_ptt=app.capturing(key_ptt),
            name_ptt=app.capturing(name_ptt),
            entry_start_ptt=entry_start_ptt,
            entry_sep_ptt=entry_sep_ptt,
            entry_end_ptt=entry_end_ptt,
            start_ptt=line_start_ptt,
            end_ptt=line_end_ptt,
        )
        for num in range(4)]

    block_ptt_ = app.capturing(block_pattern(
        symb_ptt=symb_ptt,
        key_ptt=key_ptt,
        name_ptt=name_ptt,
        entry_start_ptt=entry_start_ptt,
        entry_sep_ptt=entry_sep_ptt,
        entry_end_ptt=entry_end_ptt,
        line_start_ptt=line_start_ptt,
        line_end_ptt=line_end_ptt))

    block_ptt_ = block_ptt_ if start_ptt is None else start_ptt + block_ptt_

    if block_ptt_ is not None:
        if last:
            block_str = apf.last_capture(block_ptt_, string, case=case)
        else:
            block_str = apf.first_capture(block_ptt_, string, case=case)
    else:
        block_str = None

    if block_str is not None:
        lines = block_str.splitlines()
        nrows = len(lines)
        symbs = []
        key_mat = numpy.empty((nrows, 3), dtype=numpy.object_)
        name_mat = numpy.empty((nrows, 3), dtype=numpy.object_)
        for row_idx, line in enumerate(lines):
            ncols = min(row_idx, 3)
            caps = _cast(apf.first_capture(line_ptts_[ncols], line, case=case))
            sym = caps if ncols == 0 else caps[0]
            keys = caps[1::2]
            names = caps[2::2]

            symbs.append(sym)
            key_mat[row_idx, :ncols] = keys
            name_mat[row_idx, :ncols] = names

        symbs = tuple(symbs)
        key_mat = tuple(map(tuple, key_mat))
        name_mat = tuple(map(tuple, name_mat))
    else:
        symbs, key_mat, name_mat = None, None, None

    return symbs, key_mat, name_mat


def block_pattern(symb_ptt=par.Pattern.ATOM_SYMBOL,
                  key_ptt=KEY_PATTERN,
                  name_ptt=NAME_PATTERN,
                  entry_start_ptt=None,
                  entry_sep_ptt=ENTRY_SEP_PATTERN,
                  entry_end_ptt=None,
                  line_start_ptt=None,
                  line_end_ptt=None):
    """ Builds a pattern that can match a Z-matrix pattern from a block of
        lines (function currently assumes more than one atom).

        :param symb_ptt: matches atom symbol in first column of block
        :type symb_ptt: str
        :param key_ptt: matches key/index in columns 2, 4, 6 of block
        :type key_ptt: str
        :param name_ptt: matches z-matrix variable names in block in
            columns 3, 5, 7; can also match numbers at these positions
        :type name_ptt: str
        :param entry_start_ptt: matches before key_ptt
        :type entry_start_ptt: str
        :param entry_sep_ptt: matches between key_ptt and name_ptt
        :type entry_sep_ptt: str
        :param entry_end_ptt: matches after name_ptt
        :type entry_end_ptt: str
        :param line_start_ptt: matches at the start of a z-matrix block line
        :type line_start_ptt: str
        :param line_end_ptt: matches at the end of a z-matrix block line
        :type line_end_ptt: str
        :rtype: tuple
    """

    line_ptts = [
        line_pattern(
            num,
            symb_ptt=symb_ptt,
            key_ptt=key_ptt,
            name_ptt=name_ptt,
            entry_start_ptt=entry_start_ptt,
            entry_sep_ptt=entry_sep_ptt,
            entry_end_ptt=entry_end_ptt,
            start_ptt=line_start_ptt,
            end_ptt=line_end_ptt,
        )
        for num in range(4)]

    block_end_ptt = app.series(line_ptts[3], app.padded(app.NEWLINE))

    block_ptt = app.one_of_these([
        app.padded(app.NEWLINE).join(line_ptts[:3] + [block_end_ptt]),
        app.padded(app.NEWLINE).join(line_ptts[:3]),
        app.padded(app.NEWLINE).join(line_ptts[:2]),
        app.padded(app.NEWLINE).join(line_ptts[:1]),
    ])

    return block_ptt


def line_pattern(num,
                 symb_ptt=par.Pattern.ATOM_SYMBOL,
                 key_ptt=KEY_PATTERN,
                 name_ptt=NAME_PATTERN,
                 entry_start_ptt=None,
                 entry_sep_ptt=ENTRY_SEP_PATTERN,
                 entry_end_ptt=None,
                 start_ptt=None,
                 end_ptt=None):
    """ matrix line pattern

        :param symb_ptt: matches atom symbol in first column of block
        :type symb_ptt: str
        :param key_ptt: matches key/index in columns 2, 4, 6 of block
        :type key_ptt: str
        :param name_ptt: matches z-matrix variable names in block in
            columns 3, 5, 7; can also match numbers at these positions
        :type name_ptt: str
        :param entry_start_ptt: matches before key_ptt
        :type entry_start_ptt: str
        :param entry_sep_ptt: matches between key_ptt and name_ptt
        :type entry_sep_ptt: str
        :param entry_end_ptt: matches after name_ptt
        :type entry_end_ptt: str
        :param start_ptt: matches at the start of a z-matrix block line
        :type start_ptt: str
        :param end_ptt: matches at the end of a z-matrix block line
        :type end_ptt: str
        :rtype: tuple
    """

    assert num in range(0, 4)
    entry_ptt = entry_pattern(
        key_ptt=key_ptt,
        name_ptt=name_ptt,
        start_ptt=entry_start_ptt,
        sep_ptt=entry_sep_ptt,
        end_ptt=entry_end_ptt)

    parts = (
        [app.LINE_START] +
        ([] if start_ptt is None else [start_ptt]) +
        [symb_ptt] + num * [entry_ptt] +
        ([] if end_ptt is None else [end_ptt])
    )
    ptt = app.PADDING.join(parts)

    return ptt


def entry_pattern(key_ptt=KEY_PATTERN,
                  name_ptt=NAME_PATTERN,
                  start_ptt=None,
                  sep_ptt=ENTRY_SEP_PATTERN,
                  end_ptt=None):
    """ matrix entry pattern

        :param key_ptt: matches key/index in columns 2, 4, 6 of block
        :type key_ptt: str
        :param name_ptt: matches z-matrix variable names in block in
            columns 3, 5, 7; can also match numbers at these positions
        :type name_ptt: str
        :param start_ptt: matches before key_ptt
        :type start_ptt: str
        :param sep_ptt: matches between key_ptt and name_ptt
        :type sep_ptt: str
        :param end_ptt: matches after name_ptt
        :type end_ptt: str
        :rtype: tuple
    """

    parts = (
        ([] if start_ptt is None else [start_ptt]) +
        [key_ptt, sep_ptt, name_ptt] +
        ([] if end_ptt is None else [end_ptt])
    )
    ptt = app.PADDING.join(parts)

    return ptt
