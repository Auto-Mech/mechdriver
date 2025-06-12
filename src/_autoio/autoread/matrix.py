""" matrix parsers
"""

import numpy
from autoparse import cast as _cast
import autoparse.find as apf
import autoparse.pattern as app


VALUE_PATTERN = app.one_of_these([app.FLOAT])


def read(string,
         val_ptt=VALUE_PATTERN,
         start_ptt=None,
         block_start_ptt=None,
         line_start_ptt=None,
         last=True,
         tril=False,
         case=False):
    """ Read
    """

    mats = read_all(string,
                    val_ptt=val_ptt,
                    start_ptt=start_ptt,
                    block_start_ptt=block_start_ptt,
                    line_start_ptt=line_start_ptt,
                    tril=tril,
                    case=case)
    if mats is not None:
        mat = mats[-1] if last else mats[0]
    else:
        mat = None

    return mat


def read_all(string,
             val_ptt=VALUE_PATTERN,
             start_ptt=None,
             block_start_ptt=None,
             line_start_ptt=None,
             tril=False,
             case=False):
    """ Reads an M x N matrix from a string by capturing the matrix values,
        which can be rectangular or lower-triangular, and
        may or may not be broken into multiple blocks.

        :param val_ptt: matches numeric matrix entry
        :type val_ptt: str
        :param start_ptt: matches before start of the matrix
        :type start_ptt: str
        :param block_start_ptt: matches at start of each block within the
            matrix
        :type block_start_ptt: str
        :param line_start_ptt: matches at start of each line in the matrix or
            matrix block
        :param line_start_ptt: str
        :param last: capture the last match, instead of the first?
        :type last: bool
        :param case: make the match case-sensitive?
        :type case: bool
        :rtype: tuple(tuple(float))
    """

    line_ptt_ = line_pattern(val_ptt=val_ptt, start_ptt=line_start_ptt,
                             capture_values=True)
    block_ptt_ = block_pattern(val_ptt=val_ptt, start_ptt=block_start_ptt,
                               line_start_ptt=line_start_ptt,
                               capture_block=True)
    blocks_ptt_ = blocks_pattern(val_ptt=val_ptt, start_ptt=start_ptt,
                                 block_start_ptt=block_start_ptt,
                                 line_start_ptt=line_start_ptt,
                                 capture_blocks=True)

    blocks_str_lst = apf.all_captures(blocks_ptt_, string, case=case)
    blocks_str_lst = blocks_str_lst if blocks_str_lst is not None else ()
    # blocks_str_lst = (
    #     apf.last_capture(blocks_ptt_, string, case=case) if last else
    #     apf.first_capture(blocks_ptt_, string, case=case))

    mats = ()
    for blocks_str in blocks_str_lst:
        block_strs = apf.all_captures(block_ptt_, blocks_str, case=case)

        if block_strs is not None:
            if not tril:
                rows = numpy.concatenate(
                    [_block_rows(block_str, val_ptt, line_ptt_, case=case)
                     for block_str in block_strs], axis=1)
                mats += (_matrix(rows),)
            else:
                rows = list(_block_rows(
                    block_strs[0], val_ptt, line_ptt_, case=case))
                nrows = len(rows)
                for block_str in block_strs[1:]:
                    block_rows = _block_rows(
                        block_str, val_ptt, line_ptt_, case=case)
                    nblock_rows = len(block_rows)
                    for block_row_idx, row_idx in enumerate(
                            range(nrows-nblock_rows, nrows)):
                        rows[row_idx] += block_rows[block_row_idx]

                mats += (_symmetric_matrix_from_lower_triangle(rows),)

        else:
            mats += (None,)

    if not mats:
        mats = None

    return mats


def _matrix(rows):
    """ Format the values of matrix read from a string into a tuple-of-tuples.

        :param rows: rows of values of matrix
        :type rows: numpy.ndarray
        :rtype: tuple(tuple(float))
    """

    mat = numpy.array(rows)
    assert mat.ndim == 2

    return tuple(map(tuple, mat))


def _symmetric_matrix_from_lower_triangle(tril_rows):
    """ Build a full M x M symmetric matrix from the rows of
        a lower triangular matrix.

        :param tril_rows: rows of lower triangular matrix
        :type tril_rows: list(float)
        :rtype: tuple(tuple(float))
    """

    nrows = len(tril_rows)
    mat = numpy.zeros((nrows, nrows), dtype=numpy.object_)
    for row_idx, tril_row in enumerate(tril_rows):
        mat[row_idx, :row_idx+1] = tril_row
        mat[:row_idx+1, row_idx] = tril_row
    # make sure we ended up with a square symmetric matrix
    assert mat.ndim == 2 and mat.shape[0] == mat.shape[1]

    return tuple(map(tuple, mat))


def _block_rows(block_str, val_ptt, line_ptt_, case=False):
    """ Reads th rows of a blokc of text.

        :param block_str: string to read lines from
        :type block_str: str
        :param val_ptt: matches numeric matrix entry
        :type val_ptt: str
        :param line_ptt_: pattern for matching line of the block
        :type line_ptt_: str
        :param case: make the match case-sensitive?
        :type case: bool
        :rtype: list(float)
    """

    rows = []
    val_ptt_ = app.capturing(val_ptt)
    for line_str in apf.all_captures(line_ptt_, block_str, case=case):
        row = _cast(apf.all_captures(val_ptt_, line_str))
        rows.append(row)

    return rows


def blocks_pattern(val_ptt=VALUE_PATTERN,
                   start_ptt=None,
                   block_start_ptt=None,
                   line_start_ptt=None,
                   capture_blocks=False):
    """ Build a pattern that matches a multi-block matrix.

        :param val_ptt: matches numeric matrix entry
        :type val_ptt: str
        :param start_ptt: pattern before start of the matrix block
        :type start_ptt: str
        :param block_start_ptt: matches at start of each block of
        :type block_start_ptt: str
        :param line_start_ptt: matches at start of each line of each block
        :type line_start_ptt: str
        :param capture_blocks: add capturing pattern for the matrix block
        :type capture_blocks: bool
        :rtype: list(float)
    """

    block_ptt = block_pattern(val_ptt=val_ptt, start_ptt=block_start_ptt,
                              line_start_ptt=line_start_ptt)

    blocks_ptt_ = app.series(block_ptt, app.padded(app.NEWLINE))

    if capture_blocks:
        blocks_ptt_ = app.capturing(blocks_ptt_)

    blocks_ptt_ = blocks_ptt_ if start_ptt is None else start_ptt + blocks_ptt_

    return blocks_ptt_


def block_pattern(val_ptt=VALUE_PATTERN,
                  start_ptt=None,
                  line_start_ptt=None,
                  capture_block=False):
    """ Build a pattern that matches a block with a single matrix.

        :param val_ptt: matches numeric matrix entry
        :type val_ptt: str
        :param start_ptt: pattern before start of the matrix block
        :type start_ptt: str
        :param line_start_ptt: matches at start of each line of each block
        :type line_start_ptt: str
        :param capture_blocks: add capturing pattern for the matrix block
        :type capture_blocks: bool
        :rtype: list(float)
    """

    line_ptt = line_pattern(val_ptt=val_ptt, start_ptt=line_start_ptt)
    block_ptt_ = app.series(line_ptt, app.padded(app.NEWLINE))

    if capture_block:
        block_ptt_ = app.capturing(block_ptt_)

    block_ptt_ = block_ptt_ if start_ptt is None else start_ptt + block_ptt_

    return app.padded(block_ptt_)


def line_pattern(val_ptt=VALUE_PATTERN,
                 start_ptt=None,
                 capture_values=False):
    """ Build a pattern that matches a line from a block with a single matrix.

        :param val_ptt: matches numeric matrix entry
        :type val_ptt: str
        :param start_ptt: matches at start of the line of the block
        :type start_ptt: str
        :param capture_values: add capturing pattern for the values in the line
        :type capture_values: bool
        :rtype: list(float)
    """

    vals_ptt = app.series(val_ptt, app.LINESPACES)

    if capture_values:
        vals_ptt = app.capturing(vals_ptt)

    parts = (
        ([] if start_ptt is None else [start_ptt]) + [vals_ptt])

    ptt = app.LINE_START + app.padded(app.LINESPACES.join(parts))

    return ptt
