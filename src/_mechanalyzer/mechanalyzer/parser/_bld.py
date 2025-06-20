""" Parse mechanism build input file
"""

import ioformat
import autoparse.pattern as app


def build_input_file(bld_inp_str):
    """ Parse the input file used for build file building
    """

    # Parse the sections
    inp_block = ioformat.ptt.end_block(
        bld_inp_str, 'inp', footer='inp')
    rct1_block = ioformat.ptt.end_block(
        bld_inp_str, 'rct1', footer='rct1')
    rct2_block = ioformat.ptt.end_block(
        bld_inp_str, 'rct2', footer='rct2')
    prd_block = ioformat.ptt.end_block(
        bld_inp_str, 'prd', footer='prd')
    rcls_block = ioformat.ptt.end_block(
        bld_inp_str, 'rclasses', footer='rclasses')

    # Set the dictionary for strings
    if inp_block is not None:
        inp_dct = ioformat.ptt.keyword_dct_from_block(
            inp_block, formatvals=True)
    else:
        print('Missing inp block')
        inp_dct = None

    # Set the build series
    if all(block is not None
           for block in (rct1_block, rct2_block, rcls_block)):
        rct1_lst = ioformat.ptt.values_from_block(
            rct1_block, val_ptt=app.VARIABLE_STRING)
        rct2_lst = ioformat.ptt.values_from_block(
            rct2_block, val_ptt=app.VARIABLE_STRING)
        rcls_lst = ioformat.ptt.values_from_block(
            rcls_block, val_ptt=app.VARIABLE_STRING)

        if prd_block is not None:
            prd_lst = ioformat.ptt.values_from_block(
                prd_block, val_ptt=app.VARIABLE_STRING)
        else:
            prd_lst = None

        rseries = (
            (rct1_lst, rct2_lst, prd_lst, rcls_lst),
        )
    else:
        print('Missing rct1, rct2 or rclasses block')
        rseries = None

    return inp_dct, rseries
