""" Convert all of the input strings to data in the format that is used in the code
"""

def _conv_run(run_str):
    """ convert the run.dat input
    """

    keyword_dct = ioformat.ptt.build_keyword_dct(inp_block(run_str))

