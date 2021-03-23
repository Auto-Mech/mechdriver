""" Functions to handle reading the files of the various inputs

    Perhaps return the strings of all the sections that can be parsed out of the file
    # Read the blocks, perhaps return a dct of the strings?

    Check if the strings are none
"""

RUN_INP = 'inp/run.dat'
CSV_INP = 'inp/species.csv'
DAT_INP = 'inp/species.dat'
MECH_INP = 'inp/mechanism.dat'
SORT_INP = 'inp/sort.dat'
MODEL_INP = 'inp/models.dat'
THY_INP = 'inp/theory.dat'


def _read_run(job_path):
    """ Parse the run.dat file
    """
   
    # Read the input string
    run_str = ioformat.ptt.read_inp_str(job_path, RUN_INP, remove_comments='#')

    # Read the blocks
    inp_block = _end_block(inp_str, 'input')
    obj_block = _end_block(inp_str, 'objs')
    job_block = _end_block(inp_str, 'jobs')
    es_tsks_block = _end_block(inp_str, 'es_tsks')
    trans_tsks_block = _end_block(inp_str, 'trans_tsks')
    print_tsks_block = _end_block(inp_str, 'print_tsks')

    # Check if strings exist
    if prnt_tsks_str is None:
        print('*ERROR: No "print_tsks" section defined in run.dat')
        # sys.exit()

    return None


def _read_thy(job_path):
    """ Parse the theory.dat file
    """

    thy_str = ptt.read_inp_str(job_path, THY_INP, remove_comments='#')
    thy_blocks = _named_end_blocks(string, 'level')
    check_thy_sections_nonempty(thy_blocks)


def _read_model(job_path):
    """ Parse the models.dat file
    """

    thy_str = ptt.read_inp_str(job_path, THY_INP, remove_comments='#')
    pes_blocks = _named_end_blocks(thy_str, 'pes_model')
    spc_blocks = _named_end_blocks(thy_str, 'spc_model')


def _read_spc(job_path):
    """ a
    """
    spc_str = ioformat.ptt.read_inp_str(job_path, CSV_INP)


def _read_mech(job_path):
    """Build the PES dct
    """

    # Read the string
    mech_str = ptt.read_inp_str(job_path, MECH_INP, remove_comments='!')
    mech_info = util.read_mechanism_file(
        mech_str, mech_type, spc_dct, sort_rxns=sort_rxns)



# Patterns of use
# maybe loop over different block types for keywords
def _end_block(string, header):
    """ A pattern for a certain block
        rtype: str
    """
    return ioformat.remove_whitespace(
        apf.first_capture(ioformat.ptt.end_section(header), string))


def _named_end_blocks(string, header):
    """ A pattern for a certian block
        rtype: dict[str: str]
    """ 
    caps = apf.all_captures(ptt.end_section_wname2(header), string)
    if caps is not None:
        caps = dict(zip((cap[0] for cap in caps), (cap[1] for cap in caps)))
    return caps


def _paren_block(string, header):
    """ A patter for a certain block
    """
    return apf.first_capture(ptt.paren_section(header), string)

