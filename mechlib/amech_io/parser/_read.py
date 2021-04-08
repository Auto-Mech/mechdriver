""" Functions to handle reading the files of the various inputs

    Perhaps return the strings of all the sections that can be parsed out of
    the file
    # Read the blocks, perhaps return a dct of the strings?

    Check if the strings are none
"""

import automol
import ioformat
import mechanalyzer
from mechlib.amech_io.parser import tsks
from mechlib.amech_io import printer as ioprinter
from mechlib.amech_io.parser.mechanism import build_pes_dct
from mechlib.amech_io.parser.species import modify_spc_dct
from mechlib.amech_io.parser.species import geometry_dictionary


RUN_INP = 'inp/run.dat'
CSV_INP = 'inp/species.csv'
DAT_INP = 'inp/species.dat'
MECH_INP = 'inp/mechanism.dat'
SORT_INP = 'inp/sort.dat'
MODEL_INP = 'inp/models.dat'
THY_INP = 'inp/theory.dat'


def read_amech_inp(job_path):
    """ master read for simplicity for now
    """

    # Read method input
    thy_dct = read_thy(job_path)
    kin_mod_dct, spc_mod_dct = read_model(job_path)

    # Read the run
    a = read_run(job_path, thy_dct, kin_mod_dct, spc_mod_dct)

    # Read chemisry input
    # b = read_spc(job_path)

    return None


def read_run(job_path, thy_dct, kin_mod_dct, spc_mod_dct):
    """ Parse the run.dat file
    """

    # Read the input string
    run_str = ioformat.pathtools.read_file(
        job_path, RUN_INP, remove_comments='#', remove_whitespace=True)
    ioprinter.reading('run.dat...', newline=1)  # Add a Found <file> to msg

    # Read the main blocks
    inp_block = ioformat.ptt.end_block(run_str, 'input', footer='input')

    pes_block = ioformat.ptt.end_block(run_str, 'pes', footer='pes')
    spc_block = ioformat.ptt.end_block(run_str, 'spc', footer='spc')

    es_tsks_block = ioformat.ptt.end_block(
        run_str, 'els', footer='els')
    trans_tsks_block = ioformat.ptt.end_block(
        run_str, 'transport', footer='transport')
    therm_tsks_block = ioformat.ptt.end_block(
        run_str, 'thermo', footer='thermo')
    ktp_tsks_block = ioformat.ptt.end_block(
        run_str, 'ktp', footer='ktp')
    proc_tsks_block = ioformat.ptt.end_block(
        run_str, 'proc', footer='proc')

    # Parse information in the blocks
    inp_key_dct = _keyword_dct(inp_block, {})
    print('inp_dct\n', inp_key_dct)

    pes_idx_dct = _pes_idxs(pes_block)
    spc_idx_dct = _spc_idxs(spc_block)
    print('pes_dct\n', pes_idx_dct)
    print('spc_dct\n', spc_idx_dct)

    print('es\n', es_tsks_block)
    es_tsks_lst = tsks.es_tsk_lst(es_tsks_block, thy_dct)
    # therm_tsks_lst = tsks.trans_tsk_lst(trans_tsks_block, kin_mod_dct, spc_mod_dct)
    # ktp_tsks_lst = tsks.trans_tsk_lst(trans_tsks_block, kin_mod_dct, spc_mod_dct)
    # trans_tsks_lst = tsks.trans_tsk_lst(trans_tsks_block, thy_dct)
    # proc_tsks_lst = tsks.proc_tsk_lst(proc_tsks_block, thy_dct)

    print('es\n', es_tsks_lst)
    # print('trans\n', trans_tsks_block)
    # print('therm\n', therm_tsks_block)
    # print('ktp\n', ktp_tsks_block)
    # print('proc\n', proc_tsks_block)

    return inp_key_dct, pes_idx_dct, spc_idx_dct


def read_thy(job_path):
    """ Parse the theory.dat file
    """

    thy_str = ioformat.pathtools.read_file(
        job_path, THY_INP, remove_comments='#', remove_whitespace=True)
    ioprinter.reading('theory.dat...', newline=1)

    thy_blocks = ioformat.ptt.named_end_blocks(
        thy_str, 'level', footer='level')
    thy_dct = ioformat.ptt.keyword_dcts_from_blocks(thy_blocks)
    print(thy_dct)

    # check_dictionary2(thy_dct, THY_REQUIRED_KEYWORDS THY_SUPPORTED_KEYWORDS)

    return thy_dct


def read_model(job_path):
    """ Parse the models.dat file
    """

    mod_str = ioformat.pathtools.read_file(
        job_path, MODEL_INP, remove_comments='#', remove_whitespace=True)
    ioprinter.reading('model.dat...', newline=1)

    kin_blocks = ioformat.ptt.named_end_blocks(mod_str, 'kin')
    spc_blocks = ioformat.ptt.named_end_blocks(mod_str, 'spc')

    kin_mod_dct = _merge_glob(
        ioformat.ptt.keyword_dcts_from_blocks(kin_blocks), keep_glob=True)
    spc_mod_dct = _merge_glob(
        ioformat.ptt.keyword_dcts_from_blocks(spc_blocks), keep_glob=True)

    print('kin_mod\n', kin_mod_dct)
    print('spc_mod\n', spc_mod_dct)

    return kin_mod_dct, spc_mod_dct


def read_spc(job_path):
    """ a
    """

    # Read all of the potential species files
    spc_str = ioformat.pathtools.read_file(job_path, CSV_INP)
    ioprinter.reading('species.csv...', newline=1)

    amech_str = ioformat.pathtools.read_file(job_path, DAT_INP, remove_comments='#')
    ioprinter.reading('species.dat...', newline=1)

    geo_dct = geometry_dictionary(job_path)
    ioprinter.reading('geom.xyzs...', newline=1)

    # Build the spc dct
    spc_dct = mechanalyzer.parser.spc.build_spc_dct(spc_str, 'csv')

    amech_blocks = ioformat.ptt.named_end_blocks(amech_str, 'spc', footer='spc')
    amech_dct = _merge_glob(
        ioformat.ptt.keyword_dcts_from_blocks(amech_blocks))

    mod_spc_dct = modify_spc_dct(spc_dct, amech_dct, geo_dct)

    return spc_str


def read_mech(job_path, spc_dct, mech_type='chemkin'):
    """Build the PES dct
    """

    # Read the string
    mech_str = ioformat.pathtools.read_file(
        job_path, MECH_INP, remove_comments='!', remove_whitespace=True)
    pes_dct = build_pes_dct(
        job_path, mech_type, spc_dct, run_obj_dct, sort_rxns=False)
    # mech_info = util.read_mechanism_file(
    #     mech_str, mech_type, spc_dct, sort_rxns=False)

    return pes_dct


# Keyword dict build and check
def _keyword_dct(inp_str, def_dct):
    """ merge a dictionary from the input
        :param inp_dct: dictionary of input keywords
        :param def_dct: dictionary of default values
                        (only has ones where defaults, won't have vals for
                         keys that user must input like run_prefix)
    """
    inp_dct = ioformat.ptt.keyword_dct_from_block(inp_str)
    return automol.util.dict_.right_update(def_dct, inp_dct)


def _keyword_lst(inp_str, key_lst):
    """ build lst and check against supported lsts
        :param inp_dct: dictionary of input keywords
        :param def_dct: dictionary of default values
                        (only has ones where defaults, won't have vals for
                         keys that user must input like run_prefix)
    """
    return ioformat.ptt.idx_lst_from_line(inp_str)


def _pes_idxs(string):
    """  Build a dictionary of the PESs.
            {pes_idx: [chn_idxs]}
            breaks if pes_idx is given on two lines
    """

    run_pes = {}
    for line in string.strip().splitlines():
        [pes_nums, chn_nums] = line.split(':')
        pes_idxs = ioformat.ptt.idx_lst_from_line(pes_nums)
        chn_idxs = ioformat.ptt.idx_lst_from_line(chn_nums)
        for idx in pes_idxs:
            run_pes.update({idx: chn_idxs})

    return run_pes


def _spc_idxs(string):
    """  Build a dictionary of the PESs.
            {pes_idx: [chn_idxs]}
    """

    spc_idxs = ()
    for line in string.splitlines():
        spc_idxs += ioformat.ptt.idx_lst_from_line(line)

    return {1: spc_idxs}


def _merge_glob(dct, keep_glob=False):
    """ [change to pull out the glob dct and just use it to overwrite the spc dct
        or have parser return an additional one?

        amech dct returns nothing if only it is the global dct...and want keep_glob=false
    """


    new_dct = {}
    glob = dct.get('global', {})

    print('dct', dct)
    print('glob', glob)

    if keep_glob:
        names = tuple(x for x in dct)
    else:
        names = tuple(x for x in dct if x != 'global')
    for name in names:
        new_dct[name] = automol.util.dict_.right_update(glob, dct[name])

    return new_dct
