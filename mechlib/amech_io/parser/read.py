""" Reads all of the MechDriver input files and formats 
    the data into usable Python objects that are used by
    all subsequent drivers in the workflow.
"""

import automol
import ioformat
import mechanalyzer
from mechlib.amech_io import printer as ioprinter


from mechlib.amech_io.parser._run import pes_idxs, spc_idxs, tsk_lst
from mechlib.amech_io.parser._mech import build_pes_dct
from mechlib.amech_io.parser._spc import modify_spc_dct
from mechlib.amech_io.parser._spc import geometry_dictionary
from mechlib.amech_io.parser._keywrd import check_val_dictionary2
from mechlib.amech_io.parser._keywrd import defaults_from_val_dct
from mechlib.amech_io.parser._keywrd import THY_VAL_DCT, SPC_VAL_DCT





# Names of input file strings
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
    print('THY TEST')
    thy_dct = read_thy(job_path)
    print('thy dct\n', thy_dct)

    print('MOD TEST')
    kin_mod_dct, spc_mod_dct = read_model(job_path)
    print('kin_mod\n', kin_mod_dct)
    print('spc_mod\n', spc_mod_dct)

    # Read the run
    # print('RUN TEST')
    # inp_key_dct, pes_idx_dct, spc_idx_dct = read_run(
    #     job_path, thy_dct, kin_mod_dct, spc_mod_dct)
    # print('inp_dct\n', inp_key_dct)
    # print('pes_dct\n', pes_idx_dct)
    # print('spc_dct\n', spc_idx_dct)

    # Read chemisry input
    print('SPC TEST')
    spc_dct = read_spc(job_path)
    print('spc_dct\n', spc_dct)

    return None


def read_run(job_path, thy_dct, kin_mod_dct, spc_mod_dct):
    """ Parse the run.dat file
    """

    # Read the input string from the file
    run_str = ioformat.pathtools.read_file(
        job_path, RUN_INP, remove_comments='#', remove_whitespace=True)
    ioprinter.reading('run.dat...', newline=1)  # Add a Found <file> to msg

    # Read all of the main blocks as strings
    inp_block = ioformat.ptt.end_block(run_str, 'input', footer='input')
    pes_block = ioformat.ptt.end_block(run_str, 'pes', footer='pes')
    spc_block = ioformat.ptt.end_block(run_str, 'spc', footer='spc')
    es_block = ioformat.ptt.end_block(run_str, 'els', footer='els')
    trans_block = ioformat.ptt.end_block(run_str, 'transport', footer='transport')
    therm_tsks_block = ioformat.ptt.end_block(run_str, 'thermo', footer='thermo')
    ktp_tsks_block = ioformat.ptt.end_block(run_str, 'ktp', footer='ktp')
    proc_tsks_block = ioformat.ptt.end_block(run_str, 'proc', footer='proc')

    # Parse information in the blocks
    inp_dct = ioformat.ptt.keyword_dct_from_block(inp_block)
    # automol.util.dict_.right_update(def_dct, inp_dct)

    pes_idx_dct = _pes_idxs(pes_block)
    spc_idx_dct = _spc_idxs(spc_block)
    # hard to check these, need to crossref the mechanism

    print('es\n', es_block)
    es_tsks_lst = tsks.tsk_lst(es_block, thy_dct)
    # _new_check_dct(mod_tsk_lst, TSK_KEY_DCT, TSK_VAL_DCT, thy_dct)
    # therm_tsks_lst = tsks.tsk_lst2(trans_block, kin_mod_dct, spc_mod_dct)
    # ktp_tsks_lst = tsks.tsk_lst2(trans_block, kin_mod_dct, spc_mod_dct)
    # trans_tsks_lst = tsks.tsk_lst(trans_block, thy_dct)
    # proc_tsks_lst = tsks.tsk_lst(proc_block, thy_dct)

    print('es\n', es_tsks_lst)
    # print('trans\n', trans_tsks_lst)
    # print('therm\n', therm_tsks_lst)
    # print('ktp\n', ktp_tsks_lst)
    # print('proc\n', proc_tsks_lst)

    return inp_key_dct, pes_idx_dct, spc_idx_dct


def read_thy(job_path):
    """ Parse the theory.dat file
    """

    # Read the theory.dat input string
    thy_str = ioformat.pathtools.read_file(
        job_path, THY_INP, remove_comments='#', remove_whitespace=True)
    ioprinter.reading('theory.dat...', newline=1)

    # Format input to keyword-value theory dcts
    thy_blocks = ioformat.ptt.named_end_blocks(
        thy_str, 'level', footer='level')
    thy_dct = ioformat.ptt.keyword_dcts_from_blocks(thy_blocks)

    # Assess if the theory.dat input is valid
    for dct in thy_dct.values():
        check_val_dictionary2(dct, THY_VAL_DCT, 'thy.dat')

    return thy_dct


def read_model(job_path):
    """ Parse the models.dat file
    """

    # Read the models.dat input string
    mod_str = ioformat.pathtools.read_file(
        job_path, MODEL_INP, remove_comments='#', remove_whitespace=True)
    ioprinter.reading('model.dat...', newline=1)

    # Format the models input to the kin and spc model dcts
    kin_blocks = ioformat.ptt.named_end_blocks(mod_str, 'kin')
    spc_blocks = ioformat.ptt.named_end_blocks(mod_str, 'spc')

    kin_mod_dct = automol.util.dict_.merge_subdct(
        ioformat.ptt.keyword_dcts_from_blocks(kin_blocks), keep_subdct=True)
    spc_mod_dct = automol.util.dict_.merge_subdct(
        ioformat.ptt.keyword_dcts_from_blocks(spc_blocks), keep_subdct=True)

    # Assess if the model.dat input is valid

    return kin_mod_dct, spc_mod_dct


def read_spc(job_path):
    """ Read each of the species input files:
            (1) species.csv: CSV file with basic info like names,inchis,mults
            (2) species.dat: 
            (3) *.xyz: XYZ-files with geometries

        :param job_path: directory path where the input file(s) exist
        :type job_path: str
        :rtype dict[str: dict]
    """

    # Read the species input files: .csv, .dat and .xyz files
    spc_str = ioformat.pathtools.read_file(job_path, CSV_INP)
    ioprinter.reading('species.csv...', newline=1)

    dat_str = ioformat.pathtools.read_file(
        job_path, DAT_INP, remove_comments='#')
    ioprinter.reading('species.dat...', newline=1)

    geo_dct = geometry_dictionary(job_path)
    ioprinter.reading('geom.xyzs...', newline=1)


    # Use all species input to build the full species dct for spc and ts
    spc_dct = mechanalyzer.parser.spc.build_spc_dct(spc_str, 'csv')

    dat_blocks = ioformat.ptt.named_end_blocks(dat_str, 'spc', footer='spc')
    dat_dct = ioformat.ptt.keyword_dcts_from_blocks(dat_blocks)

    mod_spc_dct = modify_spc_dct(spc_dct, dat_dct, geo_dct)

    # Assess if the species.dat information is valid
    check_d

    return mod_spc_dct


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


# formatters, dont know where to build this
def geometry_dictionary(job_path):
    """ read in dictionary of saved geometries
    """
    geom_path = os.path.join(job_path, 'data')
    geom_dct = {}
    for dir_path, _, file_names in os.walk(geom_path):
        for file_name in file_names:
            file_path = os.path.join(dir_path, file_name)
            if file_path.endswith('.xyz'):
                xyz_str = autofile.io_.read_file(file_path)
                geo = automol.geom.from_xyz_string(xyz_str)
                ich = automol.geom.inchi(geo)
                if ich in geom_dct:
                    print('Warning: Dupilicate xyz geometry for ', ich)
                geom_dct[ich] = geo

    return geom_dct
