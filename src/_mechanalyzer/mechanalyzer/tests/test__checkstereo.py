""" test mechanalyzer.parser.sort for different mechanisms in 'data/'
    using different sorting options
"""

import os
import tempfile
import numpy as np
from ioformat import pathtools
import chemkin_io.writer
from mechanalyzer.builder import sorter
from mechanalyzer.parser import mech as mparser
from mechanalyzer.parser import ckin_ as ckin_parser
from mechanalyzer.parser import new_spc as sparser

# Set Paths to test/data directory and output directory
CWD = os.path.dirname(os.path.realpath(__file__))
TMP_OUT = tempfile.mkdtemp()

# Set types for parsing mechanisms
SPC_TYPE = 'csv'
MECH_TYPE = 'chemkin'



def test__sort_with_input():
    """ sort by using the auxlilary input files to specify parameters
        TEST FAILS. REASON IS CANON_ENANT. TO FIX
    """

    # Read the mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_C2H4_mech.dat')
    sort_path = os.path.join(CWD, 'data', 'sort.dat')

    spc_str, mech_str, sort_str = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism
    isolate_spc, sort_lst, _ = mparser.parse_sort(sort_str)
    param_dct_sort, _, cmts_dct, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst, stereo_optns=False) # stereo should be true
    index = 0



def test__readwrite_thirdbody():
    """ test mechanalyzer.parser.sort
        Checks read/write of a small set of rxns involving third bodies
        TEST FAILS. REASON IS CANON_ENANT. TO FIX
    """


    # Read the mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'NUIG_species.csv')
    mech_path = os.path.join(CWD, 'data', 'NUIG_mechred.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism by PES - No Headers Included
    isolate_spc = []
    sort_lst = ['pes', 0]

    param_dct_sort, _, _, _, _= sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst, stereo_optns=False) #stereo should be true

# Helper function


def _read_files(spc_path, mech_path, sort_path):
    """ read file names
    """

    spc_str, mech_str, sort_str = '', '', ''

    if spc_path is not None:
        with open(spc_path, encoding='utf-8') as fobj:
            spc_str = fobj.read()
    if mech_path is not None:
        with open(mech_path, encoding='utf-8') as fobj:
            mech_str = fobj.read()
    if sort_path is not None:
        with open(sort_path, encoding='utf-8') as fobj:
            sort_str = fobj.read()

    return spc_str, mech_str, sort_str


if __name__ == '__main__':
    print('tests to be fixed')
    #test__sort_with_input()   
    #test__readwrite_thirdbody()    
