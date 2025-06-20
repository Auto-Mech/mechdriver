""" replace reaction parameters in existing mechanism without rewriting the whole file
"""

import os
import ioformat
import numpy as np
from autoreact.params import RxnParams
from chemkin_io.parser.mechanism import reaction_block
from chemkin_io.parser.reaction import get_rxn_param_dct, get_rxn_cmt_dct, get_rxn_strs_dct
from chemkin_io.writer.reaction import replace_rxn_in_cki


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')

# Define stuff for testing PLOG
PLOG_DCT = {
    0.1: [[1E+15, 0.00, 28000]],
    1.0: [[1E+15, 0.00, 25000], [2E+15, 0.00, 27000]],
    100.0: [[1E+15, 0.00, 25000]]}
PLOG_PARAMS = RxnParams(plog_dct=PLOG_DCT)

keys_to_replace = [
    (('C4H2', 'C6H5'), ('C10H7',), (None,)), #dup
    (('C6H5', 'C6H4C2H'), ('C14H10',), (None,)), #simple
    (('CH3C6H4', 'C6H5C2H5'), ('C6H5C2H4C6H5', 'CH3'), (None,)), # written as irrev, but with no bw
]
params = [PLOG_PARAMS]*len(keys_to_replace)
rxn_param_dct_new = dict(zip(keys_to_replace, params))

# add random reaction to test add option:
rxn_param_dct_new[(('A','B'), ('C',), (None,))] = PLOG_PARAMS
keys0 = keys_to_replace + [(('A','B'), ('C',), (None,))]

def test_replace():
    # read existing mechanism
    ckin_str = ioformat.pathtools.read_file(DAT_PATH, 'pah_block.dat')
    rxn_block_comments = reaction_block(ckin_str, remove_comments=False)
    rxn_block_nocomments = reaction_block(ckin_str)
    rxn_param_dct_old = get_rxn_param_dct(rxn_block_nocomments, 'cal/mole', 'moles')
    
    rxn_strs_dct = get_rxn_strs_dct(rxn_block_comments)
    rxn_cmts_dct = get_rxn_cmt_dct(rxn_block_comments)

    # update rxn params for a "simple" and a "duplicate" reaction
    ckin_str_new = replace_rxn_in_cki(
        rxn_block_comments, rxn_strs_dct, rxn_param_dct_new, rxn_cmts_dct)

    # extract again the dictionaries and check that the parameters are the same
    ckin_str_new_nocmt = reaction_block('REACTIONS\n' + ckin_str_new + '\nEND\n')
    rxn_param_dct_rewritten = get_rxn_param_dct(ckin_str_new_nocmt, 'cal/mole', 'moles')
    for key, val in rxn_param_dct_rewritten.items():
        if key in keys0:
            plog_dct = val.plog
            for pressure, arr_tuples in plog_dct.items():
                for i, arr_tuple in enumerate(arr_tuples):
                    assert np.allclose(arr_tuple, rxn_param_dct_new[key].plog[pressure][i])
        else:
            plog_dct = val.plog
            arr_dct = val.arr
            if plog_dct:
                for pressure, arr_tuples in plog_dct.items():
                    for i, arr_tuple in enumerate(arr_tuples):
                        assert np.allclose(arr_tuple, rxn_param_dct_old[key].plog[pressure][i])
            if arr_dct:
                for arr_tuple in arr_dct:
                    assert np.allclose(arr_tuple, rxn_param_dct_old[key].arr[0])

if __name__ == '__main__':
    test_replace()

