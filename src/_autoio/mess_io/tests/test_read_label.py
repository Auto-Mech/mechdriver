""" test mess_io.reader.pes
"""

import os
from ioformat import pathtools
from mess_io.reader import relabel
from mess_io.reader import name_label_dct



PATH = os.path.dirname(os.path.realpath(__file__))
OPATH = os.path.join(PATH, 'data', 'out')
HOT_LOG_DBL = pathtools.read_file(OPATH, 'me_ktp_hoten_c3h7.logf')
HOT_KE_OLDFMT = pathtools.read_file(OPATH, 'ke_ped_c3h8_h.out')
LBL_DCT = {'W1': 'CH3CH2CH2', 'W2': 'CH3CHCH3', 
           'P1': 'C2H4+CH3', 'P2': 'CH3CHCH2+H'}
LBL_DCT2 = {'W0': 'W0', 'P0': 'C3H8+H',
            'P1': 'CH3CH2CH2+H2', 'P2': 'CH3CHCH3+H2'}

def test_name_label_dct():
    """ test mess_io._label.name_label_dct
    """
    lbl_dct = name_label_dct(HOT_LOG_DBL)
    assert lbl_dct == LBL_DCT
    lbl_dct2 = name_label_dct(HOT_KE_OLDFMT)
    assert lbl_dct2 == LBL_DCT2
    
def test_relabel():
    """ test mess_io.relabel
    """
    DCT0 = {(('W0',), ('P0',), (None,)): {},
            (('W0',), ('P1',), (None,)): {},
            (('W0',), ('P2',), (None,)): {},
            (('P0',), ('W0',), (None,)): {},
            (('P0',), ('P1',), (None,)): {},
            (('P0',), ('P2',), (None,)): {},
            (('P1',), ('W0',), (None,)): {},
            (('P1',), ('P0',), (None,)): {},
            (('P1',), ('P2',), (None,)): {},
            (('P2',), ('W0',), (None,)): {},
            (('P2',), ('P0',), (None,)): {},
            (('P2',), ('P1',), (None,)): {},
            }
    NEWKEYS = [(('W0',), ('C3H8', 'H'), (None,)), 
               (('W0',), ('CH3CH2CH2', 'H2'), (None,)), 
               (('W0',), ('CH3CHCH3', 'H2'), (None,)), 
               (('C3H8', 'H'), ('W0',), (None,)), 
               (('C3H8', 'H'), ('CH3CH2CH2', 'H2'), (None,)), 
               (('C3H8', 'H'), ('CH3CHCH3', 'H2'), (None,)), 
               (('CH3CH2CH2', 'H2'), ('W0',), (None,)), 
               (('CH3CH2CH2', 'H2'), ('C3H8', 'H'), (None,)), 
               (('CH3CH2CH2', 'H2'), ('CH3CHCH3', 'H2'), (None,)), 
               (('CH3CHCH3', 'H2'), ('W0',), (None,)), 
               (('CH3CHCH3', 'H2'), ('C3H8', 'H'), (None,)), 
               (('CH3CHCH3', 'H2'), ('CH3CH2CH2', 'H2'), (None,))]
    
    new_dct = relabel(DCT0, LBL_DCT2)
    assert NEWKEYS == list(new_dct.keys())
    
# def test_relabel
if __name__ == '__main__':
    test_name_label_dct()
    test_relabel()
