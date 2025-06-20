""" test mechanalyzer.parser.pes
"""

import os
import tempfile
import numpy as np
from ioformat import pathtools
from mechanalyzer.parser.pes import pes_dictionary
from mechanalyzer.parser.spc import build_spc_dct

CWD = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(CWD, 'data')
# Set types for parsing mechanisms
SPC_TYPE = 'csv'
MECH_TYPE = 'chemkin'

ref_csv_str = pathtools.read_file(DAT_PATH, 'species_forpesdct.csv')
spc_dct = build_spc_dct(ref_csv_str, 'csv')
mech_str1 = pathtools.read_file(DAT_PATH, 'mech_removecmts.dat')
mech_str2 = pathtools.read_file(DAT_PATH, 'mech_forcedcmts.dat')

results_pes_dct = {('C4H7O', 0, 0): ((0, (('C4H7ORvE4fmAB0',), ('C4H7O4H74fm1',), (None,))), 
                    (1, (('C4H7ORvE4fmAB0',), ('C4H7O-kSV4fm',), (None,))), 
                    (2, (('C4H7ORvE4fmAB0',), ('C4H6O-RvErx51', 'H-TcYTcY'), (None,))), 
                    (3, (('C4H7ORvE4fmAA0',), ('C4H7O4H74fm0',), (None,))),
                    (4, (('C4H7ORvE4fmAA0',), ('C4H7O-kSV4fm',), (None,))), 
                    (5, (('C4H7ORvE4fmAA0',), ('C4H6O-RvErx50', 'H-TcYTcY'), (None,))), 
                    (6, (('C4H7O4H74fm0',), ('C3H4OALAD-Wv9FbZ', 'CH3'), (None,))), 
                    (7, (('C4H7O4H74fm0',), ('C2H4OALD-UPQWKw', 'C2H3ALK-S58hH1'), (None,))), 
                    (8, (('C4H7O-kSV4fm',), ('C2H4OALD-UPQWKw', 'C2H3ALK-S58hH1'), (None,)))), 
     ('C4H9O2', 1, 0): ((0, (('C4H8ORvEsWvAA0', 'OH'), ('C4H7ORvE4fmAA0', 'H2O'), (None,))),), 
     ('C4H9O2', 1, 1): ((1, (('C4H8ORvEsWvAB', 'OH'), ('C4H7ORvE4fmAB0', 'H2O'), (None,))),),
     ('C4H9O3', 2, 0): ((0, (('C4H8ORvEsWvAA0', 'HO2-S580KW'), ('C4H7ORvE4fmAA0', 'H2O2-S58pAY'), (None,))),), 
     ('C5H11O', 3, 0): ((0, (('C4H8ORvEsWvAA0', 'CH3'), ('C4H7ORvE4fmAA0', 'CH4'), (None,))),), 
     ('C5H11O3', 4, 0): ((0, (('C4H8ORvEsWvAA0', 'CH3O2RO2-2LTcwB'), ('C4H7ORvE4fmAA0', 'CH4O2QOOH-2LTWKw'), (None,))),), 
     ('C5H11O2', 5, 0): ((0, (('C4H8ORvEsWvAA0', 'CH3O-S58cwB'), ('C4H7ORvE4fmAA0', 'CH4O-S58WKw'), (None,))),), 
     ('C4H8ClO', 6, 0): ((0, (('C4H8ORvEsWvAA0', 'Cl'), ('C4H7ORvE4fmAA0', 'HCl'), (None,))),)}

def test__pes_dictionary():
    pes_dct = pes_dictionary(mech_str1, spc_dct)
    print('pes_dict will be None: the ! pes.subpes.channel comments are removed')
    assert pes_dct == None
    print('pes_dict will be assigned: the # pes.subpes.channel comments are not removed')
    pes_dct = pes_dictionary(mech_str2, spc_dct)
    for key, val in pes_dct.items():
        assert val == results_pes_dct[key]

    
if __name__ == '__main__':
    test__pes_dictionary()
    #test__connected_channels_dct() #calls also find_conn_chnls
    #test__print_pes_channels()
    
    
    
