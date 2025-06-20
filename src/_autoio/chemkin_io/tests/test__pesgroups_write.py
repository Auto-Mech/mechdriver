""" Tests the writer functions in mechanism.py pertaining to species
"""

from chemkin_io.writer import pesgroups

# Define mechanism objects
PES_GROUPS_DCT_LST = [
    {'grp': 1, 'idxs': ['2:1', '2:2', '2:3', '1:1'], 
     'peds': [['C4H8-1+H=C4H71-3+H2'], ['C4H8-1+H=C4H71-4+H2'], ['C4H8-2+H=C4H71-3+H2'], []], 
     'hot': [[], [], [], ['C4H71-3', 'C4H71-4']]}, 
    {'grp': 2, 'idxs': ['3:1', '3:2', '3:3', '1:1'], 
     'peds': [['C4H8-1+O=C4H71-3+OH'], ['C4H8-1+O=C4H71-4+OH'], ['C4H8-2+O=C4H71-3+OH'], []], 
     'hot': [[], [], [], ['C4H71-3', 'C4H71-4']]}, 
    {'grp': 3, 'idxs': ['4:1', '4:2', '4:3', '1:1'], 
     'peds': [['C4H8-1+OH=C4H71-3+H2O'], ['C4H8-1+OH=C4H71-4+H2O'], ['C4H8-2+OH=C4H71-3+H2O'], []], 
     'hot': [[], [], [], ['C4H71-3', 'C4H71-4']]}, 
    {'grp': 4, 'idxs': ['5:1', '5:2', '5:3', '1:1'], 
     'peds': [['C4H8-1+CH3=C4H71-3+CH4'], ['C4H8-1+CH3=C4H71-4+CH4'], ['C4H8-2+CH3=C4H71-3+CH4'], []], 
     'hot': [[], [], [], ['C4H71-3', 'C4H71-4']]}
    ]

PES_GROUPS_STR = (
    'grp 1 \n'
    '\t idxs = [2:1, 2:2, 2:3, 1:1] \n'
    "\t peds = [['C4H8-1+H=C4H71-3+H2'], ['C4H8-1+H=C4H71-4+H2'], ['C4H8-2+H=C4H71-3+H2'], []] \n"
    "\t hot = [[], [], [], ['C4H71-3', 'C4H71-4']] \n"
    '\t modeltype = equip_phi \n'
    '\t bf_threshold = 0.1 \n'
    'end grp \n\n'
    'grp 2 \n'
    '\t idxs = [3:1, 3:2, 3:3, 1:1] \n'
    "\t peds = [['C4H8-1+O=C4H71-3+OH'], ['C4H8-1+O=C4H71-4+OH'], ['C4H8-2+O=C4H71-3+OH'], []] \n"
    "\t hot = [[], [], [], ['C4H71-3', 'C4H71-4']] \n"
    '\t modeltype = equip_phi \n'
    '\t bf_threshold = 0.1 \n'
    'end grp \n\n'
    'grp 3 \n'
    '\t idxs = [4:1, 4:2, 4:3, 1:1] \n'
    "\t peds = [['C4H8-1+OH=C4H71-3+H2O'], ['C4H8-1+OH=C4H71-4+H2O'], ['C4H8-2+OH=C4H71-3+H2O'], []] \n"
    "\t hot = [[], [], [], ['C4H71-3', 'C4H71-4']] \n"
    '\t modeltype = equip_phi \n'
    '\t bf_threshold = 0.1 \n'
    'end grp \n\n'
    'grp 4 \n'
    '\t idxs = [5:1, 5:2, 5:3, 1:1] \n'
    "\t peds = [['C4H8-1+CH3=C4H71-3+CH4'], ['C4H8-1+CH3=C4H71-4+CH4'], ['C4H8-2+CH3=C4H71-3+CH4'], []] \n"
    "\t hot = [[], [], [], ['C4H71-3', 'C4H71-4']] \n"
    '\t modeltype = equip_phi \n'
    '\t bf_threshold = 0.1 \n'
    'end grp \n\n'
)



def test_pesgroups():
    """ Tests pes groups writing
    """
    out_str = pesgroups.write_pes_groups(PES_GROUPS_DCT_LST, modeltype='equip_phi')
    assert out_str == PES_GROUPS_STR


if __name__ == '__main__':
    test_pesgroups()

