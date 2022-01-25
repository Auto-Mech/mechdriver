import automol
from mechroutines.pf import thermo
import numpy

CLOSED = ['[H][H]', 'C', 'CC(C)CCO', 'C=CCC', 'C#CC=C=C']
OPEN = ['[H]', '[CH3]', 'CC(C)CC[O]', 'CC[C](C)CO', '[CH]=CCC', 'C#C[C]=C=C']


#######################
#### Closed shell species #######
#######################

# closed_ichs = [automol.smiles.inchi(smi) for smi in CLOSED]
# print(closed_ichs)

CLOSED_ICHS = ['InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4', 'InChI=1S/C5H12O/c1-5(2)3-4-6/h5-6H,3-4H2,1-2H3', 'InChI=1S/C4H8/c1-3-4-2/h3H,1,4H2,2H3', 'InChI=1S/C5H4/c1-3-5-4-2/h1,5H,2H2']
['InChI=1S/H', 'InChI=1S/CH3/h1H3', 'InChI=1S/C5H11O/c1-5(2)3-4-6/h5H,3-4H2,1-2H3', 'InChI=1S/C5H11O/c1-3-5(2)4-6/h6H,3-4H2,1-2H3', 'InChI=1S/C4H7/c1-3-4-2/h1,3H,4H2,2H3', 'InChI=1S/C5H3/c1-3-5-4-2/h1H,2H2']


basic_refs = []
cbh0_refs = []
cbh1_refs = []
cbh2_refs = []

for ich in CLOSED_ICHS:
    basic_refs.append(thermo.heatform.get_basic(ich))
    cbh0_refs.append(thermo.heatform.get_cbhzed(ich))
    cbh1_refs.append(thermo.heatform.get_cbhone(ich))
    cbh2_refs.append(thermo.heatform.get_cbhtwo(ich))

# print(basic_refs)
# print(cbh0_refs)
# print(cbh1_refs)
# print(cbh2_refs)
 
CLOSED_BASIC = [(['InChI=1S/H2/h1H'], [1]), (['InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4'], [0.0, 1.0]), (['InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4', 'InChI=1S/H2O/h1H2'], [-5.0, 5.0, 1.0]), (['InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4'], [-4.0, 4.0]), (['InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4'], [-8.0, 5.0])]
CLOSED_CBH0 = [(['InChI=1S/H2/h1H'], [1]), (['InChI=1S/CH4/h1H4'], [1]), (['InChI=1S/CH4/h1H4', 'InChI=1S/H2O/h1H2', 'InChI=1S/H2/h1H'], [5, 1, -5.0]), (['InChI=1S/CH4/h1H4', 'InChI=1S/H2/h1H'], [4, -4.0]), (['InChI=1S/CH4/h1H4', 'InChI=1S/H2/h1H'], [5, -8.0])]
CLOSED_CBH1 = [(['InChI=1S/H2/h1H'], [1]), (['InChI=1S/CH4/h1H4'], [1]), (['InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/CH4O/c1-2/h2H,1H3', 'InChI=1S/CH4/h1H4'], [4, 1, -4]), (['InChI=1S/C2H4/c1-2/h1-2H2', 'InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/CH4/h1H4'], [1, 2, -2]), (['InChI=1S/C2H2/c1-2/h1-2H', 'InChI=1S/C2H4/c1-2/h1-2H2', 'InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/CH4/h1H4'], [1, 2, 1, -3])]
CLOSED_CBH2 = [(['InChI=1S/H2/h1H'], [1]), (['InChI=1S/CH4/h1H4'], [1]), (['InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3', 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', 'InChI=1S/C4H10/c1-4(2)3/h4H,1-3H3'], [-2, 1, 1, 1]), (['InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/C3H6/c1-3-2/h3H,1H2,2H3', 'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3'], [-1, 1, 1]), (['InChI=1S/C2H4/c1-2/h1-2H2', 'InChI=1S/C3H4/c1-3-2/h1H,2H3', 'InChI=1S/C3H4/c1-3-2/h1-2H2', 'InChI=1S/C3H6/c1-3-2/h3H,1H2,2H3', 'InChI=1S/C2H6/c1-2/h1-2H3'], [-1, 1, 1, 1, -1])]

assert(basic_refs == CLOSED_BASIC)
assert(cbh0_refs == CLOSED_CBH0)
assert(cbh1_refs == CLOSED_CBH1)
assert(cbh2_refs == CLOSED_CBH2)

#######################
##### radicals ########
#######################

# open_ichs = [automol.smiles.inchi(smi) for smi in OPEN]
# print(open_ichs)

OPEN_ICHS = ['InChI=1S/H', 'InChI=1S/CH3/h1H3', 'InChI=1S/C5H11O/c1-5(2)3-4-6/h5H,3-4H2,1-2H3', 'InChI=1S/C5H11O/c1-3-5(2)4-6/h6H,3-4H2,1-2H3', 'InChI=1S/C4H7/c1-3-4-2/h1,3H,4H2,2H3', 'InChI=1S/C5H3/c1-3-5-4-2/h1H,2H2']

basic_refs = []
cbh0_refs = []
cbh1_refs = []
cbh2_refs = []

for ich in OPEN_ICHS:
    basic_refs.append(thermo.heatform.get_basic(ich))
    cbh0_refs.append(thermo.heatform.get_cbhzed(ich))
    cbh1_refs.append(thermo.heatform.get_cbhone(ich))
    cbh2_refs.append(thermo.heatform.get_cbhtwo(ich))

# print(basic_refs)
# print(cbh0_refs)
# print(cbh1_refs)
# print(cbh2_refs)

OPEN_BASIC = [(['InChI=1S/H2/h1H'], [0.5]), (['InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4'], [-0.5, 1.0]), (['InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4', 'InChI=1S/H2O/h1H2'], [-5.5, 5.0, 1.0]), (['InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4', 'InChI=1S/H2O/h1H2'], [-5.5, 5.0, 1.0]), (['InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4'], [-4.5, 4.0]), (['InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4'], [-8.5, 5.0])]
OPEN_CBH0 = [(['InChI=1S/H'], [1]), (['InChI=1S/CH3/h1H3'], [1]), (['InChI=1S/CH4/h1H4', 'InChI=1S/HO/h1H', 'InChI=1S/H2/h1H'], [5, 1, -5.0]), (['InChI=1S/CH4/h1H4', 'InChI=1S/CH3/h1H3', 'InChI=1S/H2O/h1H2', 'InChI=1S/H2/h1H'], [4, 1, 1, -5.0]), (['InChI=1S/CH3/h1H3', 'InChI=1S/CH4/h1H4', 'InChI=1S/H2/h1H'], [1, 3, -4.0]), (['InChI=1S/CH4/h1H4', 'InChI=1S/CH3/h1H3', 'InChI=1S/H2/h1H'], [4, 1, -8.0])]
OPEN_CBH1 = [(['InChI=1S/H'], [1]), (['InChI=1S/CH3/h1H3'], [1]), (['InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/CH3O/c1-2/h1H3', 'InChI=1S/CH4/h1H4'], [4, 1, -4]), (['InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/C2H5/c1-2/h1H2,2H3', 'InChI=1S/CH4O/c1-2/h2H,1H3', 'InChI=1S/CH4/h1H4', 'InChI=1S/CH3/h1H3'], [1, 3, 1, -2, -2]), (['InChI=1S/C2H3/c1-2/h1H,2H2', 'InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/CH4/h1H4'], [1, 2, -2]), (['InChI=1S/C2H2/c1-2/h1-2H', 'InChI=1S/C2H4/c1-2/h1-2H2', 'InChI=1S/C2H5/c1-2/h1H2,2H3', 'InChI=1S/C2H3/c1-2/h1H,2H2', 'InChI=1S/CH4/h1H4', 'InChI=1S/CH3/h1H3'], [1, 1, 1, 1, -2, -1])]
OPEN_CBH2 = [(['InChI=1S/H'], [1]), (['InChI=1S/CH3/h1H3'], [1]), (['InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3', 'InChI=1S/C2H5O/c1-2-3/h2H2,1H3', 'InChI=1S/C4H10/c1-4(2)3/h4H,1-3H3'], [-2, 1, 1, 1]), (['InChI=1S/C2H5/c1-2/h1H2,2H3', 'InChI=1S/C3H7/c1-3-2/h1,3H2,2H3', 'InChI=1S/C2H5O/c1-2-3/h3H,1-2H2', 'InChI=1S/C4H9/c1-4(2)3/h1-3H3'], [-2, 1, 1, 1]), (['InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/C3H5/c1-3-2/h1,3H,2H3', 'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3'], [-1, 1, 1]), (['InChI=1S/C3H3/c1-3-2/h1H,2H2', 'InChI=1S/C3H5/c1-3-2/h1H2,2H3', 'InChI=1S/C2H5/c1-2/h1H2,2H3', 'InChI=1S/C2H3/c1-2/h1H,2H2'], [2, 1, -1, -1])]

assert(basic_refs == OPEN_BASIC)
assert(cbh0_refs == OPEN_CBH0)
assert(cbh1_refs == OPEN_CBH1)
assert(cbh2_refs == OPEN_CBH2)


#######################
### Abstraction ####
#######################

# abs_rct_smis = ['CC(C)COO', 'C[CH]C']
# abs_prd_smis = ['CC(C)CO[O]', 'CCC']
# abs_rct_ichs = [automol.smiles.inchi(smi) for smi in abs_rct_smis]
# abs_prd_ichs = [automol.smiles.inchi(smi) for smi in abs_prd_smis]
# zrxn = automol.reac.rxn_objs_from_inchi(
#     abs_rct_ichs, abs_prd_ichs, indexing='zma')[0][0]
# 
# print(zrxn)

ABS_ZRXN_STR = '''reaction class: hydrogen abstraction
forward TS atoms:
  1: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  2: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  4: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  6: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  7: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  8: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  9: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  10: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  11: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  12: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  13: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  14: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  15: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  16: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  17: {symbol: X, implicit_hydrogen_valence: 0, stereo_parity: null}
  18: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  19: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  20: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  21: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  22: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  23: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  24: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  25: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  26: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  27: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
forward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  1-4: {order: 1, stereo_parity: null}
  1-5: {order: 1, stereo_parity: null}
  2-6: {order: 1, stereo_parity: null}
  2-7: {order: 1, stereo_parity: null}
  2-8: {order: 1, stereo_parity: null}
  6-9: {order: 1, stereo_parity: null}
  6-10: {order: 1, stereo_parity: null}
  6-11: {order: 1, stereo_parity: null}
  7-12: {order: 1, stereo_parity: null}
  7-13: {order: 1, stereo_parity: null}
  7-14: {order: 1, stereo_parity: null}
  12-15: {order: 1, stereo_parity: null}
  15-16: {order: 0.9, stereo_parity: null}
  16-17: {order: 0, stereo_parity: null}
  16-18: {order: 0.1, stereo_parity: null}
  18-19: {order: 1, stereo_parity: null}
  18-20: {order: 1, stereo_parity: null}
  18-21: {order: 1, stereo_parity: null}
  19-22: {order: 1, stereo_parity: null}
  19-23: {order: 1, stereo_parity: null}
  19-24: {order: 1, stereo_parity: null}
  20-25: {order: 1, stereo_parity: null}
  20-26: {order: 1, stereo_parity: null}
  20-27: {order: 1, stereo_parity: null}
reactants keys:
- [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
- [18, 19, 20, 21, 22, 23, 24, 25, 26, 27]
backward TS atoms:
  1: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  2: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  3: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  4: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  6: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  7: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  8: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  9: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  10: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  11: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  12: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  13: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  14: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  15: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  16: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  17: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  18: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  19: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  20: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  21: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  22: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  23: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  24: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  25: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  26: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
backward TS bonds:
  1-3: {order: 1, stereo_parity: null}
  1-4: {order: 1, stereo_parity: null}
  1-5: {order: 1, stereo_parity: null}
  1-6: {order: 1, stereo_parity: null}
  2-3: {order: 1, stereo_parity: null}
  2-7: {order: 1, stereo_parity: null}
  2-8: {order: 1, stereo_parity: null}
  2-9: {order: 1, stereo_parity: null}
  3-10: {order: 1, stereo_parity: null}
  3-11: {order: 0.9, stereo_parity: null}
  11-16: {order: 0.1, stereo_parity: null}
  12-15: {order: 1, stereo_parity: null}
  12-18: {order: 1, stereo_parity: null}
  12-19: {order: 1, stereo_parity: null}
  12-20: {order: 1, stereo_parity: null}
  13-15: {order: 1, stereo_parity: null}
  13-21: {order: 1, stereo_parity: null}
  13-22: {order: 1, stereo_parity: null}
  13-23: {order: 1, stereo_parity: null}
  14-15: {order: 1, stereo_parity: null}
  14-17: {order: 1, stereo_parity: null}
  14-24: {order: 1, stereo_parity: null}
  14-25: {order: 1, stereo_parity: null}
  15-26: {order: 1, stereo_parity: null}
  16-17: {order: 1, stereo_parity: null}
products keys:
- [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
- [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]

'''

zrxn = automol.reac.from_string(ABS_ZRXN_STR)
#
##basic_refs = thermo.heatform.get_basic_ts(zrxn)
basic_refs= [], []
cbh0_refs = thermo.heatform.get_cbhzed_ts(zrxn)
cbh1_refs = thermo.heatform.get_cbhone_ts(zrxn)

print(cbh0_refs)
print(cbh1_refs)

ABS_CBH0_REFS = (['InChI=1S/CH4/h1H4', 'InChI=1S/H2O/h1H2', (('InChI=1S/CH3/h1H3', 'InChI=1S/H2O/h1H2'), ('InChI=1S/CH4/h1H4', 'InChI=1S/HO/h1H')), 'InChI=1S/H2/h1H'], [6.0, 1.0, 1.0, -7.0])
ABS_CBH1_REFS = (['InChI=1S/CH4O/c1-2/h2H,1H3', 'InChI=1S/C2H6/c1-2/h1-2H3', (('InChI=1S/CH3/h1H3', 'InChI=1S/H2O2/c1-2/h1-2H'), ('InChI=1S/CH4/h1H4', 'InChI=1S/HO2/c1-2/h1H')), (('InChI=1S/C2H5/c1-2/h1H2,2H3', 'InChI=1S/H2O/h1H2'), ('InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/HO/h1H')), 'InChI=1S/CH4/h1H4', 'InChI=1S/H2O/h1H2', (('InChI=1S/CH3/h1H3', 'InChI=1S/H2O/h1H2'), ('InChI=1S/CH4/h1H4', 'InChI=1S/HO/h1H'))], [1.0, 3.0, 1.0, 2.0, -3, -1, -2])

assert(cbh0_refs == ABS_CBH0_REFS)
assert(cbh1_refs == ABS_CBH1_REFS)

#####################
### Beta scission #####
####################

# bs_rct_smis = ['CC[CH]CO']
# bs_prd_smis = ['CCC=C', '[OH]']
# bs_rct_ichs = [automol.smiles.inchi(smi) for smi in bs_rct_smis]
# bs_prd_ichs = [automol.smiles.inchi(smi) for smi in bs_prd_smis]
# zrxn = automol.reac.rxn_objs_from_inchi(
#     bs_rct_ichs, bs_prd_ichs, indexing='zma')[0][0]
# 
# print(zrxn)

BS_ZRXN_STR = '''reaction class: beta scission
forward TS atoms:
  1: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  2: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  4: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  6: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  7: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  8: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  9: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  10: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  11: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  12: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  13: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  14: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
forward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  1-4: {order: 1, stereo_parity: null}
  1-5: {order: 1, stereo_parity: null}
  2-6: {order: 1, stereo_parity: null}
  2-7: {order: 1, stereo_parity: null}
  2-8: {order: 1, stereo_parity: null}
  6-9: {order: 1, stereo_parity: null}
  6-10: {order: 1, stereo_parity: null}
  9-11: {order: 0.9, stereo_parity: null}
  9-12: {order: 1, stereo_parity: null}
  9-13: {order: 1, stereo_parity: null}
  11-14: {order: 1, stereo_parity: null}
reactants keys:
- [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
backward TS atoms:
  1: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  2: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  3: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  4: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  6: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  7: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  8: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  9: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  10: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  11: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  12: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  13: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  14: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
backward TS bonds:
  1-3: {order: 1, stereo_parity: null}
  1-5: {order: 1, stereo_parity: null}
  1-6: {order: 1, stereo_parity: null}
  1-13: {order: 0.1, stereo_parity: null}
  2-4: {order: 1, stereo_parity: null}
  2-7: {order: 1, stereo_parity: null}
  2-8: {order: 1, stereo_parity: null}
  2-9: {order: 1, stereo_parity: null}
  3-4: {order: 1, stereo_parity: null}
  3-10: {order: 1, stereo_parity: null}
  4-11: {order: 1, stereo_parity: null}
  4-12: {order: 1, stereo_parity: null}
  13-14: {order: 1, stereo_parity: null}
products keys:
- [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
- [13, 14]
'''

zrxn = automol.reac.from_string(BS_ZRXN_STR)

cbh0_refs = thermo.heatform.get_cbhzed_ts(zrxn)
cbh1_refs = thermo.heatform.get_cbhone_ts(zrxn)

# print(cbh0_refs)
# print(cbh1_refs)

CBH0_REFS = (['InChI=1S/CH4/h1H4', (('InChI=1S/C2H5O/c1-2-3/h3H,1-2H2',), ('InChI=1S/C2H4/c1-2/h1-2H2', 'InChI=1S/HO/h1H')), 'InChI=1S/H2/h1H'], [2.0, 1.0, -2.0])
CBH1_REFS = (['InChI=1S/C2H6/c1-2/h1-2H3', (('InChI=1S/C3H7O/c1-2-3-4/h2,4H,3H2,1H3',), ('InChI=1S/C3H6/c1-3-2/h3H,1H2,2H3', 'InChI=1S/HO/h1H')), 'InChI=1S/CH4/h1H4'], [1.0, 1.0, -1])

assert(cbh0_refs == CBH0_REFS)
assert(cbh1_refs == CBH1_REFS)

