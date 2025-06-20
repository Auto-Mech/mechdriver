""" test thermfit.cbh._basic
"""

import numpy
import automol.reac
import thermfit.cbh


# Closed-Shell Species
SPC_CLOSED_ICHS = (
    'InChI=1S/H2/h1H',                                  # [HH]
    'InChI=1S/CH4/h1H4',                                # C
    'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3',                  # CCC
    'InChI=1S/C4H8/c1-3-4-2/h3H,1,4H2,2H3',             # C=CCC
    'InChI=1S/C5H4/c1-3-5-4-2/h1,5H,2H2',               # C#CC=C=C
    'InChI=1S/C5H12O/c1-5(2)3-4-6/h5-6H,3-4H2,1-2H3',   # CC(C)CCO
)
# Closed-Shell Species with less validated heteroatoms (for higher CBH-N)
SPC_CLOSED_ICHS_2 = (
    'InChI=1S/C3H9N/c1-2-3-4/h2-4H2,1H3',               # CCCN
    'InChI=1S/C3H7Cl/c1-2-3-4/h2-3H2,1H3',              # CCCCl
    'InChI=1S/C3H8S/c1-3-4-2/h3H2,1-2H3',               # CCSC
    'InChI=1S/C3H8OS/c1-5-3-2-4/h4H,2-3H2,1H3',         # OCCSC
)

# Open-Shell Species
SPC_OPEN_ICHS = (
    'InChI=1S/H',                                       # [H]
    'InChI=1S/CH3/h1H3',                                # [CH3]
    'InChI=1S/C4H7/c1-3-4-2/h1,3H,4H2,2H3',             # [CH]=CCC
    'InChI=1S/C5H3/c1-3-5-4-2/h1H,2H2',                 # C#C[C]=C=C
    'InChI=1S/C5H11O/c1-5(2)3-4-6/h5H,3-4H2,1-2H3',     # CC(C)CC[O]
    'InChI=1S/C5H11O/c1-3-5(2)4-6/h6H,3-4H2,1-2H3',     # CC[C](C)CO
)

# Transition States
HABS_ZRXN_STR = """
reaction class: hydrogen abstraction
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
"""

BS_ZRXN_STR = """
reaction class: beta scission
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
"""

HABS_ZRXN = automol.reac.from_old_string(HABS_ZRXN_STR)
BS_ZRXN = automol.reac.from_old_string(BS_ZRXN_STR)

# Set scheme variables
BASIC_SCHEME = 'basic'
CBH0_SCHEME = 'cbh0'
CBH1_SCHEME = 'cbh1'
CBH2_SCHEME = 'cbh2'


def test__closed_shell_species():
    """ test thermfit.cbh._spc
    """

    # Tests for species with bases validated up to CBH-N
    basic, cbh0, cbh1, cbh2 = (), (), (), ()
    for ich in SPC_CLOSED_ICHS:
        basic += (thermfit.cbh.species_basis(ich, BASIC_SCHEME),)
        cbh0 += (thermfit.cbh.species_basis(ich, CBH0_SCHEME),)
        cbh1 += (thermfit.cbh.species_basis(ich, CBH1_SCHEME),)
        cbh2 += (thermfit.cbh.species_basis(ich, CBH2_SCHEME),)

    ref_basic = (
        (('InChI=1S/H2/h1H',), (1,)),
        (('InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4'), (0.0, 1.0)),
        (('InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4'), (-2.0,  3.0)),
        (('InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4'), (-4.0, 4.0)),
        (('InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4'), (-8.0, 5.0)),
        (('InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4',
          'InChI=1S/H2O/h1H2'), (-5.0, 5.0, 1.0))
    )
    ref_cbh0 = (
        (('InChI=1S/H2/h1H',), (1,)),
        (('InChI=1S/CH4/h1H4',), (1,)),
        (('InChI=1S/CH4/h1H4', 'InChI=1S/H2/h1H'), (3, -2.0)),
        (('InChI=1S/CH4/h1H4', 'InChI=1S/H2/h1H'), (4, -4.0)),
        (('InChI=1S/CH4/h1H4', 'InChI=1S/H2/h1H'), (5, -8.0)),
        (('InChI=1S/CH4/h1H4', 'InChI=1S/H2O/h1H2',
          'InChI=1S/H2/h1H'), (5, 1, -5.0))
    )
    ref_cbh1 = (
        (('InChI=1S/H2/h1H',), (1,)),
        (('InChI=1S/CH4/h1H4',), (1,)),
        (('InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/CH4/h1H4'), (2, -1)),
        (('InChI=1S/C2H4/c1-2/h1-2H2', 'InChI=1S/C2H6/c1-2/h1-2H3',
          'InChI=1S/CH4/h1H4'), (1, 2, -2)),
        (('InChI=1S/C2H4/c1-2/h1-2H2', 'InChI=1S/C2H6/c1-2/h1-2H3',
          'InChI=1S/C2H2/c1-2/h1-2H', 'InChI=1S/CH4/h1H4'), (2, 1, 1, -3)),
        (('InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/CH4O/c1-2/h2H,1H3',
          'InChI=1S/CH4/h1H4'), (4, 1, -4))
    )
    ref_cbh2 = (
        (('InChI=1S/H2/h1H',), (1,)),
        (('InChI=1S/CH4/h1H4',), (1,)),
        (('InChI=1S/C3H8/c1-3-2/h3H2,1-2H3',), (1,)),
        (('InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/C3H6/c1-3-2/h3H,1H2,2H3',
          'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3'), (-1, 1, 1)),
        (('InChI=1S/C2H4/c1-2/h1-2H2', 'InChI=1S/C3H4/c1-3-2/h1H,2H3',
          'InChI=1S/C3H4/c1-3-2/h1-2H2', 'InChI=1S/C3H6/c1-3-2/h3H,1H2,2H3',
          'InChI=1S/C2H6/c1-2/h1-2H3'), (-1, 1, 1, 1, -1)),
        (('InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3',
          'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
          'InChI=1S/C4H10/c1-4(2)3/h4H,1-3H3'), (-2, 1, 1, 1))
    )

    for basis, ref_basis in zip(basic, ref_basic):
        assert basis[0] == ref_basis[0]
        assert numpy.allclose(basis[1], ref_basis[1])
    for basis, ref_basis in zip(cbh0, ref_cbh0):
        assert basis[0] == ref_basis[0]
        assert numpy.allclose(basis[1], ref_basis[1])
    for basis, ref_basis in zip(cbh1, ref_cbh1):
        assert basis[0] == ref_basis[0]
        assert numpy.allclose(basis[1], ref_basis[1])
    for basis, ref_basis in zip(cbh2, ref_cbh2):
        assert basis[0] == ref_basis[0]
        assert numpy.allclose(basis[1], ref_basis[1])

    # Test basic for unique hetereoatom species
    basic2 = ()
    for ich in SPC_CLOSED_ICHS_2:
        basic2 += (thermfit.cbh.species_basis(ich, BASIC_SCHEME),)

    ref_basic2 = (
        (('InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4',
          'InChI=1S/H3N/h1H3'), numpy.array([-3.,  3.,  1.])),
        (('InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4',
          'InChI=1S/ClH/h1H'), numpy.array([-3.,  3.,  1.])),
        (('InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4',
          'InChI=1S/O2S/c1-3-2', 'InChI=1S/H2O/h1H2'),
         numpy.array([0.,  3.,  1., -2.])),
        (('InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4',
          'InChI=1S/H2O/h1H2', 'InChI=1S/O2S/c1-3-2'),
         numpy.array([-1.,  3., -1.,  1.]))
    )

    for basis, ref_basis in zip(basic2, ref_basic2):
        assert basis[0] == ref_basis[0]
        assert numpy.allclose(basis[1], ref_basis[1])


def test__open_shell_species():
    """ test thermfit.cbh._spc
    """

    basic, cbh0, cbh1, cbh2 = (), (), (), ()
    for ich in SPC_OPEN_ICHS:
        basic += (thermfit.cbh.species_basis(ich, BASIC_SCHEME),)
        cbh0 += (thermfit.cbh.species_basis(ich, CBH0_SCHEME),)
        cbh1 += (thermfit.cbh.species_basis(ich, CBH1_SCHEME),)
        cbh2 += (thermfit.cbh.species_basis(ich, CBH2_SCHEME),)

    ref_basic = (
        (('InChI=1S/H2/h1H',), (0.5)),
        (('InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4'), (-0.5, 1.0)),
        (('InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4'), (-4.5, 4.0)),
        (('InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4'), (-8.5, 5.0)),
        (('InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4',
          'InChI=1S/H2O/h1H2'), (-5.5, 5.0, 1.0)),
        (('InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4',
          'InChI=1S/H2O/h1H2'), (-5.5, 5.0, 1.0)),
    )
    ref_cbh0 = (
        (('InChI=1S/H',), (1,)),
        (('InChI=1S/CH3/h1H3',), (1,)),
        (('InChI=1S/CH3/h1H3', 'InChI=1S/CH4/h1H4',
          'InChI=1S/H2/h1H'), (1, 3, -4.0)),
        (('InChI=1S/CH4/h1H4', 'InChI=1S/CH3/h1H3',
          'InChI=1S/H2/h1H'), (4, 1, -8.0)),
        (('InChI=1S/CH4/h1H4', 'InChI=1S/HO/h1H',
          'InChI=1S/H2/h1H'), (5, 1, -5.0)),
        (('InChI=1S/CH4/h1H4', 'InChI=1S/CH3/h1H3',
          'InChI=1S/H2O/h1H2', 'InChI=1S/H2/h1H'), (4, 1, 1, -5.0)),
    )
    ref_cbh1 = (
        (('InChI=1S/H',), (1,)),
        (('InChI=1S/CH3/h1H3',), (1,)),
        (('InChI=1S/C2H3/c1-2/h1H,2H2', 'InChI=1S/C2H6/c1-2/h1-2H3',
          'InChI=1S/CH4/h1H4'), (1, 2, -2)),
        (('InChI=1S/C2H2/c1-2/h1-2H', 'InChI=1S/C2H6/c1-2/h1-2H3',
          'InChI=1S/C2H5/c1-2/h1H2,2H3', 'InChI=1S/C2H4/c1-2/h1-2H2', 
          'InChI=1S/C2H3/c1-2/h1H,2H2', 'InChI=1S/CH4/h1H4', 
          'InChI=1S/CH3/h1H3'), 
         (1.0, 0.333333, 0.666667, 1.333333, 0.666667, -2.666667, -0.333333)),
        (('InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/CH3O/c1-2/h1H3',
          'InChI=1S/CH4/h1H4'), (4, 1, -4)),
        (('InChI=1S/C2H5/c1-2/h1H2,2H3', 'InChI=1S/C2H6/c1-2/h1-2H3',
          'InChI=1S/CH4O/c1-2/h2H,1H3', 'InChI=1S/CH4/h1H4', 'InChI=1S/CH3/h1H3'),
         (3.0, 1.0, 1.0, -2.0, -2.0)), 
    )
    ref_cbh2 = (
        (('InChI=1S/H',), (1,)),
        (('InChI=1S/CH3/h1H3',), (1,)),
        (('InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/C3H5/c1-3-2/h1,3H,2H3',
          'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3'), (-1, 1, 1)),
        (('InChI=1S/C2H2/c1-2/h1-2H', 'InChI=1S/C2H5/c1-2/h1H2,2H3',
          'InChI=1S/C3H4/c1-3-2/h1H,2H3', 'InChI=1S/C3H3/c1-3-2/h1H,2H2',
          'InChI=1S/C2H3/c1-2/h1H,2H2', 'InChI=1S/C2H4/c1-2/h1-2H2',
          'InChI=1S/C3H4/c1-3-2/h1-2H2', 'InChI=1S/C3H5/c1-3-2/h1H2,2H3', 
          'InChI=1S/C2H6/c1-2/h1-2H3'),
         (-0.333333, -0.333334, 0.666667, 1.333333, -0.333334,
          -0.666666, 0.666667, 0.333333, -0.333333)),
        (('InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3',
          'InChI=1S/C2H5O/c1-2-3/h2H2,1H3',
          'InChI=1S/C4H10/c1-4(2)3/h4H,1-3H3'), (-2, 1, 1, 1)),
        (('InChI=1S/C2H5/c1-2/h1H2,2H3', 'InChI=1S/C3H7/c1-3-2/h1,3H2,2H3',
          'InChI=1S/C2H5O/c1-2-3/h3H,1-2H2', 'InChI=1S/C4H9/c1-4(2)3/h1-3H3'),
         (-2, 1, 1, 1))
    )

    for basis, ref_basis in zip(basic, ref_basic):
        assert basis[0] == ref_basis[0]
        assert numpy.allclose(basis[1], ref_basis[1])
    for basis, ref_basis in zip(cbh0, ref_cbh0):
        assert basis[0] == ref_basis[0]
        assert numpy.allclose(basis[1], ref_basis[1])
    for basis, ref_basis in zip(cbh1, ref_cbh1):
        assert basis[0] == ref_basis[0]
        assert numpy.allclose(basis[1], ref_basis[1])
    for basis, ref_basis in zip(cbh2, ref_cbh2):
        assert basis[0] == ref_basis[0]
        assert numpy.allclose(basis[1], ref_basis[1])


def test__hydrogen_abstraction():
    """ test thermfit.cbh._ts
    """

    cbh0 = thermfit.cbh.ts_basis(HABS_ZRXN, CBH0_SCHEME)
    cbh1 = thermfit.cbh.ts_basis(HABS_ZRXN, CBH1_SCHEME)

    ref_cbh0 = (
        ['InChI=1S/CH4/h1H4',
         'InChI=1S/H2O/h1H2',
         (('InChI=1S/CH3/h1H3', 'InChI=1S/H2O/h1H2'),
          ('InChI=1S/CH4/h1H4', 'InChI=1S/HO/h1H')),
         'InChI=1S/H2/h1H'],
        [6.0, 1.0, 1.0, -7.0]
    )
    ref_cbh1 = (
        ['InChI=1S/CH4O/c1-2/h2H,1H3',
         'InChI=1S/C2H6/c1-2/h1-2H3',
         (('InChI=1S/CH3/h1H3', 'InChI=1S/H2O2/c1-2/h1-2H'),
          ('InChI=1S/CH4/h1H4', 'InChI=1S/HO2/c1-2/h1H')),
         (('InChI=1S/C2H5/c1-2/h1H2,2H3', 'InChI=1S/H2O/h1H2'),
          ('InChI=1S/C2H6/c1-2/h1-2H3', 'InChI=1S/HO/h1H')),
         'InChI=1S/CH4/h1H4',
         'InChI=1S/H2O/h1H2',
         (('InChI=1S/CH3/h1H3', 'InChI=1S/H2O/h1H2'),
          ('InChI=1S/CH4/h1H4', 'InChI=1S/HO/h1H'))],
        [1.0, 3.0, 1.0, 2.0, -3, -1, -2]
    )

    assert cbh0[0] == ref_cbh0[0]
    assert numpy.allclose(cbh0[1], ref_cbh0[1])
    assert cbh1[0] == ref_cbh1[0]
    assert numpy.allclose(cbh1[1], ref_cbh1[1])


def test__beta_scission():
    """ test thermfit.cbh._ts
    """

    cbh0 = thermfit.cbh.ts_basis(BS_ZRXN, CBH0_SCHEME)
    cbh1 = thermfit.cbh.ts_basis(BS_ZRXN, CBH1_SCHEME)

    ref_cbh0 = (
        ['InChI=1S/CH4/h1H4',
         (['InChI=1S/C2H5O/c1-2-3/h3H,1-2H2'],
          ('InChI=1S/C2H4/c1-2/h1-2H2', 'InChI=1S/HO/h1H')),
         'InChI=1S/H2/h1H'], [2.0, 1.0, -2.0])
    ref_cbh1 = (
        ['InChI=1S/C2H6/c1-2/h1-2H3',
         (['InChI=1S/C3H7O/c1-2-3-4/h2,4H,3H2,1H3'],
          ('InChI=1S/C3H6/c1-3-2/h3H,1H2,2H3', 'InChI=1S/HO/h1H')),
         'InChI=1S/CH4/h1H4'], [1.0, 1.0, -1])

    assert cbh0[0] == ref_cbh0[0]
    assert numpy.allclose(cbh0[1], ref_cbh0[1])
    assert cbh1[0] == ref_cbh1[0]
    assert numpy.allclose(cbh1[1], ref_cbh1[1])
