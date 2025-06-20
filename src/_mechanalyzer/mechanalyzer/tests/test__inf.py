""" test chemkin_io.calculator.thermo
"""

from mechanalyzer import par
from mechanalyzer import inf


# Species info
INCHI1 = 'InChI=1S/C2H5/c1-2/h1H2,2H3'
CHARGE1 = 0
MULT1 = 2

INCHI2 = 'InChI=1S/H'
CHARGE2 = 0
MULT2 = 2

INCHI3 = 'InChI=1S/C2H6/c1-2/h1-2H3'
CHARGE3 = 0
MULT3 = 1

INCHI4 = 'InChI=1S/C2H4/c1-2/h1-2H2'
CHARGE4 = 1  # Fake charge
MULT4 = 3  # Fake triplet for testing

INCHI5 = 'InChI=1S/H2/h1H'
CHARGE5 = 0
MULT5 = 1

INCHI6 = 'FakeIch6'
CHARGE6 = 0
MULT6 = 2

INCHI7 = 'FakeIch7'
CHARGE7 = 1
MULT7 = 2


SPC_DCT = {
    'C2H5': {
        par.SPC.INCHI: INCHI1,
        par.SPC.CHARGE: CHARGE1,
        par.SPC.MULT: MULT1
    },
    'H': {
        par.SPC.INCHI: INCHI2,
        par.SPC.CHARGE: CHARGE2,
        par.SPC.MULT: MULT2
    },
    'C2H6': {
        par.SPC.INCHI: INCHI3,
        par.SPC.CHARGE: CHARGE3,
        par.SPC.MULT: MULT3
    },
    'C2H4': {
        par.SPC.INCHI: INCHI4,
        par.SPC.CHARGE: CHARGE4,
        par.SPC.MULT: MULT4
    },
    'H2': {
        par.SPC.INCHI: INCHI5,
        par.SPC.CHARGE: CHARGE5,
        par.SPC.MULT: MULT5
    },
    'Fake6': {
        par.SPC.INCHI: INCHI6,
        par.SPC.CHARGE: CHARGE6,
        par.SPC.MULT: MULT6
    },
    'Fake7': {
        par.SPC.INCHI: INCHI7,
        par.SPC.CHARGE: CHARGE7,
        par.SPC.MULT: MULT7
    }
}

# Theory info
PROG1, PROG2 = 'molpro2025', 'gaussian09'
METHOD1, METHOD2 = 'ccsd(t)', 'wb97xd'
BASIS1, BASIS2 = 'cc-pvtz', '6-31G*'
ORB_RESTRICT1, ORB_RESTRICT2, ORB_RESTRICT3 = 'RR', 'UU', 'RU'

THY_DCT = {
    par.THY.PROGRAM: PROG1,
    par.THY.METHOD: METHOD1,
    par.THY.BASIS: BASIS1,
    par.THY.ORB_RESTRICT: ORB_RESTRICT1
}

# Reaction info
RXN1 = (('C2H5', 'H'), ('C2H6',))
RXN1REV = (('C2H6',), ('C2H5', 'H'))
RXN2 = (('C2H5', 'H'), ('C2H4', 'H2'))
RXN3 = (('Fake6',), ('Fake7',))


def test__spc():
    """ test mechanalyzer.inf.spc.from_data
        test mechanalyzer.inf.spc.from_dct
        test mechanalyzer.inf.spc.combine
        test mechanalyzer.inf.spc.value
    """

    # Build info objects
    spc_info1 = inf.spc.from_data(INCHI1, CHARGE1, MULT1)
    spc_info2 = inf.spc.from_data(INCHI3, CHARGE3, MULT3)
    spc_info3 = inf.spc.from_dct(SPC_DCT['H'])
    spc_info4 = inf.spc.combine(spc_info1, spc_info2)

    assert spc_info1 == (INCHI1, CHARGE1, MULT1)
    assert spc_info2 == (INCHI3, CHARGE3, MULT3)
    assert spc_info3 == (INCHI2, CHARGE2, MULT2)
    assert spc_info4 == ((INCHI1, INCHI3), 0, 2)

    # Get a value from the inf object
    ich = inf.spc.value(spc_info1, par.SPC.INCHI)

    assert ich == INCHI1


def test__thy():
    """ test mechanalyzer.inf.thy.from_dct
        test mechanalyzer.inf.thy.from_data
        test mechanalyzer.inf.thy.modify_orb_label
        test mechanalyzer.inf.thy.value
        test mechanalyzer.inf.thy.combine
    """

    # Build info objects
    thy_info1 = inf.thy.from_dct(THY_DCT)
    thy_info2 = inf.thy.from_data(PROG2, METHOD2, BASIS2, ORB_RESTRICT2)

    assert thy_info1 == (PROG1, METHOD1, BASIS1, ORB_RESTRICT1)
    assert thy_info2 == (PROG2, METHOD2, BASIS2, ORB_RESTRICT2)

    # Build modified thy objects using spins
    mod_thy_info1 = inf.thy.modify_orb_label(
            (PROG1, METHOD1, BASIS1, ORB_RESTRICT1), (INCHI3, CHARGE3, MULT3))
    mod_thy_info2 = inf.thy.modify_orb_label(
            (PROG1, METHOD1, BASIS1, ORB_RESTRICT2), (INCHI2, CHARGE2, MULT2))
    mod_thy_info3 = inf.thy.modify_orb_label(
            (PROG1, METHOD1, BASIS1, ORB_RESTRICT3), (INCHI3, CHARGE3, MULT3))
    mod_thy_info4 = inf.thy.modify_orb_label(
            (PROG1, METHOD1, BASIS1, ORB_RESTRICT3), (INCHI2, CHARGE2, MULT2))

    assert mod_thy_info1 == (PROG1, METHOD1, BASIS1, 'R')
    assert mod_thy_info2 == (PROG1, METHOD1, BASIS1, 'U')
    assert mod_thy_info3 == (PROG1, METHOD1, BASIS1, 'R')
    assert mod_thy_info4 == (PROG1, METHOD1, BASIS1, 'U')

    # Get a value from the inf object
    basis = inf.thy.value(thy_info1, par.THY.BASIS)

    assert basis == BASIS1

    # Combine theory objects, handling orb label
    comb_mod_thy_info1 = inf.thy.combine(mod_thy_info1, mod_thy_info3)
    comb_mod_thy_info2 = inf.thy.combine(mod_thy_info1, mod_thy_info2)
    comb_mod_thy_info3 = inf.thy.combine(mod_thy_info2, mod_thy_info1)

    assert comb_mod_thy_info1 == (PROG1, METHOD1, BASIS1, 'R')
    assert comb_mod_thy_info2 == (PROG1, METHOD1, BASIS1, 'U')
    assert comb_mod_thy_info3 == (PROG1, METHOD1, BASIS1, 'U')


def test__rxn():
    """ test
        test mechanalyzer.inf.rxn.from_dct
        test mechanalyzer.inf.rxn.sort
        test mechanalyzer.inf.rxn.reverse
        test mechanalyzer.inf.rxn.ts_info
        test mechanalyzer.inf.rxn.ts_chg
        test mechanalyzer.inf.rxn.ts_mult
        test mechanalyzer.inf.rxn.rgts_info
        test mechanalyzer.inf.rxn.value
    """

    # Build reaction info objects to be used below
    reacs1, prods1 = RXN1
    reacs1rev, prods1rev = RXN1REV
    reacs2, prods2 = RXN2
    reacs3, prods3 = RXN3
    rxn_info1 = inf.rxn.from_dct(reacs1, prods1, SPC_DCT, rxn_mul='low')
    rxn_info1rev = inf.rxn.from_dct(reacs1rev, prods1rev, SPC_DCT, rxn_mul='low')
    rxn_info2 = inf.rxn.from_dct(reacs2, prods2, SPC_DCT, rxn_mul='high')
    rxn_info3 = inf.rxn.from_dct(reacs3, prods3, SPC_DCT, rxn_mul='low')
    rxn_info4 = inf.rxn.from_data(
        ((INCHI1, INCHI2), (INCHI3,)), ((0, 0), (0,)), ((2, 2), (1,)), 1)

    assert rxn_info1 == (
        ((INCHI1, INCHI2), (INCHI3,)), ((0, 0), (0,)), ((2, 2), (1,)), 1
    )
    assert rxn_info2 == (
        ((INCHI1, INCHI2), (INCHI4, INCHI5)), ((0, 0), (1, 0)), ((2, 2), (3, 1)), 3
    )
    assert rxn_info3 == (
        ((INCHI6,), (INCHI7,)), ((0,), (1,)), ((2,), (2,)), 2
    )
    assert rxn_info4 == (
        ((INCHI1, INCHI2), (INCHI3,)), ((0, 0), (0,)), ((2, 2), (1,)), 1
    )

    # Rearrange the info object
    sort_rxn_info1 = inf.rxn.sort(rxn_info1)
    sort_rxn_info2 = inf.rxn.sort(rxn_info2)
    assert sort_rxn_info1 == (
        ((INCHI3,), (INCHI1, INCHI2)), ((0,), (0, 0)), ((1,), (2, 2)), 1
    )
    assert sort_rxn_info2 == (
        ((INCHI4, INCHI5), (INCHI1, INCHI2)), ((1, 0), (0, 0)), ((3, 1), (2, 2)), 3
    )

    rev_rxn_info = inf.rxn.reverse(rxn_info2)
    assert rev_rxn_info == (
        ((INCHI4, INCHI5), (INCHI1, INCHI2)), ((1, 0), (0, 0)), ((3, 1), (2, 2)), 3
    )

    # Build ts objects from the rxn_info object
    ts_info = inf.rxn.ts_info(rxn_info1)
    ts_chg = inf.rxn.ts_chg(rxn_info1)
    ts_mul = inf.rxn.ts_mult(rxn_info3)

    assert ts_info == ('', 0, 1)
    assert ts_chg == 0
    assert ts_mul == 2

    # Build reactants info objects
    rgts_info = inf.rxn.rgts_info(rxn_info2)

    assert rgts_info == (
        ((INCHI1, 0, 2), (INCHI2, 0, 2)),
        ((INCHI4, 1, 3), (INCHI5, 0, 1))
    )

    # Get a value from the inf object
    rxn_ichs = inf.rxn.value(rxn_info1, par.SPC.INCHI)

    assert rxn_ichs == ((INCHI1, INCHI2), (INCHI3,))

    # Check other values
    assert inf.rxn.radrad(rxn_info1)
    assert inf.rxn.radrad(rxn_info1rev)
    assert inf.rxn.radrad(rxn_info2)
    assert not inf.rxn.radrad(rxn_info3)
