""" test thermfit.heatform
"""

import pytest
import numpy
import thermfit.heatform


# Energies (ene+zpve [Hartrees])
DCT_H0 = {
    'InChI=1S/H2/h1H': -1.16329723,
    'InChI=1S/CH4/h1H4': -40.40245975,
    'InChI=1S/CH3/h1H3': -39.73962587,
    'InChI=1S/H2O/h1H2': -76.33504552,
    'InChI=1S/HO/h1H': -75.64954341
}
C3H7OH_H0 = -193.9754161
TS_CH4_OH_H0 = -116.0451612

# Species
C3H7OH_ICH = 'InChI=1S/C3H8O/c1-2-3-4/h4H,2-3H2,1H3'
C3H7OH_BASIS = ('InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4', 'InChI=1S/H2O/h1H2')
C3H7OH_BASIS_H0 = tuple(DCT_H0[base] for base in C3H7OH_BASIS)
C3H7OH_CFFS = numpy.array([-3., 3., 1.])

# TS
TS_CH4_OH_ICHS = (('InChI=1S/CH4/h1H4', 'InChI=1S/HO/h1H'),
                  ('InChI=1S/CH3/h1H3', 'InChI=1S/H2O/h1H2'))

# InChI not in database
NON_DB_ICH = 'InChI=1S/C8H18/c1-3-5-7-8-6-4-2/h3-8H2,1-2H3'

# Database reading parameters
TEMP1 = 0
REF_SET1 = 'ANL0'
REF_SET2 = 'ANL1'
REF_SET3 = 'Ladder'


def test__hform():
    """ test thermfit.heatform.calc_hform_0k
    """

    ref_hf0k = -0.08993069717524804

    hf0k = thermfit.heatform.calc_hform_0k(
        C3H7OH_H0, C3H7OH_BASIS_H0, C3H7OH_BASIS, C3H7OH_CFFS, REF_SET1)
    assert numpy.isclose(ref_hf0k, hf0k)

    # test reaction? (would need list of lists in calc_hform)


def test__ref_enthalpy():
    """ test thermfit.heatform.reference_enthalpy
    """

    # Species
    ref_hf0k_1 = -0.08780432505395637
    ref_hf0k_2 = -0.08785677661265853

    hf0k_1 = thermfit.heatform.reference_enthalpy(
        C3H7OH_ICH, REF_SET1, TEMP1, rxn=False)
    hf0k_2 = thermfit.heatform.reference_enthalpy(
        C3H7OH_ICH, REF_SET3, TEMP1, rxn=False)
    assert numpy.isclose(ref_hf0k_1, hf0k_1)
    assert numpy.isclose(ref_hf0k_2, hf0k_2)

    # Transition State
    ref_hf0k = -0.0033071800404707342

    hf0k = thermfit.heatform.reference_enthalpy(
        TS_CH4_OH_ICHS, REF_SET1, TEMP1, rxn=True)
    assert numpy.isclose(ref_hf0k, hf0k)

    # Missing Species
    with pytest.raises(AssertionError):
        _ = thermfit.heatform.reference_enthalpy(
            NON_DB_ICH, REF_SET1, TEMP1, rxn=False)
