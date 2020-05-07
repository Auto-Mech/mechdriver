"""
Tests calculating the 0 K heat-of-formation
"""

import os
from thermo import heatform
from thermo import util
import automol.inchi
import automol.smiles

# Inchi string for methyl nitrate (CH3ONO2)
ICH = 'InChI=1S/CH3NO3/c1-5-2(3)4/h1H3'
SMI = 'C=CC(=O)O'
ICH2 = automol.smiles.inchi(SMI)

# Thermp output file name
THERMP_OUTFILE_NAME = os.path.join(os.getcwd(), 'run', 'thermp.out')


def test__calc_hform_0k():
    """ calculates 0 K heat-of-formation for a species
    """

    # Get the molecular formula from the inchi string
    #formula = util.inchi_formula(ICH)
    formula = automol.inchi.formula(ICH)
    print('\nformula:')
    print(formula)

    # Get atom count dictionary
    #atom_dict = util.get_atom_counts_dict(formula)
    atom_dict = automol.inchi.formula_dct(ICH)
    print('\natom dict:')
    print(atom_dict)

    # Get the list of the basis
    basis = heatform.select_basis(atom_dict)
    print('\nbasis:')
    print(basis)
    
    # Get the basis list from reduced_basis
    #red_basis = heatform.select_basis(basis_ich, formula)
    #print('\nreduced basis:')
    #print(red_basis)

    # Get the coefficients for the balanced heat-of-formation eqn
    coeff = heatform.calc_coefficients(basis, atom_dict)
    print('\ncoeff:')
    print(coeff)

    # Obtain the reference energies from the database
    print('\nref e:')
    for spc in basis:
        ref_e = heatform.get_ref_h(spc, 'ATcT', 0)
        print(spc)
        print(ref_e)

    # Get the energy for the species and basis
    e_mol = -100.0
    e_basis = [-1.0, -2.0, -3.0, -4.0]
    print('\ne_mol and e_basis:')
    print(e_mol)
    print(e_basis)

    # Get the 0 K heat of formation
    hform = heatform.calc_hform_0k(e_mol, e_basis, basis, coeff, ref_set='ATcT')
    print('\nhform(0 K):')
    print(hform)


def test__read_hform_298k():
    """ reads the 298 K heat-of-formation value from thermp output
    """

    # Read the thermp output
    with open(THERMP_OUTFILE_NAME, 'r') as thermp_outfile:
        thermp_out_str = thermp_outfile.read()

    # Get the 0 K heat of formation
    hform = heatform.get_hform_298k_thermp(thermp_out_str)
    print('\nhform(298 K):')
    print(hform)


def test__cbhzed():
    """ Fragments molecule in a way that conserves each heavy-atom/heavy-atom bond
    """
    frags = heatform.cbhzed(ICH2)
    print('\nCBH0 formula: ', heatform._print_lhs_rhs(ICH2, frags))

def test__cbhone():
    """ Fragments molecule in a way that conserves each heavy-atom/heavy-atom bond
    """
    frags = heatform.cbhone(ICH2)
    print('\nCBH1 formula: ', heatform._print_lhs_rhs(ICH2, frags))

def test__cbhtwo():
    """ Fragments molecule in a way that conserves each heavy-atom/heavy-atom bond
    """
    frags = heatform.cbhtwo(ICH2)
    print('\nCBH2 formula: ', heatform._print_lhs_rhs(ICH2, frags))

if __name__ == '__main__':
    test__calc_hform_0k()
    test__read_hform_298k()
    test__cbhzed()
    test__cbhone()
    test__cbhtwo()
