"""
Computes the Heat of Formation at 0 K for a given species
"""

import os
import csv
import numpy as np
from qcelemental import constants as qcc
import autoparse.pattern as app
import autoparse.find as apf
import automol.inchi
import automol.graph
from . import util

# Conversion factors
KJ2KCAL = qcc.conversion_factor('kJ/mol', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

# Path  the database files (stored in the thermo src directory)
SRC_PATH = os.path.dirname(os.path.realpath(__file__))


def get_hform_298k_thermp(output_string):
    """
    Obtains deltaHf from thermp output
    """

    # Line pattern containing the DeltaHf value at 298 K
    dhf298_pattern = ('h298 final' +
                      app.one_or_more(app.SPACE) +
                      app.capturing(app.FLOAT))
    dhf298 = float(apf.last_capture(dhf298_pattern, output_string))

    return dhf298


def calc_hform_0k(hzero_mol, hzero_basis, basis, coeff, ref_set):
    """ calculates the heat-of-formation at 0 K
    """

    # Calculate the heat of formation
    dhzero = hzero_mol * EH2KCAL
    for i, spc in enumerate(basis):
        h_basis = get_ref_h(spc, ref_set, 0)
        if h_basis is None:
            h_basis = 0.0
        dhzero += coeff[i] * h_basis * KJ2KCAL
        dhzero -= coeff[i] * hzero_basis[i] * EH2KCAL

    return dhzero


def get_ref_h(species, ref, temp):
    """ gets a reference value
    """

    # Set path and name to thermo database file
    thermodb_name = 'thermodb_{}K.csv'.format(str(int(temp)))
    thermodb_file = os.path.join(SRC_PATH, thermodb_name)

    # Find the energy value for the given species and enery type
    h_species = None
    with open(thermodb_file, 'r') as db_file:
        reader = csv.DictReader(db_file)
        for row in reader:
            if row['inchi'] == species:
                val = row[ref]
                if val == '':
                    val = None
                h_species = float(val)

    return h_species


def select_basis(atom_dct, att=0):
    """
    Given a list of atoms, generates a list of molecules
    that is best suited to serve as a basis for those atoms

    :param atomlist: list of atoms
    :type atomlist: list
    :param att: ???
    :type att: ???

    OUPUT:
    basis    - recommended basis as a list of stoichiometries
    """

    # Determine number of basis species required
    nbasis = len(atom_dct)

    # Get a list of all the atom types in the molecule
    atoms = list(atom_dct.keys())

    # Create list of inchi keys corresponding to basis species
    basis = []
    # H2
    basis.append('InChI=1S/H2/h1H')
    # NH3
    if 'N' in atoms:
        basis.append('InChI=1S/H3N/h1H3')
    # CH4
    if 'C' in atoms:
        basis.append('InChI=1S/CH4/h1H4')
    # H2O
    if 'O' in atoms:
        basis.append('InChI=1S/H2O/h1H2')
    # SO2
    if 'S' in atoms:
        basis.append('InChI=1S/O2S/c1-3-2')
        if not 'O' in atoms:
            basis.append('InChI=1S/H2O/h1H2')

    return basis


def get_reduced_basis(basis_ich, species_formula):
    """
    Form a matrix for a given basis and atomlist
    INPUT:
    input_basis     - ich strings for set of reference molecules
    atomlist  - list of atoms (all atoms that appear
                in basis should be in atomlist)
    OUTPUT:
    mat       - matrix (length of basis by length of atomlist)
                (square if done right)
    """

    # Get the basis formulae list
    basis_formulae = [automol.inchi.formula(spc) for spc in basis_ich]

    reduced_basis = []
    for i, basis_formula in enumerate(basis_formulae):
        basis_atom_dict = util.get_atom_counts_dict(basis_formula)
        flag = True
        for key, _ in basis_atom_dict.items():
            if key not in species_formula:
                flag = False

        if flag:
            reduced_basis.append(basis_ich[i])

    return reduced_basis


def calc_coefficients(basis, mol_atom_dict):
    """
    Form a matrix for a given basis and atomlist
    INPUT:
    basis     - basis of molecules
    atomlist  - list of atoms (all atoms that appear
                in basis should be in atomlist)
    OUTPUT:
    mat       - matrix (length of basis by length of atomlist)
                (square if done right)
    """

        
    # Initialize an natoms x natoms matrix
    nbasis = len(basis)
    print('basis test:', basis) 
    basis_mat = np.zeros((nbasis, nbasis))

    # Get the basis formulae list
    basis_formulae = [automol.inchi.formula(spc) for spc in basis]
    print('basis formulae:', basis_formulae)
    #basis_atom_dict = [automol.geom.formula(automol.inchi.geom(spc) for spc in basis]
    for spc in basis_formulae:
        basis_atom_dict = util.get_atom_counts_dict(spc)
        for atom in basis_atom_dict:
            if not atom in mol_atom_dict:
                mol_atom_dict[atom] = 0
    # Set the elements of the matrix
    for i, spc in enumerate(basis_formulae):
        basis_atom_dict = util.get_atom_counts_dict(spc)
        basis_vals = []
        for key in mol_atom_dict.keys():
            if key in basis_atom_dict:
                basis_vals.append(basis_atom_dict[key])
            else:
                basis_vals.append(0)
        basis_mat[i] = basis_vals

    #  Transpose
    basis_mat = basis_mat.T

    # Form stoich vector
    stoich_vec = np.zeros(len(mol_atom_dict))
    for i, key in enumerate(mol_atom_dict.keys()):
        stoich_vec[i] = mol_atom_dict[key]

    # Solve C = M^-1 S
    basis_mat = np.linalg.inv(basis_mat)
    coeff = np.dot(basis_mat, stoich_vec)

    return coeff


def stoich(ich):
    """
    Finds the stoichiometry of a molecule
    INPUT:
    ich  -- STR inchii
    OUTPUT:
    stoich -- dictionary with key = STR atomsymbol,
                val = INT number of atomsymbol in molecule
    """

    stoich = {'H': 0}
    gra = automol.inchi.graph(ich)
    atms = automol.graph.atoms(gra)
    for atm in atms:
        stoich['H'] += atms[atm][1]
        if atms[atm][0] in stoich:
            stoich[atms[atm][0]] += 1
        else:
            stoich[atms[atm][0]] = 1
    return stoich


def cbhzed(ich):
    """
    Fragments molecule so that each heavy-atom is a seperate fragment
    INPUT:
    ich --  STR inchii name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchii name for fragments
    and value as INT their coefficient
    """

    # Graphical info about molecule
    gra = automol.inchi.graph(ich)
    rad_atms = list(automol.graph.sing_res_dom_radical_atom_keys(gra))
    atm_vals = automol.graph.atom_element_valences(gra)
    atms = automol.graph.atoms(gra)

    # Determine CBHzed fragments
    frags = {}
    for atm in atm_vals:
        if atm in rad_atms:
            atm_vals[atm] -= 1
        atm_dic = {0: (atms[atm][0], int(atm_vals[atm]), None)}
        gra = (atm_dic, {})
        frag = automol.graph.inchi(gra)
        _add2dic(frags, frag)
    print('cbhzed frags', frags)
    print(_balance_frags(ich, frags))
    return _balance_frags(ich, frags)


def cbhone(ich):
    """
    Fragments molecule in a way that conserves each heavy-atom/heavy-atom bond
    INPUT:
    ich --  STR inchii name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchii name for fragments
    and value as INT their coefficient
    """

    # Graphical info about molecule
    gra = automol.inchi.graph(ich)
    atms = automol.graph.atoms(gra)
    bnd_ords = automol.graph.one_resonance_dominant_bond_orders(gra)
    rad_atms = list(automol.graph.sing_res_dom_radical_atom_keys(gra))
    atm_vals = automol.graph.atom_element_valences(gra)
    adj_atms = automol.graph.atom_neighbor_keys(gra)

    # Determine CBHone fragments
    frags = {}
    for atm in atm_vals:
        for adj in list(adj_atms[atm]):
            if atm > adj:
                vali = atm_vals[atm]
                valj = atm_vals[adj]
                if atm in rad_atms:
                    vali -= 1
                if adj in rad_atms:
                    valj -= 1
                key = frozenset({atm, adj})
                bnd_ord = list(bnd_ords[key])[0]
                vali -= bnd_ord
                valj -= bnd_ord
                atm_dic = {0: (atms[atm][0], int(vali), None),
                           1: (atms[adj][0], int(valj), None)}
                bnd_dic = {frozenset({0, 1}): (1, None)}
                gra = (atm_dic, bnd_dic)
                frag = automol.graph.inchi(gra)
                _add2dic(frags, frag)
    frags = {k: v for k, v in frags.items() if v}

    # Balance
    balance_ = _balance(ich, frags)
    balance_ = {k: v for k, v in balance_.items() if v}

    if balance_:
        newfrags = {}
        zedfrags = cbhzed(ich)
        new = {}
        for frag in frags:
            newfrags[frag] = frags[frag]
            new = cbhzed(frag)
            for n in new:
                _add2dic(newfrags, n, - new[n] * frags[frag])
        if not frags:
            frags = cbhzed(ich)
        for frag in zedfrags:
            if frag in newfrags:
                _add2dic(newfrags, frag, zedfrags[frag])
        frags = newfrags
        frags = {k: v for k, v in frags.items() if v}
        balance_ = _balance(ich, frags)
        balance_ = {k: v for k, v in balance_.items() if v}

    if balance_:
        frags = _balance_frags(ich, frags)

    return frags


def cbhtwo(ich):
    """
    Fragments molecule for each heavy-atom to stay bonded to its adjacent atoms
    INPUT:
    ich --  STR inchii name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchii name for fragments and
    value as INT their coefficient
    """

    # Graphical info about molecule
    gra = automol.inchi.graph(ich)
    atms = automol.graph.atoms(gra)
    bnd_ords = automol.graph.one_resonance_dominant_bond_orders(gra)
    rad_atms = list(automol.graph.sing_res_dom_radical_atom_keys(gra))
    atm_vals = automol.graph.atom_element_valences(gra)
    adj_atms = automol.graph.atom_neighbor_keys(gra)

    # Determine CBHtwo fragments
    frags = {}
    for atm in atms:
        vali = atm_vals[atm]
        if atm in rad_atms:
            vali -= 1
        # First loop over all atoms of this frag to get saturation of atomi
        for adj in list(adj_atms[atm]):
            key = frozenset({atm, adj})
            bnd_ord = list(bnd_ords[key])[0]
            vali -= bnd_ord
        atm_dic = {0: (atms[atm][0], int(vali), None)}
        bnd_dic = {}
        # Then start adding bonds to the bnddic and atomdic
        j = 0
        for adj in list(adj_atms[atm]):
            j += 1
            valj = atm_vals[adj]
            if adj in rad_atms:
                valj -= 1
            key = frozenset({atm, adj})
            bnd_ord = list(bnd_ords[key])[0]
            valj -= bnd_ord
            atm_dic[j] = (atms[adj][0], int(valj), None)
            bnd_dic[frozenset({0, j})] = (1, None)
        gra = (atm_dic, bnd_dic)
        frag = automol.graph.inchi(gra)
        _add2dic(frags, frag)

    frags = {k: v for k, v in frags.items() if v}

    # Balance
    balance_ = _balance(ich, frags)
    balance_ = {k: v for k, v in balance_.items() if v}
    if balance_:
        newfrags = {}
        onefrags = cbhone(ich)
        new = {}
        for frag in frags:
            newfrags[frag] = frags[frag]
            new = cbhone(frag)
            for n in new:
                _add2dic(newfrags, n, - new[n] * frags[frag])
        if not frags:
            frags = cbhone(ich)
        for frag in onefrags:
            if frag in newfrags:
                _add2dic(newfrags, frag, onefrags[frag])
        frags = newfrags
        frags = {k: v for k, v in frags.items() if v}

        balance_ = _balance(ich, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            newfrags = {}
            zedfrags = cbhzed(ich)
            new = {}
            for frag in frags:
                newfrags[frag] = frags[frag]
                new = cbhzed(frag)
                for n in new:
                    _add2dic(newfrags, n, - new[n] * frags[frag])
            if not frags:
                frags = cbhzed(ich)
            for frag in zedfrags:
                if frag in newfrags:
                    _add2dic(newfrags, frag, zedfrags[frag])
            frags = newfrags
            frags = {k: v for k, v in frags.items() if v}

            balance_ = _balance(ich, frags)
            balance_ = {k: v for k, v in balance_.items() if v}
            if balance_:
                frags = _balance_frags(ich, frags)

    return frags


def get_basis(ich):
    formula  = automol.inchi.formula(ich)
    atm_dict = util.get_atom_counts_dict(formula)
    return select_basis(atm_dict)


def get_basic(ich):
    formula_dct  = automol.inchi.formula_dct(ich)
    spc_bas = select_basis(formula_dct)
    if len(spc_bas) == 1 and ich == spc_bas[0]:
        clist = [1]
    else:
        clist = calc_coefficients(spc_bas, formula_dct)
    return spc_bas, clist


def get_cbhzed(ich):
    frags = cbhzed(ich)
    fraglist = []
    clist = []
    for frag in frags:
        fraglist.append(frag)
        clist.append(frags[frag])
    return fraglist, clist
    #return list(cbhzed(ich).keys())


def get_cbhone(ich):
    frags = cbhone(ich)
    fraglist = []
    clist = []
    for frag in frags:
        fraglist.append(frag)
        clist.append(frags[frag])
    return fraglist, clist
    #return list(cbhone(ich).keys())


def get_cbhtwo(ich):
    frags = cbhtwo(ich)
    fraglist = []
    clist = []
    for frag in frags:
        fraglist.append(frag)
        clist.append(frags[frag])
    return fraglist, clist
    #return list(cbhtwo(ich).keys())


def _add2dic(dic, key, val=1):
    if key in dic:
        dic[key] += val
    else:
        dic[key] = val


def _lhs_rhs(frags):
    rhs = {}
    lhs = {}
    for frag in frags:
        if frags[frag] > 0:
            rhs[frag] = frags[frag]
        elif frags[frag] < 0:
            lhs[frag] = - frags[frag]
    return lhs, rhs


def _print_lhs_rhs(ich, frags):
    lhs, rhs = _lhs_rhs(frags)
    lhsprint = automol.inchi.smiles(ich)
    rhsprint = ''
    for frag in rhs:
        if rhsprint:
            rhsprint += ' +  {:.1f} {} '.format(
                rhs[frag], automol.inchi.smiles(frag))
        else:
            rhsprint = ' {:.1f} {} '.format(
                rhs[frag], automol.inchi.smiles(frag))
    for frag in lhs:
        lhsprint += ' +  {:.1f} {} '.format(
            lhs[frag], automol.inchi.smiles(frag))
    return '{} --> {}'.format(lhsprint, rhsprint)


def _balance(ich, frags):
    stoichs = {}
    for frag in frags:
        _stoich = stoich(frag)
        for atm in _stoich:
            if atm in stoichs:
                stoichs[atm] += _stoich[atm] * frags[frag]
            else:
                stoichs[atm] = _stoich[atm] * frags[frag]
    balance_ = {}
    _stoich = stoich(ich)
    for atom in _stoich:
        if atom in stoichs:
            balance_[atom] = _stoich[atom] - stoichs[atom]
        else:
            balance_[atom] = _stoich[atom]
    balance_ = {x: y for x, y in balance_.items() if y != 0}
    return balance_


def _balance_frags(ich, frags):
    balance_ = _balance(ich, frags)
    methane = automol.smiles.inchi('C')
    water = automol.smiles.inchi('O')
    ammonm = automol.smiles.inchi('N')
    hydrgn = automol.smiles.inchi('[H][H]')
    if 'C' in balance_:
        _add2dic(frags, methane, balance_['C'])
    if 'N' in balance_:
        _add2dic(frags, ammonm, balance_['N'])
    if 'O' in balance_:
        _add2dic(frags, water, balance_['O'])
    balance_ = _balance(ich, frags)
    if 'H' in balance_:
        _add2dic(frags, hydrgn, balance_['H']/2)
    return frags
