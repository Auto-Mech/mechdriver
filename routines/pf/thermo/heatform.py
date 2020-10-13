"""
Computes the Heat of Formation at 0 K for a given species
"""

import os
import csv
import numpy as np
from qcelemental import constants as qcc
import automol.inchi
import automol.graph
from . import util


# Conversion factors
KJ2KCAL = qcc.conversion_factor('kJ/mol', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

# Path  the database files (stored in the thermo src directory)
SRC_PATH = os.path.dirname(os.path.realpath(__file__))


def calc_hform_0k(hzero_mol, hzero_basis, basis, coeff, ref_set):
    """ calculates the heat-of-formation at 0 K
    """

    dhzero = hzero_mol * EH2KCAL
    for i, spc in enumerate(basis):
        ts = True
        if isinstance(spc, str):
            ts = False
        h_basis = get_ref_h(spc, ref_set, 0, ts)
        if h_basis is None:
            h_basis = 0.0
        dhzero += coeff[i] * h_basis * KJ2KCAL
        print('hzero', hzero_basis[i])
        dhzero -= coeff[i] * hzero_basis[i] * EH2KCAL

    return dhzero


def get_ref_h(species, ref, temp, ts=False):
    """ gets a reference value
    """

    # Set path and name to thermo database file
    if ts:
        thermodb_name = 'tsthermodb_{}K.csv'.format(str(int(temp)))
    else:    
        thermodb_name = 'thermodb_{}K.csv'.format(str(int(temp)))
    thermodb_file = os.path.join(SRC_PATH, thermodb_name)
    # Find the energy value for the given species and enery type
    h_species = None
    if ts:
        rcts, prds = species
        rct_str = '+'.join(rcts)
        prd_str = '+'.join(prds)
        species = '='.join([rct_str, prd_str])
    print('species is ', species)
    with open(thermodb_file, 'r') as db_file:
        reader = csv.DictReader(db_file)
        for row in reader:
            if row['inchi'] == species:
                val = row[ref]
                if val == '':
                    val = None
                h_species = float(val)
    print('h_species is ', h_species)
    return h_species


def select_basis(atom_dct):
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
    # Cl2
    if 'Cl' in atoms:
        basis.append('InChI=1S/ClH/h1H')
        # basis.append('InChI=1S/Cl2/c1-2')
    # SO2
    if 'S' in atoms:
        basis.append('InChI=1S/O2S/c1-3-2')
        if 'O' not in atoms:
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
    basis_formulae = [automol.inchi.formula_string(spc) for spc in basis_ich]

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
    basis_mat = np.zeros((nbasis, nbasis))

    # Get the basis formulae list
    basis_formulae = [automol.inchi.formula_string(spc) for spc in basis]
    # print('basis formulae:', basis_formulae)
    # basis_atom_dict = [
    # automol.geom.formula(automol.inchi.geom(spc) for spc in basis]
    for spc in basis_formulae:
        basis_atom_dict = util.get_atom_counts_dict(spc)
        for atom in basis_atom_dict:
            if atom not in mol_atom_dict:
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

def stoich_gra(gra):
    atms = automol.graph.atoms(gra)
    stoich_dct = {}
    hcount = 0
    for atm in atms:
        hcount += np.floor(atms[atm][1])
        if atms[atm][0] in stoich_dct:
            stoich_dct[atms[atm][0]] += 1
        else:
            stoich_dct[atms[atm][0]] = 1
    if not 'H' in stoich_dct:
        stoich_dct['H'] = hcount
    else:
        stoich_dct['H'] += hcount
    return stoich_dct

def stoich(ich):
    """
    Finds the stoichiometry of a molecule
    INPUT:
    ich  -- STR inchi
    OUTPUT:
    stoich -- dictionary with key = STR atomsymbol,
                val = INT number of atomsymbol in molecule
    """

    stoich_dct = {'H': 0}
    gra = automol.inchi.graph(ich)
    atms = automol.graph.atoms(gra)
    for atm in atms:
        stoich_dct['H'] += atms[atm][1]
        if atms[atm][0] in stoich_dct:
            stoich_dct[atms[atm][0]] += 1
        else:
            stoich_dct[atms[atm][0]] = 1
    return stoich_dct

def branchpoint(adj_atms_i):
    return max(1, len(adj_atms_i) - 1)

def terminalmoity(adj_atms_i):
     ret = 1
     if len(adj_atms_i) == 1:
         ret = 0     
     return ret

def cbhzed(ich, bal=True):
    """
    Fragments molecule so that each heavy-atom is a seperate fragment
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchi name for fragments
    and value as INT their coefficient
    """

    # Graphical info about molecule
    gra = automol.inchi.graph(ich)
    rad_atms = list(automol.graph.sing_res_dom_radical_atom_keys(gra))
    atm_vals = automol.graph.atom_element_valences(gra)
    atms = automol.graph.atoms(gra)
    adj_atms = automol.graph.atom_neighbor_keys(gra)

    # Determine CBHzed fragments
    frags = {}
    for atm in atm_vals:
        coeff = 1
        if not bal:
            coeff = branchpoint(adj_atms[atm]) * terminalmoity(adj_atms[atm])
        if atm in rad_atms:
            atm_vals[atm] -= 1
        atm_dic = {0: (atms[atm][0], int(atm_vals[atm]), None)}
        gra = (atm_dic, {})
        frag = automol.graph.inchi(gra)
        _add2dic(frags, frag, coeff)
    if bal:
        balance_ = _balance(ich, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            frags = _balance_frags(ich, frags)
    return frags

def _ts_graph(gra, site):
    rad_atms = list(automol.graph.sing_res_dom_radical_atom_keys(gra))
    atm_vals = automol.graph.atom_element_valences(gra)
    atms = automol.graph.atoms(gra)
    bnd_ords = automol.graph.one_resonance_dominant_bond_orders(gra)
    adj_atms = automol.graph.atom_neighbor_keys(gra)
  
    frm = [site[0], site[1]]
    brk = [site[1], site[2]]
    frm.sort()
    brk.sort()
    #update adjacent atms list to consider TS connected
    if site[1] in adj_atms[site[0]]:
        if site[1] in adj_atms[site[2]]:
            if site[1] in rad_atms:
               new_bnd_ords = bnd_ords.copy()
               #new_bnd_ords[frozenset({*frm})] = frozenset({list(bnd_ords[frozenset({*frm})])[0] - 1})
               bnd_dic = {}            
               for key in bnd_ords: 
                   bnd_dic[key] = (list(new_bnd_ords[key])[0], None)
               new_gra = (atms, bnd_dic)
               new_rad_atms = list(automol.graph.sing_res_dom_radical_atom_keys(gra))
               if site[1] not in new_rad_atms:
                   rad_atms.remove(site[1])
            bnd_ords[frozenset({*frm})] = frozenset({list(bnd_ords[frozenset({*frm})])[0] + 0.6})
            bnd_ords[frozenset({*brk})] = frozenset({list(bnd_ords[frozenset({*brk})])[0] - 0.6})
        else:
            if site[2] in rad_atms:
               new_bnd_ords = bnd_ords.copy()
               new_bnd_ords[frozenset({*frm})] = frozenset({list(bnd_ords[frozenset({*frm})])[0] - 1})
               new_bnd_ords[frozenset({*brk})] = frozenset({1})
               bnd_dic = {}            
               for key in bnd_ords: 
                   bnd_dic[key] = (list(new_bnd_ords[key])[0], None)
               new_gra = (atms, bnd_dic)
               new_rad_atms = list(automol.graph.sing_res_dom_radical_atom_keys(gra))
               if site[2] not in new_rad_atms:
                   rad_atms.remove(site[2])
            bnd_ords[frozenset({*frm})] = frozenset({list(bnd_ords[frozenset({*frm})])[0] - 0.4})
            bnd_ords[frozenset({*brk})] = frozenset({0.4})
            adj_atms[site[2]] = frozenset({site[1], *adj_atms[site[2]]})
            adj_atms[site[1]] = frozenset({site[2], *adj_atms[site[1]]})
    else:
        if site[0] in rad_atms:
           new_bnd_ords = bnd_ords.copy()
           new_bnd_ords[frozenset({*brk})] = frozenset({list(bnd_ords[frozenset({*brk})])[0] - 1})
           new_bnd_ords[frozenset({*frm})] = frozenset({1})
           bnd_dic = {}            
           for key in bnd_ords: 
               bnd_dic[key] = (list(new_bnd_ords[key])[0], None)
           new_gra = (atms, bnd_dic)
           new_rad_atms = list(automol.graph.sing_res_dom_radical_atom_keys(gra))
           if site[0] not in new_rad_atms:
               rad_atms.remove(site[0])
        bnd_ords[frozenset({*frm})] = frozenset({0.6})
        bnd_ords[frozenset({*brk})] = frozenset({list(bnd_ords[frozenset({*brk})])[0] - 0.6})
        adj_atms[site[0]] = frozenset({site[1], *adj_atms[site[0]]})
        adj_atms[site[1]] = frozenset({site[0], *adj_atms[site[1]]})
    print('ts gra', atms, bnd_ords)    
    return rad_atms, atms, bnd_ords, atm_vals

def remove_zero_order_bnds(gra):
    atms, bnds = gra
    new_bnds = {}
    for bnd in bnds:
        if bnds[bnd][0] > 0:
            new_bnds[bnd] = bnds[bnd]
    print('new bonds ', new_bnds)        
    return (atms, new_bnds)


def split_beta_gras(gras):
    rct_ichs = ['']
    prd_ichs = ['','']
    atms, bnd_ords = gras
    atms = atms.copy()
    bnd_ords = bnd_ords.copy()
    for bnd_ord in bnd_ords:
        order, tmp = bnd_ords[bnd_ord]
        if abs(np.floor(order) - (order - 0.4)) < 0.01:
            bnd_ords[bnd_ord] = (round(order + 0.6, 1), tmp)
            atmai, atmbi = bnd_ord
            if not abs(np.floor(atms[atmai][1]) - (atms[atmai][1]-0.6)) < 0.01:
                atmbi, atmai = atmai, atmbi
            atma = list(atms[atmai])
            atmb = list(atms[atmbi])
            atma[1] = round(atma[1] - 0.6, 1)
            atms[atmai] = tuple(atma)
            atms[atmbi] = tuple(atmb)
            for atmi in atms:
                if abs(np.floor(atms[atmi][1]) - (atms[atmi][1]-0.4)) < .01:
                      atm =  list(atms[atmi])
                      atm[1] = round(atm[1] - 0.4, 1)
                      atms[atmi] = tuple(atm)
                      order, tmp =  bnd_ords[frozenset({atmbi, atmi})]
                      bnd_ords[frozenset({atmbi, atmi})] = (round(order - 0.6, 1), tmp)
            rct_gra = remove_zero_order_bnds((atms, bnd_ords))
            rct_gras = automol.graph.connected_components(rct_gra)
            for idx, rgra in enumerate(rct_gras):
                if rgra:
                   rct_ichs[idx] = automol.graph.inchi(rgra)
    rct_ichs = automol.inchi.sorted_(rct_ichs)
    atms, bnd_ords = gras
    for bnd_ord in bnd_ords:
        order, tmp = bnd_ords[bnd_ord]
        if abs(np.floor(order) - (order - 0.4)) < 0.01:
            bnd_ords[bnd_ord] = (round(order - 0.4, 1), tmp)
            atmai, atmbi = bnd_ord
            if not abs(np.floor(atms[atmai][1]) - (atms[atmai][1]-0.6)) < 0.01:
                atmbi, atmai = atmai, atmbi
            atma = list(atms[atmai])
            atmb = list(atms[atmbi])
            atma[1] = round(atma[1] - 0.6, 1)
            atms[atmai] = tuple(atma)
            atms[atmbi] = tuple(atmb)
            for atmi in atms:
                if abs(np.floor(atms[atmi][1]) - (atms[atmi][1]-0.4)) < 0.01:
                      atm =  list(atms[atmi])
                      atm[1] = round(atm[1] - 0.4, 1)
                      atms[atmi] = tuple(atm)
                      order, tmp =  bnd_ords[frozenset({atmbi, atmi})]
                      bnd_ords[frozenset({atmbi, atmi})] = (round(order + 0.4, 1), tmp)
            prd_gra = remove_zero_order_bnds((atms, bnd_ords))
            prd_gras = automol.graph.connected_components(prd_gra)
            for idx, pgra in enumerate(prd_gras):
                prd_ichs[idx] = automol.graph.inchi(pgra)
    prd_ichs = automol.inchi.sorted_(prd_ichs)
    return (rct_ichs, prd_ichs)  

def split_gras(gras):
    rct_ichs = []
    prd_ichs = []
    atms, bnd_ords = gras
    atms = atms.copy()
    bnd_ords = bnd_ords.copy()
    for bnd_ord in bnd_ords:
        order, tmp = bnd_ords[bnd_ord]
        if abs(np.floor(order) - (order - 0.4)) < 0.01:
            bnd_ords[bnd_ord] = (round(order + 0.6, 1), tmp)
            atmai, atmbi = bnd_ord
            if not abs(np.floor(atms[atmai][1]) - (atms[atmai][1]-0.4)) < 0.01:
                atmbi, atmai = atmai, atmbi
            atma = list(atms[atmai])
            atmb = list(atms[atmbi])
            atma[1] = round(atma[1] - 0.6, 1)
            atms[atmai] = tuple(atma)
            atms[atmbi] = tuple(atmb)
            for atmi in atms:
                if abs(np.floor(atms[atmi][1]) - (atms[atmi][1]-0.4)) < 0.01:
                      atm =  list(atms[atmi])
                      atm[1] = round(atm[1] - 0.4, 1)
                      atms[atmi] = tuple(atm)
                      order, tmp =  bnd_ords[frozenset({atmbi, atmi})]
                      bnd_ords[frozenset({atmbi, atmi})] = (round(order - 0.6, 1), tmp)
            rct_gra = remove_zero_order_bnds((atms, bnd_ords))
            rct_gras = automol.graph.connected_components(rct_gra)
            for rgra in rct_gras:
                rct_ichs.append(automol.graph.inchi(rgra))
    if len(rct_ichs) > 1:            
        rct_ichs = automol.inchi.sorted_(rct_ichs)
    atms, bnd_ords = gras
    for bnd_ord in bnd_ords:
        order, tmp = bnd_ords[bnd_ord]
        if abs(np.floor(order) - (order - 0.4)) < 0.01:
            bnd_ords[bnd_ord] = (round(order - 0.4, 1), tmp)
            atmai, atmbi = bnd_ord
            if not abs(np.floor(atms[atmai][1]) -( atms[atmai][1]-0.4)) < 0.01:
                atmbi, atmai = atmai, atmbi
            atma = list(atms[atmai])
            atmb = list(atms[atmbi])
            atma[1] = round(atma[1] - 0.6, 1)
            atms[atmai] = tuple(atma)
            atms[atmbi] = tuple(atmb)
            for atmi in atms:
                if abs(np.floor(atms[atmi][1]) - (atms[atmi][1]-0.4)) < 0.01:
                      atm =  list(atms[atmi])
                      atm[1] = round(atm[1] - 0.4, 1)
                      atms[atmi] = tuple(atm)
                      atms[atmi] = tuple(atm)
                      order, tmp =  bnd_ords[frozenset({atmbi, atmi})]
                      bnd_ords[frozenset({atmbi, atmi})] = (round(order + 0.4, 1), tmp)
            prd_gra = remove_zero_order_bnds((atms, bnd_ords))
            prd_gras = automol.graph.connected_components(prd_gra)
            for pgra in prd_gras:
                prd_ichs.append(automol.graph.inchi(pgra))
    if len(prd_ichs) > 1:    
        prd_ichs = automol.inchi.sorted_(prd_ichs)
    return (rct_ichs, prd_ichs)   


def cbhzed_habs(gra, site, bal=True):
    """
    Fragments molecule so that each heavy-atom is a seperate fragment
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchi name for fragments
    and value as INT their coefficient
    """

    # Graphical info about molecule
    rad_atms, atms, bnd_ords, atm_vals = _ts_graph(gra, site)
    # Determine CBHzed fragments
    frags = {}
    #for atm in atm_vals:
    #    if atm in rad_atms:
    #        atm_vals[atm] -= 1
    print('geo in habs:\n', atms, bnd_ords)
    for atm in atm_vals:
        if (atms[atm][0] != 'H' or atm in site):
            if atm == site[1] or atm == site[2]:
                continue
            coeff = 1.0
            if not bal:
                coeff = branchpoint(adj_atms[atm]) * terminalmoity(adj_atms[atm])
            if atm == site[0]:
                key1 = [site[0], site[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord1 = list(bnd_ords[key])[0]
                atm_dic = {0: (atms[atm][0], atm_vals[atm]-bnd_ord1, None)}
                bnd_dct = {frozenset({0, 1}): (bnd_ord1, None)}
                key1 = [site[2], site[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord2 = list(bnd_ords[key])[0]
                atm_dic[2] = (atms[site[2]][0], atm_vals[site[2]]-bnd_ord2, None)
                atm_dic[1] = (atms[site[1]][0], atm_vals[site[1]]-bnd_ord1-bnd_ord2, None)
                bnd_dct[frozenset({1, 2})] = (bnd_ord2, None)
            else:   
                bnd_dct = {}
                atm_dic = {0: (atms[atm][0], int(atm_vals[atm]), None)}
            grai = (atm_dic, bnd_dct)
            print('frag graph grai ', grai)
            try:
                grai = automol.graph.explicit(grai)
                key = 'exp_gra'
            except:
                key = 'ts_gra'
            newname = None
            repeat = False
            for name in frags:
                if key in frags[name]:
                    if automol.graph.full_isomorphism(frags[name][key], grai):
                        newname = name
                        repeat = True
            if not repeat:
                newname = len(frags.keys())
                frags[newname] = {}
                frags[newname][key] = grai
            #frag = automol.graph.inchi(gra)
            _add2dic(frags[newname], 'coeff', coeff)
    frags = _simplify_gra_frags(frags)
    if bal:
        balance_ = _balance_ts(gra, frags)
          
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            frags = _balance_frags_ts(gra, frags)
    frags = _simplify_gra_frags(frags)
    return frags

def _simplify_gra_frags(frags):
    new_frags = {}
    for i, frag in enumerate(frags.keys()):
        if abs(frags[frag]['coeff']) > 0.01:
            new_frags[i] = frags[frag]
    return new_frags


def cbhone(ich, bal=True):
    """
    Fragments molecule in a way that conserves each heavy-atom/heavy-atom bond
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchi name for fragments
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
    if not frags:
       frags = cbhzed(ich)
    # Balance
    if bal:
        balance_ = _balance(ich, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            zedfrags = cbhzed(ich, bal=False)
            newfrags = frags.copy()
            for frag in zedfrags:
                _add2dic(newfrags, frag,  -zedfrags[frag])
            frags = {k: v for k, v in newfrags.items() if v}
        balance_ = _balance(ich, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            frags = _balance_frags(ich, frags)
    print( frags)
    return frags


def cbhtwo(ich, bal=True):
    """
    Fragments molecule for each heavy-atom to stay bonded to its adjacent atoms
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchi name for fragments and
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
        coeff = 1
        if not bal:
            coeff = branchpoint(adj_atms[atm]) * terminalmoity(adj_atms[atm])
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
        _add2dic(frags, frag, coeff)

    frags = {k: v for k, v in frags.items() if v}
    if not frags:
        frags = cbhone(frags)
    # Balance
    if bal:
        balance_ = _balance(ich, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            newfrags = frags.copy()
            onefrags = cbhone(ich, bal=False)
            for frag in onefrags:
                _add2dic(newfrags, frag, -onefrags[frag])
            frags = {k: v for k, v in newfrags.items() if v}
            balance_ = _balance(ich, frags)
            balance_ = {k: v for k, v in balance_.items() if v}
            if balance_:
                newfrags = frags.copy()
                zedfrags = cbhzed(ich, bal=False)
                for frag in zedfrags:
                    _add2dic(newfrags, frag, zedfrags[frag])
                frags = {k: v for k, v in newfrags.items() if v}
                balance_ = _balance(ich, frags)
                balance_ = {k: v for k, v in balance_.items() if v}
                if balance_:
                    frags = _balance_frags(ich, frags)

    return frags


def cbhthree(ich, bal=True):
    ''' 
    Fragments molecule to retain each heavy-atom -- heavy-atom bond, and keep the bonds of each    atm1 b1 atm2 b2 atm3 b3 atm4 b4 atm5
    would give atm1 [b1] atm2 b2 atm3, atm1 b1 atm2 [b2] atm3 b3 atm4, atm2 b2 atm3 [b3] atm4 b4 atm5, and atm3 b3 atm4 [b4] atm5
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT 
    frags -- DIC dictionary with keys as STR inchi name for fragments and value as INT their coefficient
    '''
    #Graphical info about molecule
    gra      = automol.inchi.graph(ich)
    atms     = automol.graph.atoms(gra)
    bnd_ords = automol.graph.one_resonance_dominant_bond_orders(gra)
    rad_atms = list(automol.graph.sing_res_dom_radical_atom_keys(gra))
    atm_vals = automol.graph.atom_element_valences(gra)
    adj_atms = automol.graph.atom_neighbor_keys(gra)

    #Determine CBHfour fragments
    frags = {}

    for bnd in list(bnd_ords):
        atm_dic = {}
        bnd_dic = {}
        bnd_dic[frozenset({0, 1})] = (list(bnd_ords[bnd])[0], None)
        for i, atm in enumerate(list(bnd)):
            vali = atm_vals[atm]
            if atm in rad_atms:
                vali -= 1
            for adj in list(adj_atms[atm]):
                key     = frozenset({atm, adj})
                bnd_ord = list(bnd_ords[key] )[0]
                vali   -= bnd_ord
            atm_dic[i] = (atms[atm][0], int(vali), None)
            for j, adj in enumerate(list(adj_atms[atm]), start = 1 ):
                if adj not in list(bnd):
                    valj = atm_vals[adj]
                    if adj in rad_atms:
                        valj -= 1
                    key     = frozenset({atm, adj})
                    bnd_ord = list(bnd_ords[key] )[0]
                    valj   -= bnd_ord
                    atm_dic[i*4+j+1] = (atms[adj][0], int(valj), None)
                    bnd_dic[frozenset({i,i*4+j+1})] =  (bnd_ord, None)
        gra     = (atm_dic, bnd_dic)
        frag    = automol.graph.inchi(gra)
        _add2dic(frags, frag)

    if not frags:
        frags = cbhtwo(frags)

    if bal:
        balance_ = _balance(ich, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            newfrags = frags.copy()
            twofrags = cbhtwo(ich, bal=False)
            for frag in twofrags:
                _add2dic(newfrags, frag, -twofrags[frag])
            frags = {k: v for k, v in newfrags.items() if v}
            #balance_ = _balance(ich, frags)
            #balance_ = {k: v for k, v in balance_.items() if v}
            #if balance_:
            #    newfrags = frags.copy()
            #    newerfrags = {}
            #    onefrags = cbhone(ich, bal=False)
            #    for frag in onefrags:
            #         _add2dic(newfrags, frag, - onefrags[frag])
            #    frags = {k: v for k, v in newfrags.items() if v}
            #    balance_ = _balance(ich, frags)
            #    balance_ = {k: v for k, v in balance_.items() if v}
            #    if balance_:
            #        newfrags = frags.copy()
            #        zedfrags = cbhzed(ich, bal=False)
            #        for frag in zedfrags:
            #            _add2dic(newfrags, frag, zedfrags[frag])
            #        frags = {k: v for k, v in newfrags.items() if v}
            #        balance_ = _balance(ich, frags)
            #        balance_ = {k: v for k, v in balance_.items() if v}
            #        if balance_:
            #            frags = _balance_frags(ich, frags)
    balance = _balance(ich, frags)
    balance_ = {k: v for k, v in balance_.items() if v}
    return frags

#def cbhfour(ich):
#    ''' 
#    Fragments molecule for each heavy-atom to stay bonded to its adjacent atoms, for those atoms to say bonded to their adjacent atoms atm1 b1 atm2 b2 atm3 b3 atm4 b4 atm5
#    would give [atm1] b1 atm2 b2 atm3, atm1 b1 [atm2] b2 atm3 b3 atm4,  atm1 b1 atm2 b2  [atm3] b3 atm4 b4 amt5, atm2 b2 atm3 b3 [atm4] b4 atm5, atm3 b3 atm4 b4 [atm5]
#    INPUT:
#    ich --  STR inchi name for molecule
#    OUTPUT 
#    frags -- DIC dictionary with keys as STR inchi name for fragments and value as INT their coefficient
#    '''
#    #Graphical info about molecule
#    gra      = automol.inchi.graph(ich)
#    atms     = automol.graph.atoms(gra)
#    bnd_ords = automol.graph.one_resonance_dominant_bond_orders(gra)
#    rad_atms = list(automol.graph.sing_res_dom_radical_atom_keys(gra))
#    atm_vals = automol.graph.atom_element_valences(gra)
#    adj_atms = automol.graph.atom_neighbor_keys(gra)
#
#    #Determine CBHfour fragments
#    frags = {}
#    for atmi in atms:
#        vali = atm_vals[atmi]
#        if atmi in rad_atms:
#            vali -= 1
#        #Atoms j are adj to atoms i
#        #Loop over j one time to get the valence of i
#        for atmj in list(adj_atms[atmi]):
#            key     = frozenset({atmi, atmj})
#            bnd_ord = list(bnd_ords[key] )[0]
#            vali   -= bnd_ord
#        atm_dic = {0: (atms[atmi][0], int(vali), None)}
#        bnd_dic = {}
#        #Then Loop over j a second time to get bonds to i-j, valence of j, and bonds j-k
#        for j, atmj in enumerate(list(adj_atms[atmi]), start = 1):
#            valj = atm_vals[atmj]
#            if atmj in rad_atms:
#                valj -= 1
#            for atmk in list(adj_atms[atmj]):
#            #Loop over k first to get the valence of j
#                key     = frozenset({atmj, atmk})
#                bnd_ord = list(bnd_ords[key] )[0]
#                valj   -= bnd_ord
#            key     = frozenset({atmi, atmj})
#            bnd_ord = list(bnd_ords[key] )[0]
#
#            #Add i-j bonds and j atoms
#            atm_dic[j] = (atms[atmj][0], int(valj), None)
#            bnd_dic[frozenset({0, j})] =  (1, None)
#            
#            #Loop over k to add atoms k and j-k bonds
#            for k, atmk in enumerate(list(adj_atms[atmj]), start = 1):
#               if atmk != atmi:
#                   valk = atm_vals[atmk]
#                   if atmk in rad_atms:
#                       valk -= 1
#                   key     = frozenset({atmj, atmk})
#                   bnd_ord = list(bnd_ords[key] )[0]
#                   valk   -= bnd_ord
#                   index   = k + len(list(adj_atms[atmi]))
#                   atm_dic[index] = (atms[atmk][0], int(valk), None)
#                   bnd_dic[frozenset({j, index})] =  (1, None)
#               else:
#                   k -= 1
#
#        gra     = (atm_dic, bnd_dic)
#        frag = automol.graph.inchi(gra)
#        _add2dic(frags, frag)
#    frags =  {k: v for k, v in frags.items() if v}
#    print(frags)
#    return
#    #Balance
#    return frags

def get_basis(ich):
    """ get a basis
    """

    formula = automol.inchi.formula_string(ich)
    atm_dict = util.get_atom_counts_dict(formula)
    return select_basis(atm_dict)


def get_basic(ich):
    """ get basis for basic scheme
    """
    formula_dct = automol.inchi.formula(ich)
    spc_bas = select_basis(formula_dct)
    if len(spc_bas) == 1 and ich == spc_bas[0]:
        clist = [1]
    else:
        clist = calc_coefficients(spc_bas, formula_dct)
    return spc_bas, clist


def _intersec(lst1, lst2):
   ret = None
   for atm in lst1:
       if atm in lst2:
           ret = atm
   assert ret is not None, (
       'brk_key {} and frm_key {} do not intersect'.format(lst1, lst2))
   return ret


def _xor(lst1, lst2):
   ret = None
   for atm in lst1:
       if atm not in lst2:
           ret = atm
   assert ret is not None, (
       'problem with bond_key'.format(lst1))
   return ret


def _remove_dummies(zma, frm_key, brk_key):
    """get zma and bond key idxs without dummy atoms
    """
    geo = automol.zmatrix.geometry(zma)
    dummy_idxs = automol.geom.dummy_atom_indices(geo)
    for idx in dummy_idxs:
        if frm_key:
            frm1, frm2 = frm_key
            if idx < frm1:
                frm1 -= 1
            if idx < frm2:
                frm2 -= 1
            frm_key = frozenset({frm1, frm2})    
        if brk_key:
            brk1, brk2 = brk_key
            if idx < brk1:
                brk1 -= 1
            if idx < brk2:
                brk2 -= 1
            brk_key = frozenset({brk1, brk2})   
    if dummy_idxs:        
        geo = automol.geom.without_dummy_atoms(geo)
        zma = automol.geom.zmatrix(geo, [frm_key, brk_key])
        atm_ord = automol.geom.zmatrix_atom_ordering(geo, [frm_key, brk_key])
        frm_key = frozenset({atm_ord[atm] for atm in frm_key})    
        brk_key = frozenset({atm_ord[atm] for atm in brk_key})    
    return zma, frm_key, brk_key

def _remove_frm_bnd(gra, brk_key, frm_key):
    bond_keys = automol.graph.bond_keys(gra)
    if brk_key not in bond_keys:
        gra = automol.graph.add_bonds(gra, [brk_key])
    if frm_key in bond_keys:
        gra = automol.graph.remove_bonds(gra, [frm_key])
    return gra    

def _add_appropriate_pi_bonds(gra):
    adj_atms = automol.graph.atom_neighbor_keys(gra)
    unsat_atms_dct = automol.graph.atom_unsaturated_valences(gra)
    atms, bnd_ords = gra
    brk_key = frozenset({})
    unsat_atms = []
    for atm in unsat_atms_dct:
        if unsat_atms_dct[atm] > 0:
            unsat_atms.append(atm)
    for atmi in unsat_atms:
        for atmj in unsat_atms:
            if atmi > atmj:
                print(atmi, atmj)
                if atmi in adj_atms[atmj]:
                    print('are pi bond')
                    key = [atmi, atmj]
                    key.sort()
                    key = frozenset(key)
                    brk_key = key
                    bnd, tmp = bnd_ords[key]
                    bnd_ords[key] = (bnd + 1, tmp)

    return (atms, bnd_ords), brk_key

def get_cbhzed_ts(zma, rxnclass, frm_key, brk_key):
    """ get basis for CBH0 for a TS molecule
    """
    zma, frm_key, brk_key = _remove_dummies(zma, frm_key, brk_key)
    print('bond info ', frm_key, brk_key, rxnclass)        
    gra = automol.zmatrix.connectivity_graph(zma)
    if 'addition' in rxnclass:
        gra, brk_key = _add_appropriate_pi_bonds(gra)
    print('gra in cbhzed',  gra)
    gra = _remove_frm_bnd(gra, brk_key, frm_key)
    if frm_key and brk_key:
        site = [_xor(frm_key, brk_key), _intersec(frm_key, brk_key), _xor(brk_key, frm_key)]
    elif 'beta scission' in rxnclass:
        rad_atm = list(automol.graph.sing_res_dom_radical_atom_keys(gra))[0]
        adj_atms = automol.graph.atom_neighbor_keys(gra)
        #site = [None,None,rad_atm]
        site = [rad_atm,None,None]
        for atm in brk_key:
            if rad_atm in adj_atms[atm]:
                site[1] =  atm
            else:
                #site[0] = atm
                site[2] = atm
        print('site ', site)        
        adj_atms = automol.graph.atom_neighbor_keys(gra)
    elif 'addition' in rxnclass:
        print('made it to 945')

    if 'hydrogen abstraction' in rxnclass or 'beta scission' in rxnclass or 'hydrogen migration' in rxnclass or 'addition' in rxnclass:
        frags = cbhzed_habs(gra, site)
    else:
        raise NotImplementedError
    fraglist = []
    clist = []
    for frag in frags:
        if 'exp_gra' in frags[frag]:
            fraglist.append(automol.graph.inchi(frags[frag]['exp_gra']))
            clist.append(frags[frag]['coeff'])
            print(automol.geom.string(automol.graph.geometry(frags[frag]['exp_gra'])))
        else:
            if 'beta' in rxnclass:
                fraglist.append(split_beta_gras(frags[frag]['ts_gra']))
            else:    
                fraglist.append(split_gras(frags[frag]['ts_gra']))
            clist.append(frags[frag]['coeff'])
    print('frags from cbh_ts ', fraglist)        
    return fraglist, clist
   
def get_cbhzed(ich):
    """ get basis for CBH0
    """
    frags = cbhzed(ich)
    fraglist = []
    clist = []
    for frag in frags:
        fraglist.append(frag)
        clist.append(frags[frag])
    return fraglist, clist
    # return list(cbhzed(ich).keys())


def get_cbhone(ich):
    """ get basis for CBH1
    """
    frags = cbhone(ich)
    fraglist = []
    clist = []
    for frag in frags:
        fraglist.append(frag)
        clist.append(frags[frag])
    return fraglist, clist
    # return list(cbhone(ich).keys())


def get_cbhtwo(ich):
    """ get basis for CBH2
    """
    frags = cbhtwo(ich)
    fraglist = []
    clist = []
    for frag in frags:
        fraglist.append(frag)
        clist.append(frags[frag])
    return fraglist, clist
    # return list(cbhtwo(ich).keys())

def get_cbhthree(ich):
    """ get basis for CBH3
    """
    frags = cbhthree(ich)
    fraglist = []
    clist = []
    for frag in frags:
        fraglist.append(frag)
        clist.append(frags[frag])
    return fraglist, clist

def _add2dic(dic, key, val=1):
    """ helper function to add a key to dct
    """
    if key in dic:
        dic[key] += val
    else:
        dic[key] = val


def _lhs_rhs(frags):
    """ Determine the left-hand side and right-hand side of reaction
    """
    rhs = {}
    lhs = {}
    for frag in frags:
        if frags[frag] > 0:
            rhs[frag] = frags[frag]
        elif frags[frag] < 0:
            lhs[frag] = - frags[frag]
    return lhs, rhs


def _print_lhs_rhs(ich, frags):
    """ print the fragments from each side of the reaction
    """
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
    """ balance the equation?
    """
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


def _balance_ts(gra, frags):
    """ balance the equation using graphs
    """
    stoichs = {}
    for frag in frags:
        if 'exp_gra' in frags[frag]:
            _stoich = automol.graph.formula(frags[frag]['exp_gra'])
            print('stoichs of frag', _stoich)
        elif 'ts_gra' in frags[frag]:
            _stoich = stoich_gra(frags[frag]['ts_gra'])
            print('stoichs of frag', _stoich)
        for atm in _stoich:
            if atm in stoichs:
                stoichs[atm] += _stoich[atm] * frags[frag]['coeff']
            else:
                stoichs[atm] = _stoich[atm] * frags[frag]['coeff']
    balance_ = {}
    _stoich = stoich_gra(gra)
    for atom in _stoich:
        if atom in stoichs:
            balance_[atom] = _stoich[atom] - stoichs[atom]
        else:
            balance_[atom] = _stoich[atom]
    balance_ = {x: y for x, y in balance_.items() if y != 0}
    return balance_


def _balance_frags(ich, frags):
    """ balance the equation?
    """
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


def _balance_frags_ts(gra, frags):
    """ balance the equation?
    """
    balance_ = _balance_ts(gra, frags)
    methane = automol.smiles.inchi('C')
    water = automol.smiles.inchi('O')
    ammonm = automol.smiles.inchi('N')
    hydrgn = automol.smiles.inchi('[H][H]')
    methane = automol.inchi.graph(methane)
    water = automol.inchi.graph(water)
    ammonm = automol.inchi.graph(ammonm)
    hydrgn = automol.inchi.graph(hydrgn)
    idx_dct = []
    for spc in [methane, water, ammonm, hydrgn]:
        spc = automol.graph.explicit(spc)
        found = False
        for frag in frags:
            if 'exp_gra' in frags[frag]:
                if automol.graph.full_isomorphism(frags[frag]['exp_gra'], spc):
                    idx = frag
                    found = True
                    break
        if not found:        
            idx = len(frags.keys())
            frags[idx] = {}
            frags[idx]['exp_gra'] =  spc
            frags[idx]['coeff'] = 0.0
        idx_dct.append(idx)
    if 'C' in balance_:
        _add2dic(frags[idx_dct[0]], 'coeff', balance_['C'])
    if 'N' in balance_:
        _add2dic(frags[idx_dct[1]], 'coeff', balance_['N'])
    if 'O' in balance_:
        _add2dic(frags[idx_dct[2]], 'coeff', balance_['O'])
    balance_ = _balance_ts(gra, frags)
    if 'H' in balance_:
        _add2dic(frags[idx_dct[3]], 'coeff', balance_['H']/2)
    return frags



