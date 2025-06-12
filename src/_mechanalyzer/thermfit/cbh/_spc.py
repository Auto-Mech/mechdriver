""" Species CBH schemes
    Explanation of various basis determination schemes
"""

import numpy
import automol.chi
import automol.graph
import automol.form
from thermfit.cbh import _util as util


# Main callable function
def species_basis(ich, scheme, balance=True):
    """ Get the basis for the appropriate CBH scheme

        :param ich: InChI string for spc
        :type ich: str
        :param scheme: CBH Scheme used to generate basis
        :type scheme: str
        :param balance:
        :type balance: bool
        :rtype: (tuple(str), list(float))
    """

    if scheme == 'basic':
        frag_lst, coeff_lst = basic_spc_basis(ich)
    else:
        frags = CBH_SCHEMES[scheme](ich, balance=balance)
        frag_lst, coeff_lst = (), ()
        for frag in frags:
            frag_lst += (frag,)
            coeff_lst += (frags[frag],)

    return frag_lst, coeff_lst


# Basic calculator
def basic_spc_basis(ich):
    """ Determine a basis for relative enthalpy calculations for species
        by building a list of a simple species in a inear combination which
        reproduces the number of atoms in the stoichiometry of the species.

        Basis represented by list of InChI strings and array of coefficients.

        :param ich: InChI string for spc
        :type ich: str
        :rtype: (tuple(str), numpy.ndarray)
    """

    # Get a list of all the atom types in the molecule
    fml = automol.chi.formula(ich)
    symbs = tuple(fml.keys())

    # Create list of inchi keys corresponding to basis species
    basis = ()
    # H2
    basis += ('InChI=1S/H2/h1H',)
    # CH4
    if 'C' in symbs:
        basis += ('InChI=1S/CH4/h1H4',)
    # H2O
    if 'O' in symbs:
        basis += ('InChI=1S/H2O/h1H2',)
    # NH3
    if 'N' in symbs:
        basis += ('InChI=1S/H3N/h1H3',)
    # Cl2
    if 'Cl' in symbs:
        basis += ('InChI=1S/ClH/h1H',)
        # basis += ('InChI=1S/Cl2/c1-2',)
    # SO2
    if 'S' in symbs:
        basis += ('InChI=1S/O2S/c1-3-2',)
        if 'O' not in symbs:
            basis += ('InChI=1S/H2O/h1H2',)

    # Generate the coefficients for the basis
    if len(basis) == 1 and ich == basis[0]:
        coeff_lst = (1.0,)
    else:
        coeff_lst = _coefficients(basis, fml)

    return basis, coeff_lst


def _coefficients(basis, spc_fml):
    """ Form a matrix consisting of the coefficients for the basis
        species to balance out the atoms.

        :param basis: InChI strings for basis species
        :type basis: tuple(str)
        :param spc_fml: stoichiometric formula of species basis corresponds to
        :type spc_fml: dict[str:int]
        :rtype: numpy.ndarray
    """

    # Initialize an natoms x natoms matrix
    nbasis = len(basis)
    basis_mat = numpy.zeros((nbasis, nbasis))

    # Get the basis formulae list
    basis_fml_str = [automol.chi.formula_layer(spc) for spc in basis]
    for spc in basis_fml_str:
        basis_atom_dict = automol.form.from_string(spc)
        for atom in basis_atom_dict:
            if atom not in spc_fml:
                spc_fml[atom] = 0

    # Set the elements of the matrix
    for i, spc in enumerate(basis_fml_str):
        basis_atom_dict = automol.form.from_string(spc)
        basis_vals = []
        for key in spc_fml.keys():
            if key in basis_atom_dict:
                basis_vals.append(basis_atom_dict[key])
            else:
                basis_vals.append(0)
        basis_mat[i] = basis_vals

    #  Transpose
    basis_mat = basis_mat.T

    # Form stoichometry vector vector
    stoich_vec = numpy.zeros(len(spc_fml))
    for i, key in enumerate(spc_fml.keys()):
        stoich_vec[i] = spc_fml[key]

    # Solve C = B^-1 S
    basis_mat = numpy.linalg.inv(basis_mat)
    coeff_vec = numpy.dot(basis_mat, stoich_vec)

    return coeff_vec


# Individual CBH-N calculators
def cbhzed(ich, balance=True):
    """
    Fragments molecule so that each heavy-atom is a seperate fragment
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchi name for fragments
    and value as INT their coefficient
    """

    # Graphical info about molecule
    gra = automol.chi.graph(ich)
    atms = automol.graph.atoms(gra)
    adj_atms = automol.graph.atoms_neighbor_atom_keys(gra)
    term_atms = automol.graph.terminal_atom_keys(gra)

    kek_bnd_ords = automol.graph.kekules_bond_orders(gra)
    norm_kek = 1 / len(kek_bnd_ords)

    # Determine CBHzed fragments
    frags = {}
    for bnd_ords in kek_bnd_ords:
        for atm in atms:
            grai = (
                atms.copy(),
                {key: (val, None) for (key, val) in bnd_ords.copy().items()},)
            if (atms[atm][0] == 'H' and len(atms) > 1):
                continue
            coeff = 1.0
            if not balance:
                coeff = (
                    util.branch_point(adj_atms[atm]) *
                    (atm not in term_atms)
                )
            extended_site = [atm]
            for site_atm in extended_site:
                for atm_x in adj_atms[site_atm]:
                    if atm_x not in extended_site and atms[atm_x][0] != 'H':
                        grai = util.cleave_group_and_saturate(
                            grai, bnd_ords, site_atm, atm_x)
            frag = automol.graph.chi(grai)
            util.add2dic(frags, frag, val=coeff*norm_kek)
    frags = {k: round(v, 6) for k, v in frags.items() if v}
    if balance:
        balance_ = util.balance(ich, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            frags = util.balance_frags(ich, frags)
    frags = {
        k: round(v, 6) for k, v in frags.items() if abs(round(v, 6)) != 0.0}
    return frags


def cbhone(ich, balance=True):
    """
    Fragments molecule in a way that conserves each heavy-atom/heavy-atom bond
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchi name for fragments
    and value as INT their coefficient
    """

    # Graphical info about molecule
    gra = automol.chi.graph(ich)
    atms = automol.graph.atoms(gra)
    adj_atms = automol.graph.atoms_neighbor_atom_keys(gra)

    kek_bnd_ords = automol.graph.kekules_bond_orders(gra)
    norm_kek = 1 / len(kek_bnd_ords)

    # Determine CBHone fragments
    frags = {}
    for bnd_ords in kek_bnd_ords:
        for bnd in bnd_ords:
            atma, atmb = bnd
            grai = (
                atms.copy(),
                {key: (val, None) for (key, val) in bnd_ords.copy().items()},)
            if (atms[atma][0] == 'H' or atms[atmb][0] == 'H'):
                continue
            coeff = 1.0
            extended_site = [atma, atmb]
            for site_atm in extended_site:
                for atm_x in adj_atms[site_atm]:
                    if atm_x not in extended_site and atms[atm_x][0] != 'H':
                        grai = util.cleave_group_and_saturate(
                            grai, bnd_ords, site_atm, atm_x)
            grai = automol.graph.explicit(grai)
            frag = automol.graph.chi(grai)
            util.add2dic(frags, frag, val=coeff*norm_kek)
    frags = {k: round(v, 6) for k, v in frags.items() if v}
    if not frags:
        frags = cbhzed(ich)
    # Balance
    if balance:
        balance_ = util.balance(ich, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            zedfrags = cbhzed(ich, balance=False)
            newfrags = frags.copy()
            for frag in zedfrags:
                util.add2dic(newfrags, frag, -zedfrags[frag])
            frags = {k: v for k, v in newfrags.items() if v}
            balance_ = util.balance(ich, frags)
            balance_ = {k: v for k, v in balance_.items() if v}
            assert all([v == 0 for v in balance_.values()]), \
                f"CBH0 fails to balance CBH1 for {ich} -- " + ",".join([f'{k}:{v}' for k, v in balance_.items()])

    frags = {
        k: round(v, 6) for k, v in frags.items() if abs(round(v, 6)) != 0.0}
    return frags


def cbhtwo(ich, balance=True):
    """
    Fragments molecule for each heavy-atom to stay bonded to its adjacent atoms
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchi name for fragments and
    value as INT their coefficient
    """

    # Graphical info about molecule
    gra = automol.chi.graph(ich)
    atms = automol.graph.atoms(gra)
    adj_atms = automol.graph.atoms_neighbor_atom_keys(gra)
    term_atms = automol.graph.terminal_atom_keys(gra)

    kek_bnd_ords = automol.graph.kekules_bond_orders(gra)
    norm_kek = 1 / len(kek_bnd_ords)

    # Determine CBHtwo fragments
    frags = {}
    for bnd_ords in kek_bnd_ords:
        for atm in atms:
            grai = (
                atms.copy(),
                {key: (val, None) for (key, val) in bnd_ords.copy().items()},)
            if (atms[atm][0] == 'H'):
                continue
            coeff = 1.0
            if not balance:
                coeff = (
                    util.branch_point(adj_atms[atm]) *
                    (atm not in term_atms) 
                )
            extended_site = [atm] + list(adj_atms[atm])
            #for site_atm in extended_site:
            for site_atm in list(adj_atms[atm]):
                for atm_x in adj_atms[site_atm]:
                    if atm_x != atm and atms[atm_x][0] != 'H':
                        grai = util.cleave_group_and_saturate(
                            grai, bnd_ords, site_atm, atm_x)

            frag = automol.graph.chi(grai)
            util.add2dic(frags, frag, val=coeff*norm_kek)

    frags = {k: round(v, 6) for k, v in frags.items() if v}
    if not frags:
        frags = cbhone(ich)

    # Balance
    if balance:
        balance_ = util.balance(ich, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            newfrags = frags.copy()
            onefrags = cbhone(ich, balance=False)
            for frag in onefrags:
                util.add2dic(newfrags, frag, -onefrags[frag])
            frags = {k: round(v, 6) for k, v in newfrags.items() if v}
            balance_ = util.balance(ich, frags)
            balance_ = {k: v for k, v in balance_.items() if v}
            assert all([v == 0 for v in balance_.values()]), \
                "CBH1 fails to balance CBH2 -- " + ",".join([f'{k}:{v}' for k, v in balance_.items()])
    frags = {
        k: round(v, 6) for k, v in frags.items() if abs(round(v, 6)) != 0.0}
    return frags


def cbhthree(ich, balance=True):
    """
    Fragments molecule to retain each heavy-atom -- heavy-atom bond,
    and keep the bonds of each    atm1 b1 atm2 b2 atm3 b3 atm4 b4 atm5
    would give atm1 [b1] atm2 b2 atm3, atm1 b1 atm2 [b2] atm3 b3 atm4,
    atm2 b2 atm3 [b3] atm4 b4 atm5, and atm3 b3 atm4 [b4] atm5
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchi name
    for fragments and value as INT their coefficient
    """

    # Graphical info about molecule
    gra = automol.chi.graph(ich)
    atms = automol.graph.atoms(gra)
    bnd_ords = automol.graph.kekule_bond_orders(gra)
    adj_atms = automol.graph.atoms_neighbor_atom_keys(gra)

    # Determine CBHfour fragments
    frags = {}

    for bnd in bnd_ords:
        atma, atmb = bnd
        grai = (
            atms.copy(),
            {key: (val, None) for (key, val) in bnd_ords.copy().items()},)
        if (atms[atma][0] == 'H' or atms[atmb][0] == 'H'):
            continue
        coeff = 1.0
        extended_site = [atma, atmb] + list(adj_atms[atma]) + list(adj_atms[atmb])
        for site_atm in extended_site:
            for atm_x in adj_atms[site_atm]:
                if atm_x not in extended_site and atms[atm_x][0] != 'H':
                    grai = util.cleave_group_and_saturate(
                        grai, bnd_ords, site_atm, atm_x)
        grai = automol.graph.explicit(grai)
        frag = automol.graph.chi(grai)
        util.add2dic(frags, frag, val=coeff)
    frags = {k: v for k, v in frags.items() if v}
    if not frags:
        frags = cbhtwo(ich)

    if balance:
        balance_ = util.balance(ich, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            newfrags = frags.copy()
            twofrags = cbhtwo(ich, balance=False)
            for frag in twofrags:
                util.add2dic(newfrags, frag, -twofrags[frag])
            frags = {k: v for k, v in newfrags.items() if v}
    return frags

# def cbhfour(ich):
#    """
#    Fragments molecule for each heavy-atom to stay bonded to its
#    adjacent atoms, for those atoms to say bonded to their adjacent atoms
#    atm1 b1 atm2 b2 atm3 b3 atm4 b4 atm5
#    would give [atm1] b1 atm2 b2 atm3, atm1 b1 [atm2] b2 atm3 b3 atm4,
#    atm1 b1 atm2 b2  [atm3] b3 atm4 b4 amt5, atm2 b2 atm3 b3 [atm4] b4 atm5,
#    atm3 b3 atm4 b4 [atm5]
#    INPUT:
#    ich --  STR inchi name for molecule
#    OUTPUT
#    frags -- DIC dictionary with keys as STR inchi name for fragments and
#    value as INT their coefficient
#    """
#    #Graphical info about molecule
#    gra      = automol.chi.graph(ich)
#    atms     = automol.graph.atoms(gra)
#    bnd_ords = automol.graph.kekule_bond_orders(gra)
#    rad_atms = list(automol.graph.radical_atom_keys(gra, sing_res=True))
#    atm_vals = automol.graph.atomic_valences(gra)
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
#        # Then Loop over j a second time to get bonds
#        # to i-j, valence of j, and bonds j-k
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
#        frag = automol.graph.chi(gra)
#        util.add2dic(frags, frag)
#    frags =  {k: v for k, v in frags.items() if v}
#    ioprinter.info_message(frags)
#    return
#    #Balance
#    return frags


# Dictionary of schemes
CBH_SCHEMES = {
    'cbh0': cbhzed,
    'cbh1': cbhone,
    'cbh0_0': cbhzed,
    'cbh1_0': cbhone,
    'cbh1_1': cbhone,
    'cbh2': cbhtwo,
    'cbh2_0': cbhtwo,
    'cbh2_1': cbhtwo,
    'cbh2_2': cbhtwo,
    'cbh3': cbhthree
}
