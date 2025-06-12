""" Add a species or TS to your database
usiing a log or xyz file
"""
from pathlib import Path
import sys
import autofile
import automol
from mechanalyzer.inf import thy as tinfo
from mechanalyzer.inf import rxn as rinfo
from mechanalyzer.inf import spc as sinfo
import elstruct
from automol.geom import ring_fragments_geometry as _fragment_ring_geo
from mechlib.filesys import save
import ioformat.pathtools
from mechroutines.es._routines.conformer import _saved_cnf_info
from mechroutines.es._routines.conformer import _sym_unique
from mechroutines.es._routines.conformer import _geo_unique

THEORY_DCT = {
    'lvl_wbs': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'wb97xd',
        'basis':  '6-31g*'},
    'lvl_wbm': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'wb97xd',
        'basis':  '6-31+g*'},
    'lvl_wbt': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'wb97xd',
        'basis':  'cc-pvtz'},
    'lvl_m06s': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'm062x',
        'basis':  '6-31g*'},
    'lvl_m06m': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'm062x',
        'basis':  '6-31+g*'},
    'lvl_m06t': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'm062x',
        'basis':  'cc-pvtz'},
    'lvl_b2d': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'b2plypd3',
        'basis':  'cc-pvdz'},
    'lvl_b2t': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'b2plypd3',
        'basis':  'cc-pvtz'},
    'lvl_b2q': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'b2plypd3',
        'basis':  'cc-pvqz'},
    'lvl_b3s': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'b3lyp',
        'basis':  '6-31g*'},
    'lvl_b3mg': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'b3lyp',
        'basis':  '6-311g**'},
    'lvl_b3t': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'b3lyp',
        'basis':  'cc-pvtz'},
    'cc_lvl_d': {
        'orb_res': 'RR',
        'program': 'molpro2015',
        'method':  'ccsd(t)',    
        'basis':  'cc-pvdz'},
    'cc_lvl_t': {
        'orb_res': 'RR',
        'program': 'molpro2015',
        'method':  'ccsd(t)',    
        'basis':  'cc-pvtz'},
    'cc_lvl_q': {
        'orb_res': 'RR',
        'program': 'molpro2015',
        'method':  'ccsd(t)',
        'basis':  'cc-pvqz'},
    'cc_lvl_df': {
        'orb_res': 'RR',
        'program': 'molpro2015',
        'method':  'ccsd(t)-f12',
        'basis':  'cc-pvdz-f12'},
    'cc_lvl_tf': {
        'orb_res': 'RR',
        'program': 'molpro2015',
        'method':  'ccsd(t)-f12',
        'basis':  'cc-pvtz-f12'},
    'cc_lvl_qf': {
        'orb_res': 'RR',
        'program': 'molpro2015',
        'method':  'ccsd(t)-f12',
        'basis':  'cc-pvqz-f12'},
    'mlvl_cas_dz': {
        'orb_res': 'RR',
        'program': 'molpro2015',
        'method':  'caspt2',
        'basis':  'cc-pvdz'},
    'mlvl_cas_tz': {
        'orb_res': 'RR',
        'program': 'molpro2015',
        'method':  'caspt2',
        'basis':  'cc-pvtz'}}


def parse_user_locs(insert_dct, geo, cnf_fs):
    """ read rid/cid from input dictionary 
        or generate new one, new one will be 
        consistent with rids already in FS
    """
    rid = insert_dct.get('rid', None)
    cid = insert_dct.get('cid', None)
    if rid is None:
        rid = rng_loc_for_geo(geo, cnf_fs)
    if rid is None:
        rid = autofile.schema.generate_new_ring_id()
    if cid is None:
        cid = autofile.schema.generate_new_conformer_id()
    return (rid, cid)


def parse_user_species(insert_dct):
    """ Make species info tuple (later used to make species leaf
        from insert dictionary, if info is missing break unless 
        'trusted' flag is on, which means it will be generated later
    """
    smi = insert_dct.get('smiles', None)
    ich = insert_dct.get('inchi', None)
    mult = insert_dct.get('mult', None)
    chg = insert_dct.get('charge', None)
    if ich is None and smi is None and not insert_dct['trusted']:
        print(
            'Error: user did not specify species' +
            'with an inchi or smiles in input')
        sys.exit()
    if ich is None and smi is not None:
        ich = automol.smiles.chi(smi)
    if ich is not None:
        if not automol.chi.is_complete(ich):
            ich = automol.chi.add_stereo(ich)
    if mult is None and not insert_dct['trusted']:
        print('Error: user did not specify mult in input')
        sys.exit()
    if chg is None and not insert_dct['trusted']:
        print('Error: user did not specify charge in input')
        sys.exit()
    return sinfo.from_data(ich, chg, mult)


def species_info_from_input(geo, spc_info, inp_str=None):
    """ Determine species info (inchi, charge, mult) from
        the geometry or input file TODO
    """
    spc_info = list(spc_info)
    if spc_info[0] is None:
        spc_info[0] = automol.geom.chi(geo, stereo=True)
    if spc_info[2] is None:
        #read from inp_str
        mults_allowed = automol.graph.possible_spin_multiplicities(
             automol.chi.graph(spc_info[0], stereo=False))
        spc_info[2] = mults_allowed[0]
    if spc_info[1] is None:
        #read from inp_str
        spc_info[1] = 0
    return tuple(spc_info)


def parse_user_reaction(insert_dct):
    """ read the species info for reactant and products
        from the insert dictionary
    """
    smis = insert_dct.get('smiles', None)
    ichs = insert_dct.get('inchi', None)
    mults = insert_dct.get('mult', None)
    chgs = insert_dct.get('charge', None)
    rxn_class = insert_dct.get('rxn_class', None)
    if ichs is None:
        ichs = [[], []]
        for smi in smis[0]:
            ichs[0].append(automol.smiles.chi(smi))
        for smi in smis[1]:
            ichs[1].append(automol.smiles.chi(smi))
    for idx, ich in enumerate(ichs[0]):
        if not automol.chi.is_complete(ich):
            ich = automol.chi.add_stereo(ich)
            ichs[0][idx] = ich
    for idx, ich in enumerate(ichs[1]):
        if not automol.chi.is_complete(ich):
            ich = automol.chi.add_stereo(ich)
            ichs[1][idx] = ich
    if mults is None:
        print('Error: user did not specify mults in input')
        sys.exit()
    if chgs is None:
        print('Error: user did not specify charges in input')
        sys.exit()
    ichs = list(map(list, ichs))
    mults = list(map(list, mults))
    chgs = list(map(list, chgs))
    flat_ichs = sum(ichs, [])
    if len(flat_ichs) != len(sum(mults, [])):
        print(
            'Error: number of species does not match number of mults')
        sys.exit()
    if len(flat_ichs) != len(sum(chgs, [])):
        print(
            'Error: number of species does not match number of charges')
        sys.exit()
    idx = 0
    rxn_muls = [[], []]
    rxn_chgs = [[], []]
    for ridx, ich in enumerate(ichs[0]):
        mults_allowed = automol.graph.possible_spin_multiplicities(
            automol.chi.graph(ich, stereo=False))
        if mults[0][ridx] not in mults_allowed:
            print(
                'user specified mult of {}'.format(mults[idx]) +
                'is not an allowed multiplicty for inchi {}'.format(ich))
            sys.exit()
        rxn_muls[0].append(mults[0][ridx])
        rxn_chgs[0].append(chgs[0][ridx])
        idx += 1
    for pidx, ich in enumerate(ichs[1]):
        mults_allowed = automol.graph.possible_spin_multiplicities(
            automol.chi.graph(ich, stereo=False))
        if mults[1][pidx] not in mults_allowed:
            print(
                'user specified mult of {}'.format(mults[idx]) +
                'is not an allowed multiplicty for inchi {}'.format(ich))
            sys.exit()
        rxn_muls[1].append(mults[1][pidx])
        rxn_chgs[1].append(chgs[1][ridx])
        idx += 1
    ts_mult = insert_dct['ts_mult']
    if ts_mult is None:
        print(
            'Error: user did not specify ts_mul')
        sys.exit()
    rxn_info = rinfo.sort((ichs, rxn_chgs, rxn_muls, ts_mult))
    ts_info = rinfo.ts_info(rxn_info)
    if rxn_class is None:
        print(
            'Error: user did not specify rxn_class')
        sys.exit()
    return rxn_info, ts_info, rxn_class


def parse_user_theory(insert_dct):
    """ Generate theory info tuple (later used to make theory leaf)
        from insert dictionary, they can have specified program/method/basis
        or have used a common 'theory' key (see THEORY_DCT above)
    """
    # Get input method explicitly inputted
    program = insert_dct.get('program', None)
    method = insert_dct.get('method', None)
    basis = insert_dct.get('basis', None)
    orb_res = insert_dct.get('orb_res', None)
    # Get input method from theory dictionary
    theory = insert_dct.get('theory', None)
    if theory is None:
        if program is None:
            print('Error: user did not specify program in input')
            sys.exit()
        elif method is None:
            print('Error: user did not specify method in input')
            sys.exit()
        elif basis is None:
            print('Error: user did not specify basis in input')
            sys.exit()
        elif orb_res is None:
            print('Error: user did not specify orb_res in input')
            sys.exit()
        else:
            thy_info = (program, method, basis, orb_res)
    else:
        if theory in THEORY_DCT:
            thy_info = tinfo.from_dct(THEORY_DCT[theory])
        else:
            print(
                'Error: user did not specify a theory {}'.format(theory) +
                ' that is in the THEORY_DCT' +
                'please add it to the dct in the script or use program/method/basis/orb_dct' +
                'keywords instead of theory')
            sys.exit()
    return thy_info


def create_species_filesystems(prefix, spc_info, mod_thy_info, locs=None):
    """ create SPC FS using species info, theory info,
        and rid/cid specification
    """
    # species filesystem
    spc_fs = autofile.fs.species(prefix)
    print(spc_info, spc_fs[0].path())
    spc_fs[-1].create(spc_info)
    spc_prefix = spc_fs[-1].path(spc_info)

    # theory filesystem
    thy_fs = autofile.fs.theory(spc_prefix)
    print(mod_thy_info[1:])
    thy_fs[-1].create(mod_thy_info[1:])
    thy_prefix = thy_fs[-1].path(mod_thy_info[1:])

    # conformer
    cnf_fs = autofile.fs.conformer(thy_prefix)
    if locs is not None:
        cnf_fs[-1].create(locs)
        cnf_prefix = cnf_fs[-1].path(locs)
    else:
        cnf_prefix = None

    return (
        (spc_fs, thy_fs, cnf_fs), (spc_prefix, thy_prefix, cnf_prefix))


def create_reaction_filesystems(
        prefix, rxn_info, mod_thy_info, ts_locs=None, locs=None):
    """ create RXN FS using rxn info, theory info,
        and ts_locs/rid/cid specification
    """
    # species filesystem
    rxn_fs = autofile.fs.reaction(prefix)
    sort_rxn_info = rinfo.sort(rxn_info, scheme='autofile')
    rxn_fs[-1].create(sort_rxn_info)
    rxn_prefix = rxn_fs[-1].path(sort_rxn_info)

    # theory filesystem
    thy_fs = autofile.fs.theory(rxn_prefix)
    thy_fs[-1].create(mod_thy_info[1:])
    thy_prefix = thy_fs[-1].path(mod_thy_info[1:])

    if ts_locs is None:
        ts_locs = (0,)

    ts_fs = autofile.fs.transition_state(thy_prefix)
    ts_fs[-1].create(ts_locs)
    ts_prefix = ts_fs[-1].path(ts_locs)

    # conformer
    cnf_fs = autofile.fs.conformer(ts_prefix)
    if locs is not None:
        cnf_fs[-1].create(locs)
        cnf_prefix = cnf_fs[-1].path(locs)
    else:
        cnf_prefix = None

    return (
        (rxn_fs, thy_fs, ts_fs, cnf_fs),
        (rxn_prefix, thy_prefix, ts_prefix, cnf_prefix))


def read_user_file(dct, keyword):
    """ read a file specified in the insert dictionary 
    """
    if dct[keyword] is None:
        print(
            'ERROR: No filename is specified for {}'.format(keyword) +
            'Script will exit')
        sys.exit()
    file_name = dct[keyword]
    return autofile.io_.read_file(file_name)


def read_user_filesystem(dct):
    """ break if user has not specified a location for the FS 
    """
    if dct['save_filesystem'] is None:
        print(
            'ERROR: No save_filesystem}' +
            'Script will exit')
        sys.exit()
    return dct['save_filesystem']


def choose_beta_qh_cutoff_distance(geo):
    """ Havent tested for years, originally used to automatically
        identify a heavy atom-H atom beta scission reaction 
        based on bond lengths of the geo
    """
    rqhs = [x * 0.1 for x in range(24, 38, 2)]
    double_bnd_atms = []
    h_atm = None
    conn_ts_gras = None
    good_ts_gra = None
    for rqh in rqhs:
        if len(double_bnd_atms) == 2:
            good_ts_gra = conn_ts_gras
            break
        ts_gras = automol.geom.connectivity_graph_deprecated(geo, rqq_bond_max=3.5, rqh_bond_max=rqh, rhh_bond_max=2.3)
        conn_ts_gras = automol.graph.set_stereo_from_geometry(ts_gras, geo)
        ts_gras = automol.graph.connected_components(conn_ts_gras)
        if len(ts_gras) != 2:
            continue
        unconn_bnds = []
        for ts_gra_i in ts_gras:
            vals = automol.graph.atom_unpaired_electrons(ts_gra_i)
            double_bnd_atms_i = [atm for atm in vals if vals[atm] == 1]
            if len(double_bnd_atms_i) == 2:
                double_bnd_atms = double_bnd_atms_i
            if len(double_bnd_atms_i) == 1:
                h_atm = double_bnd_atms_i[0]
            unconn_bnds.extend(list(automol.graph.bond_keys(ts_gra_i)))
    chosen_ts_gra = []
    chosen_oversaturated_atom = None
    rqhs = [x * 0.1 for x in range(26, 42, 2)]
    for rqh in rqhs:
        ts_gras = automol.geom.connectivity_graph_deprecated(geo, rqq_bond_max=3.5, rqh_bond_max=rqh, rhh_bond_max=2.3)
        ts_gras = automol.graph.set_stereo_from_geometry(ts_gras, geo)
        ts_gras = automol.graph.connected_components(ts_gras)
        if len(ts_gras) != 1:
            continue
        breaking_bond = None
        forming_bond = None
        for ts_gra_i in ts_gras:
            oversaturated_atoms = []
            conn_bnds = automol.graph.bond_keys(ts_gra_i)
            breaking_bond_lst = set(conn_bnds) - set(unconn_bnds)
            for bnd in breaking_bond_lst:
                if h_atm in bnd and any(atmi in double_bnd_atms for atmi in bnd):
                    breaking_bond = bnd
                    forming_bond = frozenset(double_bnd_atms)
            if breaking_bond is not None:
                chosen_ts_gra = automol.graph.set_bond_orders(
                    good_ts_gra, {forming_bond: 2})
                chosen_ts_gra = automol.graph.add_bonds(
                    chosen_ts_gra, [breaking_bond])
                vals = automol.graph.atom_unpaired_electrons(chosen_ts_gra, bond_order=True)
                oversaturated_atoms = [atm for atm, val in vals.items() if val < 0]
                if len(oversaturated_atoms) == 1:
                    chosen_oversaturated_atom = oversaturated_atoms[0]
                    break
    if chosen_oversaturated_atom is None:
        print('could not figure out which atom is being transfered')
        sys.exit()
    return chosen_ts_gra, breaking_bond, forming_bond


def choose_beta_cutoff_distance(geo):
    """ Havent tested for years, originally used to automatically
        identify a heavy atom-heavy atom beta scission reaction 
        based on bond lengths of the geo
    """
    rqqs = [x * 0.1 for x in range(28, 48, 2)]
    double_bnd_atms = []
    for rqq in rqqs:
        if len(double_bnd_atms) == 2:
            break
        ts_gras = automol.geom.connectivity_graph_deprecated(geo, rqq_bond_max=rqq, rqh_bond_max=2.7, rhh_bond_max=2.3)
        ts_gras = automol.graph.set_stereo_from_geometry(ts_gras, geo)
        ts_gras = automol.graph.connected_components(ts_gras)
        if len(ts_gras) != 2:
            continue
        unconn_bnds = []
        for ts_gra_i in ts_gras:
            vals = automol.graph.atom_unpaired_electrons(ts_gra_i)
            double_bnd_atms = [atm for atm in vals if vals[atm] == 1]
            unconn_bnds.extend(list(automol.graph.bond_keys(ts_gra_i)))
    chosen_ts_gra = []
    chosen_oversaturated_atom = None
    rqqs = [x * 0.1 for x in range(32, 48, 2)]
    for rqq in rqqs:
        ts_gras = automol.geom.connectivity_graph_deprecated(geo, rqq_bond_max=rqq, rqh_bond_max=2.7, rhh_bond_max=2.3)
        ts_gras = automol.graph.set_stereo_from_geometry(ts_gras, geo)
        ts_gras = automol.graph.connected_components(ts_gras)
        if len(ts_gras) != 1:
            continue
        for ts_gra_i in ts_gras:
            rad_atm = list(automol.graph.radical_atom_keys(ts_gra_i))[0]
            conn_bnds = automol.graph.bond_keys(ts_gra_i)
            breaking_bond_lst = set(conn_bnds) - set(unconn_bnds)
            if rad_atm in double_bnd_atms and len(breaking_bond_lst) < 2:
                ts_gra_i = automol.graph.set_bond_orders(
                    ts_gra_i, {frozenset(double_bnd_atms): 2})
                breaking_bond = list(breaking_bond_lst)[0]
                forming_bond = frozenset(double_bnd_atms)
            vals = automol.graph.atom_unpaired_electrons(ts_gra_i, bond_order=True)
            oversaturated_atoms = [atm for atm, val in vals.items() if val < 0]
            if len(oversaturated_atoms) == 1:
                chosen_ts_gra = ts_gras[0]
                chosen_oversaturated_atom = oversaturated_atoms[0]
                break
    if chosen_oversaturated_atom is None:
        chosen_ts_gra, breaking_bond, forming_bond = choose_beta_qh_cutoff_distance(geo)
        if breaking_bond is None:
            print('could not figure out which atom is being transfered')
            sys.exit()
    return chosen_ts_gra, breaking_bond, forming_bond


def choose_heavy_cutoff_distance(geo):
    rqqs = [x * 0.1 for x in range(32, 48, 2)]
    chosen_ts_gra = []
    chosen_oversaturated_atom = None
    for rqq in rqqs:
        ts_gras = automol.geom.connectivity_graph_deprecated(geo, rqq_bond_max=rqq, rqh_bond_max=2.7, rhh_bond_max=2.3)
        ts_gras = automol.graph.set_stereo_from_geometry(ts_gras, geo)
        ts_gras = automol.graph.connected_components(ts_gras)
        if len(ts_gras) != 1:
            continue
        for ts_gra_i in ts_gras:
            vals = automol.graph.atom_unpaired_electrons(ts_gra_i, bond_order=True)
            oversaturated_atoms = [atm for atm, val in vals.items() if val < 0]
            if len(oversaturated_atoms) == 1:
                chosen_ts_gra = ts_gras[0]
                chosen_oversaturated_atom = oversaturated_atoms[0]
                break
    if chosen_oversaturated_atom is None:
        print('could not figure out which atom is being transfered')
        sys.exit()
    return chosen_ts_gra, chosen_oversaturated_atom


def choose_cutoff_distance(geo):
    rqhs = [x * 0.1 for x in range(26, 38, 2)]
    chosen_ts_gra = []
    chosen_oversaturated_atom = None
    for rqh in rqhs:
        ts_gras = automol.geom.connectivity_graph_deprecated(geo, rqq_bond_max=3.5, rqh_bond_max=rqh, rhh_bond_max=2.3)
        ts_gras = automol.graph.set_stereo_from_geometry(ts_gras, geo)
        ts_gras = automol.graph.connected_components(ts_gras)
        if len(ts_gras) != 1:
            continue
        for ts_gra_i in ts_gras:
            vals = automol.graph.atom_unpaired_electrons(ts_gra_i, bond_order=True)
            oversaturated_atoms = [atm for atm, val in vals.items() if val < 0]
            if len(oversaturated_atoms) == 1:
                chosen_ts_gra = ts_gras[0]
                chosen_oversaturated_atom = oversaturated_atoms[0]
                break
    if chosen_oversaturated_atom is None:
        print('could not figure out which H is being transfered')
        sys.exit()
    return chosen_ts_gra, chosen_oversaturated_atom


def betasci_gra(geo):
    ts_gra, breaking_bnd, forming_bnd = choose_beta_cutoff_distance(geo)
    return ts_gra, breaking_bnd, forming_bnd


def ringformsci_gra(geo):
    ts_gra, oversat_atm = choose_heavy_cutoff_distance(geo)
    atoms_bnd = automol.graph.atoms_bond_keys(ts_gra)
    bnds = list(atoms_bnd[oversat_atm])
    bnded_atms = [list(bnd - set({oversat_atm}))[0] for bnd in bnds]

    brk_atm_idx = None
    for i in range(len(bnds)):
        gra_i = automol.graph.remove_bonds(ts_gra, (bnds[i],))
        for j in range(len(bnds)):
            if j != i:
                pth_btwn = automol.graph.shortest_path_between_atoms(
                    gra_i, bnded_atms[i], bnded_atms[j])
                if pth_btwn is None:
                    brk_atm_idx = i
    frm_atm_idx = None
    if brk_atm_idx is not None:
        len_pth_to_brk = 0
        for i in range(len(bnds)):
            if i != brk_atm_idx:
                gra_i = automol.graph.remove_bonds(ts_gra, (bnds[i],))
                pth_btwn = automol.graph.shortest_path_between_atoms(
                    gra_i, bnded_atms[i], bnded_atms[brk_atm_idx])
                if len(pth_btwn) > len_pth_to_brk:
                    print('bigger path', pth_btwn)
                    frm_atm_idx = i
                    len_pth_to_brk = len(pth_btwn)
    breaking_bond = None
    forming_bond = None
    if frm_atm_idx is not None and brk_atm_idx is not None:
        breaking_bond = bnds[brk_atm_idx]
        forming_bond = bnds[frm_atm_idx]
    return ts_gra, breaking_bond, forming_bond


def elim_gra(geo):
    ts_gra_h, oversaturated_h = choose_cutoff_distance(geo)
    ts_gra_h, _ = choose_heavy_cutoff_distance(geo)
    atoms_bnd = automol.graph.atoms_bond_keys(ts_gra_h)
    atms = automol.graph.atoms(ts_gra_h)
    bonds = atoms_bnd[oversaturated_h]
    forming_bond = None
    breaking_bond = None
    breaking_bond2 = None
    for bond in bonds:
        atma, atmb = bond
        if atma == oversaturated_h:
            atmb, atma = atma, atmb
        if atms[atma][0] == 'O':
            forming_bond = bond
            adj_bnds = atoms_bnd[atma]
            for abnd in adj_bnds:
                if abnd != bond:
                    oo_bnd = abnd
                    atmc, atmd = oo_bnd
                    if atmc == atma:
                        atmc, atmd = atmd, atmc
                    cbnds = atoms_bnd[atmc]
                    for cbnd in cbnds:
                        if cbnd != oo_bnd:
                            breaking_bond2 = cbnd
        else:
            breaking_bond = bond
    return ts_gra_h, breaking_bond, breaking_bond2, forming_bond


def h_transfer_gra(geo):
    ts_gra, oversaturated_atom = choose_cutoff_distance(geo)
    ts_gra = automol.graph.set_stereo_from_geometry(ts_gra, geo)
    atoms_bnd = automol.graph.atoms_bond_keys(ts_gra)
    bonds = atoms_bnd[oversaturated_atom]
    if len(bonds) != 2:
        print('too many bonds to transfered atom for me to figure out')
        print('I promise i will be smarter in the future')
        sys.exit()
    breaking_bond, forming_bond = bonds
    # when we move on to other reaction types we have to check for double
    # bonds when doing bond orders
    return ts_gra, breaking_bond, forming_bond


def _check_ichs_match(ts_ichs, rxn_ichs, rxn_gras): #, ts_gras, rxn_gras):
    reactant_match = False
    product_match = False
    if ts_ichs[0] == rxn_ichs[0]:
        reactant_match = True
    elif ts_ichs[0][::-1] == rxn_ichs[0]:
        ts_ichs[0] = ts_ichs[0][::-1]
        reactant_match = True
    else:
        ts_ichs = ts_ichs[::-1]
        # ts_gras = ts_gras[::-1]
        rxn_gras = rxn_gras[::-1]
        if ts_ichs[0] == rxn_ichs[0]:
            reactant_match = True
        elif ts_ichs[0][::-1] == rxn_ichs[0]:
            ts_ichs[0] = ts_ichs[0][::-1]
            reactant_match = True
    if reactant_match:
        if ts_ichs[1] == rxn_ichs[1]:
            product_match = True
        elif ts_ichs[1][::-1] == rxn_ichs[-1]:
            ts_ichs[1] = ts_ichs[1][::-1]
            product_match = True
    return reactant_match and product_match, ts_ichs, rxn_gras


def all_reaction_graphs(
        ts_gra, breaking_bond, forming_bond,
        no_form, breaking_bond2):
    """ build ts graph using reactant graph and known
        forming and breaking bonds (there are better automol functions
        for this now)
    """
    # ts_forw_gra = automol.graph.ts.graph(rct_gra, forming_bonds, breaking_bonds)
    if not no_form:
        if breaking_bond2 is None:
            forw_bnd_ord_dct = {breaking_bond: 0.9, forming_bond: 0.1}
            back_bnd_ord_dct = {breaking_bond: 0.1, forming_bond: 0.9}
        else:
            forw_bnd_ord_dct = {breaking_bond: 0.9, breaking_bond2: 0.9, forming_bond: 0.1}
            back_bnd_ord_dct = {breaking_bond: 0.1, breaking_bond2: 0.1, forming_bond: 0.9}
    else:
        forw_bnd_ord_dct = {breaking_bond: 0.9, forming_bond: 1}
        back_bnd_ord_dct = {breaking_bond: 0.1, forming_bond: 1}
    forward_gra = automol.graph.set_bond_orders(ts_gra, forw_bnd_ord_dct)
    backward_gra = automol.graph.set_bond_orders(ts_gra, back_bnd_ord_dct)
    reactant_gras = automol.graph.ts.reactants_graph(forward_gra)
    reactant_gras = automol.graph.connected_components(reactant_gras)
    product_gras = automol.graph.ts.reactants_graph(backward_gra)
    product_gras = automol.graph.connected_components(product_gras)
    ts_gras = [forward_gra, backward_gra]
    rxn_gras = [reactant_gras, product_gras]
    return ts_gras, rxn_gras


def build_zrxn_from_geo(
        geo, rxn_info, rxn_class, rct_gra, breaking_bonds, forming_bonds):
    """ build the zrxn object using known information about the reactant
        graph and the forming/breaking bonds and update with found ts geo
    """
    ts_forw_gra = automol.graph.ts.graph(rct_gra, forming_bonds, breaking_bonds)
    ts_back_gra = automol.graph.ts.reverse(ts_forw_gra)
    ts_zma, ts_zc_ = automol.geom.zmatrix_with_conversion_info(
        geo, gra=ts_forw_gra)
    ts_forw_gra = automol.graph.apply_zmatrix_conversion(ts_forw_gra, ts_zc_)
    ts_back_gra = automol.graph.apply_zmatrix_conversion(ts_back_gra, ts_zc_)
    rct_gras = automol.graph.connected_components(
        automol.graph.ts.reactants_graph(ts_forw_gra))
    prd_gras = automol.graph.connected_components(
        automol.graph.ts.products_graph(ts_forw_gra))
    rct_zc_ = []
    prd_zc_ = []
    rct_strucs = [automol.graph.geometry(gra) for gra in rct_gras]
    prd_strucs = [automol.graph.geometry(gra) for gra in prd_gras]
    rct_keys = []
    for i, gra in enumerate(rct_gras):
        keys = automol.graph.atom_keys(gra)
        key_map = {}
        rev_map = {}
        for j, key in enumerate(keys):
            key_map[key] = j
            rev_map[j] = key
        gra = automol.graph.relabel(gra, key_map)
        zma, zc_ = automol.geom.zmatrix_with_conversion_info(
            rct_strucs[i], gra=gra)
        rct_strucs[i] = zma
        rct_zc_.append(zc_)
        rct_keys.append([rev_map[zc_[key][0]] for key in sorted(zc_.keys())])
    prd_keys = []
    for i, gra in enumerate(prd_gras):
        keys = automol.graph.atom_keys(gra)
        key_map = {}
        rev_map = {}
        for j, key in enumerate(keys):
            key_map[key] = j
            rev_map[j] = key
        gra = automol.graph.relabel(gra, key_map)
        zma, zc_ = automol.geom.zmatrix_with_conversion_info(
            prd_strucs[i], gra=gra)
        prd_strucs[i] = zma
        prd_zc_.append(zc_)
        prd_keys.append([rev_map[zc_[key][0]] for key in sorted(zc_.keys())])
    # match_ich_info = _match_info(rxn_info, (ts_forw_gra, ts_back_gra))
    std_zrxn = automol.reac.from_data(
        ts_forw_gra, rct_keys, prd_keys,
        rxn_class,
        ts_struc=ts_zma, rct_strucs=rct_strucs, prd_strucs=prd_strucs, struc_typ='zmat',
        ts_zc=ts_zc_, rct_zcs=rct_zc_, prd_zcs=prd_zc_)
    ts_geo = automol.zmat.geometry(ts_zma)
    return std_zrxn, ts_zma, ts_geo, rxn_info


def _match_info(rxn_info, rxn_gras):
    rxn_ichs = [[], []]
    for i, side in enumerate(rxn_info[0]):
        for ich in side:
            rxn_ichs[i].append(ich)

    ts_ichs = [[], []]
    for rgra in rxn_gras[0]:
        try:
            rich = automol.graph.chi(rgra, stereo=True)
        except IndexError:
            rich = automol.graph.chi(rgra)
        ts_ichs[0].append(rich)
    for pgra in rxn_gras[1]:
        try:
            pich = automol.graph.chi(pgra, stereo=True)
        except IndexError:
            pich = automol.graph.chi(pgra)
        ts_ichs[1].append(pich)

    match_ich_info, ts_ichs, rxn_gras = _check_ichs_match(
        ts_ichs, rxn_ichs, rxn_gras)
    # match_ich_info, ts_ichs, ts_gras, rxn_gras = _check_ichs_match(
    #     ts_ichs, rxn_ichs, ts_gras, rxn_gras)
    if not match_ich_info:
        print('make sure these look right')
        print('my ichs  ', ts_ichs)
        print('your ichs', rxn_ichs)
        sys.exit()
    return match_ich_info


def get_zrxn(geo, rxn_info, rxn_class):
    """ Automatically determine the ts graph, and the forming and breaking
        bonds from the ts geometry and the user specified reaction class
    """
    breaking_bond2 = None
    if rxn_class in ['hydrogen abstraction', 'hydrogen migration']:
        ts_gra, breaking_bond, forming_bond = h_transfer_gra(geo)
    elif rxn_class in ['ring forming scission']:
        ts_gra, breaking_bond, forming_bond = ringformsci_gra(geo)
    elif rxn_class in ['beta scission']:
        ts_gra, breaking_bond, forming_bond = betasci_gra(geo)
    elif rxn_class in ['elimination']:
        ts_gra, breaking_bond, breaking_bond2, forming_bond = elim_gra(geo)

    ts_gras, rxn_gras = all_reaction_graphs(
        ts_gra, breaking_bond, forming_bond,
        rxn_class == 'beta scission', breaking_bond2)

    match_ich_info = _match_info(rxn_info, rxn_gras)
    if match_ich_info:
        reactant_keys = []
        for gra in rxn_gras[0]:
            reactant_keys.append(automol.graph.atom_keys(gra))
        product_keys = []
        for gra in rxn_gras[1]:
            product_keys.append(automol.graph.atom_keys(gra))
        std_rxn = automol.reac.from_forward_reverse(
            rxn_class, *ts_gras, reactant_keys, product_keys)
        std_zrxn = automol.reac.with_structures(std_rxn, "zmat")
        ts_zma = automol.reac.ts_structure(std_zrxn)
        ts_geo = automol.zmat.geometry(ts_zma)
    else:
        print(
            'The reactants and products found for the transition state' +
            ' did not match those specified in user input')
        sys.exit()
    return std_zrxn, ts_zma, ts_geo, rxn_info


def main(insert_dct):

    prefix = read_user_filesystem(insert_dct)
    # parse method from insert input file
    thy_info = parse_user_theory(insert_dct)
    prog, method, basis, _ = thy_info

    # Read in the input and output files
    # these can be filename (input_file) or strings (input_string)
    if insert_dct['input_file'] in ['string', None]:
        inp_str = insert_dct['input_string']
    else:
        inp_str = read_user_file(insert_dct, 'input_file')
    if insert_dct['output_file'] in ['string', None]:
        out_str = insert_dct['output_string']
    else:
        out_str = read_user_file(insert_dct, 'output_file')
    output_type = insert_dct['output_type']
    geo, zma, ene, hess_job = _struct_based_on_input(
        output_type, out_str, prog, method)

    # Parse out user specified save location
    if insert_dct['saddle']:
        rxn_info, spc_info, rxn_class = parse_user_reaction(insert_dct)
        # user has to have specified the reaction class or the breaking/forming bonds
        # if the later, it will be inserted as an unclassified reaction
        graph_file = insert_dct.get("graph_file")
        if graph_file:
            rxn = automol.reac.from_string(Path(graph_file).read_text())

            # Get geometry and z-matrix views
            grxn = automol.reac.with_structures(rxn, "geom")
            zrxn = automol.reac.with_structures(rxn, "zmat")

            # Add given structure and update other view
            if geo is not None:
                grxn = automol.reac.update_structures(grxn, ts_struc=geo)
                zrxn = automol.reac.with_structures(grxn, "zmat")
            elif zma is not None:
                zrxn = automol.reac.update_structures(zrxn, ts_struc=zma)
                grxn = automol.reac.with_structures(grxn, "geom")
            else:
                raise ValueError(
                    f"Must have geometry or z-matrix:\ngeo = {geo}\nzma = {zma}"
                )

            # Extract both structures, which should now be consistent
            geo = automol.reac.ts_structure(grxn)
            zma = automol.reac.ts_structure(zrxn)

            # Check that geometry is consistent
            ggra = automol.reac.ts_graph(grxn)
            assert automol.graph.geometry_matches(
                ggra, geo, stereo=False, check_ts_bonds=False
            ), f"ggra = {ggra}\ngeo = {geo}"

            # Check that z-matrix is consistent
            zgra = automol.reac.ts_graph(zrxn)
            assert automol.graph.zmatrix_matches(
                zgra, zma, stereo=False, check_ts_bonds=False
            ), f"zgra = {zgra}\nzma = {zma}"

        elif insert_dct['breaking_bonds'] is not None or insert_dct['forming_bonds'] is not None:
            zrxn, zma, geo, rxn_info = build_zrxn_from_geo(
                geo, rxn_info, rxn_class, insert_dct['rct_gra'],
                insert_dct['breaking_bonds'], insert_dct['forming_bonds'])
        else:
            zrxn, zma, geo, rxn_info = get_zrxn(geo, rxn_info, rxn_class)
    else:
        zrxn = None
        spc_info = parse_user_species(insert_dct)
        if None in spc_info:
            spc_info = species_info_from_input(geo, spc_info)
        if zma is None:
            zma = automol.geom.zmatrix(geo)
            geo = automol.zmat.geometry(zma)

    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)

    # Check that the save location matches geo information
    if not insert_dct['saddle']:
        if not species_match(geo, spc_info) and not insert_dct['trusted']:
            print(
                'I refuse to save this geometry until user specified' +
                ' info matches the info in user given output')
            sys.exit()
        # Check that the rid/cid info matches the filesystem
        fs_array, _ = create_species_filesystems(
            prefix, spc_info, mod_thy_info, locs=None)
    else:
        fs_array, _ = create_reaction_filesystems(
            prefix, rxn_info, mod_thy_info,
            ts_locs=insert_dct.get('ts_locs', None), locs=None)
    cnf_fs = fs_array[-1]
    locs = parse_user_locs(insert_dct, geo, cnf_fs)
    if not locs_match(geo, cnf_fs, locs):
        print(
            'I refuse to save this geometry until user specified' +
            ' info matches the info in user given output')
        sys.exit()

    # update info object and save conformer
    inf_obj = autofile.schema.info_objects.run(
        job=elstruct.Job.OPTIMIZATION, prog=prog, version='',
        method=method, basis=basis, status=autofile.schema.RunStatus.SUCCESS)
    inf_obj.utc_end_time = autofile.schema.utc_time()
    inf_obj.utc_start_time = autofile.schema.utc_time()
    _, saved_geos, saved_enes = _saved_cnf_info(
        cnf_fs, mod_thy_info)
    if _geo_unique(geo, ene, saved_geos, saved_enes, zrxn=zrxn):
        sym_id = _sym_unique(
            geo, ene, saved_geos, saved_enes)
        if sym_id is None:
            if cnf_fs[0].file.info.exists():
                rinf_obj = cnf_fs[0].file.info.read()
            else:
                rinf_obj = autofile.schema.info_objects.conformer_trunk(0)
                rinf_obj.nsamp = 1
            if cnf_fs[1].file.info.exists([locs[0]]):
                cinf_obj = cnf_fs[1].file.info.read([locs[0]])
                cnsampd = cinf_obj.nsamp
                cnsampd += 1
                cinf_obj.nsamp = cnsampd
            else:
                cinf_obj = autofile.schema.info_objects.conformer_branch(0)
                cinf_obj.nsamp = 1
            cnf_fs[1].create([locs[0]])
            cnf_fs[0].file.info.write(rinf_obj)
            cnf_fs[1].file.info.write(cinf_obj, [locs[0]])
            hess_ret = None
            if hess_job:
                hess_inf_obj = autofile.schema.info_objects.run(
                    job=elstruct.Job.HESSIAN, prog=prog, version='',
                    method=method, basis=basis,
                    status=autofile.schema.RunStatus.SUCCESS)
                hess_ret = (hess_inf_obj, inp_str, out_str)
            save_info = (geo, zma, ene, inf_obj, inp_str)
            save.parsed_conformer(
                save_info, cnf_fs, mod_thy_info[1:], rng_locs=[locs[0]],
                tors_locs=[locs[1]], zrxn=zrxn, hess_ret=hess_ret)
            print(
                'geometry is now saved at {}'.format(cnf_fs[-1].path(locs)))
    else:
        print('{}'.format(cnf_fs[-1].path(locs)))
        print(
            'the geometry in the output is not unique to filesystem' +
            '... not saving')


def species_match(geo, spc_info):
    match = True
    ich, _, mul = spc_info
    mults_allowed = automol.graph.possible_spin_multiplicities(
         automol.chi.graph(ich, stereo=False))
    geo_ich = automol.geom.chi(geo, stereo=True)
    if ich != geo_ich:
        print(
            'user specified inchi {}'.format(ich) +
            'does not match inchi from output {}'.format(geo_ich) +
            'which is based on geometry from output:\n' +
            '{}'.format(automol.geom.string(geo)))
        match = False
    if mul not in mults_allowed:
        print(
            'user specified mult of {}'.format(mul) +
            'is not an allowed multiplicty for inchi {}'.format(ich))
        match = False
    return match


def locs_match(geo, cnf_fs, locs):
    match = True
    rid = locs[0]
    geo_rid = rng_loc_for_geo(geo, cnf_fs)
    if geo_rid is not None:
        if geo_rid != rid:
            print(
                'Error: rid mismatch for the filesystem at' +
                ' {}'.format(cnf_fs[0].path()) +
                '\nthe expected rid for this geo is {}'.format(geo_rid) +
                '\nthe user rid in input file is {}'.format(rid))
            match = False
    return match


def rng_loc_for_geo(geo, cnf_fs):
    """ get the rid from an existing FS 
        that is consistent with the new geometry
    """
    rid = None
    frag_geo = _fragment_ring_geo(geo)
    if frag_geo is not None:
        frag_zma = automol.geom.zmatrix(frag_geo)
    checked_rids = []
    for locs in cnf_fs[-1].existing():
        current_rid, _ = locs
        if current_rid in checked_rids:
            continue
        if cnf_fs[-1].file.geometry.exists(locs):
            checked_rids.append(current_rid)
            locs_geo = cnf_fs[-1].file.geometry.read(locs)
            frag_locs_geo = _fragment_ring_geo(locs_geo)
            if frag_locs_geo is None:
                rid = locs[0]
                break
            frag_locs_zma = automol.geom.zmatrix(frag_locs_geo)
            if automol.zmat.almost_equal(
                    frag_locs_zma, frag_zma, dist_rtol=0.1, ang_atol=.4):
                rid = locs[0]
                break
    return rid


def parse_script_options(script_input_file):
    script_input = ioformat.pathtools.read_file(
        '', script_input_file, print_debug=True).splitlines()
    insert_dct = {
        'save_filesystem': None,
        'smiles': None,
        'inchi': None,
        'mult': None,
        'charge': None,
        'rid': None,
        'cid': None,
        'theory': None,
        'program': None,
        'method': None,
        'basis': None,
        'orb_res': None,
        'input_file': None,
        'output_file': None,
        'output_type': 'optimization',
        'graph_file': None,
        'ts_locs': None,
        'ts_mult': None,
        'rxn_class': None,
        'breaking_bonds': None,
        'forming_bonds': None,
        'zrxn_file': None,
        'run_path': None,
        'saddle': False,
        'trusted': False,
    }
    for i, line in enumerate(script_input):
        if len(line) < 2:
            continue
        elif '!' in line[0]:
            continue
        line = line.split('!')[0]
        if ':' not in line:
            print(
                'ERROR: line\n({}) {}\n is not parsable, '.format(i, line) +
                'script will exit until input is resolved to avoid' +
                ' filesystem contamination.' +
                'Comment lines should contain "!"' +
                'Key format should be:\n' +
                '<Keyword>: <Value>\n' +
                'Allowed keywords are:\n' +
                '{}'.format('\n'.join(list(insert_dct.keys())))
            )
            sys.exit()
        keyword, value = line.split(':')
        if keyword in insert_dct:
            if 'None' in value:
                value = None
            elif keyword in ['mult', 'charge', 'ts_mult']:
                values = []
                for val in value.split(','):
                    values.append(int(val))
                if len(values) == 1:
                    value = values[0]
                else:
                    value = values
            elif keyword in ['ts_locs']:
                value = (int(value),)
            elif keyword in ['rxn_class']:
                # strip whitespaces form either side of reaction
                # class but not in between words
                value = value.split()
                for i, val in enumerate(value):
                    value[i] = val.replace(' ', '')
                value = ' '.join(value)
            elif keyword not in ['smiles', 'inchi']:
                value = value.replace(' ', '')
            else:
                value = value.split(' = ')
                if len(value) > 1:
                    insert_dct['saddle'] = True
                    reactants, products = value
                    reactants = reactants.split(' + ')
                    products = products.split(' + ')
                    values = [[], []]
                    for reactant in reactants:
                        values[0].append(reactant.replace(' ', ''))
                    for product in products:
                        values[1].append(product.replace(' ', ''))
                    value = values
                else:
                    value = value[0].replace(' ', '')
            print(keyword, value)
            insert_dct[keyword] = value
        else:
            print(
                'ERROR: Keyword {} is not recognized'.format(keyword) +
                'script will exit until inpupt is resolved to avoid' +
                ' filesystem contamination.' +
                'Allowed keywords are:\n' +
                '{}'.format('\n'.join(list(insert_dct.keys())))
            )
            sys.exit()

    # reshape arrays
    if insert_dct['saddle']: 
        if insert_dct['mult'] is not None:
            insert_dct['mult'] = [
                [insert_dct['mult'][i] for i in range(len(reactants))],
                [insert_dct['mult'][i] for i in range(
                    len(reactants), len(reactants) + len(products))]]
        if insert_dct['charge'] is not None:
            insert_dct['charge'] = [
                [insert_dct['charge'][i] for i in range(len(reactants))],
                [insert_dct['charge'][i] for i in range(
                    len(reactants), len(reactants) + len(products))]]

    return insert_dct


def _struct_based_on_input(
        output_type, out_str, prog, method):
    """ read in geo, ene, zma, based on what user specified in
        out_str
    """
    hess_job = False
    zrxn = None
    zma = None
    if output_type in ['geo', 'geom', 'geometry', 'xyz']:
        geo = automol.geom.from_xyz_string(out_str)
        ene = float(automol.geom.xyz_string_comment(out_str))
    elif output_type in ['zma', 'zmat', 'zmatrix']:
        out_lines = out_str.splitlines()
        ene = float(out_lines[0])
        out_str = '\n'.join(out_lines[1:])
        zma = automol.zmat.from_string(out_str)
        geo = automol.zmat.geometry(zma)
    elif output_type in ['opt', 'optimization']:
        ene = elstruct.reader.energy(prog, method, out_str)
        geo = elstruct.reader.opt_geometry(prog, out_str)
        zma = elstruct.reader.opt_zmatrix(prog, out_str)
        if geo is None:
            print(
                'No geometry could be parsed from output' +
                'Check that the program matches user specied' +
                ' {}'.format(prog) + ' and method matches' +
                ' {}'.format(method))
            sys.exit()
    elif output_type in ['freqs', 'frequencies']:
        ene = elstruct.reader.energy(prog, method, out_str)
        geo = elstruct.reader.opt_geometry(prog, out_str)
        zma = elstruct.reader.opt_zmatrix(prog, out_str)
        hess_job = True
        if geo is None:
            print(
                'No geometry could be parsed from output' +
                'Check that the program matches user specied' +
                ' {}'.format(prog) + ' and method matches' +
                ' {}'.format(method))
            sys.exit()
    return geo, zma, ene, hess_job

