""" Add a species to your database
usiing a log file
"""

import sys

import autofile
import automol
from mechanalyzer.inf import thy as tinfo
from mechanalyzer.inf import rxn as rinfo
from mechanalyzer.inf import spc as sinfo
import elstruct
from automol.geom import ring_fragments_geometry as _fragment_ring_geo
from mechlib.filesys import save
from mechroutines.es._routines.conformer import _saved_cnf_info
from mechroutines.es._routines.conformer import _sym_unique
from mechroutines.es._routines.conformer import _geo_unique

THEORY_DCT = {
    'lvl_wbs': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'wb97xd',
        'basis':  '6-31g*'
    },
    'lvl_wbm': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'wb97xd',
        'basis':  '6-31+g*'
    },
    'lvl_wbt': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'wb97xd',
        'basis':  'cc-pvtz'},
    'lvl_m06s': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'm062x',
        'basis':  '6-31g*'
    },
    'lvl_m06m': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'm062x',
        'basis':  '6-31+g*'
    },
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
        'basis':  'cc-pvqz'
    },
    'lvl_b3s': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'b3lyp',
        'basis':  '6-31g*'
        },
    'lvl_b3mg': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'b3lyp',
        'basis':  '6-311g**'
    },
    'lvl_b3t': {
        'orb_res': 'RU',
        'program': 'gaussian09',
        'method':  'b3lyp',
        'basis':  'cc-pvtz'},
    'cc_lvl_d': {
        'orb_res': 'RR',
        'program': 'molpro2015',
        'method':  'ccsd(t)',    'basis':  'cc-pvdz'},
    'cc_lvl_t': {
        'orb_res': 'RR',
        'program': 'molpro2015',
        'method':  'ccsd(t)',    'basis':  'cc-pvtz'},
    'cc_lvl_q': {
        'orb_res': 'RR',
        'program': 'molpro2015',
        'method':  'ccsd(t)',    'basis':  'cc-pvqz'
    },
    'cc_lvl_df': {
        'orb_res': 'RR',
        'program': 'molpro2015',
        'method':  'ccsd(t)-f12',
        'basis':  'cc-pvdz-f12'
    },
    'cc_lvl_tf': {
        'orb_res': 'RR',
        'program': 'molpro2015',
        'method':  'ccsd(t)-f12',
        'basis':  'cc-pvtz-f12'
    },
    'cc_lvl_qf': {
        'orb_res': 'RR',
        'program': 'molpro2015',
        'method':  'ccsd(t)-f12',
        'basis':  'cc-pvqz-f12'
    },
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


def parse_user_locs(insert_dct):
    rid = insert_dct['rid']
    cid = insert_dct['cid']
    if rid is None:
        rid = autofile.schema.generate_new_ring_id()
    if cid is None:
        cid = autofile.schema.generate_new_conformer_id()
    return [rid], [cid], (rid, cid)


def parse_user_species(insert_dct):
    smi = insert_dct['smiles']
    ich = insert_dct['inchi']
    mult = insert_dct['mult']
    chg = insert_dct['charge']
    if ich is None and smi is None:
        print(
            'Error: user did not specify species' +
            'with an inchi or smiles in input')
        sys.exit()
    if ich is None:
        ich = automol.smiles.inchi(smi)
    if not automol.inchi.is_complete(ich):
        ich = automol.inchi.add_stereo(ich)
    if mult is None:
        print('Error: user did not specify mult in input')
        sys.exit()
    if chg is None:
        print('Error: user did not specify charge in input')
        sys.exit()
    return sinfo.from_data(ich, chg, mult)


def parse_user_reaction(insert_dct):
    smis = insert_dct['smiles']
    ichs = insert_dct['inchi']
    mults = insert_dct['mult']
    chgs = insert_dct['charge']
    rxn_class = insert_dct['rxn_class']
    if ichs is None:
        ichs = [[], []]
        for smi in smis[0]:
            ichs[0].append(automol.smiles.inchi(smi))
        for smi in smis[1]:
            ichs[1].append(automol.smiles.inchi(smi))
    for idx, ich in enumerate(ichs[0]):
        if not automol.inchi.is_complete(ich):
            ich = automol.inchi.add_stereo(ich)
            ichs[0][idx] = ich
    for idx, ich in enumerate(ichs[1]):
        if not automol.inchi.is_complete(ich):
            ich = automol.inchi.add_stereo(ich)
            ichs[1][idx] = ich
    if mults is None:
        print('Error: user did not specify mults in input')
        sys.exit()
    if chgs is None:
        print('Error: user did not specify charges in input')
        sys.exit()
    flat_ichs = sum(ichs, [])
    if len(flat_ichs) != len(mults):
        print(
            'Error: number of species does not match number of mults')
        sys.exit()
    if len(flat_ichs) != len(chgs):
        print(
            'Error: number of species does not match number of charges')
        sys.exit()
    idx = 0
    rxn_muls = [[], []]
    rxn_chgs = [[], []]
    for ich in ichs[0]:
        mults_allowed = automol.graph.possible_spin_multiplicities(
            automol.inchi.graph(ich, stereo=False))
        if mults[idx] not in mults_allowed:
            print(
                'user specified mult of {}'.format(mults[idx]) +
                'is not an allowed multiplicty for inchi {}'.format(ich))
            sys.exit()
        rxn_muls[0].append(mults[idx])
        rxn_chgs[0].append(chgs[idx])
        idx += 1
    for ich in ichs[1]:
        mults_allowed = automol.graph.possible_spin_multiplicities(
            automol.inchi.graph(ich, stereo=False))
        if mults[idx] not in mults_allowed:
            print(
                'user specified mult of {}'.format(mults[idx]) +
                'is not an allowed multiplicty for inchi {}'.format(ich))
            sys.exit()
        rxn_muls[1].append(mults[idx])
        rxn_chgs[1].append(chgs[idx])
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
    # Get input method explicitly inputted
    program = insert_dct['program']
    method = insert_dct['method']
    basis = insert_dct['basis']
    orb_res = insert_dct['orb_res']
    # Get input method from theory dictionary
    theory = insert_dct['theory']
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

    # species filesystem
    spc_fs = autofile.fs.species(prefix)
    spc_fs[-1].create(spc_info)
    spc_prefix = spc_fs[-1].path(spc_info)

    # theory filesystem
    thy_fs = autofile.fs.theory(spc_prefix)
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

    # species filesystem
    print('rxn_info', rxn_info)
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
    if dct[keyword] is None:
        print(
            'ERROR: No filename is specified for {}'.format(keyword) +
            'Script will exit')
        sys.exit()
    file_name = dct[keyword]
    return autofile.io_.read_file(file_name)


def read_user_filesystem(dct):
    if dct['save_filesystem'] is None:
        print(
            'ERROR: No save_filesystem}' +
            'Script will exit')
        sys.exit()
    return dct['save_filesystem']


def choose_cutoff_distance(geo):
    rqhs = [x * 0.1 for x in range(26, 38, 2)]
    chosen_ts_gra = []
    chosen_oversaturated_atom = None
    for rqh in rqhs:
        ts_gras = automol.geom.connectivity_graph(geo, rqq_bond_max=3.5, rqh_bond_max=rqh, rhh_bond_max=2.3)
        ts_gras = automol.graph.set_stereo_from_geometry(ts_gras, geo)
        ts_gras = automol.graph.connected_components(ts_gras)
        if len(ts_gras) != 1:
            continue
        for ts_gra_i in ts_gras:
            vals = automol.graph.atom_unsaturated_valences(ts_gra_i, bond_order=True)
            oversaturated_atoms = [atm for atm, val in vals.items() if val < 0]
            if len(oversaturated_atoms) == 1:
                chosen_ts_gra = ts_gras[0]
                chosen_oversaturated_atom = oversaturated_atoms[0]
                break
    if chosen_oversaturated_atom is None:
        print('could not figure out which H is being transfered')
        sys.exit()
    return chosen_ts_gra, chosen_oversaturated_atom


def get_zrxn(geo, rxn_info, rxn_class):
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
    forw_bnd_ord_dct = {breaking_bond: 0.9, forming_bond: 0.1}
    back_bnd_ord_dct = {breaking_bond: 0.1, forming_bond: 0.9}
    forward_gra = automol.graph.set_bond_orders(ts_gra, forw_bnd_ord_dct)
    backward_gra = automol.graph.set_bond_orders(ts_gra, back_bnd_ord_dct)
    reactant_gras = automol.graph.without_dummy_bonds(
        automol.graph.without_fractional_bonds(forward_gra))
    reactant_gras = automol.graph.connected_components(reactant_gras)
    product_gras = automol.graph.without_dummy_bonds(
        automol.graph.without_fractional_bonds(backward_gra))
    product_gras = automol.graph.connected_components(product_gras)
    ts_gras = [forward_gra, backward_gra]
    rxn_gras = [reactant_gras, product_gras]
    rxn_ichs = [[], []]
    for i, side in enumerate(rxn_info[0]):
        for ich in side:
            rxn_ichs[i].append(ich)
    ts_ichs = [[], []]
    for rgra in reactant_gras:
        try:
            rich = automol.graph.inchi(rgra, stereo=True)
        except IndexError:
            rich = automol.graph.inchi(rgra)
        rsmi = automol.inchi.smiles(rich)
        ts_ichs[0].append(rich)
    for pgra in product_gras:
        try:
            pich = automol.graph.inchi(pgra, stereo=True)
        except IndexError:
            pich = automol.graph.inchi(pgra)
        psmi = automol.inchi.smiles(pich)
        ts_ichs[1].append(pich)
    reactant_match = False
    product_match = False
    if ts_ichs[0] == rxn_ichs[0]:
        reactant_match = True
    elif ts_ichs[0][::-1] == rxn_ichs[0]:
        ts_ichs[0] = ts_ichs[0][::-1]
        reactant_match = True
    else:
        ts_ichs = ts_ichs[::-1]
        ts_gras = ts_gras[::-1]
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
    # force matches to be True because stereo doesn't work
    if not reactant_match or not product_match:
        print(
            'The reactants and products found for the transition state' +
            ' did not match those specified in user input',
            ts_ichs, rxn_ichs)
    reactant_match = True
    product_match = True
    # remove above lines if stereo is available for sp
    if reactant_match and product_match:
        reactant_keys = []
        for gra in rxn_gras[0]:
            reactant_keys.append(automol.graph.atom_keys(gra))
        product_keys = []
        for gra in rxn_gras[1]:
            product_keys.append(automol.graph.atom_keys(gra))
        std_rxn = automol.reac.Reaction(
            rxn_class, *ts_gras, reactant_keys, product_keys)
        ts_zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(
            std_rxn, geo)
        std_zrxn = automol.reac.relabel_for_zmatrix(
            std_rxn, zma_keys, dummy_key_dct)
        # rxn_info = (ts_ichs, *rxn_info[1:])
        ts_geo = automol.zmat.geometry(ts_zma)
    else:
        print(
            'The reactants and products found for the transition state' +
            ' did not match those specified in user input',
            ts_ichs, rxn_ichs)
        sys.exit()
    return std_zrxn, ts_zma, ts_geo, rxn_info


def main(insert_dct):

    prefix = read_user_filesystem(insert_dct)
    # parse method from insert input file
    thy_info = parse_user_theory(insert_dct)
    prog, method, basis, _ = thy_info

    hess_job = False
    zrxn = None
    zma = None
    # Read in the input and output files that we
    # Are inserting into the filesystem
    inp_str = read_user_file(insert_dct, 'input_file')
    out_str = read_user_file(insert_dct, 'output_file')
    output_type = insert_dct['output_type']
    if output_type == 'geo':
        geo = automol.geom.from_xyz_string(out_str)
        ene = float(automol.geom.xyz_string_comment(out_str))
    elif output_type == 'zma':
        out_lines = out_str.splitlines()
        ene = float(out_lines[0])
        out_str = '\n'.join(out_lines[1:])
        zma = automol.zmat.from_string(out_str)
        geo = automol.zmat.geometry(zma)
    elif output_type == 'optimization':
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
    elif output_type == 'frequencies':
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
    if zma is None:
        zma = automol.geom.zmatrix(geo)
        geo = automol.zmat.geometry(zma)

    # Parse out user specified save location
    if insert_dct['saddle']:
        rxn_info, spc_info, rxn_class = parse_user_reaction(insert_dct)
        zrxn, zma, geo, rxn_info = get_zrxn(geo, rxn_info, rxn_class)
    else:
        spc_info = parse_user_species(insert_dct)
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
    rng_locs, tors_locs, locs = parse_user_locs(insert_dct)

    # Check that the save location matches geo information
    if not insert_dct['saddle']:
        if not species_match(geo, spc_info):
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
            ts_locs=insert_dct['ts_locs'], locs=None)
    cnf_fs = fs_array[-1]
    if not locs_match(geo, cnf_fs, locs):
        print(
            'I refuse to save this geometry until user specified' +
            ' info matches the info in user given output')
        sys.exit()

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
                save_info, cnf_fs, mod_thy_info[1:], rng_locs=rng_locs,
                tors_locs=tors_locs, zrxn=zrxn, hess_ret=hess_ret)
            print(
                'geometry is now saved at {}'.format(cnf_fs[-1].path(locs)))
    else:
        print(
            'the geometry in the output is not unique to filesystem' +
            '... not saving')


def species_match(geo, spc_info):
    match = True
    ich, _, mul = spc_info
    mults_allowed = automol.graph.possible_spin_multiplicities(
         automol.inchi.graph(ich, stereo=False))
    geo_ich = automol.geom.inchi(geo, stereo=True)
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


def parse_script_input(script_input_file):
    script_input = autofile.io_.read_file(script_input_file).splitlines()
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
        'ts_locs': None,
        'ts_mult': None,
        'rxn_class': None,
        'zrxn_file': None,
        'run_path': None,
        'saddle': False,
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
    return insert_dct


if __name__ == '__main__':
    SCRIPT_INPUT_FILE = sys.argv[1]
    insert_options_dct = parse_script_input(SCRIPT_INPUT_FILE)
    main(insert_options_dct)
