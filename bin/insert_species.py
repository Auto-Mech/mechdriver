""" Add a species to your database
usiing a log file
"""

import sys
import os

import autofile
import automol
from mechanalyzer.inf import thy as tinfo
from mechanalyzer.inf import rxn as rinfo
from mechanalyzer.inf import spc as sinfo
import elstruct
import autorun
from mechroutines.es._routines.conformer import _saved_cnf_info
from mechroutines.es._routines.conformer import _sym_unique
from mechroutines.es._routines.conformer import _save_unique_conformer
from mechroutines.es._routines.conformer import _geo_unique
from mechroutines.es._routines.conformer import _fragment_ring_geo
from mechroutines.es._routines._sadpt import save_saddle_point
from mechlib.reaction.rxnid import _id_reaction

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
    return (rid, cid)


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
    zrxn_file = insert_dct['zrxn_file']
    if ichs is None:
        ichs = [[], []]
        for smi in smis[0]:
            ichs[0].append(automol.smiles.inchi(smi))
        for smi in smis[1]:
            ichs[1].append(automol.smiles.inchi(smi))
    for idx, ich in enumerate(ichs[0]):
        print('inchi', ich)
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
    if zrxn_file is not None:
        zrxn_str = autofile.io_.read_file(zrxn_file)
        zrxns = [automol.reac.from_string(zrxn_str)]
    else:
        zrxns, _ = _id_reaction(rxn_info)
    return rxn_info, ts_info, zrxns


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


def main(insert_dct):

    prefix = read_user_filesystem(insert_dct)
    # Read in the input and output files that we
    # Are inserting into the filesystem
    inp_str = read_user_file(insert_dct, 'input_file')
    out_str = read_user_file(insert_dct, 'output_file')

    # parse method from insert input file
    thy_info = parse_user_theory(insert_dct)

    # parse out geo information first, to make sure
    # user save specifications match output
    prog, method, basis, _ = thy_info
    ene = elstruct.reader.energy(prog, method, out_str)
    geo = elstruct.reader.opt_geometry(prog, out_str)
    if geo is None:
        print(
            'No geometry could be parsed from output' +
            'Check that the program matches user specied' +
            ' {}'.format(prog) + ' and method matches' +
            ' {}'.format(method))
        sys.exit()

    # Parse out user specified save location
    zrxn = None
    if insert_dct['saddle']:
        rxn_info, spc_info, zrxns = parse_user_reaction(insert_dct)
        print('TS\n', automol.graph.string(automol.geom.graph(geo)))
        for zrxn_i in zrxns:
            forw_form_key = automol.reac.forming_bond_keys(zrxn_i)
            back_form_key = automol.reac.forming_bond_keys(zrxn_i, rev=True)
            forward_gra = automol.graph.without_stereo_parities(
                automol.graph.without_dummy_bonds(
                    automol.graph.without_fractional_bonds(
                        zrxn_i.forward_ts_graph)))
            forward_gra = automol.graph.add_bonds(forward_gra, forw_form_key)
            backward_gra = automol.graph.without_stereo_parities(
                automol.graph.without_dummy_bonds(
                    automol.graph.without_fractional_bonds(
                        zrxn_i.backward_ts_graph)))
            backward_gra = automol.graph.add_bonds(backward_gra, back_form_key)
            print('forRXN', automol.graph.string(zrxn_i.forward_ts_graph))
            print('forRXN', automol.graph.string(forward_gra))
            print('bacRXN', automol.graph.string(zrxn_i.backward_ts_graph))
            print('bacRXN', automol.graph.string(backward_gra))
            if forward_gra == automol.geom.graph(geo, stereo=False):
                zrxn = zrxn_i
                zma, _, _ = automol.reac.ts_zmatrix(zrxn, geo)
            elif backward_gra == automol.geom.graph(geo, stereo=False):
                zrxn = automol.reac.reverse(zrxn_i)
                zma, _, _ = automol.reac.ts_zmatrix(zrxn, geo)
        if zrxn is None:
            print(
                'Your geometry did not match any of the attempted' +
                'zrxns, which are the following')
            for zrxn_i in zrxns:
                print(zrxns)
            sys.exit()
        hess = elstruct.reader.hessian(prog, out_str)
        if hess is None:
            print(
                'No hessian found in output, cannot save ' +
                'a transition state without a hessian')
            sys.exit()
        run_path = insert_dct['run_path']
        if run_path is None:
            run_path = os.getcwd()
        run_fs = autofile.fs.run(run_path)
        freq_run_path = run_fs[-1].path(['hessian'])
        run_fs[-1].create(['hessian'])
        script_str = autorun.SCRIPT_DCT['projrot']
        freqs, _, imags, _ = autorun.projrot.frequencies(
            script_str, freq_run_path, [geo], [[]], [hess])
        if len(imags) != 1:
            print(
                'Can only save a transition state that has a single' +
                'imaginary frequency, projrot found the following' +
                'frequencies: ' + ','.join(imags))
            sys.exit()
    else:
        spc_info = parse_user_species(insert_dct)
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
    locs = parse_user_locs(insert_dct)

    # Check that the save location matches geo information
    if not insert_dct['saddle']:
        if not species_match(geo, spc_info):
            print(
                'I refuse to save this geometry until user specified' +
                ' info matches the info in user given output')
            sys.exit()
        # Check that the rid/cid info matches the filesystem
        fs_array, prefix_array = create_species_filesystems(
            prefix, spc_info, mod_thy_info, locs=None)
    else:
        fs_array, prefix_array = create_reaction_filesystems(
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
    ret = (inf_obj, inp_str, out_str)
    saved_locs, saved_geos, saved_enes = _saved_cnf_info(
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
                cinf_obj = cnf_fs[1].file.info.read(locs[0])
                cnsampd = cinf_obj.nsamp
                cnsampd += 1
                cinf_obj.nsamp = cnsampd
            else:
                cinf_obj = autofile.schema.info_objects.conformer_branch(0)
                cinf_obj.nsamp = 1
            cnf_fs[1].create([locs[0]])
            cnf_fs[0].file.info.write(rinf_obj)
            cnf_fs[1].file.info.write(cinf_obj, [locs[0]])
            if hess is not None and zrxn is not None:
                hess_inf_obj = autofile.schema.info_objects.run(
                    job=elstruct.Job.HESSIAN, prog=prog, version='',
                    method=method, basis=basis,
                     status=autofile.schema.RunStatus.SUCCESS)
                hess_ret = (hess_inf_obj, inp_str, out_str)
                save_saddle_point(
                    zrxn, ret, hess_ret, freqs, imags,
                    mod_thy_info, {'runlvl_cnf_fs': (cnf_fs, None)}, locs,
                    zma_locs=(0,), zma=zma)
            else:
                _save_unique_conformer(
                    ret, mod_thy_info, cnf_fs, locs,
                    zrxn=zrxn, zma_locs=(0,))
            print(
                'geometry is now saved at {}'.format(prefix_array[-1]))
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
        'ts_locs': None,
        'ts_mult': None,
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
            elif keyword in ['mult', 'charge']:
                values = []
                for val in value.split(','):
                    values.append(int(val))
                if len(values) == 1:
                    value = values[0]
                else:
                    value = values
            elif keyword in ['ts_locs']:
                value = (int(value),)
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
    SCRIPT_INPUT_FILE = 'insert_options.txt'
    insert_dct = parse_script_input(SCRIPT_INPUT_FILE)
    main(insert_dct)
