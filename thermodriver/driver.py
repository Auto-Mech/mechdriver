""" drivers for thermochemistry evaluations
"""
import os
from qcelemental import constants as qcc
import automol.inchi
import automol.geom
import scripts.es
import thermo.heatform
import esdriver.driver
import autofile.fs

REF_CALLS = {"basic": "get_basic",
             "cbh0": "get_cbhzed",
             "cbh1": "get_cbhone",
             "cbh2": "get_cbhtwo"}


def run(tsk_info_lst, es_dct, spcdct, spc_queue, ref, run_prefix, save_prefix, ene_coeff=[1.],
        options=[True, True, True, False]):
    """ main driver for thermo run
    """

    # Determine options
    runes = options[0]  #run electronic structure theory (True/False)
    runmess = options[1]  #run mess (True) / just make the mess input file (False)
    runthermo = options[2]
    if not runmess:
        runthermo = False
    everypf = options[3]  #make PF for every computation (True) / only with final info (False)

    # Fix any issues in tsk_list
    tsk_info_lst = fix(tsk_info_lst)
    # Add reference molecules
    for spc in spc_queue:
        if not 'ich' in spcdct[spc]:
            spcdct[spc]['ich'] = automol.smiles.inchi(spcdct[spc]['smiles'])
    refs, msg = prepare_refs(ref, spcdct, spc_queue)
    full_queue = spc_queue + refs
    full_queue = list(dict.fromkeys(full_queue))
    print(msg)

    # prepare filesystem
    if not os.path.exists(save_prefix):
        os.makedirs(save_prefix)
    if not os.path.exists(run_prefix):
        os.makedirs(run_prefix)

    # Run ESDriver
    if runes:
        runspecies = [{'species': full_queue, 'reacs': [], 'prods': []}]
        esdriver.driver.run(tsk_info_lst, es_dct, runspecies, spcdct, run_prefix, save_prefix)

    if runmess:
        geo_lvl = ''
        harm_lvl = ''
        anharm_lvl = ''
        tors_lvl = ''
        sym_lvl = ''
        harm_lvl_ref = ''
        anharm_lvl_ref = ''
        tors_lvl_ref = ''
        sym_lvl_ref = ''

        #Get PF input header
        temp_step = 100.
        ntemps = 30
        # temp_step = TEMP_STEP
        # ntemps = NTEMPS
        global_pf_str = scripts.thermo.get_pf_header(temp_step, ntemps)

        #Gather PF model and theory level info
        spc_model = ['RIGID', 'HARM', '']
        geom = False
        hess = False
        for tsk in tsk_info_lst:
            if 'samp' in tsk[0] or 'geom' in tsk[0]:
                geo_lvl = tsk[1]
                geom = True
            if 'grad' in tsk[0] or 'hess' in tsk[0]:
                harm_lvl = tsk[1]
                harm_lvl_ref = tsk[2]
                if 'grad' in tsk[0]:
                    grad = True
                if 'hess' in tsk[0]:
                    hess = True
                if not geom:
                    ene_lvl = tsk[1]
                    geo_lvl = tsk[1]
            if 'hr' in tsk[0] or 'tau' in tsk[0]:
                tors_lvl = tsk[1]
                tors_lvl_ref = tsk[2]
                if 'md' in tsk[0]:
                    spc_model[0] = 'MDHR'
                if 'tau' in tsk[0]:
                    spc_model[0] = 'TAU'
                else:
                    spc_model[0] = '1DHR'
            if 'anharm' in tsk[0] or 'vpt2' in tsk[0]:
                anharm_lvl = tsk[1]
                anharm_lvl_ref = tsk[2]
                spc_model[1] = 'ANHARM'
                if not hess:
                    geo_lvl = tsk[1]
            if 'sym' in tsk[0]:
                sym_lvl = tsk[1]
                sym_lvl_ref = tsk[2]
                if 'samp' in tsk[0]:
                    spc_model[2] = 'SAMPLING'
                if '1DHR' in tsk[0]:
                    spc_model[2] = '1DHR'
        geo_thy_info = get_thy_info(es_dct, geo_lvl)
        harm_thy_info = get_thy_info(es_dct, harm_lvl)
        tors_thy_info = get_thy_info(es_dct, tors_lvl)
        anharm_thy_info = get_thy_info(es_dct, anharm_lvl)
        sym_thy_info = get_thy_info(es_dct, sym_lvl)
        pf_levels = [harm_thy_info, tors_thy_info, anharm_thy_info, sym_thy_info]

        harm_ref_thy_info = get_thy_info(es_dct, harm_lvl_ref)
        tors_ref_thy_info = get_thy_info(es_dct, tors_lvl_ref)
        anharm_ref_thy_info = get_thy_info(es_dct, anharm_lvl_ref)
        sym_ref_thy_info = get_thy_info(es_dct, sym_lvl_ref)
        ref_levels = [
            harm_ref_thy_info, tors_ref_thy_info, anharm_ref_thy_info, sym_ref_thy_info]

        #Collect the PF input for each species
        # Initialize the ene for each of the species
        for spc in full_queue:
            spc_info = scripts.es.get_spc_info(spcdct[spc])
            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs.leaf.create(spc_info)
            spc_save_path = spc_save_fs.leaf.path(spc_info)

            zpe, zpe_str = scripts.thermo.get_zpe(
                spc, spcdct[spc], spc_save_path, pf_levels, spc_model)
            spc_str = scripts.thermo.get_spc_input(
                spc, spcdct[spc], spc_info, spc_save_path, pf_levels, spc_model)

            spcdct[spc]['spc_info'] = spc_info
            spcdct[spc]['spc_save_path'] = spc_save_path
            spcdct[spc]['zpe'] = zpe
            spcdct[spc]['zpe_str'] = zpe_str
            spcdct[spc]['spc_str'] = spc_str
            spcdct[spc]['ene'] = 0

        # Make and Run the PF file
    if runthermo:
        for spc in spc_queue:
            spc_save_path = spcdct[spc]['spc_save_path']
            spc_str = spcdct[spc]['spc_str']
            zpe_str = spcdct[spc]['zpe_str']
            spc_info = spcdct[spc]['spc_info']
            pf_input = scripts.thermo.get_pf_input(
                spc, spc_str, global_pf_str, zpe_str)

            pf_path, nasa_path = scripts.thermo.get_thermo_paths(
                spc_save_path, spc_info, harm_thy_info)
            spcdct[spc]['pf_path'] = pf_path
            spcdct[spc]['nasa_path'] = nasa_path

            scripts.thermo.write_pf_input(pf_input, pf_path)
            scripts.thermo.run_pf(pf_path)

        # Compute Hf0K
        ene_strl = []
        ene_str = '! energy level:'
        ene_lvl = ''
        ene_idx = 0
        for tsk in tsk_info_lst:
            if 'ene' in tsk[0]:
                if ene_idx > len(ene_coeff)-1:
                    print('WARNING:',
                          'an insufficient energy coefficient list was provided')
                    break
                ene_lvl = tsk[1]
                geo_lvl = tsk[2]
                geo_thy_info = scripts.es.get_thy_info(es_dct[geo_lvl])
                ene_thy_info = scripts.es.get_thy_info(es_dct[ene_lvl])
                ene_strl.append(' {:.2f} x {}{}/{}//{}{}/{}\n'.format(
                    ene_coeff[ene_idx], ene_thy_info[3],
                    ene_thy_info[1], ene_thy_info[2],
                    geo_thy_info[3], geo_thy_info[1], geo_thy_info[2]))

                for spc in full_queue:
                    spc_info = spcdct[spc]['spc_info']
                    ene = scripts.thermo.get_electronic_energy(
                        spc_info, geo_thy_info, ene_thy_info, save_prefix)
                    print('ene test:', ene, ene_coeff[ene_idx], ene_idx)
                    spcdct[spc]['ene'] += ene*ene_coeff[ene_idx]
                ene_idx += 1
        ene_str += '!               '.join(ene_strl)
        pf_levels.append(ene_str)
        # Need to get zpe for reference molecules too
        calc_bas = True
        if isinstance(ref, list):
            spc_bas = ref
            calc_bas = False
        elif is_scheme(ref) or not ref:
            reference_function = get_function_call(ref)
        for spc in spc_queue:
            if calc_bas:
                spc_bas, clist = get_ref(spc, spcdct, reference_function)
            hf0k = scripts.thermo.get_hf0k(spc, spcdct, spc_bas, clist)
            spcdct[spc]['Hfs'] = [hf0k]

        chemkin_header_str = scripts.thermo.run_ckin_header(
            pf_levels, ref_levels, spc_model)
        chemkin_set_str = chemkin_header_str
        for spc in spc_queue:
            pf_path = spcdct[spc]['pf_path']
            nasa_path = spcdct[spc]['nasa_path']
            starting_path = scripts.thermo.go_to_path(nasa_path)
            # need to change back to starting directory after running thermp and pac99 
            # or rest of code is confused
            scripts.thermo.write_thermp_inp(spcdct[spc])
            # run thermp creats thermo and also passed back the 298 K Hf
            if spcdct[spc]['ene'] == 0.0 or spcdct[spc]['spc_str'] == '':
                print('Cannot generate thermo for species {} because information is still missing:'.format(spcdct[spc]['ich']))
                continue
            hf298k = scripts.thermo.run_thermp(pf_path, nasa_path)
            spcdct[spc]['Hfs'].append(hf298k)
            pac99_poly_str = scripts.thermo.run_pac(spcdct[spc], nasa_path)
            chemkin_poly_str = scripts.thermo.run_ckin_poly(
                spc, spcdct[spc], pac99_poly_str)
            chemkin_spc_str = chemkin_header_str + chemkin_poly_str
            chemkin_set_str += chemkin_poly_str
            scripts.thermo.go_to_path(starting_path)
            ckin_path = scripts.thermo.prepare_path(starting_path, 'ckin')
            if not os.path.exists(ckin_path):
                os.makedirs(ckin_path)
            scripts.thermo.write_nasa_file(
                spcdct[spc], ckin_path, nasa_path, chemkin_spc_str)

        with open(os.path.join(ckin_path, 'automech.ckin'), 'w') as nasa_file:
            nasa_file.write(chemkin_set_str)


def is_scheme(entry):
    """ Check whether this is a basis set scheme
    """
    calls = REF_CALLS
    return entry in calls.keys()


def get_function_call(scheme):
    """ get function call
    """
    calls = REF_CALLS
    if not scheme:
        scheme = 'basic'
    call = getattr(thermo.heatform, calls[scheme])
    return call


def get_ref(spc, spcdct, call):
    """ references for a single species
    """
    ref, clist = call(spcdct[spc]['ich'])
    return ref, clist


def get_refs(species, spcs, scheme):
    """ references for a set of species
    """
    call = get_function_call(scheme)
    ref = []
    msg = ''
    for spc in species:
        spc_ref, _ = call(spcs[spc]['ich'])
        msg += 'Species {} with basis {}\n'.format(spc, ', '.join(spc_ref))
        ref.extend(spc_ref)
    return list(dict.fromkeys(ref)), msg


def prepare_refs(refscheme, spcdct, spc_queue):
    """ add refs to species list as necessary
    """
    if isinstance(refscheme, str):
        msg = 'Determining {} reference molecules for: \n'.format(refscheme)
        if is_scheme(refscheme):
            refs, newmsg = get_refs(spc_queue, spcdct, refscheme)
            msg += newmsg
    else:
        msg = 'Reference set = {}: \n'.format(', '.join(refscheme))
        refs = refscheme
    unique_refs = []
    for ref in refs:
        needtoadd = True
        for spc in spcdct:
            if spcdct[spc]['ich'] == ref:
                needtoadd = False
                unique_refs.append(spc)
                break
        if needtoadd:
            msg += 'Adding reference species ref_{}\n'.format(ref)
            spcdct['ref_' + ref] = create_spec(ref)
            unique_refs.append('ref_' + ref)
    return unique_refs, msg


def create_spec(val, charge=0, mc_nsamp=[True, 3, 1, 3, 100, 12], hind_inc=360.):
    """ add a species to the species dictionary
    """
    spec = {}
    if isinstance(val, str):
        ich = val
        print('ich test in create_spec:', ich)
        geo = automol.inchi.geometry(ich)
        zma = automol.geom.zmatrix(geo)
        spec['zmatrix'] = zma
    else:
        geo = val
        ich = automol.geom.inchi(geo)
    form = automol.inchi.formula_dct(ich)
    rad = automol.formula.electron_count(form) % 2
    if rad:
        mult = 2
    else:
        mult = 1
    spec['ich'] = ich
    spec['chg'] = charge
    spec['mul'] = mult
    spec['mc_nsamp'] = mc_nsamp
    spec['hind_inc'] = hind_inc * qcc.conversion_factor('degree', 'radian')
    return spec


def get_thy_info(es_dct, key):
    """ setup theory info file from es dictionary
    """
    ret = []
    if key:
        ret = scripts.es.get_thy_info(es_dct[key])
    return ret


def fix(tsk_info_lst):
    """ fix the tsk_info_list to make sure the necessary jobs precede it
    """
    has_sp = False
    last_geom = []
    for tsk in tsk_info_lst:
        if tsk[0] == 'conf_energy':
            has_sp = True
        if 'find' in tsk[0] or tsk[0] == 'opt':
            last_geom = ['conf_energy']
            last_geom.extend(tsk[1:-1])
            last_geom.append(False)
    if not has_sp:
        tsk_info_lst.append(last_geom)
    return tsk_info_lst


if __name__ == "__main__":

    MSG = """
          ================================================================
          ==                        AUTOMECHANIC                        ==
          ===         Andreas Copan, Sarah Elliott, Kevin Moore,       ===
          ===     Daniel Moberg, Carlo Cavallotti, Yuri Georgievski,   ===
          ==       Ahren Jasper, Murat Keceli, Stephen Klippenstein     ==
          ================================================================
          ==                        THERMODRIVER                        ==
          ===         Sarah Elliott, Kevin Moore, Andreas Copan,       ===
          ===    Murat Keceli, Yuri Georgievski, Stephen Klippenstein   ==
          ================================================================\n"""
    print(MSG)
    # tsk_info_lst, es_dct, spcs = load_params()
    REF = 'cbh0'

    TSK_INFO_LST = [
        ['conf_samp', 'mclev', 'mclev', False],
        ['geom', 'optlev', 'optlev', False],
        ['conf_hess', 'optlev', 'optlev', False]
        # [ 'hr', 'hrlev', 'optlev', False]
        # [ 'sp', 'splev', 'optlev', False]]
        ]

    ES_DCT = {'mclev': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp', 'basis': '6-31g*',
        'ncycles': 60, 'mem': 32, 'nprocs': 8, 'econv': '1.e-8', 'gconv':
        '1.e-4'},
#   es_dct = {'mclev': {
#       'orb_res': 'RU', 'program': 'psi4', 'method': 'b3lyp', 'basis': '6-31g*',
#       'ncycles': 60, 'mem': 32, 'nprocs': 8, 'econv': '1.e-8', 'gconv': '1.e-4'},
#             'optlev': {
#                 'orb_res': 'RU', 'program': 'psi4', 'method': 'b3lyp',
#                 'basis': 'cc-pvdz', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
#                 'econv': '1.e-8', 'gconv': '1.e-4'},
#             'testlvl': {
#                 'orb_res': 'RU', 'program': 'psi4', 'method': 'b3lyp',
#                 'basis': 'cc-pvdz', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
#                 'econv': '1.e-8', 'gconv': '1.e-4'},
              'optlev': {
                  'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp',
                  'basis': 'cc-pvdz', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
                  'econv': '1.e-8', 'gconv': '1.e-4'},
              'hrlev':  {
                  'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp',
                  'basis': '6-31g*', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
                  'econv': '1.e-8', 'gconv': '1.e-4'},
              'anlev':  {
                  'orb_res': 'RU', 'program': 'psi4', 'method': 'b3lyp',
                  'basis': 'cc-pvdz', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
                  'econv': '1.e-8', 'gconv': '1.e-4'},
              '2':      {
                  'orb_res': 'RU', 'program': 'molpro', 'method': 'ccsd(t)',
                  'basis': 'cc-pvtz', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
                  'econv': '1.e-8', 'gconv': '1.e-4'},
              'splev':  {
                  'orb_res': 'RU', 'program': 'molpro', 'method': 'b3lyp',
                  'basis': 'cc-pvtz', 'ncycles': 60, 'mem': 32, 'nprocs': 8,
                  'econv': '1.e-8', 'gconv': '1.e-4'}}


    # SPCS = ['prod1']
    SPCS = ['methanol', 'ethane', 'propane']

#    SPCS = ['reac1', 'prod1', 'methyl']
    SPCDCT = {'reac1': {'chg': 0, 'mul': 1, 'geom': ' n1 \n         h2 n1 R1 \n         h3 n1 R2 h2 A2 \n         h4 n1 R3 h2 A3 h3 D3\n         R1 = 1.01899\n         R2 = 1.01899\n         A2 = 105.997\n         R3 = 1.01899\n         A3 = 105.999\n         D3 = 112.362\n  ', 'mc_nsamp': [True, 2, 1, 2, 100, 12], 'hind_inc': 6.283185307179586, 'mc_tau': {}, 'elec_levels': [[1.0, 0.0]], 'geoobj': (('N', (0.0, 0.0, 0.0)), ('H', (0.0, 0.0, 1.9256120245802724)), ('H', (0.0, 1.851044869013889, -0.5306736870295058)), ('H', (1.7118265078919384, -0.704236134131304, -0.5307383003613204))), 'ich': 'InChI=1S/H3N/h1H3'},

              'reac2':  {'chg': 0, 'mul': 1, 'geom': '  N    0.00000   0.00000   0.00000\n          H    0.00000   0.00000   1.01521\n          H    0.96607   0.00000  -0.31225\n          H   -0.42878  -0.86571  -0.31222\n  ', 'mc_nsamp': [True, 2, 1, 2, 100, 12], 'mc_tau': {}, 'hind_inc': 6.283185307179586, 'elec_levels': [[1.0, 0.0]], 'sym_factor': 1.0, 'geoobj': (('N', (0.0, 0.0, 0.0)), ('H', (0.0, 0.0, 1.9184688598260415)), ('H', (1.825607718021044, 0.0, -0.5900669826742069)), ('H', (-0.8102767680738076, -1.6359548040700964, -0.5900102908904431))), 'ich': 'InChI=1S/H3N/h1H3'},

              'prod1':  {'chg': 0, 'mul': 1, 'smiles': 'OO ', 'mc_nsamp': [False, 0, 0, 0, 0, 30], 'mc_tau': {}, 'hind_inc': 6.283185307179586, 'sym_factor': 2.0, 'ich': 'InChI=1S/H2O2/c1-2/h1-2H', 'geeobj': (('O', (-1.0723807723844054, 0.8241985074848521, 0.4600038375719564)), ('O', (1.1007348339789824, -0.8517386452175768, 0.3222368540326177)), ('H', (-2.2842812994766244, -0.3263837511162607, -0.32081457818896136)), ('H', (2.2559272378820476, 0.3539238888489811, -0.46142611341561346))), 'geoobj': (('O', (-1.144941626087497, 0.45875729409195937, -0.7208009565935234)), ('O', (1.1479480801336641, 0.32164918390178865, 0.787186322659751)), ('H', (-2.2418713776640575, -0.4244113985652187, 0.47040915503683534)), ('H', (2.238864923617891, -0.3559950794285327, -0.5367945211030657)))},

              'methyl': {'chg': 0, 'mul': 2, 'smiles': '[CH3]', 'mc_nsamp': [False, 0, 0, 0, 0, 4], 'mc_tau': {}, 'hind_inc': 6.283185307179586, 'sym_factor': 2.0, 'ich': 'InChI=1S/CH3/h1H3', 'geeobj': (('C', (6.247770443277838e-08, 1.2570703433194236e-07, -1.997542166457577e-07)), ('H', (1.9945938085890018, -0.4544422128176104, -0.007625435387349655)), ('H', (-1.3906649408743608, -1.499876364681085, 0.03744539213515074)), ('H', (-0.6039289301923418, 1.9543184517916632, -0.029819756993584162))), 'geoobj': (('C', (-6.24560014375321e-08, -2.4767479707022037e-08, -1.7813105533709988e-07)), ('H', (2.0239350412358026, -0.29715278187532795, 0.019185549336825935)), ('H', (-1.2697746247054802, -1.6027535996541415, 0.061920704516777024)), ('H', (-0.7541603540743543, 1.8999064062969138, -0.08110607572254706)))},
              'hyd':    {'chg': 0, 'mul': 2, 'smiles': '[H]', 'mc_nsamp': [False, 0, 0, 0, 0, 1], 'mc_tau': {}, 'hind_inc': 6.283185307179586, 'sym_factor': 2.0, 'ich': 'InChI=1S/H', 'geeobj': (('H', (0.0, 0.0, 0.0)),), 'geoobj': (('H', (0.0, 0.0, 0.0)),)},
              'methane': {'chg': 0, 'mul': 1, 'smiles': 'C', 'mc_nsamp': [False, 0, 0, 0, 0, 4], 'mc_tau': {}, 'hind_inc': 6.283185307179586, 'sym_factor': 2.0, 'ich': 'InChI=1S/CH4/h1H4', 'geeobj': (('C', (3.618370664195086e-07, 1.1852373479820893e-06, 9.943538704651177e-08)), ('H', (0.8896538965581517, 1.4787504779842533, 1.132139058108076)), ('H', (1.4400903250506507, -0.9770743494241632, -1.109678541281473)), ('H', (-0.9404910045941374, -1.349173540207621, 1.247045907736483)), ('H', (-1.3892535788517275, 0.8474962264101816, -1.2695065239984664))), 'geoobj': (('C', (1.7277297690425342e-08, 1.060367771821586e-08, 1.8688543989814935e-08)), ('H', (-1.1244541609856977, 0.37449815237489303, 1.6897559734348577)), ('H', (-0.8875094828643466, -1.5061736249903073, -1.0971277444851157)), ('H', (0.11048582626380793, 1.7127341328773964, -1.1464101541617095)), ('H', (1.9014778003089403, -0.5810586708656641, 0.5537819065234211)))},

              'h2':     {'chg': 0, 'mul': 1, 'smiles': '[H][H]', 'mc_nsamp': [False, 0, 0, 0, 0, 1], 'mc_tau': {}, 'hind_inc': 6.283185307179586, 'sym_factor': 2.0, 'ich': 'InChI=1S/H2/h1H', 'geeobj': (('H', (0.6693776859397009, 0.0, 0.0)), ('H', (-0.6693776859397009, 0.0, 0.0))), 'geoobj': (('H', (0.6717316790925176, 0.0, 0.0)), ('H', (-0.6717316790925176, 0.0, 0.0)))},

              'methanol': {'chg': 0, 'mul': 1, 'smiles': 'CO', 'mc_nsamp': [False, 0, 0, 0, 0, 1]},
              'ethane': {'chg': 0, 'mul': 1, 'smiles': 'CC', 'mc_nsamp': [False, 0, 0, 0, 0, 10]}, 
              'propane': {'chg': 0, 'mul': 1, 'smiles': 'CCC', 'mc_nsamp': [False, 0, 0, 0, 0, 10]} 
              }

    run(TSK_INFO_LST, ES_DCT, SPCDCT, SPCS, REF, '/lcrc/project/PACC/run', '/lcrc/project/PACC/save')
    #run(tsk_info_lst, es_dct, spcdct, spcs, ref, 'runtest', 'savetest')
